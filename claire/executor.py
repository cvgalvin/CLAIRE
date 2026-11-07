"""Run a list of shell tasks locally (in parallel) or as a SLURM/SGE array job.

Every heavy step builds a list of :class:`Task` objects. With the local backend they run immediately;
with a cluster backend a task list and array-job script are written to ``<stage>/jobs/`` (and optionally
submitted), and the step is completed later with ``claire <step> --collect``.
"""
from __future__ import annotations

import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from .project import ClaireError, Project


@dataclass
class Task:
    cmd: str            # shell command
    cwd: str            # working directory
    name: str           # short label, used for log file names


def _run_one(t: Task, logdir: Path) -> tuple[str, int]:
    log = logdir / f"{t.name}.log"
    with open(log, "w") as fh:
        fh.write(f"# cwd: {t.cwd}\n# cmd: {t.cmd}\n")
        fh.flush()
        rc = subprocess.call(t.cmd, shell=True, cwd=t.cwd, stdout=fh, stderr=subprocess.STDOUT,
                             executable="/bin/bash")
    return t.name, rc


def run_tasks(proj: Project, stage: str, tasks: list[Task], label: str = "tasks") -> bool:
    """Run tasks. Returns True if they completed now (local), False if submitted/written for a cluster."""
    if not tasks:
        return True
    cfg = proj.config["compute"]
    backend = cfg.get("backend", "local")
    sdir = proj.stage(stage)
    if backend == "local":
        logdir = sdir / "logs"
        logdir.mkdir(exist_ok=True)
        failed = []
        workers = max(1, int(cfg.get("workers", 1)))
        print(f"[{stage}] running {len(tasks)} {label} with {workers} local worker(s); logs in {logdir}")
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_run_one, t, logdir) for t in tasks]
            for i, f in enumerate(as_completed(futs), 1):
                name, rc = f.result()
                if rc != 0:
                    failed.append(name)
                if i % max(1, len(tasks) // 20) == 0 or i == len(tasks):
                    print(f"  {i}/{len(tasks)} done", flush=True)
        if failed:
            print(f"[{stage}] WARNING: {len(failed)} task(s) failed, e.g. {failed[:5]} (see {logdir})")
        return True
    if backend not in ("slurm", "sge"):
        raise ClaireError(f"Unknown compute.backend '{backend}' (use local, slurm or sge).")
    jobs = sdir / "jobs"
    (jobs / "logs").mkdir(parents=True, exist_ok=True)
    tfile = jobs / f"{label}.tasks"
    tfile.write_text("".join(f"cd {shlex.quote(t.cwd)} && {t.cmd}\n" for t in tasks))
    setup = cfg.get("setup", "") or ""
    n = len(tasks)
    if backend == "slurm":
        s = cfg.get("slurm", {})
        hdr = [f"#SBATCH --job-name=claire_{stage}", f"#SBATCH --array=1-{n}", f"#SBATCH --time={s.get('time', '24:00:00')}",
               f"#SBATCH --mem={s.get('mem', '8G')}", f"#SBATCH --output={jobs}/logs/%a.log"]
        if s.get("partition"):
            hdr.append(f"#SBATCH --partition={s['partition']}")
        if s.get("extra"):
            hdr.append(s["extra"])
        idx = "${SLURM_ARRAY_TASK_ID}"
        submit = ["sbatch"]
    else:
        s = cfg.get("sge", {})
        hdr = [f"#$ -N claire_{stage}", f"#$ -t 1-{n}", "#$ -cwd", f"#$ -l mem_free={s.get('mem', '8G')}",
               f"#$ -l h_rt={s.get('time', '24:00:00')}", f"#$ -o {jobs}/logs", "#$ -j y"]
        if s.get("extra"):
            hdr.append(s["extra"])
        idx = "${SGE_TASK_ID}"
        submit = ["qsub"]
    script = jobs / f"{label}.sh"
    script.write_text("#!/bin/bash\n" + "\n".join(hdr) + f"\n{setup}\n"
                      f'CMD=$(sed -n "{idx}p" {tfile})\necho "$CMD"\neval "$CMD"\n')
    script.chmod(0o755)
    if cfg.get("submit"):
        out = subprocess.run(submit + [str(script)], capture_output=True, text=True)
        print(out.stdout.strip() or out.stderr.strip())
        print(f"[{stage}] submitted {n} {label}. When finished, run `claire {stage} --collect`.")
    else:
        print(f"[{stage}] wrote {script} ({n} array tasks). Submit with `{submit[0]} {script}`, "
              f"then run `claire {stage} --collect` when the jobs finish.")
    return False
