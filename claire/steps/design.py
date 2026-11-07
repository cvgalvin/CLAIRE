"""Steps 6-12: matching, burial filter, enzyme design, HBRefine, scoring, ProteinMPNN and FastDesign.

Every step reads the previous stage's manifest and writes its own, so the chain runs without paths:
    claire match -> burial -> enzdes -> hbrefine -> score -> mpnn -> fastdesign
Rosetta flags are those used in the original workflow.py.
"""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import numpy as np

from .. import tools
from ..executor import Task, run_tasks
from ..project import ClaireError, Project
from .common import apply_filters, link_or_copy, match_decoy, pdbs, read_scores

CHUNK = 50  # structures per scoring task


# ============================================================ helpers
def motif_residues(name: str) -> list[int]:
    """Motif residue numbers from a matcher-style name 'UM_<n>_<D45N84L97>_...' (as in match_resfile2.py)."""
    try:
        code = Path(name).name.split("_")[2]
    except IndexError:
        return []
    nums, cur = [], ""
    for ch in code:
        if ch.isdigit():
            cur += ch
        else:
            if cur:
                nums.append(int(cur))
            cur = ""
    if cur:
        nums.append(int(cur))
    return nums


def _atoms(pdb: Path):
    for l in Path(pdb).read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            yield l


def binding_site(pdb: Path, lig: str, cutoff: float) -> set[tuple[str, int]]:
    lig_xyz, prot = [], {}
    for l in _atoms(pdb):
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if l[17:20].strip().upper() == lig.upper():
            lig_xyz.append(xyz)
        elif l.startswith("ATOM") and l[76:78].strip() != "H" and not l[12:16].strip().startswith("H"):
            prot.setdefault((l[21], int(l[22:26])), []).append(xyz)
    L = np.array(lig_xyz)
    if len(L) == 0:
        raise ClaireError(f"Ligand {lig} not found in {pdb}")
    return {k for k, v in prot.items() if np.min(np.linalg.norm(np.array(v)[:, None] - L[None], axis=2)) <= cutoff}


def _score_tasks(proj: Project, stage: str, inputs: list[Path], xml: str, extra: list[str], tag: str) -> list[Task]:
    d = proj.stage(stage)
    (d / "lists").mkdir(exist_ok=True)
    (d / "scores").mkdir(exist_ok=True)
    tasks = []
    for i in range(0, len(inputs), CHUNK):
        k = i // CHUNK + 1
        lst = d / "lists" / f"{tag}_{k}.list"
        lst.write_text("\n".join(str(p) for p in inputs[i:i + CHUNK]) + "\n")
        cmd = [tools.rosetta(proj, "rosetta_scripts"), "-l", str(lst), "-parser:protocol", tools.xml(xml),
               "-parser:view", "-run:preserve_header", "-in:file::extra_res_fa", str(proj.params()),
               "-load_PDB_components", "False", "-scorefile_format", "json",
               "-out:file:score_only", str(d / "scores" / f"{tag}_{k}.json"), *extra]
        tasks.append(Task(" ".join(cmd), str(d), f"{tag}_{k}"))
    return tasks


def _collect_scores(d: Path, tag: str) -> list[dict]:
    return read_scores(sorted((d / "scores").glob(f"{tag}_*.json")))


def _by_stem(paths: list[Path]) -> dict:
    return {p.stem: p for p in paths}


# ============================================================ 6. match
def auto_pos(pdb: Path, cutoff: float) -> list[int]:
    """Heuristic matcher positions: residues whose CA->CB vector points into the protein core within `cutoff`."""
    ca, cb = {}, {}
    for l in _atoms(pdb):
        if not l.startswith("ATOM"):
            continue
        r = int(l[22:26])
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if l[12:16].strip() == "CA":
            ca[r] = xyz
        elif l[12:16].strip() == "CB":
            cb[r] = xyz
    centre = np.mean(list(ca.values()), axis=0)
    pos = []
    for r, a in ca.items():
        b = cb.get(r)
        if b is None:
            continue
        v, w = b - a, centre - a
        if np.linalg.norm(b - centre) <= cutoff and np.dot(v, w) / (np.linalg.norm(v) * np.linalg.norm(w)) > 0.3:
            pos.append(r)
    return sorted(pos)


def match(proj: Project, args) -> None:
    d = proj.stage("match")
    if getattr(args, "collect", False):
        return _match_collect(proj, d)
    sdir = proj.path(proj.config["scaffolds"]["dir"])
    scaffolds = pdbs(sdir)
    if not scaffolds:
        raise ClaireError(f"No scaffold PDBs in {sdir}")
    proj.manifest("motifs")
    mdir = proj.stage("motifs", False)
    sets = [s for s in proj.config["match"]["motif_sets"] if list((mdir / s).glob("*.cst"))]
    if not sets:
        raise ClaireError("No motif constraint files found; run `claire motifs`.")
    posdir = d / "pos"
    posdir.mkdir(exist_ok=True)
    n_auto = 0
    for s in scaffolds:
        user = s.with_suffix(".pos")
        dst = posdir / f"{s.stem}.pos"
        if user.exists():
            shutil.copy(user, dst)
        elif not dst.exists():
            pos = auto_pos(s, float(proj.config["scaffolds"]["auto_pos_cb_cutoff"]))
            if not pos:
                raise ClaireError(f"No core-facing positions found for {s.name}; provide {s.stem}.pos or increase "
                                  "scaffolds.auto_pos_cb_cutoff.")
            dst.write_text(" ".join(map(str, pos)) + "\n")
            n_auto += 1
    if n_auto:
        print(f"[match] generated {n_auto} .pos files heuristically (core-facing residues); "
              "provide <scaffold>.pos files to control matcher positions.")
    rosetta_match = tools.rosetta(proj, "match")
    tasks = []
    for st in sets:
        csts = sorted((mdir / st).glob("*.cst"))
        for s in scaffolds:
            out = d / st / s.stem
            out.mkdir(parents=True, exist_ok=True)
            loop = " ".join(str(c) for c in csts)
            cmd = (f"for c in {loop}; do n=$(basename $c .cst); mkdir -p $n && cd $n && "
                   f"{rosetta_match} -s {s} -extra_res_fa {proj.params()} -match:geometric_constraint_file $c "
                   f"-match::scaffold_active_site_residues {posdir / (s.stem + '.pos')} -match:output_format PDB "
                   f"-match:match_grouper SameSequenceGrouper -match:consolidate_matches -match:output_matches_per_group 1 "
                   f"-use_input_sc -ex1 -ex2 -extrachi_cutoff 0 -enumerate_ligand_rotamers false "
                   f"-match::lig_name {proj.lig} -load_PDB_components false > match.log 2>&1; cd ..; done")
            tasks.append(Task(cmd, str(out), f"match_{st}_{s.stem}"))
    if run_tasks(proj, "match", tasks, label="match"):
        _match_collect(proj, d)


def _rosetta_errors(logs: list[Path], limit: int = 5) -> dict:
    """Unique Rosetta error messages found in log files -> count."""
    found = {}
    for log in logs:
        try:
            lines = log.read_text(errors="replace").splitlines()
        except OSError:
            continue
        for i, l in enumerate(lines):
            if ("ERROR" in l or "Assertion" in l or "terminate called" in l or "Library not loaded" in l
                    or "cannot open shared object" in l or "image not found" in l or "Traceback" in l
                    or "Error:" in l or "No such file" in l or "command not found" in l):
                msg = " ".join(x.strip() for x in lines[i:i + 3] if x.strip())[:300]
                found[msg] = found.get(msg, 0) + 1
                break
    return dict(sorted(found.items(), key=lambda kv: -kv[1])[:limit])


def explain_empty(proj: Project, stage: str, what: str) -> None:
    """Called when a step produced nothing: print errors found in its logs and where to look."""
    d = proj.stage(stage, False)
    logs = sorted(list(d.glob("logs/*.log")) + list(d.glob("jobs/logs/*")) + list(d.glob("**/*.resfile.log"))
                  + list(d.glob("**/match.log")))
    errs = _rosetta_errors(logs, limit=4)
    print(f"[{stage}] no {what} were produced.")
    if errs:
        print("  Errors found in the logs:")
        for msg, c in errs.items():
            print(f"   ({c}x) {msg}")
    else:
        tails = [l for l in logs if l.stat().st_size > 0][:1]
        for l in tails:
            last = [x for x in l.read_text(errors="replace").splitlines() if x.strip()][-8:]
            print(f"  No explicit error found; last lines of {l.relative_to(proj.root)}:")
            for x in last:
                print(f"    {x[:200]}")
    print(f"  Logs: {d.relative_to(proj.root)}/logs/ (local) or jobs/logs/ (cluster)")


def _match_collect(proj: Project, d: Path) -> None:
    logs = sorted(d.glob("*/*/*/match.log"))
    errs = _rosetta_errors(logs)
    if any("Library not loaded" in m or "shared object" in m for m in errs):
        print("[match] Rosetta could not load its shared libraries. This usually means the Rosetta folder was moved or "
              "renamed after compiling; symlink the original path to the new location (see `claire check`).")
    if errs:
        n = sum(errs.values())
        print(f"[match] Rosetta reported errors in {n} of {len(logs)} matcher runs. Most common:")
        for msg, c in errs.items():
            print(f"  ({c}x) {msg}")
        print("  Full logs: 06_matches/<motif set>/<scaffold>/<motif>/match.log")
    mdir = proj.stage("motifs", False)
    records = []
    for p in sorted(d.glob("*/*/*/UM_*.pdb")):
        st, cst = p.parts[-4], p.parts[-2]
        records.append({"pdb": str(p.relative_to(proj.root)), "cst": str((mdir / st / f"{cst}.cst").relative_to(proj.root)),
                        "motif_set": st})
    proj.write_manifest("match", {"matches": records, "n_outputs": len(records)})
    if not records:
        explain_empty(proj, "match", "matches")
        return
    print(f"[match] {len(records)} matches. Next: `claire burial`.")


# ============================================================ 7. burial
def burial(proj: Project, args) -> None:
    d = proj.stage("burial")
    recs = proj.manifest("match")["matches"]
    if not getattr(args, "collect", False):
        inputs = [proj.root / r["pdb"] for r in recs]
        if not run_tasks(proj, "burial", _score_tasks(proj, "burial", inputs, "match_analysis_standard.xml", [], "burial"),
                         label="burial scoring"):
            return
    scores = _collect_scores(d, "burial")
    if not scores:
        explain_empty(proj, "burial", "burial scores")
    passed = apply_filters(scores, [["dsasa", ">=", proj.config["burial"]["dsasa_min"]]], "burial")
    stems = _by_stem([proj.root / r["pdb"] for r in recs])
    keep = {match_decoy(x["decoy"], stems).stem for x in passed if match_decoy(x["decoy"], stems) is not None}
    out = d / "buried"
    out.mkdir(exist_ok=True)
    kept = []
    for r in recs:
        p = proj.root / r["pdb"]
        stem = f"{p.stem}__{r['motif_set']}__{Path(r['cst']).stem}"
        if p.stem in keep:
            # rename so that the motif residues stay at name.split('_')[2] and the source is unique
            link_or_copy(p, out / f"{stem}.pdb")
            shutil.copy(proj.root / r["cst"], out / f"{stem}.cst")
            kept.append(stem)
    proj.write_manifest("burial", {"designs": kept, "dir": "07_buried/buried", "n_outputs": len(kept)})
    print(f"[burial] {len(kept)} well-buried matches in {out}. Next: `claire enzdes`.")


# ============================================================ 8. enzyme design
def enzdes(proj: Project, args) -> None:
    """Binding-site design, then (enzdes.filter) fast-metric filters from the Methods before HBRefine."""
    d = proj.stage("enzdes")
    st_file = d / "state.json"
    st = json.loads(st_file.read_text()) if (getattr(args, "collect", False) and st_file.exists()) else {"phase": "design"}
    cfg = proj.config["enzdes"]
    while True:
        if st["phase"] == "design":
            src = proj.stage("burial", False) / "buried"
            inp = d / "inputs"
            inp.mkdir(exist_ok=True)
            (d / "designs").mkdir(exist_ok=True)
            tasks = []
            for stem in proj.manifest("burial")["designs"]:
                link_or_copy(src / f"{stem}.pdb", inp / f"{stem}.pdb")
                shutil.copy(src / f"{stem}.cst", inp / f"{stem}.cst")
                cmd = (f"{tools.python(proj)} {tools.legacy('match_resfile2.py')} {stem}.pdb {proj.params()} {cfg['design_shell']} no "
                       f"> {stem}.resfile.log 2>&1 && "
                       f"{tools.rosetta(proj, 'enzyme_design')} -s {stem}.pdb -in:file::extra_res_fa {proj.params()} "
                       f"-resfile {stem}.resfile -ex1 -ex2 -extrachi_cutoff 0 -enzdes:cst_opt -enzdes:bb_min -enzdes:chi_min "
                       f"-enzdes:cst_design -enzdes:design_min_cycles 3 -enzdes:lig_packer_weight 2.5 -enzdes:cst_min "
                       f"-out:nstruct {cfg['nstruct']} -load_PDB_components False -enzdes:cstfile {stem}.cst "
                       f"-out:path:all {d / 'designs'} -out:path:pdb {d / 'designs'} -overwrite")
                tasks.append(Task(cmd, str(inp), f"enzdes_{stem}"))
            ok = run_tasks(proj, "enzdes", tasks, label="enzyme design")
            st["phase"] = "score" if cfg.get("filter", True) else "done"
        elif st["phase"] == "score":
            ok = run_tasks(proj, "enzdes", _score_tasks(proj, "enzdes", pdbs(d / "designs"), "edanalysis_std_sfn.xml",
                                                        ["-corrections::beta_nov16"], "enzdes"), label="fast metrics")
            st["phase"] = "done"
        else:
            # enzyme_design names outputs <input>__DE_<n>.pdb; collect them from the output folder and, if Rosetta
            # ignored the output path, from the inputs folder (input matches are never touched)
            for p in (d / "inputs").glob("*__DE_*.pdb"):
                shutil.move(str(p), d / "designs" / p.name)
            designs = sorted((d / "designs").glob("*__DE_*.pdb"))
            out = d / "passed"
            shutil.rmtree(out, ignore_errors=True)
            out.mkdir()
            if cfg.get("filter", True):
                stems = _by_stem(designs)
                kept = apply_filters(_collect_scores(d, "enzdes"), proj.config["filters"]["enzdes"], "enzdes")
                designs = [p for p in (match_decoy(x["decoy"], stems) for x in kept) if p is not None]
            for p in designs:
                link_or_copy(p, out / p.name)
            proj.write_manifest("enzdes", {"dir": "08_enzdes/passed", "designs": sorted(p.name for p in designs),
                                           "n_outputs": len(designs)})
            if not designs:
                explain_empty(proj, "enzdes", "enzyme_design outputs")
            print(f"[enzdes] {len(designs)} designs passed to HBRefine. Next: `claire hbrefine`.")
            st_file.write_text(json.dumps({"phase": "complete"}))
            return
        st_file.write_text(json.dumps(st))
        if not ok:
            return


# ============================================================ 9. HBRefine
def hbrefine(proj: Project, args) -> None:
    d = proj.stage("hbrefine")
    raw = d / "raw"
    if not getattr(args, "collect", False):
        src = proj.stage("enzdes", False) / "passed"
        raw.mkdir(exist_ok=True)
        tasks = [Task(f"{tools.python(proj)} {tools.legacy('hbrefine.py')} {proj.params()} {src / n} {raw / Path(n).stem}",
                      str(d), f"hbrefine_{Path(n).stem}") for n in proj.manifest("enzdes")["designs"]]
        if not run_tasks(proj, "hbrefine", tasks, label="HBRefine"):
            return
    out = d / "refined"
    out.mkdir(exist_ok=True)
    kept = []
    # hbrefine.py writes <name>_unrefined.pdb when the ligand had no unsatisfied polar atoms to begin with, and
    # <name>_<n>mutations_fixed.pdb when refinement satisfied them; both go forward (unsats/ubp/bothproblems do not)
    n_fixed = n_ok = 0
    for p in list(raw.glob("*/*_fixed.pdb")) + list(raw.glob("*/*_unrefined.pdb")):
        link_or_copy(p, out / p.name)
        kept.append(p.name)
        n_fixed += p.name.endswith("_fixed.pdb")
        n_ok += p.name.endswith("_unrefined.pdb")
    proj.write_manifest("hbrefine", {"dir": "09_hbrefine/refined", "designs": sorted(kept), "n_outputs": len(kept),
                                     "already_satisfied": n_ok, "fixed_by_hbrefine": n_fixed})
    if not list(raw.glob("*/*.pdb")):
        explain_empty(proj, "hbrefine", "HBRefine outputs")
    print(f"[hbrefine] {len(kept)} designs go forward: {n_ok} had no unsatisfied ligand polar atoms, "
          f"{n_fixed} were fixed by HBRefine. Next: `claire score`.")


# ============================================================ 10. scoring + binding filters
FULL_FLAGS = ("-ignore_unrecognized_res -run:nblist_autoupdate -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ "
              "-packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 "
              "-corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 "
              "-corrections:score:lj_hbond_OH_donor_dis=2.6 -corrections:beta_nov16 -enzdes:bb_min_allowed_dev=0.05 "
              "-enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10").split()


def score(proj: Project, args) -> None:
    d = proj.stage("score")
    src = proj.stage("hbrefine", False) / "refined"
    designs = [src / n for n in proj.manifest("hbrefine")["designs"]]
    st_file = d / "state.json"
    st = json.loads(st_file.read_text()) if (getattr(args, "collect", False) and st_file.exists()) else {"phase": "fast"}
    full = bool(proj.config["score"]["full_metrics"])
    while True:
        if st["phase"] == "fast":
            ok = run_tasks(proj, "score", _score_tasks(proj, "score", designs, "edanalysis_std_sfn.xml",
                                                       ["-corrections::beta_nov16"], "fast"), label="fast metrics")
            st["phase"] = "full" if full else "filter"
        elif st["phase"] == "full":
            fast = _collect_scores(d, "fast")
            pre = apply_filters(fast, proj.config["filters"]["binding"], "binding (fast metrics)", skip_missing=True)
            stems = _by_stem(designs)
            todo = [p for p in (match_decoy(x["decoy"], stems) for x in pre) if p is not None]
            ok = run_tasks(proj, "score", _score_tasks(proj, "score", todo, "ligbinderanalysis_jump1.xml", FULL_FLAGS, "full"),
                           label="full metrics")
            st["phase"] = "filter"
        else:
            scores = _collect_scores(d, "full") if full else _collect_scores(d, "fast")
            if full:   # merge fast metrics into full records
                fast = {s["decoy"]: s for s in _collect_scores(d, "fast")}
                scores = [{**fast.get(s["decoy"], {}), **s} for s in scores]
            kept = apply_filters(scores, proj.config["filters"]["binding"], "binding", skip_missing=not full)
            out = d / "passed"
            out.mkdir(exist_ok=True)
            names = []
            stems = _by_stem(designs)
            for s in kept:
                p = match_decoy(s["decoy"], stems)
                if p is not None:
                    link_or_copy(p, out / p.name)
                    names.append(p.name)
            (d / "passed_scores.json").write_text(json.dumps(kept, indent=0))
            proj.write_manifest("score", {"dir": "10_scored/passed", "designs": sorted(names), "n_outputs": len(names)})
            print(f"[score] {len(names)} designs pass binding filters. Next: `claire mpnn`.")
            st["phase"] = "done"
            ok = True
        st_file.write_text(json.dumps(st))
        if not ok or st["phase"] == "done":
            return


# ============================================================ 11. ProteinMPNN + sequence-profile resfiles
def _chain_residues(pdb: Path, lig: str):
    seen = []
    for l in _atoms(pdb):
        if l.startswith("ATOM") and l[17:20].strip().upper() != lig.upper():
            k = (l[21], int(l[22:26]))
            if not seen or seen[-1] != k:
                seen.append(k)
    return seen


def mpnn(proj: Project, args) -> None:
    d = proj.stage("mpnn")
    src = proj.stage("score", False) / "passed"
    designs = [src / n for n in proj.manifest("score")["designs"]]
    cfg = dict(proj.config["mpnn"])
    if getattr(args, "cutoff", None) is not None:
        cfg["binding_site_cutoff"] = args.cutoff
    print(f"[mpnn] fixing motif residues and residues within {cfg['binding_site_cutoff']} A of the ligand")
    inp, out = d / "inputs", d / "output"
    if not getattr(args, "collect", False):
        shutil.rmtree(inp, ignore_errors=True)
        inp.mkdir()
        out.mkdir(exist_ok=True)
        fixed = {}
        for p in designs:
            link_or_copy(p, inp / p.name)
            site = binding_site(p, proj.lig, float(cfg["binding_site_cutoff"]))
            motif = set(motif_residues(p.name))
            res = _chain_residues(p, proj.lig)
            chains = {}
            for i, (ch, num) in enumerate(res):
                chains.setdefault(ch, [])
            for ch in chains:
                idx = [i + 1 for i, k in enumerate([r for r in res if r[0] == ch]) if k in site or k[1] in motif]
                chains[ch] = idx                       # ProteinMPNN fixed positions are 1-based per chain
            fixed[p.stem] = chains
        (d / "fixed_pdbs.jsonl").write_text(json.dumps(fixed) + "\n")
        kind, mp = tools.mpnn(proj)
        py = tools.python(proj)
        if kind == "proteinmpnn":
            if cfg.get("model_type", "protein_mpnn") != "protein_mpnn":
                raise ClaireError("mpnn.model_type ligand_mpnn needs a LigandMPNN checkout (tools.proteinmpnn_dir).")
            helper = mp / "helper_scripts" / "parse_multiple_chains.py"
            cmd = (f"{py} {helper} --input_path={inp} --output_path={d / 'parsed_pdbs.jsonl'} && "
                   f"{py} {mp / 'protein_mpnn_run.py'} --jsonl_path {d / 'parsed_pdbs.jsonl'} --out_folder {out} "
                   f"--fixed_positions_jsonl {d / 'fixed_pdbs.jsonl'} --num_seq_per_target {cfg['num_seq']} "
                   f"--sampling_temp \"{cfg['temperature']}\" --batch_size 1")
        else:
            model = cfg.get("model_type", "protein_mpnn")
            flag, ck = tools.ligandmpnn_checkpoint(proj, mp, model)
            # LigandMPNN takes fixed residues as "A14 A17 ..." (chain + PDB residue number) per input file
            pdb_multi, fixed_multi = {}, {}
            for p in designs:
                key = str(inp / p.name)
                pdb_multi[key] = ""
                res = [r for r in _chain_residues(p, proj.lig)]
                tok = []
                for ch, idx in fixed[p.stem].items():
                    chain_res = [r for r in res if r[0] == ch]
                    tok += [f"{ch}{chain_res[i - 1][1]}" for i in idx]
                fixed_multi[key] = " ".join(tok)
            (d / "pdb_paths.json").write_text(json.dumps(pdb_multi, indent=0))
            (d / "fixed_residues.json").write_text(json.dumps(fixed_multi, indent=0))
            n = int(cfg["num_seq"])
            batch = min(10, n)
            cmd = (f"{py} {mp / 'run.py'} --model_type {model} --{flag} {ck} "
                   f"--pdb_path_multi {d / 'pdb_paths.json'} --fixed_residues_multi {d / 'fixed_residues.json'} "
                   f"--out_folder {out} --temperature {cfg['temperature']} --batch_size {batch} "
                   f"--number_of_batches {-(-n // batch)}")
        if not run_tasks(proj, "mpnn", [Task(cmd, str(mp), "mpnn")], label=f"{kind} ({cfg.get('model_type', 'protein_mpnn')})"):
            return
    # sequence profile -> resfile allowing only residues sampled by MPNN at designable positions
    fixed = json.loads((d / "fixed_pdbs.jsonl").read_text())
    rdir = d / "resfiles"
    rdir.mkdir(exist_ok=True)
    made = []
    for p in designs:
        fa = out / "seqs" / f"{p.stem}.fa"
        if not fa.exists():
            continue
        seqs = [l.strip().replace(":", "").replace("/", "") for l in fa.read_text().splitlines()
                if l.strip() and not l.startswith(">")][1:]                                     # first record = input
        res = [r for r in _chain_residues(p, proj.lig)]
        lines = ["NATAA", "start"]
        for ch, fx in fixed[p.stem].items():
            chain_res = [r for r in res if r[0] == ch]
            for i, (c, num) in enumerate(chain_res, 1):
                if i in fx:
                    continue
                allowed = sorted({s[i - 1] for s in seqs if len(s) >= i and s[i - 1].isalpha()})
                if allowed:
                    lines.append(f"{num} {c} PIKAA {''.join(allowed)}")
        (rdir / f"{p.stem}.resfile").write_text("\n".join(lines) + "\n")
        made.append(p.name)
    proj.write_manifest("mpnn", {"designs": made, "resfiles": "11_mpnn/resfiles", "n_outputs": len(made)})
    print(f"[mpnn] sequence-profile resfiles for {len(made)} designs. Next: `claire fastdesign`.")


# ============================================================ 12. FastDesign with MPNN profile
def fastdesign(proj: Project, args) -> None:
    d = proj.stage("fastdesign")
    src = proj.stage("score", False) / "passed"
    rdir = proj.stage("mpnn", False) / "resfiles"
    if not getattr(args, "collect", False):
        (d / "designs").mkdir(exist_ok=True)
        (d / "scores").mkdir(exist_ok=True)
        tasks = []
        for n in proj.manifest("mpnn")["designs"]:
            stem = Path(n).stem
            cmd = [tools.rosetta(proj, "rosetta_scripts"), "-s", str(src / n), "-parser:protocol", tools.xml("fd_mpnn.xml"),
                   "-parser:view", "-run:preserve_header", "-in:file::extra_res_fa", str(proj.params()),
                   "-resfile", str(rdir / f"{stem}.resfile"), "-ex1", "-ex2", "-extrachi_cutoff", "0", "-use_input_sc",
                   "-out:nstruct", str(proj.config["fastdesign"]["nstruct"]), "-beta_nov16", "-load_PDB_components", "False",
                   "-out:path:all", str(d / "designs"), "-scorefile_format", "json",
                   "-out:file:scorefile", str(d / "scores" / f"fd_{stem}.json")]
            tasks.append(Task(" ".join(cmd), str(d), f"fd_{stem}"))
        if not run_tasks(proj, "fastdesign", tasks, label="FastDesign"):
            return
    scores = read_scores(sorted((d / "scores").glob("fd_*.json")))
    if not scores:
        explain_empty(proj, "fastdesign", "FastDesign scores")
    kept = apply_filters(scores, proj.config["filters"]["fastdesign"], "fastdesign")
    out = d / "passed"
    out.mkdir(exist_ok=True)
    names = []
    for s in kept:
        p = d / "designs" / (s["decoy"] if s["decoy"].endswith(".pdb") else s["decoy"] + ".pdb")
        if p.exists():
            link_or_copy(p, out / p.name)
            names.append(p.name)
    (d / "passed_scores.json").write_text(json.dumps(kept, indent=0))
    proj.write_manifest("fastdesign", {"dir": "12_fastdesign/passed", "designs": sorted(names), "n_outputs": len(names)})
    print(f"[fastdesign] {len(names)} designs pass stability/binding filters. Next: `claire predict`.")
