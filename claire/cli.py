"""`claire` command-line interface: one subcommand per workflow step."""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .project import STAGES, ClaireError, Project, init_project
from .steps import cluster, contacts, design, fragments, motifs, predict, target

STEPS = {
    "target": (target.run, "Rosetta params + model for the ligand (01_target)"),
    "fragments": (fragments.run, "build fragments and find PDB entries containing them (02_fragments)"),
    "contacts": (contacts.run, "download structures, extract/filter contacts, split by residue (03_contacts)"),
    "cluster": (cluster.run, "spatially cluster contacts; screen H-bonding clusters (04_clusters)"),
    "motifs": (motifs.run, "assemble polar + nonpolar binding motifs as match constraint files (05_motifs)"),
    "match": (design.match, "RosettaMatch motifs into scaffolds (06_matches)"),
    "burial": (design.burial, "keep matches with buried ligand (07_buried)"),
    "enzdes": (design.enzdes, "binding-site design with Rosetta enzyme_design (08_enzdes)"),
    "hbrefine": (design.hbrefine, "hydrogen-bond refinement of binding sites (09_hbrefine)"),
    "score": (design.score, "binding metrics and filters (10_scored)"),
    "mpnn": (design.mpnn, "ProteinMPNN outside the binding site -> sequence-profile resfiles (11_mpnn)"),
    "fastdesign": (design.fastdesign, "FastDesign with MPNN profile + stability filters (12_fastdesign)"),
    "predict": (predict.predict, "ColabFold single-sequence prediction (13_colabfold)"),
    "final": (predict.final, "select designs by pLDDT and Ca RMSD (14_final)"),
}
ORDER = list(STEPS)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="claire", description="CLAIRE: combinatorial design of small-molecule binding proteins")
    ap.add_argument("--version", action="version", version=f"claire {__version__}")
    ap.add_argument("-p", "--project", help="project directory (default: search upwards from the current directory)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("init", help="create a new project directory with an annotated claire.yaml")
    s.add_argument("directory")
    s.add_argument("--ligand", help="ligand model (.sdf), copied to inputs/ligand.sdf")
    s.add_argument("--name", help="3-letter ligand residue name")
    s.add_argument("--scaffolds", help="folder of scaffold PDBs (and optional .pos files), copied to inputs/scaffolds")

    sub.add_parser("status", help="show which steps are complete")
    sub.add_parser("check", help="check that Python packages and external tools are installed and configured")
    r = sub.add_parser("run", help="run several steps in order (local backend)")
    r.add_argument("--from", dest="start", default="target", choices=ORDER)
    r.add_argument("--to", dest="end", default="final", choices=ORDER)

    for name, (_, helptext) in STEPS.items():
        s = sub.add_parser(name, help=helptext)
        s.add_argument("--collect", action="store_true", help="finish the step after cluster jobs have completed")
        if name == "motifs":
            s.add_argument("--fast", action="store_true",
                           help="use motifs.fast_preset (wide tolerances, ideal values only) for quick test runs")
        if name == "mpnn":
            s.add_argument("--cutoff", type=float, help="binding-site distance (A) kept fixed; overrides mpnn.binding_site_cutoff")
        if name == "fragments":
            s.add_argument("--list-atoms", action="store_true", help="print ligand atom names to help define fragments")

    args = ap.parse_args(argv)
    try:
        if args.cmd == "init":
            cfg = init_project(args.directory, args.ligand, args.name, args.scaffolds)
            print(f"Created {cfg}\nEdit it (ligand name, tool paths, fragments, polar atoms), then run `claire target`.")
            return 0
        proj = Project(args.project) if args.project else Project.find()
        if args.cmd == "check":
            from .check import run_check
            return run_check(proj)
        if args.cmd == "status":
            print(f"Project: {proj.root}  (ligand {proj.lig}, backend {proj.config['compute']['backend']})")
            for key, folder, st in proj.status():
                print(f"  {key:<11} {folder:<15} {st}")
            return 0
        if args.cmd == "run":
            if proj.config["compute"]["backend"] != "local":
                raise ClaireError("`claire run` chains steps on the local backend only; "
                                  "on a cluster run each step and its --collect in turn.")
            steps = ORDER[ORDER.index(args.start):ORDER.index(args.end) + 1]
            for st in steps:
                print(f"\n=== {st} ===")
                STEPS[st][0](proj, argparse.Namespace(collect=False, list_atoms=False))
            return 0
        STEPS[args.cmd][0](proj, args)
        return 0
    except ClaireError as e:
        print(f"claire: error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
