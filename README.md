# CLAIRE

**C**ombinatoria**L** **A**ssembly with **I**ntegrated **RE**finement: design of small-molecule binding proteins by
combinatorially placing statistically enriched protein–ligand interaction motifs into libraries of de novo scaffolds,
followed by hydrogen-bond refinement, ProteinMPNN stabilization and AlphaFold2-based filtering.

> Galvin CV, Guo AB, Chen MN, Alfonso IL, Grisingher D, Kretschmer S, Kranthi D, Kelly MJS, Kortemme T.
> *A Hybrid Physics–Deep Learning Framework for Combinatorial De Novo Design of Small-Molecule Binding Proteins.*
> bioRxiv (2026). https://doi.org/10.64898/2026.04.12.717919

This package turns the original research scripts into a command-line tool. You provide a ligand model and a scaffold
library; each workflow step is one command, every step finds its inputs in the previous step's folder, and the whole
project is configured from a single `claire.yaml`.

---

## Installation

```bash
conda create -n claire python=3.10 -c conda-forge openbabel
conda activate claire
pip install -e .           # from a clone of this repository
```

`pip install` installs every Python dependency (numpy, scipy, pandas, matplotlib, PyYAML, RDKit, ProDy, Biopython).
The following programs cannot be installed through pip, because they need separate licensing, compilation or GPU
setup. Install them and point `claire.yaml` at them:

| Program | Needed for | How to install |
|---|---|---|
| [Rosetta](https://rosettacommons.org) | target, match, burial, enzdes, score, fastdesign | licensed download; set `tools.rosetta_bin` (or `$ROSETTA_BIN`) |
| [PyRosetta](https://www.pyrosetta.org) | cluster (H-bond screen), motifs, enzdes, hbrefine | licensed; `pip install pyrosetta-installer` then `python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"` |
| [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) or [LigandMPNN](https://github.com/dauparas/LigandMPNN) | mpnn | `git clone`; set `tools.proteinmpnn_dir` (auto-detected; needs PyTorch). For LigandMPNN run `bash get_model_params.sh "./model_params"` once. `mpnn.model_type: ligand_mpnn` uses the ligand-aware model |
| [ColabFold](https://github.com/sokrypton/ColabFold) | predict | e.g. LocalColabFold; `colabfold_batch` on PATH (GPU recommended) |
| Open Babel | fragments (optional) | `conda install -c conda-forge openbabel`; RDKit is used if absent |

Run `claire check` inside a project to see what is found and which steps each missing piece affects.

## Quick start

```bash
claire init my_target --ligand ligand.sdf --name PRG --scaffolds path/to/scaffolds/
cd my_target
#   edit claire.yaml: tool paths, compute backend
claire target                     # params file + model         -> 01_target
claire fragments --list-atoms     # print ligand atom names
#   edit claire.yaml: fragments, hbond_fragments, polar_atoms
claire fragments                  # PDB ligands/entries containing each fragment -> 02_fragments
claire contacts                   # download structures, extract + filter contacts -> 03_contacts
claire cluster                    # cluster contacts per residue type            -> 04_clusters
claire motifs                     # polar + nonpolar motif constraint files      -> 05_motifs
claire match                      # RosettaMatch into scaffolds                  -> 06_matches
claire burial                     # keep matches with a buried ligand            -> 07_buried
claire enzdes                     # binding-site design                          -> 08_enzdes
claire hbrefine                   # hydrogen-bond refinement                     -> 09_hbrefine
claire score                      # binding metrics + filters                    -> 10_scored
claire mpnn                       # ProteinMPNN sequence profiles                -> 11_mpnn
claire fastdesign                 # FastDesign with MPNN profile + filters       -> 12_fastdesign
claire predict                    # ColabFold (single sequence)                  -> 13_colabfold
claire final                      # pLDDT / Ca RMSD selection                    -> 14_final
claire status                     # what has been run
claire check                      # which dependencies are installed
```

On a workstation, `claire run --from match --to final` chains steps.

## Inputs

* **Ligand model** (`ligand.sdf`): protonated for pH 7.4 and minimized (e.g. MMFF94 in Avogadro). For multiple
  conformers set `ligand.conformers` > 1 (RDKit conformer generation).
* **Scaffold library** (`inputs/scaffolds/*.pdb`): e.g. the 1,816 LUCS NTF2 scaffolds from the paper
  ([Dryad](https://doi.org/10.5061/dryad.7pvmcvf8c)). An optional `<scaffold>.pos` file lists matcher positions;
  without one, core-facing residues are chosen heuristically (`scaffolds.auto_pos_cb_cutoff`). Providing `.pos`
  files is recommended.
* **Fragments** (`fragments:` in `claire.yaml`): ligand substructures, as lists of atom names, used to mine protein
  contacts from the PDB. One entry per fragment:
  ```yaml
  fragments:
    - name: ringA
      atoms: [C1, C2, C3, C4, C5, C10]
    - name: acetyl
      atoms: [C17, C20, O20, C21]
  ```
  (shorthand: `fragments: [[C1, C2, C3, C4, C5, C10], [C17, C20, O20, C21]]`). Keep fragments small: each is
  searched for in all PDB ligands, and a large fragment only matches close analogues of your target, giving few
  contacts. Add hydrogens where a substructure should terminate (e.g. on a carboxylate oxygen so esters
  are not matched). Nonpolar fragments are usually enough; list H-bonding fragments under `hbond_fragments` to keep
  only clusters that genuinely hydrogen bond.
* **Polar atoms** (`polar_atoms:`): each ligand polar atom that should receive an explicit hydrogen bond, with two
  bonded atoms defining the geometry, donor/acceptor role and hybridization.

## Structures and the PDB-REDO requirement

The original scripts needed a full local PDB-REDO mirror. This package instead:

* downloads only the entries containing your fragments, per entry, from PDB-REDO into a shared cache
  (`~/.cache/claire`, override with `$CLAIRE_CACHE`), reused across projects;
* falls back to the RCSB deposited model if an entry is not in PDB-REDO, or PDB-REDO is unreachable
  (`sources.fallback_to_rcsb`); `03_contacts/structure_sources.json` records where each structure came from;
* uses a local mirror directly if you have one (`sources.local_pdb_redo`).

Fragment → ligand → entry lookup uses the RCSB search API by default (current PDB). Set
`sources.ligand_search: bundled` with `sources.bundled_db` pointing to the unzipped `db/` folder from the original
repository to reproduce the published (2022) searches exactly. Either way, every candidate ligand is re-checked with
the same RDKit substructure match as the original code.

## Matching tolerances

Motif constraint files carry, for each interaction, an ideal geometry, a tolerance, and a number of values sampled
around the ideal by RosettaMatch. By default these are exactly the values from the original `workflow.py` (they differ
by H-bond case). Override any of them in `claire.yaml`:

```yaml
motifs:
  constraints:
    hbond:    {distance_tolerance: 0.4, angle_samples: 2}
    nonpolar: {torsion_samples: 1}
```

Keys: `distance_tolerance`, `angle_A_tolerance`, `angle_B_tolerance`, `torsion_A_tolerance`, `torsion_AB_tolerance`,
`torsion_B_tolerance`, `distance_samples`, `angle_samples`, `torsion_samples`. Torsions left free in the defaults
(tolerance 360) stay free. `claire motifs --fast` applies `motifs.fast_preset` instead (wide tolerances, ideal values
only), which makes matching much faster but finds fewer matches; useful for test runs. The settings used are recorded
in `05_motifs/constraint_settings.json`.

## Clusters (SLURM / SGE)

```yaml
compute:
  backend: slurm        # or sge
  submit: true          # submit automatically; false = only write the job script
  setup: "source ~/miniconda3/bin/activate claire"
  slurm: {partition: normal, time: "24:00:00", mem: 8G}
```

Heavy steps write `NN_stage/jobs/<name>.sh` (an array job over `<name>.tasks`) and exit. When the jobs finish,
`claire <step> --collect` gathers the results (multi-phase steps such as `contacts` and `score` submit their next
phase on `--collect`). Logs are in `NN_stage/jobs/logs/` (cluster) or `NN_stage/logs/` (local).

## Project layout

```
my_target/
  claire.yaml
  inputs/ligand.sdf, inputs/scaffolds/
  01_target/ ... 14_final/      each with manifest.json describing its outputs
  14_final/final_designs.csv    sequences, best ColabFold model, pLDDT, Ca RMSD
  14_final/designs/             design models + matching ColabFold predictions
```

## Defaults and differences from the original scripts

Thresholds and Rosetta flags are the original values, all editable in `claire.yaml`
(`filters.binding`, `filters.fastdesign`, `burial`, `final`, `motifs`). Scientific code from the original
repository (contact extraction/filtering, clustering, HBRefine, resfile generation, conformer generation) is vendored
in `claire/legacy/` with only path handling and pandas 2 compatibility changed; logic that lived inline in
`workflow.py` (motif assembly, H-bond cluster screen, score filtering, AF2 analysis) is ported into `claire/steps/`.
Changes worth knowing:

* **ProteinMPNN design region**: the motif residues and all residues within `mpnn.binding_site_cutoff` (10 Å) of
  the ligand are fixed; ProteinMPNN redesigns everything outside the binding-site neighborhood.
* **MPNN profile resfiles** allow, at each designable position, the residues sampled by ProteinMPNN; binding-site
  positions are NATAA.
* **Final selection** computes Ca RMSD with a numpy superposition instead of PyRosetta.
* **HBRefine outputs** carried forward are designs whose ligand polar atoms were already satisfied (`*_unrefined.pdb`)
  and designs HBRefine fixed (`*_fixed.pdb`); designs it could not fix (`_unsats`, `_ubp`, `_bothproblems`) are dropped.

## Module map (original → package)

| Original | Package |
|---|---|
| `main/new_target.py`, molfile_to_params | `claire init`, `claire target` |
| `main/Fragment_Search.py` + `db/` | `claire fragments` (RCSB API or bundled db) |
| `extract_align_contacts.py`, `process_fragment_search.py`, `setup_filter_contacts.py`, `filter_contacts_array.py`, `consolidate_filtered_contacts.py`, `process_filtered_contacts2.py` | `claire contacts` |
| `spatial_cluster.py`, workflow.py §4 | `claire cluster` |
| workflow.py §5 | `claire motifs` |
| workflow.py §6–15, `match_resfile2.py`, `hbrefine.py`, XMLs, `run_colabfold.py` | `claire match` … `claire final` |

## License

Free for non-commercial research and teaching; see [LICENSE](LICENSE). Commercial use requires a separate license
(UCSF Innovation Ventures). Rosetta, PyRosetta, ProteinMPNN and ColabFold are not distributed with this package and
are covered by their own licenses.

## Citation

Galvin CV, Guo AB, Chen MN, Alfonso IL, Grisingher D, Kretschmer S, Kranthi D, Kelly MJS, Kortemme T.
A Hybrid Physics–Deep Learning Framework for Combinatorial De Novo Design of Small-Molecule Binding Proteins.
bioRxiv (2026). https://doi.org/10.64898/2026.04.12.717919
