Combinatorial Assembly with Integrated REfinement for de novo Design of Small Molecule Binding Proteins
---------------
This repository contains code necessary to recapitulate the workflow described in [this paper:](https://doi.org/10.64898/2026.04.12.717919) 
	
	Connor V. Galvin, Amy B. Guo, Maple N. Chen, Isabella L. Alfonso, Dominic Grisingher, Simon Kretschmer, Divya Kranthi, Mark J. S. Kelly, Tanja Kortemme. A Hybrid Physics-Deep Learning Framework for Combinatorial De Novo Design of Small-Molecule Binding Proteins

where it was successfully applied to design binders for two target molecules, achieving high experimental success rates in both cases (4/26 tested designs bound target), sufficient to bypass the need for large-scale experimental screening. 

Provided are scripts and instructions to perform motif library generation for an arbitrary small molecule target, along with matching, optimization, and evaluation of resultant designs. Excluded, however, is the generation of a protein scaffold library against which to screen binding motifs, for which the [LUCS](https://github.com/Kortemme-Lab/loop_helix_loop_reshaping) software was used in the referenced publication. The 1,816 LUCS-generated NTF2 fold proteins used therein is available as [supplementary data](https://doi.org/10.5061/dryad.7pvmcvf8c), under "LUCS_NTF2_scaffold_library".

'workflow.py' provides a complete walkthrough of the workflow, with example terminal commands for executing every step as well as descriptions of inputs/outputs for each. Many of the scripts were produced according to the system specifications at the time of their application, and thus reference paths and distributed computing options that individual users must change to align with their system, but note is made where this is the case and the relevant lines are specified. 

Necessary Software/Modules for 1:1 Execution of Workflow
---------------

prody <br />
babel <br />
openbabel <br />
numpy <br />
scipy <br />
io <br />
biopython <br />
shutil <br />
collections <br />
pandas <br />
amber antechamber <br />
proteinMPNN <br />
colabfold <br />
Rosetta/pyrosetta <br />


Overview of Workflow
---------------

1. Prepare target model/params

2. Define fragments/perform fragment search

3. Scrape pdb 4 contacts + filter 

4. Cluster contacts

5. Generate motifs

6. RosettaMatch

7. Filter matches for ligand burial(+residue preferences)

8. EnzymeDesign application on match binding sites

9. HBond Refinement

10. Fast scoring - binding metrics (eg hbonds+burial)

11. Additional scoring binding metrics (preorganization, molecular contact surface, etc.)

12. ProteinMPNN on non binding site positions

13. Rosetta FastDesign with PMPNN sequence profile + score stability metrics 

14. Colabfold

15. Score all metrics + final filtering of designs


Description of Contents
---------------

#db contains lists for matching user defined ligand substructures with existing ligands in pdb, then finding pdb id codes that those ligands are a part of to scrape contacts from 

#workflow.py script includes information on general workflow with example code and flags used for submitting rosetta applications and xml scripts 

#main contains the code for the motif library generation, ie scraping pdb for protein-ligand contacts with user defined ligand substructures, transforming relative to target ligand, filtering for quality+redundancy, and clustering the resultant interaction pools 

#tools contains most useful design evaluation xmls and refinement scripts (hbrefine.py, fd_mpnn.xml)

#extra_tools contains additional scripts that may be useful, eg

	-genpot folder contains params generation, enzyme design application + 	metric analysis, fastdesign (with and without 3bop), and coupledmoves 		application when using rosetta generic potential scorefunction 

	-other_metric_eval contains rosetta xmls for design evaluation for different cases 

	-cm_standard_jump1.xml = rosetta xml for running coupledmoves protocol 	with standard all atom scorefunction

	-FD_XXXXX.xml scripts run rosetta fastdesign application with different options 

	-motif_to_cst.py this script takes a pdb file as input and outputs a rosetta match cst file that describes the protein-ligand interction geometries in the input pdb 

	-rfdaa_example.sh contains example command for running rosetta fold diffusion all atom and subsequently ligandmpnn on outputs 
