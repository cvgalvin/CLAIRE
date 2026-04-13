Necessary Software
---------------
prody

babel

openbabel

numpy

scipy

io

biopython

shutil

collections

pandas

amber antechamber

proteinMPNN 

colabfold

rosetta/pyrosetta


Overview of Workflow
---------------

1.prepare target model/params
2.define fragments/perform fragment search
3.scrape pdb 4 contacts + filter 
4.cluster contacts
5.generate motifs
6.match
7.filter matches for ligand burial(+residue preferences)
8.enzyme design application on matches
9.hbond refinement
10.fast scoring binding metrics (eg hbonds+burial)
11.additional scoring binding metrics (preorganization, molecular contact surface, etc.)
12.proteinmpnn on non binding site positions
13.fastdesign with mpnn seq profile - score stability metrics 
14.colabfold
15.score all metrics repeat - final filtering of designs

##Description of Contents

#db contains lists for matching user defined ligand substructures with existing ligands in pdb, then finding pdb id codes that those ligands are a part of to scrape contacts from 

#workflow.py script includes information on dependencies, general workflow, some example code 
and options used for submitting rosetta applications and xml scripts 

#main contains the code for the primary workflow of scraping pdb for protein-ligand contents with user defined ligand substructures, transforming relative to target ligand, filtering for quality+redundancy, and clustering the resultant interaction pools 

#tools contains most useful design evaluation xmls and also refinement scripts (hbrefine.py, fd_mpnn.xml)

#extra_tools contains additional scripts that are potentially useful 

	-genpot folder contains params generation, enzyme design application + 	metric analysis, fastdesign (with and without 3bop), and coupledmoves 		application when using rosetta generic potential scorefunction 

	-other_metric_eval contains rosetta xml for design evaluation for different cases 

	-cm_standard_jump1.xml = rosetta xml for running coupledmoves protocol 	with standard all atom scorefunction

	-FD_XXXXX.xml scripts run rosetta fastdesign application with different options 

	-motif_to_cst.py was used in previous iterations of workflow where monte carlo methods were used for motif generation, they firsts made pdb files of motif residues with the target ligand. this script takes a pdb file as input and outputs a rosetta match cst file that describes the protein ligand interctions in the input pdb 

	-rfdaa_example.sh contains example command for running rosetta fold diffusion all atom and subsequently ligandmpnn on outputs 
