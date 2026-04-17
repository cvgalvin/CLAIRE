#############################################
# 1.prepare target model/params
#############################################

#create a folder for target with the proper directory structure
time python main/new_target.py target_3_letter_code
#eg 
time python main/new_target.py 38e

'''
^^^THIS FOLDER IS NOW OUR 'TARGET PARENT DIRECTORY'
'''


'''
download sdf model of target, rcsb and pubchem are good sources
if ONLY 1 CONFORMER
AVOGADRO open target sdf 
    hydrogens for ph 7.4
    optimize energy with mm94 forcefield until dE=0


move target sdf file to target parent directory


for multiple conformers
AVOGADRO
    hydrogens for ph 7.4
    optimize energy with mm94 ff
    Extensions --> Molecular Mechanics --> Setup Force Field...
        mm94, default
    Extensions --> Molecular Mechanics --> Conformer Search...

alternatively use included conformer generation script:
'''
time python tools/gencon.py input_file.sdf number_desired_conformers max_attempts prune_rms_threshold
#eg
time python tools/gencon.py 38e.sdf 50 1000 0.1



'''
create rosetta params file (will also output pdb of target molecule)
'''
#in terminal, directory where you have target sdf file:
compounds=("38e") #your ligand 3 letter code 
for I in ${compounds[@]}; do
    ~/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n $I -p $I ${I}.sdf
done
'''
 NOTE:  for eg a carboxylate you want to run the fragment search (next step) with the hydrogen on the carboxylate O
        even though it is not protonated at physiological pH
        otherwise it will take contacts for esters 
        that is, hydrogen serves as sort of terminal cap for the fragment search, so add it wherever you dont 
        want to include contacts with continued alkyl chains
        when you then wanna do further design and everything,
        you use a separate model of the target without the hydrogen 
        so here you would make two versions of params and pdb model, one with one without hydrogen on relevant atoms




move params and pdb to Inputs/Rosetta_Inputs:
'''
compounds=("38e")#your ligand 3 letter code 
for I in ${compounds[@]}; do
    mv ${I}.params Inputs/Rosetta_Inputs/${I}.params
    mv ${I}_0001.pdb Inputs/Rosetta_Inputs/${I}_0001.pdb
done








#############################################
# 2.define frags/fragment search
#############################################
'''
fragment search

make a text file called ligand_3_letter_code_frags.txt
                eg 38efrags.txt

in it specify the atom names composing fragments between brackets,
eg for a benzene ring and a dihydroxyl group this file might look like:

------------------
<
C1
C2
C3
C4
C5
C6
>

<
C1
C2
O1
O2
H1
H2
>
------------------
these atom names should correspond to those in the pdb file for your target molecule 
that was generated along with its params file 

you can open that target molecule pdb in chimera and hover over atoms to get their
names real quick, or also click the atoms in pymol to get names

YOU USUALLY ONLY NEED TO DO NONPOLAR FRAGMENTS, since we manually
make hbonds with polar fragments, unless there is a 'special' polar fragment where you
want discrete contacts or hbonds. eg guanidinium really wants a bidentate hbond with
a carboxylate, so you can do frag search for that and then use derived contacts
to specify specific rotamers that make those interactions. you may also wanna do that
for a carboxylate, which prefers bidentate with arg. 
or two hydroxyl groups next to each other you could 
try to scrape contacts that hbond with both simultaneously.


NO PROBLEM WITH BEING REDUNDANT WITH FRAGMENTS.
ie you could do benzene ring alone, but then also use the same ring with the substituents
and see if you can find any contacts for that.


now open /main/Fragment_Search.py
and change the paths on lines 12 and 21 to reflect where your db folder is 
then in target parent directory, with this frags.txt file:
'''
time python /main/Fragment_Search.py ${I}frags.txt Inputs/Rosetta_Inputs/${I}_0001.pdb
# where ${I} is your target 3 letter code 
#and move fragment pdb and mol files to fragment inputs folder
mv *.pdb Inputs/Fragment_Inputs/
mv *.mol Inputs/Fragment_Inputs/
'''
FROM THIS WE GET output: PDBSEARCHRESULTS.JSON

'''
















####################################
# 3.scrape pdb 4 contacts
####################################
'''
we next run main/extract_align_contacts.py on the json output from the previous step 

USES LOCAL COPY OF PDBREDO DATABASE, MUST MANUALLY CHANGE PATH IN SCRIPT TO LOCATION IN YOUR SYSTEM
'''
time ipython main/extract_align_contacts.py  pdb_search_results_json
'''
however, if there are many outputs, that we want to split for distributed computing 
we can first use the process_fragment_search.py script
this will split the json file into multiple outputs 
the final block in the script is for outputting shellfile (run_align.sh) 
for my personal account on ucsf wynton cluster 
that references all of these file output locations so i can easily submit job, 
but you can change to reflect your account/system
'''
time python ~/main/process_fragment_search.py
'''
NOTE THE TEXT OUTPUT IN TERMINAL FOR THIS SCRIPT, THIS IS NUMBER OF JOBS
WE SUBMIT ON CLUSTER, eg:
getting smiles...
474
then ->
'''
qsub -cwd -t 1-474 run_align.sh
'''
THE OUTPUT OF THIS IS DIRECTORY 
Transformed_Aligned_pdbs
WHICH CONTAINS SUBDIRECTORIES FOR EVERY FRAGMENT, AND IN THESE IS THE RAW PDB FILES 
FOR THE PDB SCRAPED CONTACTS WITH THE RELEVANT LIGAND FRAGMENT





NEXT WE FILTER CONTACTS FOR REDUNDANCY/QUALITY, and consolidate them into one pdb file (fuzzball)
THIS IS CURRENTLY SETUP ONLY TO RUN ON UCSF WYNTON AS ARRAY JOB 
BUT THE SCRIPT :

filter_contacts_array.py 

can be modified to work with your system 


FIRST WE 
SET UP SHELLFILE FOR SUBMITTING FILTER JOB
in target parent directory:
'''
#!/bin/bash
time python main/setup_filter_contacts.py
'''

TERMINAL OUTPUT OF THIS SCRIPT IS THE NUMBER YOU NEED FOR THE NEXT ARRAY JOB
go to cluster_output folder and do
cat run_setupfilter.sh.o.......
eg:
38e
147

THEN WE RUN THE GENERATED SHELLFILE WITH APPROPRIATE INPUT -t 1-X
'''
qsub -cwd -t 1-147 run_filter_array.sh
'''

THIS IS WILL CREATE A FRAGMENT_X_CONTACT_FUZZBALL.PDB IN EACH TRANSFORMED_ALIGNED_PDBS/FRAGMENT_X



CONSOLIDATE FILTERED CONTACTS
'''
#!/bin/bash
time python /main/consolidate_filtered_contacts.py
'''
THIS REORDERS ALL RESIDUES IN FUZZBALLS FROM 1 




NEXT WE FURTHER PROCESS THESE FILTERED CONTACTS - do for each fragment
TO MAKE FUZZBALLS FOR EACH INDIVIDUAL RESIDUE TYPE
MAKE SURE YOU CHANGE PYTHON SCRIPT COMMAND below to relevant ligand 3 letter code, show example with 38e
bash:
'''
compounds=("Fragment_1" "Fragment_2" "Fragment_3")
for I in ${compounds[@]}; do
    time python3 /main/process_filtered_contacts2.py 38e/Transformed_Aligned_PDBs/${I}/${I}_contact_fuzzball.pdb 38e/Inputs/Rosetta_Inputs/38e_0001.pdb 38e $I
done
'''
WE MOVE ALL PROCESSED CONTACTS TO ONE DIRECTORY, THE OUTPUT FOLDERS FROM LAST STEP ALL 
HAVE _residue_contacts SUFFIX 
'''
mkdir residue_contact_fuzzballs
mv *_residue_contacts residue_contact_fuzzballs

##################################
# 4.cluster contacts
##################################

'''

CLUSTER ALL
    CHANGE THE PATHS HERE AT cmd='....'
    AND CHANGE CONDA STUFF IF NECESSARY

we use :


main/spatial_cluster.py 

script, below example for how I did it  
to split up as array job for every residue type fuzzball for each fragment
but must be reworked for whatever system you are using
'''
import os
dirs=[i for i in os.listdir() if os.path.isdir(i)==True]
#
residues=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP','ILE','LEU','MET','PHE','VAL','GLY']
#
c=1
for fragdir in dirs:
    os.chdir(fragdir)
    l=[]
    for i in residues:
        for x in os.listdir():
            if x[-3:]=='pdb':
                if i in x:
                    l.append(os.path.join(os.getcwd(),x))
    jname='sc'+str(c)+'.sh'
    f=open(jname,'w')
    f.write('#!/bin/bash')
    f.write('\n')
    f.write('\n')
    f.write('tasks=(0\n')
    for match in l[:-1]:
        f.write('       '+str(match)+'\n')
    f.write('       '+l[-1]+')')
    f.write('\n')
    ########################################
    ########################################
    cmd='time python3 ~/main/spatial_cluster.py ${tasks[$SGE_TASK_ID]} 38e/Inputs/Rosetta_Inputs/38e_0001.pdb'
    ########################################
    ########################################
    f.write(cmd)
    f.write('\nqstat -j "$JOB_ID"')
    f.close()
    print(len(l))
    os.system('qsub -cwd -t 1-17 -l mem_free=16G '+jname)
    os.chdir('..')
    c+=1



'''
okay now we have all of our clustered contacts, we can look at a few things 

for our example 38e here, with the fragments we defined:
1 benzene
2 imidazole
3 benzene with fluorines


we might want to look at hbonding stats for fragment (2)
so in fragment 2 directory:
'''
import os
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP'] #
rescounts={}
for i in os.listdir():
    if os.path.isdir(i)==True:
        rn=i.split('_')[3]
        if rn in polar_residues:
            os.chdir(i)
            l=[i for i in os.listdir() if i[-3:]=='pdb']
            count=0
            for clust in l:
                cc=clust.split('_')[-1].split('.')[0]
                ncc=int(cc)
                count+=ncc
            os.chdir('..')
            rescounts[rn]=count

'''
In [2]: rescounts
Out[2]:
{'LYS': 632,
 'ASP': 572,
 'GLU': 369,
 'ARG': 1066,
 'GLN': 461,
 'HIS': 634,
 'ASN': 1509,
 'SER': 968,
 'THR': 1123,
 'TRP': 652,
 'TYR': 833}
'''
#plot
import matplotlib.pyplot as plt
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP'] #no his
tc=0
for key in d.keys():
    tc+=d[key]
freqs=[]
for i in polar_residues:
    try:
        freqs.append((i,d[i]/float(tc)))
    except:
        pass
freqs=sorted(freqs, reverse=True, key=lambda nmem: nmem[1])
labs=[]
vals=[]
for a,b in freqs:
    labs.append(a)
    vals.append(b)

fig = plt.figure()
ax = fig.add_subplot()
ax.bar(labs,vals,color='orange',width=0.1)
ax.set_title('Hbond Residue Preferences')
ax.set_ylabel('Frequency')
ax.set_xlabel('Residue Type')
plt.savefig('respref.pdf')
plt.clf()
plt.close()

'''
CHECK WITHIN THESE CLUSTERS WHETHER OR NOT THERE ACTUALLY
IS A HYDROGEN BOND? THEN GET SAME STATS BUT FOR JUST THOSE EXAMPLES?
    also move them to their own directory if we want to use in motif csts

USING ROSETTA HERE
'''
import os
from pyrosetta import *
init('-pack_missing_sidechains False')
#
sf=get_fa_scorefxn()

##########################################
##########################################
##########################################
ligparams=['Inputs/Rosetta_Inputs/38e.params']
##########################################
##########################################
##########################################

clusts_data={}
#
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP']
#
hbond_clusts=[]
for resdir in os.listdir():
    if os.path.isdir(resdir)==True:
        rn=resdir.split('_')[3]
        if rn in polar_residues:
            os.chdir(resdir)
            clusts=[i for i in os.listdir() if i[-3:]=='pdb']
            clust_data=[]
            for clust in clusts:
                nhb=0
                clust_pop=clust.split('_')[-1].split('.')[0]
                nclust_pop=int(clust_pop)
                p = Pose()
                generate_nonstandard_residue_set(p,ligparams)
                pose_from_file(p, clust)
                res_scores=[]
                #
                ligand_resnum=p.total_residue()
                for res_ind in range(1,ligand_resnum):
                    contact_pose=Pose()
                    ligand_pose=p.residue(ligand_resnum).clone()
                    res_pose=p.residue(res_ind).clone()
                    contact_pose.append_residue_by_jump(res_pose, 1)
                    contact_pose.append_residue_by_jump(ligand_pose, 1)
                    res_resname=contact_pose.residue(1).name()[0:3]
                    lig_resname=contact_pose.residue(2).name()[0:3]
                    if res_resname==lig_resname:
                        print('uh oh')
                        continue
                    else:
                        pass
                    hbond_set = rosetta.core.scoring.hbonds.HBondSet()
                    contact_pose.update_residue_neighbors()
                    rosetta.core.scoring.hbonds.fill_hbond_set(contact_pose, False, hbond_set)
                    if hbond_set.nhbonds()>0: #IF YOU WANT DOUBLE HBONDS, CHANGE TO >1
                        nhb+=1
                    if nhb>=(0.1*nclust_pop):
                        hbond_clusts.append(os.path.join(resdir,clust))
                        break
            os.chdir('..')


rescounts={}
polar_residues=['SER','THR','GLN','ASN','TYR','TRP','GLY','HIS','ARG','LYS','GLU','ASP']
for i in polar_residues:
    c=0
    for j in hbond_clusts:
        rn=j.split('_')[3]
        if rn==i:
            cp=j.split('_')[-1].split('.')[0]
            cpn=int(cp)
            c+=cpn
    rescounts[i]=c

os.makedirs('hbond_clusts',exist_ok=True)
for clust in hbond_clusts:
    pop=clust.split('_')[-1].split('.')[0]
    npop=int(pop)
    if npop>=15:
        os.system('cp '+clust+' hbond_clusts/'+clust.split('/')[-1])

'''

{'SER': 37,
 'THR': 32,
 'GLN': 81,
 'ASN': 507,
 'TYR': 29,
 'TRP': 12,
 'GLY': 0,
 'HIS': 55,
 'ARG': 72,
 'LYS': 78,
 'GLU': 6,
 'ASP': 4}

lets look at benzene and plot score vs clust freq
WITH JUST BENZEE FRAG SO NO CLASEHES
    in fragment 1 directory:
'''
#

import os
dirs=[i for i in os.listdir() if os.path.isdir(i)==True]
#
from pyrosetta import *
init('-pack_missing_sidechains False')
#
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep,fa_elec, fa_sol,hbond_sc, fa_elec, hbond_bb_sc,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
#etc. ... all terms 
#
ligparams=['Inputs/Rosetta_Inputs/38e.params']
#
clusts_data={}
#
for resdir in dirs:
    os.chdir(resdir)
    clusts=[i for i in os.listdir() if i[-3:]=='pdb']
    clust_data=[]
    for clust in clusts:
        clust_pop=clust.split('_')[-1].split('.')[0]
        nclust_pop=int(clust_pop)
        p = Pose()
        generate_nonstandard_residue_set(p,ligparams)
        pose_from_file(p, clust)
        res_scores=[]
        #
        ligand_resnum=p.total_residue()
        for res_ind in range(1,ligand_resnum):
            contact_pose=Pose()
            ligand_pose=p.residue(ligand_resnum).clone()
            res_pose=p.residue(res_ind).clone()
            contact_pose.append_residue_by_jump(res_pose, 1)
            contact_pose.append_residue_by_jump(ligand_pose, 1)
            res_resname=contact_pose.residue(1).name()[0:3]
            lig_resname=contact_pose.residue(2).name()[0:3]
            if res_resname==lig_resname:
                print('uh oh')
                continue
            else:
                pass
            # min_mover.apply(contact_pose)
            score_this_one=sf(contact_pose)
            res_scores.append(score_this_one)
        bestscore=min(res_scores)
        clust_data.append((nclust_pop,bestscore))
    clusts_data[resdir.split('_')[3]]=clust_data
    os.chdir('..')

import json
json.dump(clusts_data,open('clustsdata.json','w'))

'''
AFTER IT FINISHES WE PLOT:
'''
import json
with open('clustsdata.json','r') as f:
    clusts_data=json.load(f)

#plot score vs freq across all
x=[]
y=[]
for key in clusts_data.keys():
    for i in clusts_data[key]:
        if -5<i[1]<5:
            x.append(i[1])
            y.append(i[0])

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot()
#
ax.scatter(x,y)
#
ax.set_title('Rosetta Energy Score vs. Cluster Population')
ax.set_ylabel('Population')
ax.set_xlabel('Rosetta Energy Score (REU)')
plt.savefig('scorevspop.pdf')
plt.clf()

'''
and we can see that score does not correlate with frequency, so we definitely want to 
bias motif generation towards most highly populated clusters
'''































########################################
# 5.generate motifs
########################################



'''
OKAY - MOTIF GENERATION WITH MANUAL GENERAL HBOND CONSTRAINTS
'''
#input ligand name and number of unique polar csts you want
lig_resname='prg' #progesterone
#sp2 yes/no/s=sp.donor 
#is the ligand polar atom an hbond donor? y or n 
#.a1a2a3
relevant_atoms={'1':['y','n','O2','C20','C19'],
                '2':['y','n','O1','C17','C21'],
                '3':['y','n','O2','C20','C19']}
#define wheteher or not charged residues are allowed to be used with frag
charged_allowed={'1':'yes',
                '2':'yes',
                '3':'yes'}
#########
allowed_problematic_motifs=15 #problematic motif=has problematic res, still only 1 allowed
#new directory where motifs will be output
polar_cst_paths='hb_motifs'
n_desired_motifs=100



#define polar residues that can be used in csts (donate hb, no asp/glu)
#and problematic residues that will be limited to 1 per cst
q_polar_residues=['LYS','ARG','SER','THR','TYR','TRP','GLY'] #includes charged
polar_residues=['SER','THR','TYR','TRP','GLY'] #no his
problematic_residues=['LYS','ARG','GLN','ASN']
acceptor_residues=['GLU','ASP','GLY']
acceptor_res=[]
######
import random
import os



#3 to 1 letter code dict
one_lcd={'LYS':'K',
         'ARG':'R',
         'SER':'S',
         'THR':'T',
         'GLN':'Q',
         'ASN':'N',
         'HIS':'H',
         'TYR':'Y',
         'TRP':'W',
         'GLY':'G',
         'GLU':'E',
         'ASP':'D'}

#pass as residue type either Hpol for all sc or HNbb for bb gly
def generate_single_constraint_block_base(distanceAB,angleA,angleB,torsionA,torsionAB,torsionB,
                                          residue_resname,ligand_resname,ligand_atoms,residue_atoms,
                                          distance_tolerance=0.15, angle_A_tolerance=15, angle_B_tolerance=15,
                                          torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                          torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                          distance_constraint_sample_number=0,isbb=False):
    """
    Generate a single constraint block for one residue-ligand interaction
    :return:
    """

    residue_tag = 'residue3'
    if isbb==False:
        constraint_block = [
            '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(ligand_atoms)),
            '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
            '  TEMPLATE::   ATOM_MAP: 2 atom_type: {}'.format(residue_atoms),
            '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, residue_resname),
            '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
                distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
            '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(angleA), angle_A_tolerance, 100, angle_constraint_sample_number),
            '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(angleB), angle_B_tolerance, 100, angle_constraint_sample_number),
            '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  180.00  {3:3}'.format(
                float(torsionA), torsion_A_tolerance, 100,
                torsion_constraint_sample_number),
            '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(torsionB), torsion_B_tolerance, 100,
                torsion_constraint_sample_number),
            '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(torsionAB), torsion_AB_tolerance, 100,
                torsion_constraint_sample_number)]
    else:
        constraint_block = [
            '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(ligand_atoms)),
            '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
            '  TEMPLATE::   ATOM_MAP: 2 atom_type: {}'.format(residue_atoms),
            '  TEMPLATE::   ATOM_MAP: 2 is_backbone',
            '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, residue_resname),
            '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
                distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
            '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(angleA), angle_A_tolerance, 100, angle_constraint_sample_number),
            '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(angleB), angle_B_tolerance, 100, angle_constraint_sample_number),
            '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  180.00  {3:3}'.format(
                float(torsionA), torsion_A_tolerance, 100,
                torsion_constraint_sample_number),
            '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(torsionB), torsion_B_tolerance, 100,
                torsion_constraint_sample_number),
            '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
                float(torsionAB), torsion_AB_tolerance, 100,
                torsion_constraint_sample_number)]




    return constraint_block




#for each fragment, create a cst block
motifs=[]
nproblematic=0
ntries=0
while len(motifs)<n_desired_motifs:
    problem_res=0
    ngly=0
    current_motif=[]
    residues_used=[]
    for frag in relevant_atoms.keys():
        isdonor=relevant_atoms[frag][1]
        sp2_acc=relevant_atoms[frag][0]
        lig_atoms=relevant_atoms[frag][2:]
        if isdonor=='n':
            if charged_allowed[frag]=='yes':
                res_resname=random.choice(q_polar_residues)
            else:
                res_resname=random.choice(polar_residues)
            residues_used.append(res_resname)
            if res_resname=='GLY':
                res_atoms='HNbb'
                ngly+=1
                if sp2_acc=='y':
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=120,angleB=165,torsionA=180,torsionAB=0,torsionB=0,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15,
                                                                   torsion_A_tolerance=20, torsion_AB_tolerance=360, torsion_B_tolerance=360,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=3,isbb=True)
                elif sp2_acc=='s':
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=180,angleB=165,torsionA=180,torsionAB=0,torsionB=0,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15,
                                                                   torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=360,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=3,isbb=True)
                else:
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=120,angleB=180,torsionA=180,torsionAB=0,torsionB=0,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15,
                                                                   torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=360,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=5,isbb=True)
                current_motif.append(cstblock)
            else:
                res_atoms='Hpol'
                if res_resname in problematic_residues:
                    problem_res+=1
                if sp2_acc=='y':
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=120,angleB=165,torsionA=180,torsionAB=0,torsionB=0,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15,
                                                                   torsion_A_tolerance=20, torsion_AB_tolerance=360, torsion_B_tolerance=360,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=3,isbb=False)
                elif sp2_acc=='s':
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=180,angleB=165,torsionA=180,torsionAB=0,torsionB=0,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15,
                                                                   torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=360,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=3,isbb=False)
                    
                else:
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=120,angleB=180,torsionA=180,torsionAB=0,torsionB=0,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=10, angle_B_tolerance=15,
                                                                   torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=360,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=5,isbb=False)
                current_motif.append(cstblock)
        elif isdonor=='y':
            res_resname=random.choice(acceptor_residues)
            residues_used.append(res_resname)
            if res_resname=='GLY':
                res_atoms='OCbb'
                ngly+=1
                if sp2_acc=='y':
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=165,angleB=120,torsionA=0,torsionAB=0,torsionB=180,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=15, angle_B_tolerance=10,
                                                                   torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=20,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=3,isbb=True)
                else:
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=165,angleB=120,torsionA=0,torsionAB=0,torsionB=180,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=15, angle_B_tolerance=10,
                                                                   torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=20,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=3,isbb=True)
                current_motif.append(cstblock)
            else:
                res_atoms='OOC'
                if res_resname in problematic_residues:
                    problem_res+=1
                if sp2_acc=='y':
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=165,angleB=120,torsionA=0,torsionAB=0,torsionB=180,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=15, angle_B_tolerance=10,
                                                                   torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=20,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=3,isbb=False)
                else:
                    cstblock=generate_single_constraint_block_base(distanceAB=1.9,angleA=165,angleB=120,torsionA=0,torsionAB=0,torsionB=180,
                                                                   residue_resname=res_resname,ligand_resname=lig_resname,ligand_atoms=lig_atoms,residue_atoms=res_atoms,
                                                                   distance_tolerance=0.5, angle_A_tolerance=15, angle_B_tolerance=10,
                                                                   torsion_A_tolerance=360, torsion_AB_tolerance=360, torsion_B_tolerance=20,
                                                                   torsion_constraint_sample_number=5, angle_constraint_sample_number=1,
                                                                   distance_constraint_sample_number=3,isbb=False)
                current_motif.append(cstblock)
    if problem_res==1:
        if nproblematic<=allowed_problematic_motifs:
            if ngly<=2:
                if len(list(set(residues_used)))>=len(list(relevant_atoms.keys()))-3:
                    if current_motif not in motifs:
                        nproblematic+=1
                        motifs.append(current_motif)
    elif problem_res==0:
        if ngly<=2:
            if len(list(set(residues_used)))>=len(list(relevant_atoms.keys()))-3:
                if current_motif not in motifs:
                    motifs.append(current_motif)
    elif problem_res>1:
        pass
    ntries+=1
    if ntries>=10000:
        break



#output the csts
os.makedirs(polar_cst_paths,exist_ok=True)
c=1
for cstf in motifs:
    cstfn='cst_'+str(c)+'.cst'
    of=open(os.path.join(polar_cst_paths,cstfn),'w')
    for block in cstf:
        of.write('CST::BEGIN\n')
        of.write('\n'.join(block))
        of.write('\n')
        of.write('CST::END\n')
    of.close()
    c+=1





'''
ADD EXPLICIT CSTS TO EXISTING POLAR CSTS
'''

##################################################
##################################################
#gotta give paths to relevant files/directories
#the main dirs for fragment clusters for all fragments we want to add interaxns to
#ligand params 
#output directory 
#the hbbond motifs we've just made
##################################################
##################################################
fragment_fuzzballs_dirs=['38e/residue_contact_fuzzballs/38e_Fragment_3_residue_contacts']
ligparams=['38e/Inputs/Rosetta_Inputs/38e.params']
motifsoutdir='npmotifpdbs'
polar_cst_paths='38e/motifs'
############
##################################################
##################################################
##################################################
#rosetta score thresholds for accepting each number of added cluster intneracxns
single_contact_threshold=-2.0
double_contact_score_threshold=-5.0
triple_contact_score_threshold=-9.0
##################################################
##################################################
##################################################


#
import os
import numpy as np
import math
from pyrosetta import *
init('-pack_missing_sidechains False')
import random
#####
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep,fa_elec#, fa_sol,hbond_sc, fa_elec, hbond_bb_sc,lk_ball_iso
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
# sf.set_weight(fa_sol, 1)
# sf.set_weight(hbond_sc, 1)
# sf.set_weight(fa_elec, 1)
# sf.set_weight(hbond_bb_sc,1)
min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
mm4060 = pyrosetta.rosetta.core.kinematics.MoveMap()
mm4060.set_bb(False)
mm4060.set_chi(1,True)
# mm4060.set_chi(1,True)
# mm4060.set_jump(1,True)
min_mover.movemap(mm4060)
min_mover.score_function(sf)
# min_mover.max_iter(1)
min_mover.tolerance(1e-6)

#we make sure we are in the right working directory
if not os.getcwd()==polar_cst_paths:
    os.chdir(polar_cst_paths)
#make the nop motif pdb output directory
os.makedirs(motifsoutdir,exist_ok=True)





'''
only use nonpolar res
'''
npres=['THR','VAL','LEU','ILE','PHE','TYR','TRP','MET']
#get the paths of the np residue fuzzball pdbs to open
clust_paths=[]
for fragment_fuzzballs_dir in fragment_fuzzballs_dirs:
    for i in os.listdir(fragment_fuzzballs_dir):
        if os.path.isdir(os.path.join(fragment_fuzzballs_dir,i))==True:
            for clust in os.listdir(os.path.join(fragment_fuzzballs_dir,i)):
                if clust[-3:]=='pdb':
                    resn=clust.split('_')[1]
                    if resn in npres:
                        clust_paths.append(os.path.join(fragment_fuzzballs_dir,i,clust))


'''
use all res
'''
#get the paths of the residue fuzzball pdbs to open
# clust_paths=[]
# for fragment_fuzzballs_dir in fragment_fuzzballs_dirs:
#     for i in os.listdir(fragment_fuzzballs_dir):
#         if os.path.isdir(os.path.join(fragment_fuzzballs_dir,i))==True:
#             for clust in os.listdir(os.path.join(fragment_fuzzballs_dir,i)):
#                 if clust[-3:]=='pdb':
#                     clust_paths.append(os.path.join(fragment_fuzzballs_dir,i,clust))
'''
use only hbond res for fragment 2
delete path to frag2 clusters directory
in fragment_fuzzballs_dirs list and do this instead
'''
hbond_clusts_path='38e/residue_contact_fuzzballs/38e_Fragment_2_residue_contacts/hbond_clusts'
for clust in os.listdir(hbond_clusts_path):
    if clust[-3:]=='pdb':
        clust_paths.append(os.path.join(hbond_clusts_path,clust))


#get the populations
clust_paths_pops=[]
for clustpath in clust_paths:
    pop=int(clustpath.split('/')[-1].split('.')[0].split('_')[3])
    clust_paths_pops.append((clustpath,pop))


#sort by population across all fragments
nps=sorted(clust_paths_pops, reverse=True, key=lambda nmem: nmem[1])


#functions to get distance between 2 points
def displace(p1,p2):
    x = p1[0] - p2[0]
    y = p1[1] - p2[1]
    z = p1[2] - p2[2]
    return (x,y,z)
def norm(x):
    return math.sqrt(sum(i**2 for i in x))
def dist(p1,p2):
    v = displace(p1,p2)
    return norm(v)
#function to get angle between 3 points
def getAngle(p1, p2, p3):
    a=np.array(p1)
    b=np.array(p2)
    c=np.array(p3)
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)
#function to get dihedral from 4 points
def dihedral(i1,i2,i3,i4):
    p0=np.array(i1)
    p1=np.array(i2)
    p2=np.array(i3)
    p3=np.array(i4)
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b0xb1 = np.cross(b0, b1)
    b1xb2 = np.cross(b2, b1)
    b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)
    y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
    x = np.dot(b0xb1, b1xb2)
    return np.degrees(np.arctan2(y, x))

#this function i will use just to add np cst at end
def generate_single_constraint_block_base(distanceAB,angleA,angleB,torsionA,torsionAB,torsionB,
                                          residue_resname,ligand_resname,ligand_atoms,residue_atoms,
                                          distance_tolerance=0.25, angle_A_tolerance=15, angle_B_tolerance=15,
                                          torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                          torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                          distance_constraint_sample_number=2):
    """
    Generate a single constraint block for one residue-ligand interaction
    :return:
    """

    residue_tag = 'residue3'


    constraint_block = [
        '  TEMPLATE::   ATOM_MAP: 1 atom_name: {}'.format(' '.join(ligand_atoms)),
        '  TEMPLATE::   ATOM_MAP: 1 residue3: {}\n'.format(ligand_resname),
        '  TEMPLATE::   ATOM_MAP: 2 atom_name: {}'.format(' '.join(residue_atoms)),
        '  TEMPLATE::   ATOM_MAP: 2 {}: {}\n'.format(residue_tag, residue_resname),
        '  CONSTRAINT:: distanceAB: {0:7.2f} {1:6.2f} {2:6.2f}       0  {3:3}'.format(
            distanceAB, distance_tolerance, 100, distance_constraint_sample_number),
        '  CONSTRAINT::    angle_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(angleA), angle_A_tolerance, 100, angle_constraint_sample_number),
        '  CONSTRAINT::    angle_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(angleB), angle_B_tolerance, 100, angle_constraint_sample_number),
        '  CONSTRAINT::  torsion_A: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(torsionA), torsion_A_tolerance, 100,
            torsion_constraint_sample_number),
        '  CONSTRAINT::  torsion_B: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(torsionB), torsion_B_tolerance, 100,
            torsion_constraint_sample_number),
        '  CONSTRAINT:: torsion_AB: {0:7.2f} {1:6.2f} {2:6.2f}  360.00  {3:3}'.format(
            float(torsionAB), torsion_AB_tolerance, 100,
            torsion_constraint_sample_number)]
    return constraint_block





#look at top x clusters and score (check for clashes), retain for motifs if meet score threshold 
np_cst_blocks=[]
accepted_clusts=[]
ac=0
pdboutpaths=[]
clustnames_blocks=[]
for clust,freq in nps[:50]: ####################### the number here is the top X most populated clusters you want to include 
    p = Pose()
    generate_nonstandard_residue_set(p,ligparams)
    pose_from_file(p, clust)
    ligand_resnum=p.total_residue()
    clust_cont_scores=[]
    for indd in range(1,p.total_residue()):
        contact_pose=Pose()
        ligand_pose=p.residue(ligand_resnum).clone()
        res_pose=p.residue(indd).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        res_resname=contact_pose.residue(1).name()[0:3]
        lig_resname=contact_pose.residue(2).name()[0:3]
        if res_resname==lig_resname:
            print('uh oh')
            continue
        else:
            pass
        # min_mover.apply(contact_pose)
        score_this_one=sf(contact_pose)
        clust_cont_scores.append((indd,score_this_one))
    sccs=sorted(clust_cont_scores, key=lambda nmem: nmem[1])
    if sccs[0][1]<=single_contact_threshold:
        accepted_clusts.append((clust,freq,sccs[0][1]))
        ac+=1
        contact_pose=Pose()
        ligand_pose=p.residue(ligand_resnum).clone()
        res_pose=p.residue(sccs[0][0]).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        pdboutpath=os.path.join(motifsoutdir,'npmotif_'+str(ac)+'.pdb')
        contact_pose.dump_pdb(pdboutpath)
        pdboutpaths.append(pdboutpath)
        res_resname=contact_pose.residue(1).name()[0:3]
        lig_resname=contact_pose.residue(2).name()[0:3]
        #get x,y,z coordinates for every atom in residue and ligand
        ligand_atoms_xyz={}#atomindex=(x,y,z,index)
        residue_atoms_xyz={}
        n_residue_atoms=contact_pose.residue(1).natoms()
        n_ligand_atoms=contact_pose.residue(2).natoms()
        for k in range(1,n_ligand_atoms):
            x,y,z=contact_pose.residue(2).atom(k).xyz()
            ligand_atoms_xyz[(contact_pose.residue(2).atom_name(k)).strip()]=(x,y,z,k)
        for j in range(1,n_residue_atoms):
            x,y,z=contact_pose.residue(1).atom(j).xyz()
            residue_atoms_xyz[(contact_pose.residue(1).atom_name(j)).strip()]=(x,y,z,j)
        #find 2 atoms with shortest distance, will define atom1 for each res in constraint block
        distances=[]
        for key in ligand_atoms_xyz.keys():
            p1=ligand_atoms_xyz[key][:3]
            index1=ligand_atoms_xyz[key][3]
            for key2 in residue_atoms_xyz.keys():
                p2=residue_atoms_xyz[key2][:3]
                index2=residue_atoms_xyz[key2][3]
                d=dist(p1,p2)
                distances.append((key,key2,d,index1,index2))
        sd=sorted(distances, key=lambda x: x[2])
        #remove hydrogens
        for r in range(10):#not clear to me why but need to repeat to be thorough
            for i in sd:
                if 'H' in i[0] or 'H' in i[1]:
                    sd.remove(i)
        ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[0]
        if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
            try:
                ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[1]
                if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                    try:
                        ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[2]
                        if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                            try:
                                ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[3]
                                if residue_atom1=='N' or residue_atom1=='C' or residue_atom1=='CA' or residue_atom1=='O':
                                    try:
                                        ligand_atom1,residue_atom1,distance_AB,indexlig,indexres=sd[4]
                                    except:
                                        pass
                            except:
                                pass
                    except:
                        pass
            except:
                pass
        #now find base atoms for res and lig, will be atoms 2 and 3 for each
        ligatomsequence=[]
        bondedto1=[]
        for ia in range(1,n_ligand_atoms):
            atom1 = AtomID(indexlig, 2)
            if ia!=indexlig:
                atom2 = AtomID(ia, 2)
                if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                    bondedto1.append(ia)
        for x in bondedto1:
            for ia2 in range(1,n_ligand_atoms):
                atom1 = AtomID(x, 2)
                if ia2!=x and ia2!=indexlig:
                    atom2 = AtomID(ia2, 2)
                    if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                        if 'H' in (contact_pose.residue(2).atom_name(ia2)).strip():
                            continue
                        else:
                            ligatomsequence.append((indexlig,x,ia2))
        resatomsequence=[]
        bondedto1r=[]
        for ia in range(1,n_residue_atoms):
            atom1 = AtomID(indexres, 1)
            if ia!=indexres:
                atom2 = AtomID(ia, 1)
                if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                    bondedto1r.append(ia)
        for x in bondedto1r:
            for ia2 in range(1,n_residue_atoms):
                atom1 = AtomID(x, 1)
                if ia2!=x and ia2!=indexres:
                    atom2 = AtomID(ia2, 1)
                    if contact_pose.conformation().is_bonded(atom1,atom2)==True:
                        if 'H' in (contact_pose.residue(1).atom_name(ia2)).strip():
                            continue
                        else:
                            resatomsequence.append((indexres,x,ia2))
    #
        if len(ligatomsequence)>0:
            ligbase1=ligatomsequence[0][1]
            ligbase2=ligatomsequence[0][2]
            ligand_atom2=(contact_pose.residue(2).atom_name(ligbase1)).strip()
            ligand_atom3=(contact_pose.residue(2).atom_name(ligbase2)).strip()
        else:
            continue
        if len(resatomsequence)>0:
            resbase1=resatomsequence[0][1]
            resbase2=resatomsequence[0][2]
            residue_atom2=(contact_pose.residue(1).atom_name(resbase1)).strip()
            residue_atom3=(contact_pose.residue(1).atom_name(resbase2)).strip()
        else:
            continue
        #fix oxt to O in res if present
        if residue_atom1=='OXT':
            residue_atom1='O'
        if residue_atom2=='OXT':
            residue_atom2='O'
        if residue_atom3=='OXT':
            residue_atom3='O'
        #save res and ligand atoms
        lig_atoms=[ligand_atom1,ligand_atom2,ligand_atom3]
        res_atoms=[residue_atom1,residue_atom2,residue_atom3]
        if len(set(res_atoms))<3:
            continue
        if len(set(lig_atoms))<3:
            continue
        res_atom_coords=[]
        lig_atom_coords=[]
        x,y,z=contact_pose.residue(1).atom(indexres).xyz()
        res_atom_coords.append((x,y,z))
        x,y,z=contact_pose.residue(1).atom(resbase1).xyz()
        res_atom_coords.append((x,y,z))
        x,y,z=contact_pose.residue(1).atom(resbase2).xyz()
        res_atom_coords.append((x,y,z))
        x,y,z=contact_pose.residue(2).atom(indexlig).xyz()
        lig_atom_coords.append((x,y,z))
        x,y,z=contact_pose.residue(2).atom(ligbase1).xyz()
        lig_atom_coords.append((x,y,z))
        x,y,z=contact_pose.residue(2).atom(ligbase2).xyz()
        lig_atom_coords.append((x,y,z))
        #okay getting angles
        #RES1 IN CONSTRAINT FILE IS GONNA BE LIGAND
        angle_A=getAngle(lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0])
        angle_B=getAngle(lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1])
        #finally, getting dihedrals
        torsion_A=dihedral(lig_atom_coords[2],lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0])
        torsion_AB=dihedral(lig_atom_coords[1],lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1])
        torsion_B=dihedral(lig_atom_coords[0],res_atom_coords[0],res_atom_coords[1],res_atom_coords[2])
        #get constraint block
        cstblock=generate_single_constraint_block_base(distance_AB,angle_A,angle_B,torsion_A,torsion_AB,torsion_B,
                                                       res_resname,lig_resname,lig_atoms,res_atoms,
                                                       distance_tolerance=0.25, angle_A_tolerance=15, angle_B_tolerance=15,
                                                       torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                                       torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                       distance_constraint_sample_number=2)
        for wp in range(freq):
            np_cst_blocks.append(cstblock)
        clustnames_blocks.append((pdboutpath,cstblock))



'''
np_cst_blocks list now holds the clusters we will add to motifs 





ADD 1 EXPLICIT CST
'''

if os.path.exists('1np'):
    os.system('rm -r 1np')
    os.mkdir('1np')
else:
    os.mkdir('1np')
os.chdir('1np')
#to add 1 np cst to csts
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(1):          #this is how many unique hybrid cst i want to make per polar cst
            rblock=random.choice(np_cst_blocks)
            if rblock not in bs: #check that unique 
                bs.append(rblock)
                f=open(os.path.join(polar_cst_paths,cst),'r')
                lines=[line for line in f.readlines()]
                f.close()
                ofname='hybrid_'+str(c)+'.cst'
                of=open(ofname,'w')
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblock))
                of.write('\n')
                of.write('CST::END\n')
                for line in lines:
                    of.write(line)
                of.close()
                print(c)
                c+=1
            else: 
                rblock=random.choice(np_cst_blocks)
                if rblock not in bs:
                    bs.append(rblock)
                    f=open(os.path.join(polar_cst_paths,cst),'r')
                    lines=[line for line in f.readlines()]
                    f.close()
                    ofname='hybrid_'+str(c)+'.cst'
                    of=open(ofname,'w')
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblock))
                    of.write('\n')
                    of.write('CST::END\n')
                    for line in lines:
                        of.write(line)
                    of.close()
                    print(c)
                    c+=1
                else:
                    pass


os.chdir('..')
















#to make 2 residue motifs, we gotta score all 2 res comobos
#with double interaxn threshold we defined 
#can take a few minutes 
clust_cont_scores=[]
for cind,clust in enumerate(pdboutpaths):
    p = Pose()
    generate_nonstandard_residue_set(p,ligparams)
    pose_from_file(p, clust)
    for clust2 in pdboutpaths[cind+1:]:
        p2 = Pose()
        generate_nonstandard_residue_set(p2,ligparams)
        pose_from_file(p2, clust2)
        #
        contact_pose=Pose()
        res_pose=p.residue(1).clone()
        res_pose2=p2.residue(1).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(res_pose2, 1)
        ligand_pose=p.residue(2).clone()
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        score_this_one=sf(contact_pose)
        clust_cont_scores.append((clust,clust2,score_this_one))


sccs=sorted(clust_cont_scores, key=lambda nmem: nmem[2])
cstblockss=[]
for x in sccs:
    if x[2]<=double_contact_score_threshold:
        cstblocks=[]
        for y in clustnames_blocks:
            if x[0]==y[0]:
                cstblocks.append(y[1])
            elif x[1]==y[0]:
                cstblocks.append(y[1])
        cstblockss.append((cstblocks[0],cstblocks[1]))

#cstblockss
'''
len(cstblockss)
'''

os.mkdir('2np')
os.chdir('2np')
c=1
for cst in os.listdir(polar_cst_paths):
    if cst[-3:]=='cst':
        bs=[]
        for i in range(2):          #this is how many unique hybrid cst i want to make per polar cst
            rblock=random.choice(cstblockss)
            if rblock not in bs:
                bs.append(rblock)
                f=open(os.path.join(polar_cst_paths,cst),'r')
                lines=[line for line in f.readlines()]
                f.close()
                ofname='hybrid_'+str(c)+'.cst'
                of=open(ofname,'w')
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblock[0]))
                of.write('\n')
                of.write('CST::END\n')
                of.write('\n')
                of.write('CST::BEGIN\n')
                of.write('\n'.join(rblock[1]))
                of.write('\n')
                of.write('CST::END\n')
                for line in lines:
                    of.write(line)
                of.close()
                print(c)
                c+=1
            else:
                rblock=random.choice(cstblockss)
                if rblock not in bs:
                    bs.append(rblock)
                    f=open(os.path.join(polar_cst_paths,cst),'r')
                    lines=[line for line in f.readlines()]
                    f.close()
                    ofname='hybrid_'+str(c)+'.cst'
                    of=open(ofname,'w')
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblock[0]))
                    of.write('\n')
                    of.write('CST::END\n')
                    of.write('\n')
                    of.write('CST::BEGIN\n')
                    of.write('\n'.join(rblock[1]))
                    of.write('\n')
                    of.write('CST::END\n')
                    for line in lines:
                        of.write(line)
                    of.close()
                    print(c)
                    c+=1
                else:
                    pass







os.chdir('..')


















###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
# 6. MATCH
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################

            matcher_cmd = ['~/main/source/bin/match.linuxgccrelease',
                           '-s', '${tasks[$SGE_TASK_ID]}.pdb',  ###########
                           '-extra_res_fa', paramspath,
                           '-match:geometric_constraint_file', cst_path,
                           '-match::scaffold_active_site_residues', '${tasks[$SGE_TASK_ID]}.pos',   #########
                           '-match:output_format', 'PDB',
                           '-match:match_grouper', 'SameSequenceGrouper',
                           '-match:consolidate_matches',
                           '-match:output_matches_per_group', '1',
                           '-use_input_sc',
                           '-ex1', '-ex2','-extrachi_cutoff 0',
                           '-enumerate_ligand_rotamers', 'false',
                           '-match::lig_name', lig_name,
                           '-load_PDB_components', 'false']
   



####
# 7.filter for burial
###
'''
we apply /tools/match_analysis_standard.xml script to score with rosetta and evaluate ligand burial 
in matches directory:
'''

        cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
             '-s',matchpath,         ##
             '-parser:protocol','tools/match_analysis_standard.xml',
             '-parser:view','-run:preserve_header',
             '-in:file::extra_res_fa',prm1,
             '-load_PDB_components','False',
             '-scorefile_format','json',
             '-out:file:score_only',os.path.join(outputdir,scorefilename_id)]
        cmdd=' '.join(cmd)


'''

NOW FILTER THE WELL BURIED ONES


'''
#analysis of scores from json file
import os
import json
jsonoutputdir=os.getcwd()
sfs=[os.path.join(jsonoutputdir,i) for i in os.listdir(jsonoutputdir) if i[-4:]=='json']
c=1
scores=[]
for sf in sfs:
    f=open(sf,'r')
    lines=[line for line in f.readlines()]
    f.close()
    for line in lines:
        scores.append(json.loads(line))
    print(str(c))
    c+=1
terms=list(scores[0].keys())


def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores
###############
###############
###############
f1=return_filtered(scores,'dsasa','>',0.7)
print(len(scores))
print(len(f1))
###############
###############
###############

import os
filtered_strc=[]
for d in f1:
    filtered_strc.append(d['decoy'])













































#################################################################################
#################################################################################
# 8.enzyme design application
#################################################################################
#################################################################################

#resfile generation
'''
we use  tools/match_resfile2.py script 

assuumes matches have typical naming scheme UM_X_D45N84L97

takes command line arguments INPUTMATCH.PDB TARGETROSETTAPARAMS.PARAMS DISTANCE ISGENPOT
distance is the distance from ligand that design is allowed (ie if 5.0, all residues within 5 angstroms of ligand designable)
isgenpot specifies whether we are using standard params or genpot params file, takes 'genpot' (if yes) or 'no'

makes match resfile with code notaa cp (for designable residues, allow design to anything except cysteine/proline)
'''
#################################################################################
#################################################################################
#################################################################################
#################################################################################

cmdd='time python /tools/match_resfile2.py '+matchpath+' '+prm1+' 5.0 no'



# # ROSETTA ENZYME DESIGN APPLICATION
    cmd=['~/main/source/bin/enzyme_design.linuxgccrelease',
         '-s','${tasks[$SGE_TASK_ID]}.pdb',         ##
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-ex1','-ex2','-extrachi_cutoff','0',
         '-enzdes:cst_opt','-enzdes:bb_min',
         '-enzdes:chi_min',
         '-enzdes:cst_design','-enzdes:design_min_cycles 3', '-enzdes:lig_packer_weight 2.5',
         '-enzdes:cst_min',
         '-out:nstruct',ndesignseachmatch,
         # '-score::weights beta_genpot', #if we are using gen potential
         # '-corrections::gen_potential',
         '-load_PDB_components','False',
         '-enzdes:cstfile','${tasks[$SGE_TASK_ID]}.cst',
         '-out:prefix',output_prefix]





















#########
# 9.hbrefine
########

'''
tools/hbrefine.py script 

tries to build hbonds with unsat ligand polar atoms 
and to remove unsat polar amino acids in binding site

use:
'''
time python tools/hbrefine.py paramspath.params pdbpath.pdb outputdirectorypath

















#


'''
scoring designs
~/BSFF/tools/edanalysis_std_sfn.xml script 

can use extratools/edanalysis.xml with genpot

calculates fast metrics that we can run relatively quickly on a lot of designs (>1 million if needed)

'''
 cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
             '-s',matchpath,
             # '-parser:protocol','tools/edanalysis.xml',
             '-parser:protocol','tools/edanalysis_std_sfn.xml',
             '-parser:view','-run:preserve_header',
             '-in:file::extra_res_fa',prm1,
             # '-score::weights beta_genpot',
             # '-corrections::gen_potential',
             '-corrections::beta_nov16',
             '-load_PDB_components','False',
             '-scorefile_format','json',
             '-out:file:score_only',scorefile_name]


'''
now on smaller filtered set of designs 
we can calculate more involved metrics 
that take longer per structure (around an hour per input is typical)

tools/ligbinderanalysis_jump1.xml script

this requires adding a few options/flags 
'''

'main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:s='+matchpath+' -parser:protocol tools/ligbinderanalysis_jump1.xml -in:file::extra_res_fa '+params+' -ignore_unrecognized_res -load_PDB_components False -scorefile_format json -out:file:score_only '+scorefile_name+' -run:nblist_autoupdate -parser:view -jd2:ntrials=1 -packing:no_optH=false -packing:flip_HNQ -packing:extrachi_cutoff=1 -packing:use_input_sc -packing:linmem_ig=10 -packing:ex1 -packing:ex2 -corrections:score:no_his_his_pairE -corrections:score:lj_hbond_hdis=1.75 -corrections:score:lj_hbond_OH_donor_dis=2.6 -corrections:beta_nov16 -enzdes:bb_min_allowed_dev=0.05 -enzdes:detect_design_interface -enzdes:cut1=4 -enzdes:cut2=6 -enzdes:cut3=8 -enzdes:cut4=10')





#analysis of scores from json file
import os
import json
jsonoutputdir=os.getcwd()
sfs=[os.path.join(jsonoutputdir,i) for i in os.listdir(jsonoutputdir) if i[-4:]=='json']
c=1
scores=[]
for sf in sfs:
    f=open(sf,'r')
    lines=[line for line in f.readlines()]
    f.close()
    for line in lines:
        scores.append(json.loads(line))
    print(str(c))
    c+=1
terms=list(scores[0].keys())



def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores

#
f1=return_filtered(scores,'buns2interface','<',0.0)
f2=return_filtered(f1,'boltz','<',-0.3)
f3=return_filtered(f2,'ddg','<',-20.0)
f4=return_filtered(f3,'contact_molsurf','>',150)
f5=return_filtered(f4,'hbtolig','>',2.0)
f6=return_filtered(f5,'shape_comp','>',0.65)
f7=return_filtered(f6,'lighphobesasa','<',20.0)
f8=return_filtered(f7,'ligoversat','<',0)

























'''
RUNNING MPNN ON MATCHES TO GENERATE SEQUENCE PROFILE
design at all positions not in binding site first shell

EXAMPLE OF FIXED POSITIONS DICT
{"3HTN": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": []}, "4YOW": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": [], "D": [], "E": [], "F": []}}

'''
# get binding site residues to hold constant

prm1='/Inputs/Rosetta_Inputs/dog.params'
ofjsonname='mpnn_params.json'
################
import os
from pyrosetta import*
import json
init('-load_PDB_components False')
##########
sf=get_fa_scorefxn()

filters_xml = f'''
                <FILTERS>
                    <DSasa name="dsasa" confidence="0.0"/>
                </FILTERS>'''
dsasa_filter = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(filters_xml).get_filter("dsasa")
#
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
tfdata={}

for pdb in pdbs:
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    p.update_residue_neighbors()
    dsasa_val=dsasa_filter.report_sm(p)
    if dsasa_val>0.8:
        tofreeze=[]
        #######################################################
        #######################################################
        #######################################################
        #######################################################
        #######################################################
        #######################################################
        sss=pdb.split('_')[2]
        ss=''
        for i in sss:
            if i.isdigit()==True:
                ss+=str(i)
            else:
                ss+='_'
        ds=ss.split('_')
        dss=[i for i in ds if i.isdigit()==True]
        for i in dss:
            tofreeze.append(int(i))
        ####################
        ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
        neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
        neighborhood_selector_bool = neighborhood_selector.apply(p)
        neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
        first_shell_res=list(neighborhood_residues_resnums)
        ###############
        for resnum in range(1,p.total_residue()):
            if resnum not in first_shell_res:
                if resnum not in tofreeze:
                    tofreeze.append(resnum)
        tfdata[pdb]=tofreeze

print(len(list(tfdata.keys())))
json.dump(tfdata,open(ofjsonname,'w'))
###############
###############
###############
###############
###############
###############
###############




#make the fixed position dictionaries and put them in the
#proper directories
#then
#make the shellfile scripts for mpnn jobs:
'''
{"3HTN": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": []}, "4YOW": {"A": [1, 2, 3, 4, 5, 6, 7, 8, 23, 25], "C": [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 40], "B": [], "D": [], "E": [], "F": []}}
'''

    of.write('#!/bin/bash')
    of.write('\n')
    of.write('\n')
    of.write('folder_with_pdbs="'+id+'"')
    of.write('\n')
    of.write('output_dir="'+od+'"')
    of.write('\n')
    of.write('if [ ! -d $output_dir ]')
    of.write('\n')
    of.write('then')
    of.write('\n')
    of.write('    mkdir -p $output_dir')
    of.write('\n')
    of.write('fi')
    of.write('\n')
    of.write('path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"')
    of.write('\n')
    of.write('path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"')
    of.write('\n')
    of.write('python ProteinMPNN-main/vanilla_proteinmpnn/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains')
    of.write('\n')
    cmdl=['python',
    'ProteinMPNN-main/vanilla_proteinmpnn/protein_mpnn_run.py',
    '--jsonl_path $path_for_parsed_chains',
    '--out_folder $output_dir',
    '--fixed_positions_jsonl $path_for_fixed_positions',
    '--num_seq_per_target 500',
    '--sampling_temp "0.15"',
    '--batch_size 1']
    cmd=' '.join(cmdl)




#use results to make resfiles that allow design to 
#residues seen in mpnn sequence profile 




'''
Rosetta
FASTDESIGN with mpnn sequence profile
some filters but no ddg/boltz

~tools/fd_mpnn.xml script 
'''
    cmd=['~/main/source/bin/rosetta_scripts.default.linuxgccrelease',
         '-s','${tasks[$SGE_TASK_ID]}.pdb',         ##
         '-parser:protocol','tools/fd_mpnn.xml script ',
         '-parser:view','-run:preserve_header',
         '-in:file::extra_res_fa',prm1,
         '-resfile','${tasks[$SGE_TASK_ID]}.resfile',  ##
         '-ex1','-ex2','-extrachi_cutoff','0','-use_input_sc',
         '-out:nstruct',ndesignseachmatch,
         '-beta_nov16',
         '-load_PDB_components','False',
         '-out:path:all',outputdirectory,
         '-out:prefix',output_prefix,
         '-scorefile_format','json',
         '-out:file:scorefile',scorefile_name]







'''

SCORES

'''


#analysis of scores from json file
import os
import json
jsonoutputdir=os.getcwd()
sfs=[os.path.join(jsonoutputdir,i) for i in os.listdir(jsonoutputdir) if i[-4:]=='json']
c=1
scores=[]
for sf in sfs:
    f=open(sf,'r')
    lines=[line for line in f.readlines()]
    f.close()
    for line in lines:
        scores.append(json.loads(line))
    print(str(c))
    c+=1
terms=list(scores[0].keys())





def return_filtered(scores,term,condition,threshold):
    filtered_scores=[]
    for d in scores:
        try:
            termval=float(d[term])
            if condition=='<':
                if termval<=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
            elif condition=='>':
                if termval>=float(threshold):
                    filtered_scores.append(d)
                else:
                    pass
        except:
            print(d)
    return filtered_scores

#
f1=return_filtered(scores,'buns2interface','<',1.0)
f2=return_filtered(f1,'contact_molsurf','>',165)
f3=return_filtered(f2,'hbtolig','>',2.0)
f4=return_filtered(f3,'shape_comp','>',0.65)
f5=return_filtered(f4,'lighphobesasa','<',40.0)
f6=return_filtered(f5,'packstat','>',0.6)
f7=return_filtered(f6,'buns_bb_heavy','<',5)
f8=return_filtered(f7,'buns_sc_heavy','<',1)
f9=return_filtered(f8,'ligoversat','<',0)
f10=return_filtered(f9,'oversat','<',0)
f11=return_filtered(f10,'exphyd','<',1200)
f12=return_filtered(f11,'cav','<',140)
































'''
SUBMIT TO COLABFOLD


time python ~/BSFF/tools/run_colabfold.py input_dir output_dir
'''



















'''
ANALYZING COLABFOLD RESULTS ON FILTERED DESIGNS


compare design model and predicted structure using rosetta to grab rmsd 
'''

import os
import numpy as np
import json
from pyrosetta import *
init('-ignore_unrecognized_res')

directory_prefix='filtdesigns'



#####specifying fastas directory and rosetta design model directory 
allfastasdir='bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered/cf/fastas'
nonaf2desdir='bestmatches/enzdes/filtered/analysis/run/filtered2/mpnn/resfiles/fd_mpnn/filtered2/refined/filtered'
#####


#
prediction_pdb_paths=[]
input_directories=[os.path.join(allfastasdir,i) for i in os.listdir(allfastasdir) if os.path.isdir(os.path.join(allfastasdir,i))==True and i[:len(directory_prefix)]==directory_prefix]
for id in input_directories:
    outputdirname=os.path.join(id,id.split('/')[-1]+'_output')
    for i in os.listdir(outputdirname):
        if i[-3:]=='pdb':
            prediction_pdb_paths.append(os.path.join(outputdirname,i))
#
data_allstrc={}
seen_names=[]
ogds=[i for i in os.listdir(nonaf2desdir) if i[-3:]=='pdb']
for p in prediction_pdb_paths:
    # print(str(p))
    f=open(p,'r')
    fuzzball_lines=[line for line in f.readlines() if line[0:4]=='ATOM' or line[:6]=='HETATM']
    f.close()
    #first gotta get indices of unique residues in fuzzball_lines
    fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
    starts=[];lasts=[]
    for index, line in enumerate(fuzzball_lines): #collecting start/end indices for unique res
            try:
                resnum=int(fuzzball_lines[index][22:29].strip())
            except:
                print(index)
            # print((fuzzball_lines[index][22:29].strip()))
            resname=line[17:20]
            try:
                lastresnum=int(fuzzball_lines[index-1][22:29].strip())
                lastresname=fuzzball_lines[index-1][17:20]
                if resnum!=lastresnum or resname!=lastresname:
                    start=index
                    starts.append(start)
            except:
                start=index
                starts.append(start)
            try:
                nextresname=fuzzball_lines[index+1][17:20]
                next_resnum=int(fuzzball_lines[index+1][22:29].strip())
                if resnum!=next_resnum or resname!=nextresname:
                    last=index+1
                    lasts.append(last)
            except:
                last=len(fuzzball_lines)
                lasts.append(last)
    for index,start in enumerate(starts): #put the indices together for each res
        fuzzball_residue_indices.append((start,lasts[index]))
    fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices
    fuzzball_residue_indices=sorted(fuzzball_residue_indices, key=lambda first: first[0])
    #
    #
    #
    currstrc_lddts=[]
    for rin in fuzzball_residue_indices:
        rlddt=fuzzball_lines[rin[0]][60:66].strip(' ')
        currstrc_lddts.append(float(rlddt))
    currstrc_plddt=np.mean(currstrc_lddts)
    ############
    currstrc_name='_'.join(p.split('/')[-1].split('_')[:-5])###########################################
    csnogds=currstrc_name+'.pdb'
    if csnogds in ogds:
        p11=pose_from_pdb(p)
        p22=pose_from_pdb(os.path.join(nonaf2desdir,csnogds))
        carmsd=pyrosetta.rosetta.core.scoring.CA_rmsd(p11,p22)
    ##########
    if currstrc_name not in seen_names:
         seen_names.append(currstrc_name)
         data_allstrc[currstrc_name]={}
    else:
        pass
    #########
    currstrc_model=str(p.split('/')[-1].split('_')[-3])###########################################
    ##################
    data_allstrc[currstrc_name][currstrc_model]=[currstrc_plddt,currstrc_lddts,carmsd]

#output the data so its easy to load later and compare with designs
##################################################################
##################################################################
##################################################################
json.dump(data_allstrc,open('af2_data.json','w'))
##################################################################
##################################################################
##################################################################

import numpy as np
#id designs with plddt over threshold
#use best vals:

plddt_threshold=88.0
carmsd_threshold=1.5
aplddts=[]
accepted={}
for key in data_allstrc.keys():
    plddts=[]
    carmsds=[]
    for k2 in data_allstrc[key]:
        cplddt=data_allstrc[key][k2][0]
        carmsd=data_allstrc[key][k2][2]
        if cplddt>=plddt_threshold:
            if carmsd<=carmsd_threshold:
                if key not in list(accepted.keys()):
                    accepted[key]=[cplddt,carmsd]
                    plddts.append((cplddt,key,k2,carmsd))
                else:
                    if cplddt>accepted[key][0] and carmsd<accepted[key][1]:
                        accepted[key]=[cplddt,carmsd]
                        plddts.append((cplddt,key,k2,carmsd))
    try:
        aplddts.append(plddts[-1])
    except:
        pass

print(len(aplddts))




#















'''
MIGHT THEN WANT TO REPEAT SCORING ON AF2 FILTERED OUTPUTS TO ENSURE RETENTION OF KEY FEATURES 
'''
