

import json
import os
from prody import *
from rdkit import Chem


pdb_search_results_json='PDB_search_results.json'
f=open(pdb_search_results_json,'r')
pdb_search_results_json_data=json.load(f)
f.close()

d_ligands={}
for current_fragment in list(pdb_search_results_json_data.keys()):
    l=[]
    for frag_containing_pdb in pdb_search_results_json_data[current_fragment]['Ligands']:
        l.append(frag_containing_pdb)
    d_ligands[current_fragment]=l

d_ligand_smiles={}
print('getting smiles...')
for current_fragment in d_ligands:
    frag_dict={}
    for pdb_ligand in d_ligands[current_fragment]:
        # print(pdb_ligand)
        try:
            pdb_ligand_smiles_canon = pdb_search_results_json_data[current_fragment]['SMILES'][pdb_ligand]
            frag_dict[pdb_ligand]=pdb_ligand_smiles_canon
        except:
            print('cannot get smiles of ligand in pdb '+pdb_ligand)
            continue
    d_ligand_smiles[current_fragment]=frag_dict



d={}
for current_fragment in list(pdb_search_results_json_data.keys()):
    l=[]
    for frag_containing_pdb in pdb_search_results_json_data[current_fragment]['PDBs']:
        l.append(frag_containing_pdb)
    d[current_fragment]=l

output_dicts=[]
#####################################################################################################################################
max_strc_per_job=150 #######change to your desired value, number of inputs per node
#####################################################################################################################################
for key in d.keys():
    n_strc_current_frag=len(d[key])
    if n_strc_current_frag<=max_strc_per_job:
        output_dict={}
        sub_dict={}
        sub_dict['PDBs']=d[key]
        sub_dict['Ligands']=d_ligands[key]
        sub_dict['SMILES']=d_ligand_smiles[key]
        output_dict[key]=sub_dict
        output_dicts.append(output_dict)
    else:
        remainder=n_strc_current_frag%max_strc_per_job
        n_divisions=int(n_strc_current_frag/max_strc_per_job)
        for i in range(0,(n_divisions)):
            output_dict={}
            sub_dict={}
            strc_selection=list(d[key][(i*max_strc_per_job):((i+1)*max_strc_per_job)])
            sub_dict['PDBs']=strc_selection
            sub_dict['Ligands']=d_ligands[key]
            sub_dict['SMILES']=d_ligand_smiles[key]
            output_dict[key]=sub_dict
            output_dicts.append(output_dict)
        output_dict={}
        sub_dict={}
        strc_selection=list(d[key][(n_strc_current_frag-remainder):(n_strc_current_frag)])
        sub_dict['PDBs']=strc_selection
        sub_dict['Ligands']=d_ligands[key]
        sub_dict['SMILES']=d_ligand_smiles[key]
        output_dict[key]=sub_dict
        output_dicts.append(output_dict)

os.makedirs('Alignment_Inputs',exist_ok=True)
count=1
input_filesnames=[]
for dict in output_dicts:
    filename='Alignment_Inputs/PDB_search_results_'+str(count)+'.json'
    with open(filename, 'w') as fp:
        json.dump(dict, fp)
    input_filesnames.append(filename)
    count+=1



###################designed for my person account on wynton ucsf cluster, outputs shell files
###################change for your account/cpu system 
##################################################################################################################
job_shellfile='run_align.sh'
jsf=open(job_shellfile,'w')
jsf.write('#!/bin/bash')
jsf.write('\n#$ -S /bin/bash')
jsf.write('\nsource ~/anaconda3/etc/profile.d/conda.sh')
jsf.write('\nconda activate pyr37')
jsf.write('\ntrap "conda deactivate" EXIT')
jsf.write('\ntasks=(0\n')
for filename in input_filesnames[:-1]:
    jsf.write('       '+str(filename)+'\n')
jsf.write('       '+input_filesnames[-1]+')')
jsf.write('\npython ~/BSFF/extract_align_contacts.py ${tasks[$SGE_TASK_ID]}')
jsf.close()

print(str(len(input_filesnames)))
###############################################################################################




