#time python ~/bsff_scripts/setup_filter_array_job.py

import sys
import os
from shutil import copy2

parent_dir=os.getcwd()
tapdbs=os.path.join(parent_dir,'Transformed_Aligned_PDBs')
fragment_dirs=[i for i in os.listdir(tapdbs) if os.path.isdir(os.path.join(tapdbs,i))==True]

max_strc_per_job=500
input_directory_paths=[]
for fragment_dir in fragment_dirs:
    print(fragment_dir)
    #get the names of the pdb files in working directory
    pdbfiles=[file for file in os.listdir(os.path.join(tapdbs,fragment_dir)) if file[-3:]=='pdb']
    n_pdbs=len(pdbfiles)
    if n_pdbs<=max_strc_per_job:
        input_directory_paths.append(os.path.join(tapdbs,fragment_dir))
    else:
        n_divisions=int(n_pdbs/max_strc_per_job)
        remainder=n_pdbs%max_strc_per_job
        c=0
        for i in range(n_divisions):
            os.makedirs(os.path.join(tapdbs,fragment_dir,'filter_input',str(i+1)),exist_ok=True)
            input_directory_paths.append(os.path.join(tapdbs,fragment_dir,'filter_input',str(i+1)))
            pdbs_to_copy=pdbfiles[(i*max_strc_per_job):((i+1)*max_strc_per_job)]
            for x in pdbs_to_copy:
                copy2(os.path.join(tapdbs,fragment_dir,x),os.path.join(tapdbs,fragment_dir,'filter_input',str(i+1),x))
            c=i
        os.makedirs(os.path.join(tapdbs,fragment_dir,'filter_input',str(c+2)),exist_ok=True)
        input_directory_paths.append(os.path.join(tapdbs,fragment_dir,'filter_input',str(c+2)))
        pdbs_to_copy=pdbfiles[(n_pdbs-remainder):(n_pdbs)]
        for x in pdbs_to_copy:
            copy2(os.path.join(tapdbs,fragment_dir,x),os.path.join(tapdbs,fragment_dir,'filter_input',str(c+2),x))



job_shellfile='run_filter_array.sh'
jsf=open(job_shellfile,'w')
jsf.write('#!/bin/bash')
jsf.write('\nsource ~/anaconda3/etc/profile.d/conda.sh')
jsf.write('\nconda activate pyr37')
jsf.write('\ntasks=(0\n')
for filename in input_directory_paths[:-1]:
    jsf.write('       '+str(filename)+'\n')
jsf.write('       '+input_directory_paths[-1]+')')
jsf.write('\npython ~/BSFF/filter_contacts_array.py ${tasks[$SGE_TASK_ID]}')
jsf.close()

print('\n\n\n\n\n'+str(len(input_directory_paths)))
print(str(parent_dir))

'''
2575
591
97
286
13

x=0
for path in input_directory_paths:
    print(len(os.listdir(path)))
    x+=len(os.listdir(path))

print(x)

for fragment_dir in fragment_dirs:
    path=os.path.join(tapdbs,fragment_dir,'filter_input')
    try:
        os.system('rm -r '+path)
    except:
        pass

37
'''
