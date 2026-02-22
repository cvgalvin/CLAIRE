#time ipython ~/desktop/tools/fragment_quality.py inputs_json wt.pdb wt9mers jsonoutputname des.pdb
#script assumes 200 fragments in frag libraries
#inputs_json = dictionary structure_path:designable positions
#reports wildtype value - design value

import sys
import json
# import os
import numpy as np
import pandas as pd
from pyrosetta import *
init('-ignore_unrecognized_res')
#from  pyrosetta.rosetta.core.scoring import *

#manually input designed positions
#store pdb numbering of binding site residues for query protein
designable_positions_json=sys.argv[1]
outputjsonname=sys.argv[4]
designs=[sys.argv[5]]

#load ninemers
with open(designable_positions_json,'r') as f:
    desinfo=json.load(f)
fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
fragset9.read_fragment_file(list(desinfo.keys())[0])
#load design pdb and params file
designed_positions=desinfo[(list(desinfo.keys())[0])]


#load wildtype scaffold pdb and 9mers
wt_ninemers=pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
wt_ninemers.read_fragment_file(sys.argv[3])
wt_scaffolds=[sys.argv[2]]  #this doesnt need to be a list im just lazy

#open text file to write results
outfilename='frag_filters.txt'
of=open(outfilename,'a')



#to store stuff for pandas df
resultslist=[]
resultsdict={}
resultsdict['SCORE:']='SCORE:'
#get the rmsd of design position containing fragments in wildtype structure so
#that they may be compared with designs
#mostly just copying code from below (i wrote the block for design evaluation first)
#and changing names to wt 9mers/pdb
wt_values_for_comparison=[] #(n_good_frags/total_frags,coverage,(pos,avg_rmsd_pos),avg_rmsd_all_pos)
for wt_scaffold in wt_scaffolds:
    of.write('\n'+'WILD TYPE SCAFFOLD'+'\n')
    of.write('\n'+str(wt_scaffold[:-4].upper())+'\n')
    p=pose_from_pdb(wt_scaffold)
    #run fragqual filter on the whole ptn for n good frags + fraction of ptn covered by good frags
    calc9=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(wt_ninemers)
    calc9.begin(1)
    calc9.end(p.total_residue()-1)
    calc9.rmsd_cutoff(1.0)
    fragscore9=calc9.get('num_goodfrag',p)
    total_frags=wt_ninemers.size() #total number of fragments
    of.write('The number of good fragments is: '+str(fragscore9)+' out of: '+str(total_frags)+' total fragments')
    wt_values_for_comparison.append(float(fragscore9)/total_frags)
    ###############coverage
    calc9r=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(wt_ninemers)
    calc9r.begin(1)
    calc9r.end(p.total_residue()-1)
    calc9r.ratio_cutoff(0.3)
    fragscore9r=calc9r.get('coverage',p)
    of.write('\nThe fraction of positions covered by good fragments is : '+str(fragscore9r)+'\n')
    wt_values_for_comparison.append(float(fragscore9r))
    ############### design position specific metrics for 9mers
    of.write('\nDesign position specific metrics: '+'\n')
    fr=pyrosetta.rosetta.core.fragment.FragmentRmsd(wt_ninemers)
    designed_position_fragments=[]
    for i in designed_positions:
        for k in range(i-8,i+1):
            if k>0 and k<=p.total_residue()-9:
                designed_position_fragments.append(k)
    designed_position_fragments=list(set(designed_position_fragments))
    #first calculate rmsd for all unique fragments containg des pos, to avoid redundant calculation
    all_rmsd_data={}
    print('\n'+'\n'+'\n'+'\n'+'calculating fragment rmsd values for designed positions for wildtype'.upper())
    for fragment_number in designed_position_fragments:
        design_fraglib_rmsd=[]
        for x in range(1,201):
            curr_rmsd=fr.rmsd(fragment_number,x,p)
            design_fraglib_rmsd.append(curr_rmsd)
            print(str(curr_rmsd))
        all_rmsd_data[fragment_number]=design_fraglib_rmsd
    #now go thru each design position and get the rmsd values for frags its a part of
    for i in designed_positions:
        frags=[]
        pos_all_rmsd=[]
        for k in range(i-8,i+1):
            if k>0 and k<=p.total_residue()-9:
                frags.append(k)
        for k in frags:
            for x in all_rmsd_data[k]:
                pos_all_rmsd.append(x)
        position_avg_rmsd=np.mean(pos_all_rmsd)
        position_max_rmsd=np.max(pos_all_rmsd)
        position_min_rmsd=np.min(pos_all_rmsd)
        position_med_rmsd=np.median(pos_all_rmsd)
        of.write('\nFor Design Position: '+str(i))
        of.write('\nThe Average RMSD is: '+str(position_avg_rmsd))
        of.write('\nThe Median RMSD is: '+str(position_med_rmsd))
        of.write('\nThe Min RMSD is: '+str(position_min_rmsd))
        of.write('\nThe Max RMSD is: '+str(position_max_rmsd)+'\n')
        wt_values_for_comparison.append((i,position_avg_rmsd))
    #now get the average rmsd for all fragments that contain a design position
    all_design_position_means=[]
    for key in all_rmsd_data.keys():
        for value in all_rmsd_data[key]:
            all_design_position_means.append(np.mean(value))
    of.write('\nThe Average RMSD for all fragments containing a designed position is: '+str(np.mean(all_design_position_means)))
    wt_values_for_comparison.append(float(np.mean(all_design_position_means)))


#evaluate designs and compare with wt
for design in designs:
    des_values_for_comparison=[]
    of.write('\n'+'\n'+'DESIGNS'+'\n')
    of.write('\n'+str(design[:-4].upper())+'\n')
    p=pose_from_pdb(design)
    #run fragqual filter on the whole ptn for n good frags + fraction of ptn covered by good frags
    calc9=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
    calc9.begin(1)
    calc9.end(p.total_residue()-1)
    calc9.rmsd_cutoff(1.0)
    fragscore9=calc9.get('num_goodfrag',p)
    total_frags=fragset9.size() #total number of fragments
    of.write('The number of good fragments is: '+str(fragscore9)+' out of: '+str(total_frags)+' total fragments')
    des_values_for_comparison.append(float(fragscore9)/total_frags)
    resultsdict['design_fragqual_good/total']=float(fragscore9)/total_frags
    ###############coverage
    calc9r=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
    calc9r.begin(1)
    calc9r.end(p.total_residue()-1)
    calc9r.ratio_cutoff(0.3)
    fragscore9r=calc9r.get('coverage',p)
    of.write('\nThe fraction of positions covered by good fragments is : '+str(fragscore9r)+'\n')
    des_values_for_comparison.append(float(fragscore9r))
    resultsdict['design_fragqual_coverage']=float(fragscore9r)
    ############### design position specific metrics for 9mers
    of.write('\nDesign position specific metrics: '+'\n')
    fr=pyrosetta.rosetta.core.fragment.FragmentRmsd(fragset9)
    designed_position_fragments=[]
    for i in designed_positions:
        for k in range(i-8,i+1):
            if k>0 and k<=p.total_residue()-9:
                designed_position_fragments.append(k)
    designed_position_fragments=list(set(designed_position_fragments))
    #first calculate rmsd for all unique fragments containg des pos, to avoid redundant calculation
    all_rmsd_data={}
    print('\n'+'\n'+'\n'+'\n'+'calculating fragment rmsd values for designed positions'.upper())
    for fragment_number in designed_position_fragments:
        design_fraglib_rmsd=[]
        for x in range(1,201):
            curr_rmsd=fr.rmsd(fragment_number,x,p)
            design_fraglib_rmsd.append(curr_rmsd)
            print(str(curr_rmsd))
        all_rmsd_data[fragment_number]=design_fraglib_rmsd
    #now go thru each design position and get the rmsd values for frags its a part of
    for i in designed_positions:
        frags=[]
        pos_all_rmsd=[]
        for k in range(i-8,i+1):
            if k>0 and k<=p.total_residue()-9:
                frags.append(k)
        for k in frags:
            for x in all_rmsd_data[k]:
                pos_all_rmsd.append(x)
        position_avg_rmsd=np.mean(pos_all_rmsd)
        position_max_rmsd=np.max(pos_all_rmsd)
        position_min_rmsd=np.min(pos_all_rmsd)
        position_med_rmsd=np.median(pos_all_rmsd)
        of.write('\nFor Design Position: '+str(i))
        of.write('\nThe Average RMSD is: '+str(position_avg_rmsd))
        of.write('\nThe Median RMSD is: '+str(position_med_rmsd))
        of.write('\nThe Min RMSD is: '+str(position_min_rmsd))
        of.write('\nThe Max RMSD is: '+str(position_max_rmsd)+'\n')
        des_values_for_comparison.append((i,position_avg_rmsd))
    #now get the average rmsd for all fragments that contain a design position
    all_design_position_means=[]
    for key in all_rmsd_data.keys():
        for value in all_rmsd_data[key]:
            all_design_position_means.append(np.mean(value))
    of.write('\nThe Average RMSD for all fragments containing a designed position is: '+str(np.mean(all_design_position_means)))
    des_values_for_comparison.append(float(np.mean(all_design_position_means)))
    resultsdict['design_avg_rmsd_des_pos']=float(np.mean(all_design_position_means))
    #now comparing the design with the wildtype scaffold
    of.write('\n'+'\n'+'COMPARISON OF DESIGN '+str(design[:-4].upper())+' WITH WILDTYPE (wt_value - design_value)'+'\n')
    of.write('\n'+'n_good_frags/total_frags, coverage, (pos,pos,avg_rmsd_pos), avg_rmsd_all_pos')
    comparisonlist=[]
    ccc=1
    for i,e in enumerate(wt_values_for_comparison):
        if type(e)==float:
            diff=e-des_values_for_comparison[i]
            of.write('\n'+str(diff))
            comparisonlist.append(diff)
        elif type(e)==tuple:
            diff=e[1]-des_values_for_comparison[i][1]
            of.write('\n'+str(e[0])+','+str(des_values_for_comparison[i][0])+': '+str(diff))
            resultsdict['wt_diff_'+str(ccc)]=diff
            ccc+=1
    resultsdict['wt_diff_fragqual']=comparisonlist[0]
    resultsdict['wt_diff_coverage']=comparisonlist[1]
    resultsdict['wt_diff_rmsd']=comparisonlist[2]

of.close()

resultsdict['description']=sys.argv[5]#always this way in scorefiles
resultslist.append(resultsdict)
df=pd.DataFrame(resultslist)
dfname='fragment_filters.sc'
if os.path.exists(dfname):
    df.to_csv(dfname, header=None,index=False, sep='\t', mode='a')
else:
    dffnew=open(dfname,'w')
    dffnew.write('SEQUENCE:\n')
    dffnew.close()
    df.to_csv(dfname, index=False, sep='\t', mode='a')



json.dump(resultsdict,open(outputjsonname,'w'))
