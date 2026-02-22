#time ipython fragqual_filter.py des9mers des.pdb outputjsonname
#script assumes 200 fragments in frag libraries for each segment

#example command:

import sys
# import os
import numpy as np
import pandas as pd
from pyrosetta import *
init('-ignore_unrecognized_res')
# init( "-ignore_zero_occupancy false" )
#from  pyrosetta.rosetta.core.scoring import *


#load ninemers
fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
fragset9.read_fragment_file(sys.argv[1])
#load design pdb and params file
designs=[sys.argv[2]]
# paramsfile=[sys.argv[3]]
outputjsonname=sys.argv[3]

#open text file to write results
outfilename='frag_filters.txt'
of=open(outfilename,'a')

#to store stuff for pandas df
resultslist=[]
resultsdict={}
resultsdict['SCORE:']='SCORE:'

#evaluate designs
for design in designs:
    des_values_for_comparison=[]
    of.write('\n'+'\n'+'DESIGNS'+'\n')
    of.write('\n'+str(design[:-4].upper())+'\n')
    p=pose_from_pdb(design)
    # p=Pose()
    # generate_nonstandard_residue_set(p,paramsfile)
    # pose_from_file(p, design)
    #run fragqual filter on the whole ptn for n good frags + fraction of ptn covered by good frags
    calc9=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
    calc9.begin(1)
    calc9.end(p.total_residue())
    calc9.rmsd_cutoff(1.0)
    fragscore9=calc9.get('num_goodfrag',p)
    total_frags=fragset9.size() #total number of fragments
    of.write('The number of good fragments is: '+str(fragscore9)+' out of: '+str(total_frags)+' total fragments')
    des_values_for_comparison.append(float(fragscore9)/total_frags)
    resultsdict['design_fragqual_good/total']=float(fragscore9)/total_frags
    ###############coverage
    calc9r=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
    calc9r.begin(1)
    calc9r.end(p.total_residue())
    calc9r.ratio_cutoff(0.3)
    fragscore9r=calc9r.get('coverage',p)
    of.write('\nThe fraction of positions covered by good fragments is : '+str(fragscore9r)+'\n')
    des_values_for_comparison.append(float(fragscore9r))
    resultsdict['design_fragqual_coverage']=float(fragscore9r)
    #
    ############### average rmsd all fragments
    fr=pyrosetta.rosetta.core.fragment.FragmentRmsd(fragset9)
    #calculate rmsds
    all_rmsd_data={}
    minimum_rmsd_values=[]
    print('\n\n\n\n\ncalculating fragment rmsd values...\n\n\n\n\n')
    for fragment_number in range(1,p.total_residue()-8):
    # for fragment_number in range(1,2):
        try:
            design_fraglib_rmsd=[]#list of rmsd values for 200 frags for each frag number
            for x in range(1,201):
                curr_rmsd=fr.rmsd(fragment_number,x,p)
                design_fraglib_rmsd.append(curr_rmsd)
                # print(str(curr_rmsd))
            all_rmsd_data[fragment_number]=design_fraglib_rmsd
            minimum_rmsd_values.append(min(design_fraglib_rmsd))
        except:
            print(str(p.total_residue()))
            print(fragment_number)
    #now get the average rmsd for all fragments
    all_means=[]
    for key in all_rmsd_data.keys():
        for value in all_rmsd_data[key]:
            all_means.append(np.mean(value))
    of.write('\nThe Average RMSD for all fragments is: '+str(np.mean(all_means)))
    resultsdict['design_avg_rmsd']=float(np.mean(all_means))
    of.write('\n'+str(minimum_rmsd_values))
    resultsdict['max_min_rmsd']=float(max(minimum_rmsd_values))


of.close()

resultsdict['description']=sys.argv[2][:-4]
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

import json
json.dump(resultsdict,open(outputjsonname,'w'))
