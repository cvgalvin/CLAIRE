#ex: time ipython ~/match_resfile.py 1F4P_0001.pdb FMN.params 4.0 genpot/no
#assumes match pdb has motif res in name eg UM_1_W36W79Q43....pdb

import math
import sys
from pyrosetta import *
if sys.argv[4]=='genpot':
	init('-gen_potential -ignore_unrecognized_res -load_PDB_components False')
else:
	init('-ignore_unrecognized_res -load_PDB_components False')

#load pose
lig=[sys.argv[2]]
p=Pose()
generate_nonstandard_residue_set(p,lig)
matchname=sys.argv[1]
pose_from_file(p, matchname)

motifres=matchname.split('_')[2]
resind=[]
for e,i in enumerate(motifres):
	try:
		int(i)
	except:
		resind.append(e)
cstres=[]
for i,e in enumerate(resind[:-1]):
	resnum=motifres[e+1:resind[i+1]]
	cstres.append(int(resnum))
lastresnum=motifres[resind[-1]+1:]
cstres.append(int(lastresnum))


#identify the ligand residue number
lig_resnum=p.total_residue()
# for i in range(p.total_residue()):
# 	resnum=i+1
# 	if p.residue(resnum).is_ligand()==True:
# 		lig_resnum+=resnum

#functions to get distance between atoms
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

##########
cut1=float(sys.argv[3])
designable=[]
repackable=[]
native_sc=[]
designable_aro=[]
aro_res=['TYR','PHE','TRP']
##############
if lig_resnum==0: #check that ligand was identified
	print('COULD NOT IDENTIFY LIGAND')
	sys.exit()
else: #iterate over all ligand-protein atom atom distances and store res that meet cutoffs
	for v in p.residue(lig_resnum).atoms():
		v1=[float(str(v).split(',')[0]),float(str(v).split(',')[1]),float(str(v).split(',')[2])]
		for i in range(1,p.total_residue()):
			if i in cstres:
				repackable.append(i)
			else:
				if i!=lig_resnum:
					for vv in p.residue(i).atoms():
						v2=[float(str(vv).split(',')[0]),float(str(vv).split(',')[1]),float(str(vv).split(',')[2])]
						d=dist(v1,v2)
						if d <=cut1:
							if i not in cstres:
								designable.append(i)
								# resname=p.residue(i).name()[:3]
								# if resname in aro_res:
								# 	designable_aro.append(i)
								# else:
								# 	designable.append(i)
						# elif cut1<d<=cut2:
						# 	if i not in cstres:
						# 		repackable.append(i)
						# elif cut2<d:
						# 	native_sc.append(i)
# designable_aro=set(designable_aro)
designable=list(set(designable))


print(len(designable))
#
#get the sasa of first shell res to make sure they're facing binding site
for residd in designable:
    individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    individual_res_selector.set_index(residd)
    individual_res_selector.apply(p)
    sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
    sasa_metric.set_residue_selector(individual_res_selector)
    this_res_sasa = sasa_metric.calculate(p)
    if this_res_sasa>=12.:
        designable.remove(residd)
#this is one of those stupid instances wbere i have to filter a list twice
#for no apparent reason, pisses me off -.-
for residd in designable:
    individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    individual_res_selector.set_index(residd)
    individual_res_selector.apply(p)
    sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
    sasa_metric.set_residue_selector(individual_res_selector)
    this_res_sasa = sasa_metric.calculate(p)
    if this_res_sasa>=12.:
        designable.remove(residd)
print(len(designable))


#okay now find repackable
for resnum in designable:
	for v in p.residue(resnum).atoms():
		v1=[float(str(v).split(',')[0]),float(str(v).split(',')[1]),float(str(v).split(',')[2])]
		for i in range(1,p.total_residue()):
			if i!=lig_resnum:
				if i not in cstres:
					for vv in p.residue(i).atoms():
						v2=[float(str(vv).split(',')[0]),float(str(vv).split(',')[1]),float(str(vv).split(',')[2])]
						d=dist(v1,v2)
						if d<=cut1:
							repackable.append(i)
repackable=set(repackable)
for i in range(1,p.total_residue()):
	if i not in cstres:
		if i!=lig_resnum:
			if i not in designable:
				if i not in repackable:
					native_sc.append(i)
native_sc=set(native_sc)















#make sure no residues in multiple lists
for i in designable:
	if i in repackable:
		repackable.remove(i)
	if i in native_sc:
		native_sc.remove(i)
for i in repackable:
	if i in native_sc:
		native_sc.remove(i)



ofilename=matchname.split('.')[0]+'.resfile'
ofile=open(ofilename,'w')
#write header
ofile.write('USE_INPUT_SC')
ofile.write('\nstart'+'\n')

for i in range(1,p.total_residue()+1):
	if i in designable:
		s=str(p.pdb_info().pose2pdb(i))+'NOTAA CP'
		ofile.write(s+'\n')
	# elif i in designable_aro:
	# 	s=str(p.pdb_info().pose2pdb(i))+'PIKAA FWY'
	# 	ofile.write(s+'\n')
	elif i in repackable:
		s=str(p.pdb_info().pose2pdb(i))+'NATAA'
		ofile.write(s+'\n')
	elif i in native_sc:
		s=str(p.pdb_info().pose2pdb(i))+'NATRO'
		ofile.write(s+'\n')
	elif i==lig_resnum:
		s=str(p.pdb_info().pose2pdb(i))+'NATAA'
		ofile.write(s+'\n')

ofile.close()

print('designable:')
print(str(len(designable)))
print('repackable:')
print(str(len(repackable)))
