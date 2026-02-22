'''

use:
time python paramspath.params pdbpath.pdb outputdirectorypath


'''
import os
import sys

prm1=sys.argv[1]
pdb=sys.argv[2]
curr_soln_dirname=sys.argv[3]

#
#

from pyrosetta import *
from pyrosetta.toolbox import mutate_residue


init('-load_PDB_components False -ex1 -ex2')

#
sf=get_fa_scorefxn()

#
polar_res=['S','T','Y','K','W'] #res to mutate to looking for new hbonds
polar_residues=['SER','THR','ARG','LYS','ASN','GLN','ASP','GLU'] #res to find in bs and try mutate to np if not hbonding
npres=['L','I','V','F','M'] #npres to mutate to

#
#
#
noligunsats=[]
nouseless_buried_polar=[]
if not os.path.exists(curr_soln_dirname):
    os.makedirs(curr_soln_dirname,exist_ok=True)


#load pose
lig=[prm1]
p=Pose()
generate_nonstandard_residue_set(p,lig)
pose_from_file(p, pdb)
p.update_residue_neighbors()


#identify first shell residues
ligand_residue_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector('X')
#the residue selector is weird, so we use value 10 and cast a broad net to begin 
neighborhood_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10., False)
neighborhood_selector_bool = neighborhood_selector.apply(p)
neighborhood_residues_resnums = pyrosetta.rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
first_shell_res=list(neighborhood_residues_resnums)



# get the sasa of first shell res to make sure they're facing binding site
print(first_shell_res)
for residd in first_shell_res:
    individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    individual_res_selector.set_index(residd)
    individual_res_selector.apply(p)
    sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
    sasa_metric.set_residue_selector(individual_res_selector)
    this_res_sasa = sasa_metric.calculate(p)
    if this_res_sasa>=20.: #using a value of 20 but this can be played with 
        first_shell_res.remove(residd)
print(first_shell_res)
#need to double filter on sasa for reasons i dont understand 
for residd in first_shell_res:
    individual_res_selector=pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    individual_res_selector.set_index(residd)
    individual_res_selector.apply(p)
    sasa_metric = pyrosetta.rosetta.core.simple_metrics.metrics.SasaMetric()
    sasa_metric.set_residue_selector(individual_res_selector)
    this_res_sasa = sasa_metric.calculate(p)
    if this_res_sasa>=20.:
        first_shell_res.remove(residd)
print(first_shell_res)






#get atom indices of ligand polar atoms in network pose
lig_pol_atoms=[]
for i in range(1,len(p.residue(p.total_residue()).atoms())+1):
    s=str(p.residue(p.total_residue()).atom_type(i))
    ligatom=(s.split('\n'))[0].split('Atom Type:')[1].strip()   #####MAY NEED OTHER ROSETTA ATOM TYPES FOR OTHER CASES BUT THESE COVER ONES I USED 
    if ligatom=='OOC' or ligatom=='OH' or ligatom=='Nbb' or ligatom=='NLys':
        lig_pol_atoms.append(i)
        print(ligatom)






#find whether all ligand polar atoms hbond
hbond_set = rosetta.core.scoring.hbonds.HBondSet()
p.update_residue_neighbors()
rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbond_set)
hbond_data=[]
residues_hb_with_lig=[]
s=p.sequence()
if hbond_set.nhbonds()>0:
    for hbond_index in range(1,hbond_set.nhbonds()+1):
        drip=hbond_set.hbond(hbond_index).don_res_is_protein()
        arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
        if drip==True and arip==True:
            continue
        donres_ind=int(hbond_set.hbond(hbond_index).don_res())
        accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
        acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
        donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
        don_atom_index=int(p.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
        if drip==True and accres_ind==p.total_residue():
            hbond_data.append(acc_atom_index)
            residues_hb_with_lig.append(donres_ind)
        elif drip==False and donres_ind==p.total_residue():
            hbond_data.append(don_atom_index)
            residues_hb_with_lig.append(accres_ind)
unsat_lig_atoms=[]
for atind in lig_pol_atoms:
    if atind not in hbond_data:
        unsat_lig_atoms.append(atind)


#
useless_buried_polar=[]
for resnum in first_shell_res:
    restype=p.residue(resnum).name()[:3]
    if restype in polar_residues:
        if resnum not in residues_hb_with_lig:
            important_hb=0
            for resnum2 in residues_hb_with_lig:
                try:
                    w=p.energies().energy_graph().find_energy_edge(resnum,resnum2)
                    w.fill_energy_map()
                    hbsc=w[rosetta.core.scoring.hbond_sc]
                    if hbsc<=-0.5:
                        important_hb+=1
                except:
                    pass
            if important_hb>0:
                pass
            else:
                useless_buried_polar.append(resnum)



if len(unsat_lig_atoms)==0:
    os.system('cp '+pdb+' '+curr_soln_dirname+'/'+os.path.basename(pdb).split('.')[0]+'_unrefined.pdb')
    print('\n\n\nNO PROBLEMS\n\n\n')
    sys.exit()
else:
    pass



############################################################
############################################################
############################################################
############################################################
############################################################
ulp=0
ghb=0
gm=0
nmut=0

'''
tests=p.sequence()
'''
#mutate first shell residues not already hbonding with lig to polar, minimize,
#see if theres hbond now
solncount=0
if len(unsat_lig_atoms)==0:
    noligunsats.append(pdb)
else:
    for xx in range(len(unsat_lig_atoms)+1):
        good_hbond_mutations=[]
        for resid in first_shell_res:
            if resid not in residues_hb_with_lig:
                for polrestype in polar_res:
                    test_pose=p.clone()
                    start_energy=sf(test_pose)
                    mutate_residue(test_pose,resid,polrestype)
                    task_pack = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(test_pose)
                    task_pack.restrict_to_repacking()
                    task_pack.temporarily_fix_everything()
                    for fsr in first_shell_res:
                        if fsr not in residues_hb_with_lig:
                            task_pack.temporarily_set_pack_residue(fsr, True)
                    pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(sf, task_pack)
                    pack_mover.apply(test_pose)
                    new_energy=sf(test_pose)
                    w=test_pose.energies().energy_graph().find_energy_edge(resid,test_pose.total_residue())
                    w.fill_energy_map()
                    hbsc=w[rosetta.core.scoring.hbond_sc]
                    if hbsc<=-0.5: 
                        solncount+=1
                        good_hbond_mutations.append((new_energy,p.residue(resid).name(),resid,polrestype))
        if len(good_hbond_mutations)>0:
            print(good_hbond_mutations)
            ghb+=1
            good_hbond_mutations=sorted(good_hbond_mutations, key=lambda first: first[0])
            best_hbond_mut=good_hbond_mutations[0]
            best_hbond_mut_resid=best_hbond_mut[2]
            best_hbond_mut_restype=best_hbond_mut[3]
            test_pose=p.clone()
            mutate_residue(test_pose,best_hbond_mut_resid,best_hbond_mut_restype)
            task_pack = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(test_pose)
            task_pack.restrict_to_repacking()
            task_pack.temporarily_fix_everything()
            for fsr in first_shell_res:
                if fsr not in residues_hb_with_lig:
                    task_pack.temporarily_set_pack_residue(fsr, True)
            pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(sf, task_pack)
            pack_mover.apply(test_pose)     
            p=test_pose.clone()
            nmut+=1
            residues_hb_with_lig.append(best_hbond_mut_resid)
        else:
            continue




############################################################
############################################################
############################################################
############################################################
############################################################
#identify polar residues in binding site that are not hydrogen
#bonding with the ligand
#mutate them to np
useless_buried_polar=[]
for resnum in first_shell_res:
    restype=p.residue(resnum).name()[:3]
    if restype in polar_residues:
        if resnum not in residues_hb_with_lig:
            important_hb=0
            for resnum2 in residues_hb_with_lig:
                try:
                    w=p.energies().energy_graph().find_energy_edge(resnum,resnum2)
                    w.fill_energy_map()
                    hbsc=w[rosetta.core.scoring.hbond_sc]
                    if hbsc<=-0.5:
                        important_hb+=1
                except:
                    pass
            if important_hb>0:
                continue
            else:
                useless_buried_polar.append(resnum)
solncount2=0
if len(useless_buried_polar)==0:
    nouseless_buried_polar.append(pdb)
else:
    for resid in useless_buried_polar:
        good_mutations=[]
        for nprestype in npres:
            test_pose=p.clone()
            start_energy=sf(test_pose)
            mutate_residue(test_pose,resid,nprestype)
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            mm.set_bb(resid,True)
            mm.set_chi(resid,True) ## For side chain minimization
            minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            minmover.movemap(mm)
            minmover.score_function(sf)
            minmover.apply(test_pose)
            new_energy=sf(test_pose)
            print('\n\n\n')
            print(start_energy)
            print(new_energy)
            print('\n\n\n')
            if new_energy<=start_energy+2.:
                solncount2+=1
                good_mutations.append((new_energy,p.residue(resid).name(),resid,nprestype))
        if len(good_mutations)>0:
            gm+=1
            good_mutations=sorted(good_mutations, key=lambda first: first[0])
            best_mut=good_mutations[0]
            best_mut_resid=best_mut[2]
            best_mut_restype=best_mut[3]
            test_pose=p.clone()
            mutate_residue(test_pose,best_mut_resid,best_mut_restype)
            mm = MoveMap()
            mm.set_bb(False)
            mm.set_chi(False)
            mm.set_bb(resid,True)
            mm.set_chi(resid,True) ## For side chain minimization
            minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
            minmover.movemap(mm)
            minmover.score_function(sf)
            minmover.apply(test_pose)
            p=test_pose.clone()
            nmut+=1
        else:
            continue


#find whether all ligand polar atoms hbond
hbond_set = rosetta.core.scoring.hbonds.HBondSet()
p.update_residue_neighbors()
rosetta.core.scoring.hbonds.fill_hbond_set(p, False, hbond_set)
hbond_data=[]
residues_hb_with_lig=[]
s=p.sequence()
if hbond_set.nhbonds()>0:
    for hbond_index in range(1,hbond_set.nhbonds()+1):
        drip=hbond_set.hbond(hbond_index).don_res_is_protein()
        arip=hbond_set.hbond(hbond_index).acc_res_is_protein()
        if drip==True and arip==True:
            continue
        donres_ind=int(hbond_set.hbond(hbond_index).don_res())
        accres_ind=int(hbond_set.hbond(hbond_index).acc_res())
        acc_atom_index=int(hbond_set.hbond(hbond_index).acc_atm())
        donh_atom_index=int(hbond_set.hbond(hbond_index).don_hatm())
        don_atom_index=int(p.residue(donres_ind).first_adjacent_heavy_atom(donh_atom_index))
        if drip==True and accres_ind==p.total_residue():
            hbond_data.append(acc_atom_index)
            residues_hb_with_lig.append(donres_ind)
        elif drip==False and donres_ind==p.total_residue():
            hbond_data.append(don_atom_index)
            residues_hb_with_lig.append(accres_ind)
unsat_lig_atoms=[]
for atind in lig_pol_atoms:
    if atind not in hbond_data:
        unsat_lig_atoms.append(atind)


useless_buried_polar=[]
for resnum in first_shell_res:
    restype=p.residue(resnum).name()[:3]
    if restype in polar_residues:
        if resnum not in residues_hb_with_lig:
            important_hb=0
            for resnum2 in residues_hb_with_lig:
                try:
                    w=p.energies().energy_graph().find_energy_edge(resnum,resnum2)
                    w.fill_energy_map()
                    hbsc=w[rosetta.core.scoring.hbond_sc]
                    if hbsc<=-0.5:
                        important_hb+=1
                except:
                    pass
            if important_hb>0:
                continue
            else:
                useless_buried_polar.append(resnum)


new_name=os.path.join(curr_soln_dirname,os.path.basename(pdb).strip('.pdb')+'_'+str(nmut))+'mutations'


if len(unsat_lig_atoms)==0:
    p.dump_pdb(new_name+'_fixed.pdb')
    sys.exit()

if solncount2>=1 or solncount>=1:
    if len(unsat_lig_atoms)>0 and len(useless_buried_polar)==0:
        p.dump_pdb(new_name+'_unsats.pdb')
    elif len(unsat_lig_atoms)==0 and len(useless_buried_polar)>0:
        p.dump_pdb(new_name+'_ubp.pdb')
    elif len(unsat_lig_atoms)>0 and len(useless_buried_polar)>0:
        p.dump_pdb(new_name+'_bothproblems.pdb')
else:
    print('no solutions')

