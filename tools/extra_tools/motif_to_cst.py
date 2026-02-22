#time ipython motif_to_cst.py motif.pdb ligand.params cstfile.cst
#ensure ligand is last residue in pdb file
import sys
import math
import numpy as np
import os
from pyrosetta import *
init('-load_PDB_components False')

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

def generate_single_constraint_block_base(distanceAB,angleA,angleB,torsionA,torsionAB,torsionB,
                                          residue_resname,ligand_resname,ligand_atoms,residue_atoms,
                                          distance_tolerance=0.5, angle_A_tolerance=5, angle_B_tolerance=5,
                                          torsion_A_tolerance=5, torsion_AB_tolerance=5, torsion_B_tolerance=5,
                                          torsion_constraint_sample_number=1, angle_constraint_sample_number=1,
                                          distance_constraint_sample_number=0, greasy_sampling=True):
    """
    Generate a single constraint block for one residue-ligand interaction
    :return:
    """

    residue_tag = 'residue3'
    # Increase tolerance/sampling by 5/1 for greasy residues (ACFILMVWY)
    if greasy_sampling and residue_resname in ['ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR']:
        angle_A_tolerance += 10
        angle_B_tolerance += 10
        torsion_A_tolerance += 10
        torsion_AB_tolerance += 10
        torsion_B_tolerance += 10
        torsion_constraint_sample_number += 1
        angle_constraint_sample_number += 1

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


#load the motif
frag = [sys.argv[2]] #ligand params
p = Pose()
generate_nonstandard_residue_set(p,frag)
pose_from_file(p, sys.argv[1])
ligand_resnum=p.total_residue()


'''
load the binding site in terminal (comment out previous block to use this way)
'''
# frag = ['DEX.params'] #ligand params
# p = Pose()
# generate_nonstandard_residue_set(p,frag)
# pose_from_file(p, '3mnp_motif.pdb')
# ligand_resnum=p.total_residue()

#open cst file to write
of=open(sys.argv[3],'w')
# of=open('csttest.cst','w')
#go through bs res and write cst blocks
for i in range(p.total_residue()-1):
    resnum=i+1
    #make new pose for ligand and each res
    contact_pose=Pose()
    ligand_pose=p.residue(p.total_residue()).clone()
    res_pose=p.residue(resnum).clone()
    contact_pose.append_residue_by_jump(res_pose, 1)
    contact_pose.append_residue_by_jump(ligand_pose, 1)
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
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    if len(resatomsequence)>0:
        resbase1=resatomsequence[0][1]
        resbase2=resatomsequence[0][2]
        residue_atom2=(contact_pose.residue(1).atom_name(resbase1)).strip()
        residue_atom3=(contact_pose.residue(1).atom_name(resbase2)).strip()
    else:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
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
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nREDUNDANT RES ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
    if len(set(lig_atoms))<3:
        of.close()
        os.system('rm '+sys.argv[3])
        print('\n\n\n\n\nREDUNDANT LIG ATOMS, ABORTED CST GENERATION FOR '+str(sys.argv[1])+'\n\n\n\n')
        sys.exit()
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
                                                   distance_tolerance=0.25, angle_A_tolerance=10, angle_B_tolerance=10,
                                                   torsion_A_tolerance=15, torsion_AB_tolerance=15, torsion_B_tolerance=15,
                                                   torsion_constraint_sample_number=2, angle_constraint_sample_number=2,
                                                   distance_constraint_sample_number=2, greasy_sampling=False)
    of.write('CST::BEGIN\n')
    of.write('\n'.join(cstblock))
    of.write('\n')
    of.write('CST::END\n')
    of.write('\n')
of.close()


