
import json
import sys
import os
from prody import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
import numpy as np
from io import StringIO
from Bio import pairwise2 as pw2

#suppress rdkit warnings
RDLogger.DisableLog('rdApp.*')

#james' function to calculate the transformation matrix to align reference
#onto target
def calculate_transformation_matrix(reference_conformer_atoms, target_conformer_atoms):
    """
    Caluclate transformation matrix by hand...
    When applied to target, transformation matrix transforms target onto reference.

    :param reference_conformer_atoms: target (in ProDy speak)
    :param target_conformer_atoms: mobile (in ProDy speak)
    :return:
    """
    # Calculate centroid for selected atoms
    reference_centroid = np.sum(reference_conformer_atoms, axis=0) / len(reference_conformer_atoms)
    target_centroid = np.sum(target_conformer_atoms, axis=0) / len(target_conformer_atoms)

    # Generate centered coordinate sets
    reference_centered = np.asarray(
        [atom - reference_centroid for atom in reference_conformer_atoms])
    target_centered = np.asarray([atom - target_centroid for atom in target_conformer_atoms])

    # Calculate Covariance Matrix
    covariance_matrix = np.dot(reference_centered.T, target_centered)

    # Singular Value Decomposition of covariance matrix
    V, s, Ut = np.linalg.svd(covariance_matrix)

    # Calculate determinant of Ut.T and V.T
    det = np.linalg.det(np.dot(Ut.T, V.T))
    det_matrix = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, det]])

    # Calculate Rotation Matrix
    rotation_matrix = np.dot(np.dot(Ut.T, det_matrix), V.T)

    # Translation vector
    translation_vector = np.dot(rotation_matrix, -(reference_centroid[np.newaxis].T)) + target_centroid[
        np.newaxis].T

    return prody.Transformation(rotation_matrix, np.ravel(translation_vector))

#let's use the og bsff search json as a starting point still,
#gives the pdb ligands which contain the fragment as a substructure
#and the pdbids of strc bound to these ligands
pdb_search_results_json=sys.argv[1]
f=open(pdb_search_results_json,'r')
pdb_search_results_json_data=json.load(f)
f.close()

#
os.makedirs('Transformed_Aligned_PDBs',exist_ok=True)
n_contact_pdbs_output={}
for current_fragment in list(pdb_search_results_json_data.keys()):
    output_count=0
    pdb_paths=[]
    relevant_ligands=list(pdb_search_results_json_data[current_fragment]['Ligands'])
    # os.makedirs('pdb_bank/'+current_fragment,exist_ok=True)
    os.makedirs('Transformed_Aligned_PDBs/'+current_fragment,exist_ok=True)
    print('\nretrieving local pdbs for '+current_fragment+'; from https://pdb-redo.eu')
    '''
    getting pdb redo strc from local copy of database
    '''
    for frag_containing_pdb in pdb_search_results_json_data[current_fragment]['PDBs']:
        pdb_redo_strc_path=r'/wynton/group/databases/PDB-REDO/'+frag_containing_pdb.lower()[1:3]+'/'+frag_containing_pdb.lower()+'/'+frag_containing_pdb.lower()+'_final.pdb'
        # print(pdb_redo_strc_path)
        if os.path.exists(pdb_redo_strc_path):
            pdb_paths.append(pdb_redo_strc_path)
            print(frag_containing_pdb)
        else:
            print('\ncould not find local copy of pdb: '+frag_containing_pdb+'\n')
            continue
    '''
    going through pdb strc,finding which relevant ligand/s they have,
    identifying fragment substructure in ligand
    align pdb to fragment (put pdb in frag coord frame)
    check for redundancy between aligned contacts
    write pdb files of nonredundant contacts
    '''
    fragment_mol='Inputs/Fragment_Inputs/'+current_fragment+'.mol' #mol file needs to be in working directory
    m=Chem.AllChem.MolFromMolFile(fragment_mol,removeHs=False)
    # m=Chem.AllChem.MolFromMolFile(fragment_mol)
    fragment_coords=[]
    for i in range(0, m.GetNumAtoms()):
        fpos = m.GetConformer().GetAtomPosition(i)
        fragment_coords.append((fpos.x,fpos.y,fpos.z))#m.GetAtomWithIdx(i).GetSymbol(),
    print('\nprocessing pdbs...')
    if m.GetNumAtoms()<=5:
        uselow=1
    else:
        uselow=0
    for pdb_to_process in pdb_paths:
        try:
            pdb_object = prody.parsePDB(pdb_to_process)
            pdb_lig=pdb_object.select('not water hetero')
            len(set(list(pdb_lig.getResnames())))
        except:
            print('\nproblem processing pdb '+pdb_to_process)
            continue
        for resname in set(list(pdb_lig.getResnames())):
            if resname in relevant_ligands:
                print('\n'+pdb_to_process+': '+resname)
                try:
                    pdb_ligand_smiles_canon=pdb_search_results_json_data[current_fragment]['SMILES'][resname]
                    pdbm_ideal=Chem.MolFromSmiles(pdb_ligand_smiles_canon)
                except:
                    print('cannot get smiles of ligand in pdb '+pdb_to_process)
                    continue
                rel_lig_pdb=pdb_lig.select('resname '+resname)
                hv = rel_lig_pdb.getHierView()
                instances=1
                good_instances=0
                aligned_pdbs=[]
                for chain in list(hv):#iterating different chains in pdb file containing target ligand
                    for residue in chain:#going through each instance of ligand
                        try:
                            output = StringIO()
                            writePDBStream(output, residue)
                            pdb_string = output.getvalue()
                            pdbm=AllChem.MolFromPDBBlock(pdb_string)
                            pdbm_with_bonds=AllChem.AssignBondOrdersFromTemplate(pdbm_ideal, pdbm)
                            if int(pdbm_with_bonds.GetNumAtoms())<6:
                                continue
                            else:
                                pass
                        except:
                            print('could not assign bond orders to '+resname)
                            continue
                        try:
                            pdbm_with_bonds_withHs=Chem.AddHs(pdbm_with_bonds,addCoords=True)
                        except:
                            print('could not assign hydrogens to '+resname)
                            continue
                        mappings=pdbm_with_bonds_withHs.GetSubstructMatches(m)
                        # mappings=pdbm_with_bonds.GetSubstructMatches(m)
                        if len(mappings)<1:
                            print('no mappings;fragment atoms could not be properly mapped for instance '+str(instances))
                            instances+=1
                            continue
                        else:
                            mapping_count=0
                            for mapping in mappings: #going through each mapping for given instance of ligand
                                if len(mapping)!=len(fragment_coords):
                                    print('could not map all fragment atoms; instance '+str(instances)+'; mapping '+str(mappings.index(mapping)))
                                    continue
                                else:
                                    mapping_count+=1
                                    # print(mapping)
                                    # print(fragment_coords)
                                    substructure_atom_coordinates=[]
                                    for atom_idx in mapping:
                                        pos=pdbm_with_bonds_withHs.GetConformer().GetAtomPosition(atom_idx)
                                        # pos=pdbm_with_bonds.GetConformer().GetAtomPosition(atom_idx)
                                        substructure_atom_coordinates.append((pos.x,pos.y,pos.z))
                                    transform_matrix=calculate_transformation_matrix(np.asarray(substructure_atom_coordinates),np.asarray(fragment_coords))
                                    # print(fragment_coords)
                                    transformed_substrc=prody.applyTransformation(transform_matrix,np.asarray(substructure_atom_coordinates))
                                    rmsd = prody.calcRMSD(np.asarray(fragment_coords), transformed_substrc)
                                    print('rmsd instance '+str(instances)+' mapping '+str(mapping_count)+': '+str(rmsd))
                                    if uselow==0:
                                        if rmsd<=1.0:
                                            good_instances+=1
                                            l=[atom for atom in residue]
                                            l2=[atom for idx,atom in enumerate(l) if idx in mapping]
                                            substructure_atom_serials=[str(atom.getSerial()) for atom in l2]
                                            temp_pdb_object = prody.parsePDB(pdb_to_process)
                                            target_shell = temp_pdb_object.select('serial {0} or (protein and within 12 of (serial {0}) and not resnum {1})'.format(' '.join(substructure_atom_serials), residue.getResnum()))
                                            transformed_pdb=prody.applyTransformation(transform_matrix, target_shell)
                                            aligned_pdb_string = StringIO()
                                            writePDBStream(aligned_pdb_string, transformed_pdb)
                                            transformed_pdb_name='Transformed_Aligned_PDBs/'+current_fragment+'/'+pdb_to_process[-14:-10]+'_'+resname+'_'+str(instances)+'_'+str(mapping_count)+'.pdb'
                                            aligned_pdbs.append((transformed_pdb_name,aligned_pdb_string))
                                        else:
                                            print('rejected: RMSD too high ')
                                    elif uselow==1:
                                        if rmsd<=0.5:
                                            good_instances+=1
                                            l=[atom for atom in residue]
                                            l2=[atom for idx,atom in enumerate(l) if idx in mapping]
                                            substructure_atom_serials=[str(atom.getSerial()) for atom in l2]
                                            temp_pdb_object = prody.parsePDB(pdb_to_process)
                                            target_shell = temp_pdb_object.select('serial {0} or (protein and within 12 of (serial {0}) and not resnum {1})'.format(' '.join(substructure_atom_serials), residue.getResnum()))
                                            transformed_pdb=prody.applyTransformation(transform_matrix, target_shell)
                                            aligned_pdb_string = StringIO()
                                            writePDBStream(aligned_pdb_string, transformed_pdb)
                                            transformed_pdb_name='Transformed_Aligned_PDBs/'+current_fragment+'/'+pdb_to_process[-14:-10]+'_'+resname+'_'+str(instances)+'_'+str(mapping_count)+'.pdb'
                                            aligned_pdbs.append((transformed_pdb_name,aligned_pdb_string))
                                        else:
                                            print('rejected: RMSD too high ')
                            instances+=1
                print(str(instances-1)+' instances of '+resname+' in structure, '+str(good_instances)+' low RMSD mappings')
                '''
                check structures for redundancy using first atom count comparison
                if <100 atom difference between two structures, do sequence comparison
                >75% seq ID considered redundant and one is deleted
                '''
                if len(aligned_pdbs)>0:
                    print('checking for redundancy...')
                    for i,aligned_pdb_1 in enumerate(aligned_pdbs[:-1]):
                        try:
                            aln_pdb_1=prody.parsePDBStream(StringIO(aligned_pdb_1[1].getvalue()))
                        except:
                            aligned_pdbs.remove(aligned_pdb_1)
                            print('\nproblem loading structure '+aligned_pdb_1[0])
                            continue
                        aln1_atoms=aln_pdb_1.numAtoms()
                        if i<len(aligned_pdbs)-1:
                            for z,aligned_pdb_2 in enumerate(aligned_pdbs[i+1:]):
                                try:
                                    aln_pdb_2=prody.parsePDBStream(StringIO(aligned_pdb_2[1].getvalue()))
                                except:
                                    print('\nproblem loading structure '+aligned_pdb_2[0])
                                    aligned_pdbs.remove(aligned_pdb_2)
                                    continue
                                aln2_atoms=aln_pdb_2.numAtoms()
                                if abs(aln2_atoms-aln1_atoms)>100:
                                    continue
                                else:
                                    aln1_hv = aln_pdb_1.getHierView()
                                    aln2_hv = aln_pdb_2.getHierView()
                                    aln1_sequence = [chain.getSequence() for chain in aln1_hv]
                                    aln2_sequence = [chain.getSequence() for chain in aln2_hv]
                                    # for chain in aln1_hv:
                                    #     aln1_sequence = [chain.getSequence() for chain in aln1_hv]
                                    # for chain in aln2_hv:
                                    #     aln2_sequence = [chain.getSequence() for chain in aln2_hv]
                                    try:
                                        global_align = pw2.align.globalms(''.join(aln1_sequence), ''.join(aln2_sequence),1, 0, -1.0, -0.1)
                                        matches=global_align[0]
                                        print(matches)
                                        seq_length = min(len(''.join(aln1_sequence)), len(''.join(aln2_sequence)))
                                        # print(seq_length)
                                        nmatches = float(matches[2])
                                        # print(matches)
                                        percent_match = (nmatches/seq_length) * 100
                                        # print(aligned_pdb_1[0])
                                        # print(aligned_pdb_2[0])
                                        print('percent match: '+str(percent_match))
                                        if percent_match>70.0:
                                            print('eliminating redundant structure '+aligned_pdb_2[0])
                                            aligned_pdbs.remove(aligned_pdb_2)
                                        else:
                                            continue
                                    except:
                                        print('problem aligning sequences:'+aligned_pdb_1[0]+':'+aligned_pdb_2[0])
                                        continue
                    '''
                    write nonredundant structures to pdb file
                    '''
                    for aligned_pdb in aligned_pdbs:
                        writePDB(aligned_pdb[0],prody.parsePDBStream(StringIO(aligned_pdb[1].getvalue())))
                        output_count+=1
                    print(str(len(aligned_pdbs))+' unique structures written to PDB')
                else:
                    continue
    n_contact_pdbs_output[current_fragment]=output_count

print('\nContact PDBs output:')
print(n_contact_pdbs_output)
