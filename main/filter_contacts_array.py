#intended to be run on the pdb files generated in the fragment_x directories
#within Transformed_Aligned_pdbs dir from the 'align'command in BSFF protocol
#returns a single pdb file containing all residues with at least one atom
#within X angstroms of the ligand fragment, where x=4 for nonpolar and 3.5 A for polar res
#excludes cys, pro, ala
#excludes high b factor or zero occupancy

#USAGE: time ipython path/to/filter_contacts_array.py directory

import sys
import math
import os
import pandas as pd
import collections
from scipy.spatial import distance


#filter out bad residues
def filters(fuzzball_lines,fuzzball_residue_indices):
    startnres=len(fuzzball_residue_indices)
    for a,b in fuzzball_residue_indices: #filter res without at least 4 atoms
        l1=b-a
        if l1<4:
            fuzzball_residue_indices.remove((a,b))
            continue
    #now extend to requirements for individual residues to have complete sidechains
        resname=fuzzball_lines[a][17:20]
        if resname=='CYS' or resname=='ALA' or resname=='PRO': #remove all cys,ala,pro
            fuzzball_residue_indices.remove((a,b))
            continue
        rlines=fuzzball_lines[a:b]
        res_dists=[]
        for line in rlines:
            atom_coords=(float(line[31:39].strip()),float(line[39:47].strip()),float(line[47:54].strip()))
            for frag_atom in fragment_coords:
                d=distance.euclidean(atom_coords,frag_atom)
                res_dists.append(d)
            if line[11:16].split()[0][0]=='H' or line[11:16].split()[0]=='OXT':
                rlines.remove(line)
        if resname=='ASN' or resname=='ARG' or resname=='TYR' or resname=='LYS' or resname=='HIS' or resname=='GLU' or resname=='ASP' or resname=='GLN' or resname=='THR' or resname=='SER':
            if min(res_dists)<=3.5:
                pass
            else:
                fuzzball_residue_indices.remove((a,b))
                continue
        else:
            if min(res_dists)<=4.0:
                pass
            else:
                fuzzball_residue_indices.remove((a,b))
                continue
        l=len(rlines) #only looking at non hydrogens
        if resname=='SER':
            if l<6:
                fuzzball_residue_indices.remove((a,b))
                continue
        elif resname=='VAL' or resname=='THR':
            if l<7:
                fuzzball_residue_indices.remove((a,b))
                continue
        elif resname=='ILE' or resname=='LEU' or resname=='MET' or resname=='ASN' or resname=='ASP':
            if l<8:
                fuzzball_residue_indices.remove((a,b))
                continue
        elif resname=='GLN' or resname=='LYS' or resname=='GLU':
            if l<9:
                fuzzball_residue_indices.remove((a,b))
                continue
        elif resname=='HIS':
            if l<10:
                fuzzball_residue_indices.remove((a,b))
                continue
        elif resname=='PHE' or resname=='ARG':
            if l<11:
                fuzzball_residue_indices.remove((a,b))
                continue
        elif resname=='TYR':
            if l<12:
                fuzzball_residue_indices.remove((a,b))
                continue
        elif resname=='TRP':
            if l<14:
                fuzzball_residue_indices.remove((a,b))
                continue
        for line in rlines:
            occupancy=line[54:60]
            temp_factor=line[60:66]
            if float(occupancy)<0.5 or float(temp_factor)>40.0:
                fuzzball_residue_indices.remove((a,b))
                break
    end_nres=len(fuzzball_residue_indices)
    return startnres, end_nres


parent_dir=sys.argv[1]

pdbfiles=[file for file in os.listdir(parent_dir) if file[-3:]=='pdb']
#storing ligand coordinates from one of the pdb files in list fragment_coords
#storing full ligand lines in list liglines
#storing contact residue lines in list fuzzball_lines
liglines=[]
fragment_coords=[]
f=open(os.path.join(parent_dir,pdbfiles[0]),'r') #get lig lines
for line in f.readlines():
    if line[0:6]=='HETATM':
        liglines.append(line)
        x=float(line[31:38]) ;y=float(line[38:46]) ;z=float(line[46:54])
        fragment_coords.append((x,y,z))
f.close()
fuzzball_lines=[]
for file in pdbfiles: #
    f=open(os.path.join(parent_dir,file),'r')
    atomlines=[line for line in f.readlines() if line[0:4]=='ATOM'];f.close()
    for i in atomlines:
        fuzzball_lines.append(i)
#first gotta get indices of unique residues in fuzzball_lines
fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
starts=[];lasts=[]
for index, line in enumerate(fuzzball_lines): #collecting start/end indices for unique res
        resnum=int(fuzzball_lines[index][22:26])
        resname=line[17:20]
        try:
            lastresnum=int(fuzzball_lines[index-1][22:26])
            lastresname=fuzzball_lines[index-1][17:20]
            if resnum!=lastresnum or resname!=lastresname:
                start=index
                starts.append(start)
        except:
            start=index
            starts.append(start)
        try:
            nextresname=fuzzball_lines[index+1][17:20]
            next_resnum=int(fuzzball_lines[index+1][22:26])
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
og_n_res=str(len(fuzzball_residue_indices))
#polar residues (DEHKNQRSTY) >3.5Å
#for some reason the filters remove more and more residues
#upon additional application, so i am applying them over and over until the
#number of remaining res converges
print('starting to filter contacts...')
for i in range(100):
    s,n = filters(fuzzball_lines,fuzzball_residue_indices)
    print(s,n)
    if s==n:
        break
#reorder the residue numbers of the filtered contact residues because
#it is useful to have a unique identifier for residues
#(before doing this they have their original residue numbers from source pdb files)
clean_fuzzball_lines=[] #the renumbered lines
for line in liglines:
    newline=line[0:23]+'1'+'     '+line[29:]
    clean_fuzzball_lines.append(newline)
new_residue_number=2
for (a,b) in fuzzball_residue_indices: #go through each residue and edit line w new resnum
    current_residue_lines=[i for i in fuzzball_lines[a:b]]
    for line in current_residue_lines:
        if new_residue_number<10:
            newline=line[0:22]+' '+str(new_residue_number)+'     '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9<new_residue_number<100:
            newline=line[0:22]+' '+str(new_residue_number)+'    '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99<new_residue_number<1000:
            newline=line[0:22]+' '+str(new_residue_number)+'   '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 999<new_residue_number<10000:
            newline=line[0:22]+' '+str(new_residue_number)+'  '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9999<new_residue_number<100000:
            newline=line[0:22]+' '+str(new_residue_number)+' '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99999<new_residue_number:
            newline=line[0:22]+' '+str(new_residue_number)+line[29:]
            clean_fuzzball_lines.append(newline)
        else:
            print(str(new_residue_number))
    new_residue_number+=1
#print this stuff cus why not
quickstats=open(os.path.join(parent_dir,'contactsfound.txt'),'a')
n_files=str(len(os.listdir(os.path.join(parent_dir))))
quickstats.write(n_files)
quickstats.write('\n'+str(len(fuzzball_residue_indices)))
quickstats.write('\n'+og_n_res)
quickstats.close()
#write the output pdb fuzzball
ofile=open(os.path.join(parent_dir,'contact_fuzzball.pdb'),'w')
for line in clean_fuzzball_lines:
    ofile.write(line)
ofile.close()
