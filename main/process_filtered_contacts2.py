import os
import sys
from collections import defaultdict
#
fogpath=sys.argv[1]
ligand_path=sys.argv[2]
target_name=sys.argv[3]
frag_name=sys.argv[4]
resultdir=('_').join([target_name,frag_name,'residue_contacts'])
os.makedirs(str(resultdir),exist_ok=True)
#
fog=open(fogpath,'r')
lines=[line for line in fog.readlines()]
fog.close()
all_fuzzball_lines=lines
#first gotta get indices of unique residues in all_fuzzball_lines
fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
starts=[];lasts=[]
for index, line in enumerate(all_fuzzball_lines): #collecting start/end indices for unique res
        resnum=int(all_fuzzball_lines[index][22:26])
        resname=line[17:20]
        try:
            lastresnum=int(all_fuzzball_lines[index-1][22:26])
            lastresname=all_fuzzball_lines[index-1][17:20]
            if resnum!=lastresnum or resname!=lastresname:
                start=index
                starts.append(start)
        except:
            start=index
            starts.append(start)
        try:
            nextresname=all_fuzzball_lines[index+1][17:20]
            next_resnum=int(all_fuzzball_lines[index+1][22:26])
            if resnum!=next_resnum or resname!=nextresname:
                last=index+1
                lasts.append(last)
        except:
            last=len(all_fuzzball_lines)
            lasts.append(last)
for index,start in enumerate(starts): #put the indices together for each res
    fuzzball_residue_indices.append((start,lasts[index]))
fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices
fuzzball_residue_indices=sorted(fuzzball_residue_indices, key=lambda first: first[0])
#
aas=['LYS','ARG','ASP','GLU','SER','THR','GLN','ASN','HIS','TYR','TRP','GLY','TYR','TRP','ALA','LEU','ILE','VAL','PHE','MET'] #ONLY TRUE hbonders AND GLY
#separate into lists for each res type
reslines=defaultdict(list)
for (a,b) in fuzzball_residue_indices:
    resname=all_fuzzball_lines[a][17:20]
    reslines[resname].append((a,b))
#load the full ligand lines to add to fuzzball outputs
#wanna put it as last line in fuzzball files
#
f=open(ligand_path,'r')
liglines=[line for line in f.readlines()]
f.close()
#now lets build the fuzzballs for each res
new_residue_number=1
clean_fuzzball_lines=[]
#
for aa in aas:
    new_residue_number=1
    clean_fuzzball_lines=[]
    for (a,b) in reslines[aa]: #go through each residue and edit line w new resnum
        current_residue_lines=[i for i in all_fuzzball_lines[a:b]]
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
            else:
                print(str(new_residue_number))
        new_residue_number+=1
    #now add the ligand lines and output the polar cluster
    if len(clean_fuzzball_lines)>0:
        for line in liglines:
            clean_fuzzball_lines.append(line)
        offbname=('_').join([target_name,frag_name,aa,'contact_fuzzball.pdb'])
        of=open(os.path.join(resultdir,offbname),'w')
        for line in clean_fuzzball_lines:
            of.write(line)
        of.close()
    else:
        continue
