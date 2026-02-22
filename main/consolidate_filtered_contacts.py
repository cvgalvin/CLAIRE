
import os

parent_dir=os.getcwd()
tapdbs=os.path.join(parent_dir,'Transformed_Aligned_PDBs')
fragment_dirs=[i for i in os.listdir(tapdbs) if os.path.isdir(os.path.join(tapdbs,i))==True]

contacts_d={}
for fragment_dir in fragment_dirs:
    print(fragment_dir)
    all_fuzzball_lines=[]
    all_contacts_found=[]
    try:
        f=open(os.path.join(tapdbs,fragment_dir,'contact_fuzzball.pdb'),'r')
        lines=[line for line in f.readlines()]
        f.close()
        for line in lines:
            all_fuzzball_lines.append(line)
        f2=open(os.path.join(tapdbs,fragment_dir,'contactsfound.txt'),'r')
        cflines=[line for line in f2.readlines()]
        all_contacts_found.append((cflines[0],cflines[1],cflines[2]))
        f2.close()
    except:
        dirs=[d for d in os.listdir(os.path.join(tapdbs,fragment_dir,'filter_input')) if os.path.isdir(os.path.join(tapdbs,fragment_dir,'filter_input',d))==True]
        for dir in dirs:
            try:
                f=open(os.path.join(tapdbs,fragment_dir,'filter_input',dir,'contact_fuzzball.pdb'),'r')
                lines=[line for line in f.readlines()]
                f.close()
                for line in lines:
                    all_fuzzball_lines.append(line)
                f2=open(os.path.join(tapdbs,fragment_dir,'filter_input',dir,'contactsfound.txt'),'r')
                cflines=[line for line in f2.readlines()]
                all_contacts_found.append((cflines[0],cflines[1],cflines[2]))
                f2.close()
            except:
                print('problem '+str(fragment_dir)+'; '+str(dir))
    #get ligand lines and erase repeat instances of ligand
    liglines=[]
    for line in all_fuzzball_lines:
        if line[:6]=='HETATM':
            liglines.append(line)
        else:
            break
    all_fuzzball_lines=[line for line in all_fuzzball_lines if line[:4]=='ATOM']
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
    clean_fuzzball_lines=[] #the renumbered lines
    for line in liglines:
        newline=line[0:23]+'1'+'     '+line[29:]
        clean_fuzzball_lines.append(newline)
    new_residue_number=2
    for (a,b) in fuzzball_residue_indices: #go through each residue and edit line w new resnum
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
            elif 99999<new_residue_number:
                newline=line[0:22]+' '+str(new_residue_number)+line[29:]
                clean_fuzzball_lines.append(newline)
            else:
                print(str(new_residue_number))
        new_residue_number+=1
    of=open(os.path.join(tapdbs,fragment_dir,fragment_dir+'_contact_fuzzball.pdb'),'w')
    for line in clean_fuzzball_lines:
        of.write(line)
    of.close()
    contacts_d[fragment_dir]=all_contacts_found

of2=open(os.path.join(os.getcwd(),os.getcwd()+'_contactsfound.txt'),'w')
for key in contacts_d.keys():
    allcontacts_current_frag=contacts_d[key]
    nfiles=0
    nfilt=0
    nog=0
    for entry in allcontacts_current_frag:
        nfiles+=int(entry[0])
        nfilt+=int(entry[1])
        nog+=int(entry[2])
    of2.write('\n'+key)
    of2.write('\nn_files: '+str(nfiles))
    of2.write('\nnfilt: '+str(nfilt))
    of2.write('\nnog: '+str(nog))
of2.close()


for fragment_dir in fragment_dirs:
    try:
        os.system('rm -r '+os.path.join(tapdbs,fragment_dir,'filter_input'))
    except:
        pass
