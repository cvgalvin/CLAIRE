
#
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from collections import defaultdict
import pandas as pd

def pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:4]=='ATOM':
            d={}
            d['recordname']=line[0:6]
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=line[20:22]
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']='\n'
            l.append(d)
    df=pd.DataFrame(l)
    return df

#convert dataframe to pdb file
def df_to_pdb(dataframe,ofile):
    of=open(ofile,'w')
    for i in range(dataframe.shape[0]):
        line=''.join(dataframe.iloc[i])
        of.write(line)
    of.close()

#convert ligand pdb to dataframe
def lig_pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:6]=='HETATM':
            d={}
            d['recordname']='ATOM'+'  '
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=' X'
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']='\n'
            l.append(d)
        else:
            break
    df=pd.DataFrame(l)
    return df

#return the rmsd between two sets of cartesian coordinates
def calc_rmsd(v1,v2):
    diff = np.array(v1) - np.array(v2)
    N = len(v1)
    return np.sqrt((diff * diff).sum() / N)

#
rfb=sys.argv[1]
ligand_path=sys.argv[2]
#


#
pdfname=rfb.split('.')[0]+'_clustered.pdf'
pdf = PdfPages(pdfname)

df=pdb_to_df(rfb)
print('df loaded')
#get coordinates of residues
residue_coords_dict=defaultdict(list)
for i in range(df.shape[0]):
        print(str(i))
        resnum=df.iloc[i][6]
        atom_name=str(df.iloc[i][2].split()[0])
        ptn_atom_coords=( float(df.iloc[i][8]), float(df.iloc[i][9]), float(df.iloc[i][10]) )
        if atom_name[0]!='H' and atom_name!='OXT' and atom_name!='C' and atom_name!='CA'and atom_name!='O' and atom_name!='N':
            contact_atom_info=(atom_name,float(ptn_atom_coords[0]),float(ptn_atom_coords[1]),float(ptn_atom_coords[2]))
            if contact_atom_info not in residue_coords_dict[resnum]:
                residue_coords_dict[resnum].append(contact_atom_info)
#create clusters based off of sidechain rmsd for sc contacts
print('\nstarting to cluster contacts...')
print(str(rfb))
cluster_dict=defaultdict(list)
already_clustered=[]
for residue_ in residue_coords_dict.keys():
    print(residue_)
    if residue_ not in already_clustered:
        already_clustered.append(residue_)
        current_cluster=[]
        for residue_2 in residue_coords_dict.keys():
            if residue_!=residue_2 and residue_2 not in already_clustered:
                vs1=[];vs2=[]
                for a,b,c,d in residue_coords_dict[residue_]:
                    for e,f,g,h in residue_coords_dict[residue_2]:
                        if a==e:
                            vs1.append((b,c,d))
                            vs2.append((f,g,h))
                        else:
                            pass
                x=len(vs1)
                y=len(vs2)
                if x!=y:
                    print('different number of atoms between residues '+str(residue_)+' and '+str(residue_2))
                    print('cannot calculate rmsd :(')
                else:
                    try:
                        rmsd=calc_rmsd(vs1,vs2)
                        if rmsd<=0.5:
                            already_clustered.append(residue_2)
                            cluster_dict[residue_].append(residue_2)
                        else:
                            pass
                    except:
                        pass
#okay so next thing is to create pdbs of clusters
##############################
#compile and organize cluster populations
cluster_populations=[]
for key in cluster_dict.keys():
    n_members=len(cluster_dict[key])+1
    cluster_populations.append((key,int(n_members)))
cluster_populations=sorted(cluster_populations, reverse=True, key=lambda nmem: nmem[1])
#plot cluster populations
ncm=[y for x,y in cluster_populations]
clabs=[i+1 for i in range(len(cluster_populations))]
plt.bar(clabs, ncm)
plt.xlabel('Cluster ID')
plt.ylabel('Number of Members')
plt.title('Cluster Population Distribution')
pdf.savefig()
plt.clf()
#make a nice home for stats and cluster pdbs
cluster_results_path=rfb.split('.')[0]+'_clustered'
if len(cluster_populations)>0:
    os.makedirs(cluster_results_path,exist_ok=True)
    n_contacts_check=[]
    fullligdf=lig_pdb_to_df(ligand_path)
    #export cluster fuzzballs with full ligand
    rns=[]
    for i in df['resnum']:
        rns.append(int(i))
    for a,b in cluster_populations:
        cluster_df=pd.DataFrame()
        if b>=0:
            n_contacts_check.append(b)
            resnumbers=[int(i) for i in cluster_dict[a]]
            resnumbers.append(int(a))
            for resnum in resnumbers:
                rows=[i for i,e in enumerate(rns) if e == resnum]
                res_df=df.iloc[min(rows):max(rows)+1]
                cluster_df=cluster_df.append(res_df,ignore_index=True)
            cluster_df=cluster_df.append(fullligdf,ignore_index=True)
            resnamee=str(res_df.iloc[0][4])
            ofilename=str(a.strip())+'_'+resnamee+'_cluster_'+str(b)+'.pdb'
            df_to_pdb(cluster_df,os.path.join(cluster_results_path,ofilename))
        else:
            continue
    n_contacts_check_sum=sum(n_contacts_check)
    print('there are '+str(n_contacts_check_sum)+' contacts after clustering')
else:
    print('there are '+str(len(cluster_populations))+' clusters')



pdf.close()
