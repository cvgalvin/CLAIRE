import sys
import os
#
import json
from rdkit import Chem

fragdefs=sys.argv[1]
ligpdb=sys.argv[2]


##################################################
pdb_ligands='db/PDB_small_molecule_SMILES_Inchi.json' #add path
##################################################



f=open(pdb_ligands,'r')
pdb_ligands_json_data=json.load(f)
f.close()
##################################################
f=open('db/pdb_small_molecules_ligandexpo.txt','r') #add path
##################################################



lines=[line.strip('\n') for line in f.readlines()]
f.close()
ligtopdb={}
for line in lines:
    ll=line.replace('\t',' ')
    l=ll.split(' ')
    l.remove('')
    ligtopdb[l[0]]=l[1:]

with open(fragdefs,'r') as f:
    lines=[line for line in f.readlines()]
with open(ligpdb,'r') as f:
    liglines=[line for line in f.readlines()]

starts=[]
ends=[]
for i,line in enumerate(lines):
    if line.strip('\n')=='<':
        starts.append(i)
    elif line.strip('\n')=='>':
        ends.append(i)
    else:
        pass
fragindices=[]
for i,e in enumerate(starts):
    fragindices.append((e+1,ends[i]))

fragmolfiles=[]
for i,e in enumerate(fragindices):
    fragatoms=[]
    for line in lines[e[0]:e[1]]:
        fragatoms.append(line.strip('\n').strip())
    newfragpdblines=[]
    for ligline in liglines:
        atomname=ligline[12:16]
        for atom in fragatoms:
            if atomname.strip()==atom:
                newfragpdblines.append(ligline)
    ofname='Fragment_'+str(i+1)+'.pdb'
    if len(newfragpdblines)==len(fragatoms):
        with open(ofname,'w') as of:
            for z in newfragpdblines:
                of.write(z)
        os.system('obabel -i pdb '+ofname+' -o mol -O '+ofname[:-3]+'mol')
    else:
        print('Problem converting mol: fragment '+str(i))
    fragmolfiles.append(ofname[:-3]+'mol')

print('all fragments identified and converted to pdb and mol files')
print(fragmolfiles)



#
smilesresults={}
# inchiresults={}
for fragmentmol in fragmolfiles:
    c=1
    #
    smilesligands=[]
    smilespdbs=[]
    #
    # inchiligands=[]
    # inchipdbs=[]
    #get smiles string of fragment
    m = Chem.MolFromMolFile(fragmentmol)
    frag_smiles=Chem.MolToSmiles(m)
    for ligand in list(pdb_ligands_json_data.keys()):
        checkligandsmiles=pdb_ligands_json_data[ligand]['SMILES_ccd']
        pdbm_ideal=Chem.MolFromSmiles(checkligandsmiles)
        # pdbm_with_bonds_withHs=Chem.AddHs(pdbm_with_bonds,addCoords=True)
        mappings=pdbm_ideal.HasSubstructMatch(m)
        # print(c)
        c+=1
        if mappings==True:
            if ligand not in smilesligands:
                smilesligands.append(ligand)
        else:
            pass
    td={}
    td['Ligands']=smilesligands
    listofpdbs=[]
    smilesd={}
    for tlig in smilesligands:
        for x in ligtopdb[tlig]:
            if x not in listofpdbs:
                listofpdbs.append(x)
        smilesd[tlig]=pdb_ligands_json_data[tlig]['SMILES_ccd']
    td['PDBs']=listofpdbs
    td['SMILES']=smilesd
    smilesresults[fragmentmol.split('.')[0]]=td
    print('\n\ndone')
    print(str(fragmentmol))




with open("PDB_search_results.json", "w") as outfile:
    json.dump(smilesresults, outfile)
