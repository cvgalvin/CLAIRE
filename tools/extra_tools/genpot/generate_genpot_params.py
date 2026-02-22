import os
import sys
#first I need to get gp params files and replace lig coords in pdb with mol2params output pdb
ligand_charge=sys.argv[3]
ligname=sys.argv[2]
cleandirname='genpot'
os.makedirs(cleandirname,exist_ok=True)
#####################
initial_match=sys.argv[1]
fname=initial_match.split('.')[0]
fname2=os.path.join(cleandirname,initial_match.split('.')[0]+'_lig.pdb')
f=open(initial_match,'r')
lines=[]
alines=[]
for line in f.readlines():
    if line[0:6]=='HETATM':
        lines.append(line)
    elif line[0:4]=='ATOM':
        alines.append(line)
    elif line[:6]=='REMARK':
        alines.append(line)
f.close()
newligpdb=open(fname2,'w')
for line in lines:
    newligpdb.write(line)
newligpdb.close()
os.chdir(cleandirname)
os.system('obabel -i pdb '+fname+'_lig.pdb'+' -o mol2 -O '+fname+'_lig.mol2')
os.system('obabel -i mol2 '+fname+'_lig.mol2'+' -o mol2 -O '+fname+'_ligH.mol2 -p 7.4')
os.system('obabel -i mol2 '+fname+'_ligH.mol2'+' -o pdb -O '+fname+'_ligH.pdb')
os.system('time antechamber -i '+fname+'_ligH.pdb'+' -fi pdb -o '+fname+'_ligH_bcc.mol2'+' -fo mol2 -c bcc -nc '+ligand_charge)
##########
####################
##########
##########CHANGE TO YOUR PATH TO ROSETTA
##########
os.system('~/main/source/scripts/python/public/generic_potential/mol2genparams.py --nm '+ligname+' -s '+fname+'_ligH_bcc.mol2')
##########
##########
##########
##########
##########
##########
##########
newlig=open(fname+'_ligH_bcc_0001.pdb','r')
newliglines=[line for line in newlig.readlines() if line[:6]=='HETATM']
newpdb=open(os.path.join(fname+'_clean.pdb'),'w')
for line in alines:
    newpdb.write(line)
newpdb.write('TER\n')
for line in newliglines:
    newpdb.write(line)
newpdb.write('TER\n')
newpdb.close()
os.chdir('..')
