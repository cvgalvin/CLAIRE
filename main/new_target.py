#usage:
#time python new_target.py target_name(3 letter code)

import os
import sys

target_name=sys.argv[1]

os.makedirs(target_name,exist_ok=True)
os.makedirs(os.path.join(target_name,'Inputs'),exist_ok=True)
os.makedirs(os.path.join(target_name,'Inputs','Fragment_Inputs'),exist_ok=True)
os.makedirs(os.path.join(target_name,'Inputs','Rosetta_Inputs'),exist_ok=True)
