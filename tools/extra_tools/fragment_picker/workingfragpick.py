'''
so this just uses xingjies software packages
and the default paths already in his make_fragments.pl script
the only extra thing i had to do was give permission to execute the files
via chmod ugo+x makefrags.pl

the code below is from xingjies runforwardfolding script, it is simply
a function to run the fragment picking protocol
though I believe I could also just run the makefrags perl script itself
'''
#
#
#
import os
import re
import subprocess
import shutil
import sys

fasta_inp=sys.argv[1]
frag_out=sys.argv[2]

if not os.path.exists(frag_out):
    os.makedirs(frag_out,exist_ok=True)



def generate_fragments(input_fasta, data_path):
    '''Generate fragments and save the fragments into the data_path.
    Return the job id.
    '''
    cwd = os.getcwd()
    os.chdir(data_path)

    output_path = 'fragments'
    if os.path.exists('fragments'):
      shutil.rmtree(output_path)

    os.makedirs('fragments', exist_ok=True)
    os.chdir('fragments')

    proc = subprocess.Popen(['qsub',
                            os.path.join('/wynton/home/kortemme/cgalvin/BSFF/tools/fragment_picker', 'fragment_picking_scripts/submit_frag_picking_Wynton.sh'),
                            os.path.join('/wynton/home/kortemme/cgalvin/BSFF/tools/fragment_picker', 'fragment_picking_scripts/makefragsxj.pl'),
                           input_fasta],
                           stdout=subprocess.PIPE)

    if proc.stdout:
      stdout = proc.stdout.read().decode()
      print(stdout)
    if proc.stderr:
      print(proc.stderr.read().decode())

    for line in stdout.split('\n'):
        m = re.match(r'Your job (\d+) \("submit_frag_picking_Wynton.sh"\) has been submitted', line)
        if m:
          break

    os.chdir(cwd)

    return m.group(1)


generate_fragments(fasta_inp,frag_out)


