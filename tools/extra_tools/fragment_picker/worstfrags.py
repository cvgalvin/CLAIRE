'''
conda install h5py


'''


import os
from pyrosetta import *
init('-indexed_structure_store:fragment_store /wynton/home/kortemme/lhong/software/ss_grouped_vall_all.h5')
prm1='ed1UM_4_M59T56Y42T63_1_relaxed_relaxed_5tpj_293640_design_3_unrelaxed_model_1_rank_1_0001_design_2_unrelaxed_model_2_rank_1_0001_design_6_unrelaxed_model_2_rank_1_0001_hybrid_100_1_0001_clean__DE_1_oglig_0001.params'


filters_xml = f'''
                <FILTERS>
                    <worst9mer name="worst_9mer" confidence="0" rmsd_lookup_threshold="0.01" report_mean_median="true"/>
                </FILTERS>'''
w9f = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(filters_xml).get_filter('worst_9mer')
#
pdbs=[i for i in os.listdir() if i[-3:]=='pdb']
#
for pdb in pdbs[:1]:
    lig=[prm1]
    p=Pose()
    generate_nonstandard_residue_set(p,lig)
    pose_from_file(p, pdb)
    p.update_residue_neighbors()
    w9f.report_sm(p)
    w9f.score(p)




##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
import os
def generate_fragments_klab(input_pdb, data_path):
    '''Generate fragments and save the fragments into the data_path.
    Return the job id.
    '''
    cwd = os.getcwd()
    os.chdir(data_path)

    output_path = 'fragments'
    if os.path.exists('fragments'):
      shutil.rmtree(output_path)

    proc = subprocess.Popen(['klab_generate_fragments',
                           '-m', '100',
                           '-x', '10',
                           '-r', '48',
                           '-d', output_path,
                           input_pdb],
                           stdout=subprocess.PIPE)

    if proc.stdout:
      stdout = proc.stdout.read().decode()
      print(stdout)
    if proc.stderr:
      print(proc.stderr.read().decode())

    for line in stdout.split('\n'):
        m = re.match(r"Fragment generation jobs started with job ID (\d+).", line)
        if m:
          break

    os.chdir(cwd)

    return m.group(1)

generate_fragments_klab('relaxed_relaxed_99934_design_1_unrelaxed_model_3_rank_1_0001_design_5_unrelaxed_model_1_rank_1_0001_0001.pdb', 'testdp')
