# time python run_colabfold.py input_dir output_dir

import sys


input_dir = sys.argv[1]
result_dir = sys.argv[2]

# number of models to use
msa_mode = 'single_sequence'
# msa_mode = "MMseqs2 (UniRef only)"
num_models = 5
num_recycles = 12
stop_at_score = 100
use_custom_msa = False
use_amber = False
use_templates = False
do_not_overwrite_results = True
zip_results = True



from colabfold.batch import get_queries, run
from colabfold.download import default_data_dir
from colabfold.utils import setup_logging
from pathlib import Path

# For some reason we need that to get pdbfixer to import
# if use_amber and '/usr/local/lib/python3.7/site-packages/' not in sys.path:
#     sys.path.insert(0, '/usr/local/lib/python3.7/site-packages/')

setup_logging(Path(result_dir).joinpath('log.txt'))

queries, is_complex = get_queries(input_dir)

run(queries=queries,
    result_dir=result_dir,
    use_templates=use_templates,
    use_amber=use_amber,
    msa_mode=msa_mode,
    # model_type='auto',
    num_models=num_models,
    num_recycles=num_recycles,
    model_order=[3, 4, 5, 1, 2],
    is_complex=is_complex,
    data_dir=default_data_dir,
    keep_existing_results=do_not_overwrite_results,
    rank_mode='auto',
    pair_mode='unpaired+paired',
    stop_at_score=stop_at_score)
