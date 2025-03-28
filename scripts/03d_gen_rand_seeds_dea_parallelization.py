import numpy as np
import pandas as pd
import anndata as ad
import os
import sys
import pickle
import argparse
import yaml
from utils.get_rand_seed import get_rand_seed

# This script generates a set of random seeds, one
# for each potential group and fold specified for DEA,
# to use for initializing the global numpy random number
# generator for each of those groups.
# This step is meant to ensure repeatability in the
# parallelized DEA runs for each fold-removed subset
# of each partition in the expression data.

# The random state generator from after cell typing is
# used to initialize this random seed generation so
# that DEA can be started for some partitions before
# setup for all partitions has been completed.
# Since these integers will themselves be used as
# random seeds, the relative randomness of DEA and
# fold IDs should be preserved.

# Script setup

# 00. Create an argparse object for reading input arguments
parser = argparse.ArgumentParser(description='K-fold cross-validated differential expression' + \
        ' analysis using diffexpr.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
args = parser.parse_args()

# 01. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01a.iii. Post-cell-typing input expression data directory
cell_typed_data_dir = qc_root_dir + f'/' + cfg.get('input_data_dir')
# 01a.iv. Post-cell-typing input expression data file name
cell_typed_input_data_fn = cfg.get('input_data_fn')

# 01b. DEA partitioning parameters
# 01b.i. DEA partition variable
grouping = cfg.get('dea_partition_variable')
# 01b.ii. The number of folds to use
# for k-fold cross-validation
N_FOLDS_KFCV = cfg.get('dea_n_folds_kfcv')

# 01c. Random state parameters
# 01c.i. Input random state directory
input_rand_state_dir = qc_root_dir + f'/' + cfg.get('input_rng_state_dir')
# 01c.ii. Input random state file name
rand_state_input_fn = cfg.get('rand_state_input_fn')
# 01c.iii. Output random state directory
output_rand_state_dir = qc_root_dir + f'/' + cfg.get('output_rng_state_dir')
if not os.path.isdir(output_rand_state_dir):
    os.system(f'mkdir -p {output_rand_state_dir}')
# 01c.iv. Output random state file name
rand_state_output_fn = cfg.get('rand_state_output_fn__gen_dea_seeds')
# 01c.v. DEA random seed directory 
rand_seed_output_dir = output_rand_state_dir + f'/' + cfg.get('rng_dea_seeds_dir')
# 01c.v. DEA random seed output file name tag
rand_seed_fn_tag = cfg.get('rng_dea_seed_file_tag')

# 02. Make DEA random seed output directory if needed
if not os.path.isdir(rand_seed_output_dir):
    os.system(f'mkdir -p {rand_seed_output_dir}')

# 03. Import expression data
ahc_fn = f'{cell_typed_data_dir}/{cell_typed_input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 04. Initialize a random number generator and
# assign the input random state to it
random_state_full_fn = f'{input_rand_state_dir}/{rand_state_input_fn}'
numpy_random_state = np.random.RandomState()
working_np_state = None
with open(random_state_full_fn,'rb') as f:
    working_np_state = pickle.load(f)
numpy_random_state.set_state(working_np_state)

# 04. Get the groups of interest
groups = list(
        np.unique(
            adatas_human_qc.obs[grouping].values.tolist()
            )
        )

# 05. Delete the AnnData object to save space
del adatas_human_qc

# 06. Spawn a file for each group and fold containing a
# random seed
for group in groups:
    for fold in np.arange(1,N_FOLDS_KFCV+1):
        rs_curr = get_rand_seed(numpy_random_state)
        output_full_fn = f'{rand_seed_output_dir}/dea_rand_seed__{rand_seed_fn_tag}__{group}__{fold}.txt'
        with open(output_full_fn,'w') as f:
            f.write(f'{rs_curr}\n')

# 07. Save the current state of the random
# number generator as a branch, in case
# additional rounds of random seed generation
# are desired
numpy_random_state_dict = numpy_random_state.get_state(
        legacy=False
        )
numpy_random_state_full_output_fn = f'{output_rand_state_dir}/{rand_state_output_fn}'
with open(numpy_random_state_full_output_fn,'wb') as outfile:
    pickle.dump(numpy_random_state_dict,
            outfile)

sys.exit()

