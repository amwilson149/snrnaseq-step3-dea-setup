#!/bin/bash

set -ev

# 01. Set up environment
exec_dir=$( pwd )
cd "${exec_dir}"
sud_dea_dir="${exec_dir}/scripts"

# 02. Set up config files
# 02a. Specify config file path
cfg="${exec_dir}/configs/config_DEA_setup_example_SUD_DEA.yaml"
echo "${cfg}"
# 02b. Add root directory to config
# file if it not specified
if ! $( grep -q 'root_dir' ${cfg} ); then
	echo "Initializing config file with current directory as root directory"
	echo "root_dir: '${exec_dir}'" >> ${cfg}
else
	echo "Config file already contains root directory; proceeding"
fi

# 03. Specify cell type(s) to set up for DEA (comma-separated list)
group="dopaminergic_neuron"

# 04. Run DEA setup
# 04a. Reformat dataset for input to DESeq2
echo -e "Reformatting expression data (into count and metadata files)...."
python "${sud_dea_dir}/03a_reformat_anndata_for_dea.py" --config-yaml-path ${cfg}
echo -e "\n"

# 04b. Partition DESeq2-formatted files for specified cell types into data
# subsets for k-fold cross-validation
echo -e "Generating data subsets for k-fold cross-validation...."
python "${sud_dea_dir}/03b_set_up_dea_run_kfcv_and_rs.py" --config-yaml-path ${cfg} --group-val-to-process-override "${group}" --new-writes-only False --turn-off-pseudocount False

# 04c. Run file setup checks
echo -e "Running file setup checks...."
python "${sud_dea_dir}/03c_check_dea_setup.py" --config-yaml-path ${cfg}

# 04d. Run random seed generation for individual cell types
# to create reproducible DEA runs for each
echo -e "Generating random seeds for individual cell types to facilitate reproducible DEA runs"
python "${sud_dea_dir}/03d_gen_rand_seeds_dea_parallelization.py" --config-yaml-path ${cfg}


