#!/bin/bash
#BSUB -P acc_motor # project name
#BSUB -q premium # queue name ('premium' is standard on Minerva)
#BSUB -n 10 # number of tasks in a parallel job (also submits as a parallel job)
#BSUB -R span[hosts=1] # resource requirements
#BSUB -R rusage[mem=20000] # resource requirements
#BSUB -W 10:00 # job runtime limit (HH:MM)
#BSUB -J /sc/arion/projects/motor/WILSOA28/demuxlet_analysis/snrnaseq__03_dea_setup/log_files/run__03_dea_setup_example__0
#BSUB -o /sc/arion/projects/motor/WILSOA28/demuxlet_analysis/snrnaseq__03_dea_setup/log_files/run__03_dea_setup_example__0.o
#BSUB -e /sc/arion/projects/motor/WILSOA28/demuxlet_analysis/snrnaseq__03_dea_setup/log_files/run__03_dea_setup_example__0.e
#BSUB -L /bin/bash

set -ev

# 01. Set up environment
# 01a. Activate conda environment
ml anaconda3/2020.11
ml -python
source /hpc/packages/minerva-centos7/anaconda3/2020.11/etc/profile.d/conda.sh
conda activate CO-deg-setup-env

# 01b. Get root directory
exec_dir=$( pwd )
cd "${exec_dir}"
# 01c. Define code path
sud_dea_dir="${exec_dir}/scripts"

# 02. Set up config files
# 02a. Specify config file path
cfg="${exec_dir}/configs/config_DEA_setup_example_SUD_DEA.yaml"
echo "${cfg}"
# 02b. Add root directory to config
# file if it does not exist in there
# yet (meaning the script hasn't been
# run before)
if ! $( grep -q 'root_dir' ${cfg} ); then
	echo "Initializing config file with current directory as root directory"
	echo "root_dir: '${exec_dir}'" >> ${cfg}
else
	echo "Config file already contains root directory; proceeding"
fi

# 03. Specify cell type upon which to run DESeq2 (for SUD DEA example)
group="dopaminergic_neuron"

# 04. Set up and run DEA on dopaminergic neuron expression dataset as an example
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


