import numpy as np
import pandas as pd
import anndata as ad
import os
import sys
import argparse
import yaml

# This script puts QCed, cell-typed expression data and
# clinical variables into the required format for running
# differential expression analysis (DEA) with DESeq2
# (through the Python wrapper program diffexpr).

# Script setup

# 00. Create an argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Reformatting of expression data and annotations' + \
        ' for differential expression analysis using diffexpr (a python wrapper of DESeq2).')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
args = parser.parse_args()

# 01. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01a.ii. Input expression data directory
input_data_dir = qc_root_dir + f'/' + cfg.get('input_data_dir')
# 01a.iii. Input expression data file name
input_data_fn = cfg.get('input_data_fn')
# 01a.iv. Differential expression analysis directory
dea_data_dir = qc_root_dir + f'/' + cfg.get('dea_data_dir')
if not os.path.isdir(dea_data_dir):
    os.system(f'mkdir -p {dea_data_dir}')
# 01a.v. DESeq2-formatted data output directory
dea_reformatting_output_dir = dea_data_dir + f'/' + cfg.get('dea_deseq_formatted_data_dir')
if not os.path.isdir(dea_reformatting_output_dir):
    os.system(f'mkdir -p {dea_reformatting_output_dir}')

# 01b. DEA design parameters
# 01b.i. Cell type label to use for dataset partitioning
partition_variable = cfg.get('dea_partition_variable')
# 01b.ii. List of covariates
all_covariates = cfg.get('dea_covariate_list')
# 01b.iii. Test condition name
tc_name = cfg.get('dea_test_condition')
if tc_name not in all_covariates:
    all_covariates.append(tc_name)
# 01b.iv. Dictionaries with instructions
# for generating the test condition of interest
# 01b.iv.A. Pre-build instructions
pre_build_instructions = cfg.get('dea_test_pre_build_instructions')
# 01b.iv.B. Build instructions
build_instructions = cfg.get('dea_test_cond_build_instructions')

# Build input DESeq2-formatted expression data
# and metadata for each group in the specified
# partitioning

# 02. Import expression data
ahc_fn = f'{input_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 03. Split expression data into counts and metadata to match DESeq2
# input data format
# Two components are required:
# (1) the counts array as a DataFrame with columns 1. gene_name 2.,...,N. obs_ID
# (i.e. the .X object transposed to be counts x observations, and with gene/row
# names also as the first column)
# (2) Barcode covariate annotations as a metadata DataFrame,
# with both the index and first column having observation IDs,
# and with all relevant metadata (e.g. partition group, covariates, and any
# other metadata of interest)

# 03a. Define a helper function to generate DESeq2 counts and metadata DataFrames
def generate_DESEQ2_dataframes(adata,
        partition_var,
        covariate_col_names):
    partition_values = np.unique(adata.obs[partition_var].values).tolist()
    print(f'Partitioning anndata by {partition_var}.')
    for p in partition_values:
        print(f'/tProcessing partition {p}.')
        adata_part = adata[adata.obs[partition_var] == p, :].copy()
        ahc_X = adata_part.X.todense().copy()
        # 03a.i. Take the count matrix transpose for the count DataFrame
        ahc_X_deseq2 = ahc_X.transpose()
        obs_barcode_list = adata_part.obs.index.values.tolist()
        count_df = pd.DataFrame(ahc_X_deseq2,
                columns=obs_barcode_list)
        # 03a.ii. Add a 'gene_name' column to the count DataFrame
        count_df['gene_name'] = adata_part.var_names.values.tolist()
        cols = count_df.columns.tolist()
        new_col_order = cols[-1:] + cols[:-1]
        count_df = count_df[new_col_order]
        # 03a.iii. Generate the metadata DataFrame
        metadata_col_names = [
                'patient_ID',
                'pool_name',
                'mito_frac',
                'n_umi',
                'n_gene'
                ] + [partition_var] + covariate_col_names
        metadata_df = adata_part.obs[metadata_col_names].copy()
        count_fn = f'{dea_reformatting_output_dir}/ahc_preproc_count_df__{partition_var}_{p}.csv'
        metadata_fn = f'{dea_reformatting_output_dir}/ahc_preproc_metadata_df__{partition_var}_{p}.csv'
        count_df.to_csv(count_fn,index=True)
        metadata_df.to_csv(metadata_fn,index=True)
        del count_df, metadata_df

# 03b. Split the expression data by partition variable (cell type) value
# and produce count and metadata DataFrames for each
generate_DESEQ2_dataframes(adata = adatas_human_qc,
        partition_var = partition_variable,
        covariate_col_names = all_covariates)


sys.exit()
