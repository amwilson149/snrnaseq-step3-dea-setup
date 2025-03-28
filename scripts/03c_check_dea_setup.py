import numpy as np
import pandas as pd
import os
import sys
import argparse
import yaml

# This script checks that the fold-removed metadata 
# and count files were appropriately written, by 
# determining
# (1) whether the number of rows in the counts DataFrame
#     matches the number of columns in the corresponding
#     metadata DataFrame
# (2) whether the sum of barcodes across all the fold-removed
#     cases match the set of barcodes in the original, not-fold
#     removed version of that partition's count/metadata data

# It then writes a text file summarizing whether each diffexpr-ready
# partition's files passed or failed each of these checks.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='K-fold cross-validated differential expression' + \
        ' analysis using diffexpr.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
parser.add_argument(f'--groups-to-check-override',type=str,required=False)
args = parser.parse_args()

# 01. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01a.ii. Differential expression analysis directory
dea_data_dir = qc_root_dir + f'/' + cfg.get('dea_data_dir')
# 01a.iii. Diffexpr-formatted file directory
dea_diffexpr_dir = dea_data_dir + f'/' + cfg.get('dea_deseq_formatted_data_dir')
# 01a.iv. Directory with diffexpr-ready files
dea_diffexpr_out_dir = dea_diffexpr_dir + f'/' + cfg.get('dea_ready_deseq_formatted_data_dir')
# 01a.v. Validity check summary output file name root
output_summary_fn_root = cfg.get('dea_pass_fail_check_fn_root')

# 02. Get the list of groups in the diffexpr-ready file directory
groups_to_check = [
        _ 
        for _ in 
        os.listdir(dea_diffexpr_out_dir)
        if os.path.isdir(f'{dea_diffexpr_out_dir}/{_}')
        ]
if args.groups_to_check_override:
    groups_to_check = args.groups_to_check_override.split('[')[
            -1].split(']')[
                    0].split(',')

groups_str = ', '.join(groups_to_check)
print(f'Groups that will be checked: {groups_str}')

# 03. Pull the DESeq2-formatted files for each group,
# verifying that they have been fully written
# 03a. Pull the list of DESeq2 formatted files that
# preceded k-fold cross-validation setup
presetup_files = [
        _ for _ in
        os.listdir(dea_diffexpr_dir)
        if os.path.isfile(f'{dea_diffexpr_dir}/{_}')
        ]

# 03b. Check the k-fold cross-validation files
# for groups in the above DESeq2-formatted file list
# (this step restricts checks to avoid any spurious
# test files)
groups_checked = []
checks_df = pd.DataFrame()
c_idx = 0
for group in groups_to_check:
    print(f'Processing {group}')
    # 03b.i. Pull the current group/cell type's
    # pre-setup files
    m_barcodes_pre = []
    c_barcodes_pre = []
    precursor_files = [
            _ for _ in 
            presetup_files
            if group in _
            ]
    if len(precursor_files) > 0:
        groups_checked.append(group)
        checks_df.at[c_idx,'group'] = group
        # 03b.ii. Get count and metadata barcodes
        # in these pre-setup files
        mfn_pre = [_ for _ in precursor_files if 'metadata' in _][0]
        cfn_pre = [_ for _ in precursor_files if 'count' in _][0]
        mfn_pre_full = f'{dea_diffexpr_dir}/{mfn_pre}'
        cfn_pre_full = f'{dea_diffexpr_dir}/{cfn_pre}'
        mdf_pre = pd.read_csv(mfn_pre_full,index_col=0)
        cdf_pre = pd.read_csv(cfn_pre_full,index_col=0)
        m_barcodes_pre = mdf_pre.index.values.tolist()
        c_barcodes_pre = [
                _ for _ in cdf_pre.columns.values.tolist()
                if _!='gene_name'
                ]
        # 03b.iii. Delete pre-setup files to save space
        del mdf_pre, cdf_pre
    # 03b.iv. Iterate through post-setup (k-fold subsetted)
    # files to pull barcodes in each
    postsetup_dir = f'{dea_diffexpr_out_dir}/{group}'
    postsetup_files = [
            _ for _ in os.listdir(postsetup_dir)
            if '.csv' in _           
            ]
    fold_tags = list(
            np.unique(
                [
                    _.split('__')[-1]
                    for _ in postsetup_files
                    ]
                )
            )
    m_barcodes_post = []
    c_barcodes_post = []
    for fold_tag in fold_tags:
        print(f'Processing fold {fold_tag}')
        mfn_post = [
                _ for _ in postsetup_files
                if (
                    ('metadata' in _)
                    &
                    (fold_tag in _)
                    )
                ][0]
        cfn_post = [
                _ for _ in postsetup_files
                if (
                    ('counts' in _)
                    &
                    (fold_tag in _)
                    )
                ][0]
        mfn_post_full = f'{postsetup_dir}/{mfn_post}'
        cfn_post_full = f'{postsetup_dir}/{cfn_post}'
        mdf_post_curr = pd.read_csv(mfn_post_full,index_col=0)
        cdf_post_curr = pd.read_csv(cfn_post_full,index_col=0)
        mdf_bc_post_curr = mdf_post_curr.index.values.tolist()
        cdf_bc_post_curr = [
                _ for _ in cdf_post_curr.columns.values.tolist()
                if _!='gene_name'
                ]
        print(f'\tMetadata rows equal to counts columns: {mdf_bc_post_curr==cdf_bc_post_curr}')
        m_barcodes_post.extend(mdf_bc_post_curr)
        c_barcodes_post.extend(cdf_bc_post_curr)
        # 03b.v. Record whether the current fold-removed
        # metadata and counts barcodes match
        checks_df.at[c_idx,f'fold_{fold_tag}_eq'] = (
                'P' if mdf_bc_post_curr==cdf_bc_post_curr
                else 'F'
                )
    print(f'\n')
    # 03b.vi. Record whether everything is repeated across
    # folds n_folds - 1 times
    mbc_post_u, n_mbc_post_u = np.unique(
            m_barcodes_post,
            return_counts=True
            )
    cbc_post_u, n_cbc_post_u = np.unique(
            c_barcodes_post,
            return_counts=True
            )
    nu_mbs_post_u = list(np.unique(n_mbc_post_u))
    cu_mbs_post_u = list(np.unique(n_cbc_post_u))
    print(f'\tN repeats observed per barcode:\n' + \
            f'\t\tmetadata: {nu_mbs_post_u}\n' + \
            f'\t\tcounts: {cu_mbs_post_u}')
    checks_df.at[c_idx,f'm_n_repeats_okay'] = (
            'P' if nu_mbs_post_u==[len(fold_tags)-1]
            else 'F'
            )
    checks_df.at[c_idx,f'c_n_repeats_okay'] = (
            'P' if cu_mbs_post_u==[len(fold_tags)-1]
            else 'F'
            )
    # 03b.vii. Record whether the total number
    # of barcodes covered across folds in c and m
    # matches the precursor files
    print(f'\n')
    print(f'\tN barcodes pre, total:\n' + \
            f'\t\tmetadata: {len(m_barcodes_pre)}\n' + \
            f'\t\tcounts: {len(c_barcodes_pre)}')
    print(f'\tN barcodes post, total:\n' + \
            f'\t\tmetadata: {len(mbc_post_u)}\n' + \
            f'\t\tcounts: {len(cbc_post_u)}')
    m_post_sort = sorted(list(mbc_post_u))
    m_pre_sort = sorted(m_barcodes_pre)
    checks_df.at[c_idx,f'm_bc_covered'] = (
            'P' if m_post_sort==m_pre_sort
            else 'F'
            )
    c_post_sort = sorted(list(cbc_post_u))
    c_pre_sort = sorted(list(c_barcodes_pre))
    checks_df.at[c_idx,f'c_bc_covered'] = (
            'P' if c_post_sort==c_pre_sort
            else 'F'
            )
    print(f'\n\n')
    c_idx += 1

# 04. Write the summary table to file
groups_checked_str = '_'.join(groups_checked)
output_full_fn = dea_diffexpr_out_dir + f'/{output_summary_fn_root}__{groups_checked_str}.csv'
checks_df.to_csv(
        output_full_fn,
        index=True
        )

sys.exit()













