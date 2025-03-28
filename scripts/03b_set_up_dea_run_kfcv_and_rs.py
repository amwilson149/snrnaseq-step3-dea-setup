import numpy as np
import pandas as pd
import math
import os
import sys
import time
import pickle
import argparse
import yaml
from utils.consts import *
from utils.get_rand_seed import *

print('Packages imported.')

# This script sets up differential expression analysis
# (DEA) for the specified expression dataset, by
# (1) assigning fold indices for k-fold cross-validation partitioning,
# (2) composing and validating the specified design function,
# (3) determining whether a sufficient number of barcodes
#     exists to meet statistical power requirements (determined
#     in another script),
# (4) writing the barcode-to-fold map and design function
#     for each data group being analyzed,
# (5) writing a list of the test conditions and groups that
#     are powered for analysis.

# Script setup

# 00. Create argparse objects for reading input arguments
# 00a. Configuration file
parser = argparse.ArgumentParser(description='K-fold cross-validated differential expression' + \
        ' analysis using diffexpr.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
# 00b. Groups to process for DEA
parser.add_argument(f'--group-val-to-process-override',type=str,required=False)
# 00c. Flag specifying whether to only construct files that
# have not yet been written (to handle the case of a long-duration run
# that is interrupted)
parser.add_argument(f'--new-writes-only',type=str,required=False)
# 00d. Flag specifying whether to add a pseudocount or not
parser.add_argument(f'--turn-off-pseudocount',type=str,required=False)
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
# 01a.iii. Diffexpr-formatted input file directory
dea_diffexpr_dir = dea_data_dir + f'/' + cfg.get('dea_deseq_formatted_data_dir')
# 01a.iv. Output directory for post-setup, diffexpr-ready files
dea_diffexpr_out_dir = dea_diffexpr_dir + f'/' + cfg.get('dea_ready_deseq_formatted_data_dir')
# 01a.v. File name tag for post-setup, diffexpr-ready files
dea_diffexpr_out_fn_tag = cfg.get('dea_ready_deseq_formatted_data_fn_tag')

# 01b. Differential expression analysis parameters
# 01b.i. Test type
TEST_TYPE = cfg.get('dea_test_type')
# 01b.ii. Minimum number of barcodes per group
# required to permit analysis
# Note: If this value is specified to be 'None',
# this code uses as a minimum the greater of the
# following numbers:
# (1) the minimum value needed to determine
#     a solution for all parameters (the number
#     of design function covariates plus a constant)
# (2) the rough sample size needed to meet the
#     specified DEA test's assumptions or a user-specified
#     minimum number of barcodes (specified
#     in the configuration file)
N_BARCODES_MIN_DESEQ2 = cfg.get('dea_min_n_barcodes')
if N_BARCODES_MIN_DESEQ2 == 'None':
    N_BARCODES_MIN_DESEQ2 = None
# 01b.iii. Minimum allowed sample size (number of barcodes)
MIN_SAMPLE_SIZE = cfg.get('dea_min_sample_size')
# 01b.iii. Number of folds to use for k-fold cross-validation
N_FOLDS_KFCV = cfg.get('dea_n_folds_kfcv')
# 01b.iii.A. Optional dictionary identifying k-fold
# splitting for any exception groups
kfcv_exceptions_dict = {}
try:
    kfcv_exceptions_dict = cfg.get('dea_n_fold_kfcv_exceptions_dict')
except:
    kfcv_exceptions_dict = {}
if kfcv_exceptions_dict is None:
    kfcv_exceptions_dict = {}

# 01b.iv. List of non-covariate columns
non_covariates = cfg.get('dea_non_covariates')
# 01b.v. Grouping variable used to split dataset
# for DEA (e.g. cell type)
grouping = cfg.get('dea_partition_variable')
# 01b.v.A. If the grouping variable is a mapping,
# generate it for the dataset
if type(grouping)==dict:
    grouping = grouping['new_col_name']
non_covariates += [grouping]
# 01b.vi. Optional override list of group values
# to process (if not None, overrides the default
# of processing all groups)
grouping_value_override = args.group_val_to_process_override
if grouping_value_override == 'None':
    grouping_value_override = None
# 01b.vi. Additional covariates to ignore for DEA
covariates_to_ignore = cfg.get('dea_covariates_to_ignore')
if covariates_to_ignore == 'None':
    covariates_to_ignore = None
non_covariates += covariates_to_ignore
# 01b.vii. DEA test condition
test_condition = cfg.get('dea_test_condition')
# 01b.viii. List of test condition levels
poss_tc_levels = cfg.get('dea_test_cond_levels')
# 01b.ix. List of test condition level comparisons
condition_comparison_list = cfg.get('dea_cond_comparison_list')
# 01b.x. Test condition build instructions
tc_build_instruction_dict = cfg.get('dea_test_cond_build_instructions')

# 01c. Random state information
# 01c.i. Input random state directory
input_rand_state_dir = qc_root_dir + f'/' + cfg.get('input_rng_state_dir')
# 01c.ii. Output random state directory
output_rand_state_dir = qc_root_dir + f'/' + cfg.get('output_rng_state_dir')
if not os.path.isdir(output_rand_state_dir):
    os.system(f'mkdir -p {output_rand_state_dir}')
# 01c.iii. Input random state file name
rand_state_input_fn = cfg.get('rand_state_input_fn')
# 01c.iv. Output random state file name
rand_state_output_fn = cfg.get('rand_state_output_fn__kfcv_setup')

# 01d. K-fold cross-validation
# summary information
nfolds_dir = qc_root_dir + f'/' + cfg.get('kfcv_n_fold_log_dir')
if not os.path.isdir(nfolds_dir):
    os.system(f'mkdir -p {nfolds_dir}')
nfolds_out_fn = cfg.get('n_folds_used_fn__03_dea__01')

# 01e. Flag specifying whether to only write non-existing files
new_writes_only = args.new_writes_only
if not new_writes_only:
    new_writes_only = False
else:
    new_writes_only = (
            True if new_writes_only.lower() in ['t','true']
            else False
            )

# 01f.  Optional flag to turn off pseudocount addition
# (for an alternative form of DEA)
pseudocount_modifier = args.turn_off_pseudocount
ADD_PSEUDOCOUNT = True
if pseudocount_modifier.lower()=='true':
    ADD_PSEUDOCOUNT = False
    print(f'Preparing expression data for DEA WITHOUT adding a pseudocount.')

# 02. Make the diffexpr-ready file directory if it does not exist
if not os.path.isdir(dea_diffexpr_out_dir):
    os.system(f'mkdir -p {dea_diffexpr_out_dir}')

# 03. Pull DESeq2-reformatted data files for the specified groups
dea_diffexpr_contents = os.listdir(dea_diffexpr_dir)
group_files = [_ for _ in dea_diffexpr_contents
        if (
            (os.path.isfile(os.path.join(dea_diffexpr_dir,_))
                & ('ahc_preproc' in _)
                )
            )
        ]

# 04. Identify groups for which all reformatted data needed
# for setup are present
grouping_values = []
if grouping_value_override:
    grouping_values = [grouping_value_override]
else:
    grouping_values_raw, gv_counts = np.unique(
            [_.split(f'{grouping}_')[-1].split('.')[0]
                for _ in group_files
                ],
            return_counts=True
            )
    grouping_values = [_ for i,_ in enumerate(grouping_values_raw)
            if gv_counts[i] == 2]

print(f'Group values that will be set up for DEA:',
        f'{grouping_values}\n')

# 05. Initialize a random number generator and
# assign the input random state to it
random_state_full_fn = f'{input_rand_state_dir}/{rand_state_input_fn}'
numpy_random_state = np.random.RandomState()
working_np_state = None
with open(random_state_full_fn,'rb') as f:
    working_np_state = pickle.load(f)
numpy_random_state.set_state(working_np_state)

# 06. Define infrastructure for DEA setup
# 06a. Function to standardize numeric covariates
def standardize_data(df_col):
    raw_vals_list = df_col.copy().values.tolist()
    mean_val = np.mean(raw_vals_list)
    std_val = np.std(raw_vals_list)
    standardized_values = [(_ - mean_val)/std_val
            for _ in raw_vals_list]
    return standardized_values

# 06b. Function to determine whether any
# categorical variables in a design function
# are equivalent
def test_covariate_equivalence(df,
        cov_col_list):
    # 06b.i. Binarize each covariate, always
    # starting with zero as the first element
    binary_cols_dict = {}
    for col in cov_col_list:
        vals_raw = df[col].values.tolist()
        vals_bin = [0
                if _ == vals_raw[0]
                else 1
                for _ in vals_raw]
        binary_cols_dict[col] = vals_bin
    # 06b.ii. Compare all possible covariate
    # pairs
    equivalent_factors = []
    for i in range(len(cov_col_list)-1):
        col1 = cov_col_list[i]
        col2 = cov_col_list[i+1]
        b1 = binary_cols_dict[col1]
        b2 = binary_cols_dict[col2]
        if b1 == b2:
            equivalent_factors.append(
                    [col1,col2])
    return equivalent_factors

# 06c. Function to build the DESeq2 design formula
def build_design_formula(design_form_ordered_terms,
        test_term,
        equiv):
    perfect_tc_cov = False
    if len(equiv) > 0:
        for eqpair in equiv:
            print(f'Terms {eqpair[0]} and {eqpair[1]} encode equivalent information.')
            # 06c.i. If either term is equivalent to the test condition,
            # print a flag
            if test_term in eqpair:
                t_eq = [_ for _ in eqpair if _!=test_term]
                print(f'{t_eq} is equivalent to the test condition!')
                perfect_tc_cov = True
            # 06c.ii. Check whether both terms are still in the design formula
            # if so, remove the later term; if not, do nothing
            if all(
                    [
                        _ in design_form_ordered_terms
                        for _ in
                        eqpair
                        ]
                    ):
                print(f'Removing {eqpair[1]} from design formula.')
                design_form_ordered_terms.remove(eqpair[1])
    # 06c.iii. Build design formula with covariate degeneracies
    # removed (unless they are equivalent to the test condition)
    design_formula = '~ ' + ' + '.join(design_form_ordered_terms)
    n_terms = len(design_form_ordered_terms)
    return design_formula, n_terms, perfect_tc_cov

# 06d. Initialize list of groups ready for DEA
diffexpr_ready_groups = []

# 07. Iterate over groups and perform DEA setup
for g in grouping_values:
    print(f'Performing DEA setup for group {g}....\n')
    # 07a. Import counts and metadata
    print(f'Importing count and metadata DataFrames in DESeq2 format for {grouping} {g}...')
    counts_fn = f'{dea_diffexpr_dir}/ahc_preproc_count_df__{grouping}_{g}.csv'
    metadata_fn = f'{dea_diffexpr_dir}/ahc_preproc_metadata_df__{grouping}_{g}.csv'
    counts_df = pd.read_csv(counts_fn,index_col=0)
    metadata_df = pd.read_csv(metadata_fn,index_col=0)

    # 07b. Assign barcodes folds for k-fold cross-validation,
    # using either the number of folds specified in the
    # configuration file, or for groups in the dictionary
    # of k-fold cross-validation exceptions, the specified
    # number of folds and sampling framework
    # 07b.i. Get kfcv partitioning parameters
    # 07b.i.A. Set default partitioning parameters
    n_folds_to_use = N_FOLDS_KFCV
    n_folds_removed_per_partition = 1
    n_partitions = N_FOLDS_KFCV
    # 07b.i.B. Set flag allowing partitions
    # to have overlapping sets of folds by
    # default (which is appropriate when
    # doing leave-one-out kfcv or otherwise
    # when only a small number of folds is
    # removed per partition, for example when
    # the total number of observations is on
    # the smaller side; this approach allows
    # much of the data to still be used when
    # controlling for idiosyncracies in the
    # data using kfcv)
    use_non_overlapping_kept_folds = False
    # 07b.i.C. Adjust kfcv parameters
    # if the current group is listed
    # as a kfcv exception
    if g in kfcv_exceptions_dict.keys():
        n_folds_to_use = kfcv_exceptions_dict[g]['n_folds_total']
        n_folds_removed_per_partition = kfcv_exceptions_dict[g]['n_folds_removed_per_partition']
        n_partitions = kfcv_exceptions_dict[g]['n_partitions']
    # 07b.i.D. Compute the number of folds
    # to be kept in each partition
    n_folds_kept_per_partition = n_folds_to_use - n_folds_removed_per_partition
    # 07b.i.E. Require the fold sets in 
    # each partition to be non-overlapping
    # if the kfcv setup is such that the
    # total number of folds required to build
    # all partitions for kfcv is less than or
    # equal to the number of folds available
    # This approach allows for more of the
    # data to be used in kfcv when sampling
    # of the dataset is sparser
    if n_folds_kept_per_partition*n_partitions <= n_folds_to_use:
        use_non_overlapping_kept_folds = True
    # 07b.ii. Determine the folds that will be
    # used for each data partition given the
    # parameters above
    all_fold_ids = np.arange(1,n_folds_to_use+1)
    fold_ids_left = all_fold_ids
    folds_removed_per_partition_dict = {}
    print(f'Creating {n_partitions} partitions for k-fold cross validation,\n' + \
            f'each with {n_folds_removed_per_partition} fold' + \
            ('s ' if (n_folds_removed_per_partition>1) else ' ')  + \
            f'removed and {n_folds_kept_per_partition} fold' + \
            ('s ' if (n_folds_kept_per_partition>1) else ' ') + \
            f'retained (from {n_folds_to_use} folds total)'
            )
    # 07b.ii.A. Initialize fold set variables
    folds_to_keep_curr = []
    folds_to_remove_curr = []
    chosen_fold_set_list = []
    # 07b.ii.B. Evolve the numpy random state
    # by generating a number of random values
    # that is specific to the ASCII number sequence
    # for the group's name
    # This step prevents the same sequence
    # of fold sets from being generated
    # for every grouping if this code
    # is run in parallel for a collection
    # of groupings
    grouping_ascii_vals = [
            ord(char)
            for char in
            g
            ]
    for iter_curr in range(len(grouping_ascii_vals)):
        dummy = numpy_random_state.choice(
                grouping_ascii_vals,
                size=1,
                replace=True
                )
    del dummy
    if use_non_overlapping_kept_folds==True:
        # 07b.ii.B. Build partitions to have
        # mutually exclusive sets of folds
        print(f'\tForcing kept fold sets to be non-overlapping')
        for partition_idx in range(1,n_partitions+1):
            print(f'\tGenerating partition {partition_idx}')
            folds_to_keep_curr = numpy_random_state.choice(
                fold_ids_left,
                size=n_folds_kept_per_partition,
                replace=False
                )
            # 07b.ii.C. Print summary of kept and removed
            # folds for the current partition
            fk_str = [f'{_}' for _ in folds_to_keep_curr]
            fkcurr_str = ', '.join(fk_str)
            folds_to_remove_curr = list(
                    set(all_fold_ids) - set(folds_to_keep_curr)
                    )
            fr_str = [f'{_}' for _ in folds_to_remove_curr] 
            ftrc_str = '_'.join(fr_str)
            frcurr_str = ', '.join(fr_str)
            print(f'\t\tFolds kept:\n\t\t\t{fkcurr_str}')
            print(f'\t\tFolds removed to create current partition:\n\t\t\t{frcurr_str}')
            # 07b.ii.D. Add information about removed folds
            # to the kfcv partition metadata dictionary
            # Note: partitions are defined and named in
            # terms of the removed folds to be consistent
            # with the leave-one-out kfcv naming convention,
            # in which the dataset with fold 1 removed
            # would be given the suffix '_1'.
            folds_removed_per_partition_dict[ftrc_str] = folds_to_remove_curr
            # 07b.ii.E. Remove currently selected kept folds
            # from the list of possible choices so that
            # subsequent partitions cannot have those
            # folds kept
            fold_ids_left = [
                    _ for _ in
                    fold_ids_left
                    if _ not in
                    folds_to_keep_curr
                    ]
            print(f'\n')
    else:
        # 07b.ii.F. Build partitions to allow kept folds
        # to be partially overlapping, but to maximize
        # data usage during kfcv by requiring the sets
        # of removed folds to be non-overlapping
        print(f'Forcing removed fold sets to be non-repeating')
        # 07b.ii.G. Check that the number of partitions
        # to be built is fewer than the maximum number
        # of possible samplings, and if so, force the
        # number of partitions to be smaller than
        # that. Otherwise, the fold computation would
        # be unnecessarily repetitive.
        # The maximum number of different partitions
        # is the number of ways to choose the number
        # of removed folds from the total number of folds.
        max_n_different_partitions = math.comb(
                n_folds_to_use,
                n_folds_removed_per_partition
                )
        if n_partitions > max_n_different_partitions:
            print(f'\tThe number of partitions specified is greater than\n' + \
                    f'\tthe number of unique partitions that can be created.\n' + \
                    f'\tSetting the number of partitions to {max_n_different_partitions}' + \
                    f' and proceeding.'
                    )
            n_partitions = max_n_different_partitions
        for partition_idx in range(1,n_partitions+1):
            print(f'\tGenerating partition {partition_idx}')
            folds_to_remove_curr = numpy_random_state.choice(
                    fold_ids_left,
                    size=n_folds_removed_per_partition,
                    replace=False
                    )
            print(f'\t\t\tfolds to remove current: {folds_to_remove_curr}')
            print(f'\t\t\tchosen fold set list: {chosen_fold_set_list}')
            while any(
                    [
                        folds_to_remove_curr==_
                        for _ in
                        chosen_fold_set_list
                        ]
                    ):
                folds_to_remove_curr = numpy_random_state.choice(
                        fold_ids_left,
                        size=n_folds_removed_per_partition,
                        replace=False
                        )
            # 07b.ii.H. Print summary information about kept and removed folds
            folds_to_keep_curr = list(
                    set(all_fold_ids) - set(folds_to_remove_curr)
                    )
            fk_str = [f'{_}' for _ in folds_to_keep_curr]
            fkcurr_str = ', '.join(fk_str)
            fr_str = [f'{_}' for _ in folds_to_remove_curr]
            ftrc_str = '_'.join(fr_str)
            frcurr_str = ', '.join(fr_str)
            print(f'\t\tFolds kept:\n\t\t\t{fkcurr_str}')
            print(f'\t\tFolds removed to create current partition:\n\t\t\t{frcurr_str}')
            # 07b.ii.I. Add information about the current
            # partition to the partition metadata dictionary
            folds_removed_per_partition_dict[ftrc_str] = folds_to_remove_curr
            # 07b.ii.J. In the special case of leave-one-out
            # kfcv, where one fold is removed at a time,
            # speed up the building of subsequent partitions
            # by removing the selected fold from the list
            # of possible options
            if n_folds_removed_per_partition==1:
                fold_ids_left = [
                        _ for _ in
                        fold_ids_left
                        if _ not in
                        folds_to_remove_curr
                        ]
            # 07b.ii.K. Add the removed fold set for the
            # current partition to the list of previously
            # chosen fold sets
            chosen_fold_set_list.append(folds_to_remove_curr)
            print(f'\n')
            print(f'\n')

    print(f'\tSummary of partitions:')
    for pname,pinfo in folds_removed_per_partition_dict.items():
        print(f'\t{pname}:\t{pinfo}')

    # 07b.iii. Assign all possible fold IDs to
    # barcodes for the current group
    barcodes = metadata_df.index.values.tolist()
    barcode_fold_df = pd.DataFrame(
            index = barcodes,
            columns = ['fold_ID']
            )
    barcode_fold_df['fold_ID'] = list(
            numpy_random_state.randint(
                low=1,
                high=n_folds_to_use+1,
                size=len(barcode_fold_df)
                )
            )

    # 07c. Preprocess counts and metadata
    # 07c.i. If the test_condition (e.g. 'hiv_sud_status')
    # is not in the list of covariate columns,
    # build it and remove the component covariates
    # 07c.ii. Get list of component names and how to
    # map their values
    components_dict = tc_build_instruction_dict['components']
    if test_condition not in metadata_df.columns.values.tolist():
        print(f'\n\tGenerating test condition values')
        for row_idx in barcodes:
            row = metadata_df.loc[row_idx]
            new_val_list = []
            # 07c.iii. Get the character that should
            # be used to join the order list of mapped
            # values from the component covariates
            join_char = tc_build_instruction_dict['join_char']
            for component,comp_map in components_dict.items():
                comp_curr = row[component]
                comp_curr_mapped = comp_map[comp_curr]
                new_val_list.append(comp_curr_mapped)
            metadata_df.at[row_idx,test_condition] = join_char.join(new_val_list)
    # 07c.iv. Remove the component covariates from the metadata DataFrame,
    # should they still exist, regardless of whether the test
    # condition was already in the DataFrame
    metadata_df.drop(
            columns=list(
                components_dict.keys()
                ),
            inplace=True
            )

    print(f'\n\tProcessing covariates')
    # 07c.v. Get the list of covariates to process
    # from the metadata DataFrame
    covariates_to_process = [
            _ for _ in
            metadata_df.columns.values.tolist()
            if _ not in non_covariates
            ]
    
    # 07c.vi. Separate numeric and categorical covariates
    categorical_covariates = [
            _ for _ in
            covariates_to_process
            if metadata_df[_].dtype=='object'
            ]
    numeric_covariates = [
            _ for _ in
            covariates_to_process
            if _ not in categorical_covariates
            ]
        
    # 07c.vii. Set categorical types as dtype 'category'
    # to override any legacy import as dtype 'object'
    for categorical_covariate in categorical_covariates:
        metadata_df[categorical_covariate] = metadata_df[
                categorical_covariate
                ].astype('category')
        
    # 07c.viii. Check that the columns of count matrices and
    # rows of metadata matrices are in the same order
    # and have the same names
    # Note: in the diffexpr version of DESeq2, the first column
    # of the counts matrix contains the gene names, so for
    # our comparison we need to ignore it
    print(f'\tChecking that count and metadata observations are in the same order')
    count_cols = counts_df.columns.tolist()[1:]
    metadata_rows = metadata_df.index.tolist()
    records_match = (count_cols == metadata_rows)
    if records_match:
        print(f'\t\tCount and metadata observations are in the same order. Proceeding...')
    else:
        print(f'\t\tCount and metadata observations are not in the same order.\n\t\t' +
                f'Reordering count columns to match metadata rows...')
        # Reorder columns and make sure the first one
        # still contains the gene names
        counts_df = counts_df[['gene_name'] + metadata_rows]
        del count_cols, metadata_rows, records_match

    # 07d. For each fold-removal case in the partition, 
    # (1) Check whether there are any equivalent pairs of categorical covariates.
    #     If yes for any fold-removed sample, remove the degeneracy from all
    #     fold-removed design formulas for consistency.
    # (2) Check whether there are any non-varying covariates; if yes for any
    #     fold-removed set, remove the covariate for all of them
    # (3) Check whether there are sufficient sample sizes for each comparison.
    eq_pairs_across_folds = []
    nonvarying_covs_across_folds = []
    tc_ss_map_dict = {}
    print(f'\n\tBuilding counts and metadata files for each partition')
    #for fold_removed in np.arange(1,N_FOLDS_KFCV+1):
    for fold_removed_name,fold_removed in folds_removed_per_partition_dict.items():
        print(f'\n\t\tBuilding partition {fold_removed_name}')
        # 07d.i. Get barcodes not in the fold being removed
        bc_idxs_curr = barcode_fold_df.loc[
                ~barcode_fold_df['fold_ID'].isin(fold_removed)
                ].index.values.tolist()
        # 07d.ii. Get the fold-removed metadata DataFrame
        metadata_df_curr = metadata_df.loc[bc_idxs_curr].copy()
        # 07d.iii. Find equivalent covariate pairs
        eq_pairs_curr = test_covariate_equivalence(
                df = metadata_df_curr,
                cov_col_list = categorical_covariates 
                )
        # 07d.iv. Find nonvarying covariates
        for cov_curr in (categorical_covariates + numeric_covariates):
            if len(
                    np.unique(
                        metadata_df_curr[cov_curr].values.tolist()
                        )
                    ) == 1:
                nonvarying_covs_across_folds.append(cov_curr)
        eq_pairs_across_folds.extend(eq_pairs_curr)
        # 07d.v. Check the test sample sizes
        tc_ss_curr = {}
        for tcl in poss_tc_levels:
            ss_curr = len(
                    metadata_df_curr.loc[
                        metadata_df_curr[test_condition] == tcl
                        ]
                    )
            tc_ss_curr[tcl] = ss_curr
        tc_ss_map_dict[fold_removed_name] = tc_ss_curr
    # 07d.vi. Get unique equivalent covariate pairs
    # across folds
    eq_pairs_across_folds = list(
            np.unique(
                eq_pairs_across_folds
                )
            )
    # 07d.vii. Get the unique nonvarying covariates
    # across folds
    nonvarying_covs_across_folds = list(
            np.unique(
                nonvarying_covs_across_folds
                )
            )

    # 08. Print equivalent covariates, nonvarying covariates,
    # the sizes of each fold-removal sample's test levels, and
    # move to the next partition.
    eq_cov_str = '\n\t'.join(eq_pairs_across_folds)
    nv_cov_str = '\n\t'.join(nonvarying_covs_across_folds)
    print(f'Design formula information for group {g}:')
    print(f'\tPairs of equivalent covariates:\n{eq_cov_str}')
    print(f'\tNon-varying covariates:\n{nv_cov_str}')
    print(f'\tTest level sizes per fold-removed sample:\n\t{tc_ss_map_dict}')
    
    # 09. Build the design function for the DESeq objects
    covariates_to_use = [
            _ for _ in covariates_to_process
            if _ not in nonvarying_covs_across_folds]
    formula,n_terms,perfect_tc_cov = build_design_formula(
            covariates_to_use,
            test_condition,
            eq_pairs_across_folds
            )
    print(f'\tDesign formula: {formula}')

    # 10. Check whether there are enough
    # barcodes to power each DEA test
    N_BARCODES_MIN_DESEQ2 = np.max(
            [
                n_terms,
                MIN_SAMPLE_SIZE
                ]
            )
    print(f'\tNumber of barcodes per arm needed to power DEA: {N_BARCODES_MIN_DESEQ2}.')

    tc_ss_map_pf_dict = {} 
    powered_conditions_per_fold_rem_sample = {}
    powered_conditions_across_folds = []
    powered_condition_levels_across_folds = []
    for foldrem,fr_tc_dict in tc_ss_map_dict.items():
        print(f'\t\t{foldrem}:')
        tc_ss_map_pf_curr = {}
        for frtcl,frtcl_ss in fr_tc_dict.items():
            pf_frtcl_ss = ('PASS' if frtcl_ss > N_BARCODES_MIN_DESEQ2 else 'FAIL')
            print(f'\t\t\t{frtcl}: {pf_frtcl_ss}')
            tc_ss_map_pf_curr[frtcl] = pf_frtcl_ss
        tc_ss_map_pf_dict[foldrem] = tc_ss_map_pf_curr
        powered_cond_curr = []
        for cond_curr in condition_comparison_list:
            if (
                    (tc_ss_map_pf_curr[cond_curr[1]] == 'PASS') 
                    & 
                    (tc_ss_map_pf_curr[cond_curr[2]] == 'PASS')
                    ):
                powered_cond_curr.append(cond_curr)
        powered_conditions_per_fold_rem_sample[foldrem] = powered_cond_curr
        powered_conditions_across_folds.extend(powered_cond_curr)

    print(f'\tPowered comparisons per fold: {powered_conditions_per_fold_rem_sample}')
    print(f'\n')

    # 11. Get the list of comparisons that are powered across all folds
    join_str = '__'
    powered_conditions_across_folds_str = [
            join_str.join(_)
            for _ in
            powered_conditions_across_folds
            ]
    powered_conditions_u, powered_conditions_count = np.unique(
            powered_conditions_across_folds_str,
            return_counts=True
            )
    powered_conditions_all_folds = [
            _ for i,_ in enumerate(powered_conditions_u)
            if powered_conditions_count[i] == N_FOLDS_KFCV
            ]
    pcaf_str = '\n'.join(powered_conditions_all_folds)
    print(f'\tConditions powered for analysis across all folds: ' + \
            f'{pcaf_str}')
    for cond_curr in condition_comparison_list:
        check_str = join_str.join(cond_curr)
        if check_str in powered_conditions_all_folds:
            powered_condition_levels_across_folds.extend(cond_curr[1:])
    powered_condition_levels_across_folds = list(
            set(
                powered_condition_levels_across_folds
                )
            )

    # 12. If the partition is powered for analysis, reformat the 
    # count and metadata matrices, and save these, the barcode map,
    # and the powered condition list to a file
    if len(powered_conditions_all_folds) > 0:
        # 12a. Add the current partition to the
        # list of those validated for DEA
        diffexpr_ready_groups.append(g)

        # 12b. Set up an output directory for diffexpr-ready
        # count and metadata files
        out_dir_curr = dea_diffexpr_out_dir + f'/{g}'
        if not os.path.isdir(out_dir_curr):
            os.system(f'mkdir -p {out_dir_curr}')

        # 12c. Split the counts and metadata by fold
        # and save each as a separate pair of files
        #for fold_removed in np.arange(1,N_FOLDS_KFCV+1):
        for fold_removed_name,fold_removed in folds_removed_per_partition_dict.items():
            # 12c.i. Set up count and metadata output
            # file names for this fold-removed case
            mdf_curr_fn = f'{out_dir_curr}/metadata_std_{g}{dea_diffexpr_out_fn_tag}__{fold_removed_name}.csv'
            cdf_curr_fn = f'{out_dir_curr}/counts_1p_{g}{dea_diffexpr_out_fn_tag}__{fold_removed_name}.csv'
            # 12c.ii. If the 'new_writes_only' flag is set,
            # only proceed if both output files do not exist
            to_proceed = False
            if new_writes_only:
                if not (
                        (os.path.isfile(mdf_curr_fn))
                        &
                        (os.path.isfile(cdf_curr_fn))
                        ):
                    to_proceed = True
            else:
                to_proceed = True
            if to_proceed:
                bc_idxs_curr = barcode_fold_df.loc[
                        ~barcode_fold_df['fold_ID'].isin(fold_removed)
                        ].index.values.tolist()
                # 12c.iii. Get the metadata DataFrame for this fold
                metadata_df_curr = metadata_df.loc[bc_idxs_curr].copy()
                # 12c.iv. Standardize numeric covariates in this subset
                # of the metadata DataFrame
                num_cols = [
                        _ for _ in metadata_df_curr.columns.values.tolist()
                        if _ in numeric_covariates]
                for col_curr in num_cols:
                    metadata_df_curr[col_curr] = standardize_data(
                            metadata_df_curr[col_curr]
                                )
                # 12c.v. Save the metadata DataFrame
                metadata_df_curr.to_csv(
                        mdf_curr_fn,
                        index=True
                        )
                # 12c.vi. Delete the metadata DataFrame to save space
                del metadata_df_curr
                # 12c.vii. Get the counts DataFrame for this fold removal
                count_cols_curr = ['gene_name'] + bc_idxs_curr
                counts_df_curr = counts_df[count_cols_curr].copy()
                # 12c.viii. Add a pseudocount of 1 to all data if
                # specified (to avoid issues with zero counts)
                if ADD_PSEUDOCOUNT==True:
                    for col_curr in bc_idxs_curr:
                        counts_df_curr[col_curr] = counts_df_curr[col_curr] + 1
                # 12c.ix. Save current counts DataFrame to file
                counts_df_curr.to_csv(
                        cdf_curr_fn,
                        index=True
                        )
                # 12c.vii. Delete the counts DataFrame to save space
                del cdf_curr_fn

        # 12d. Write the design formula for this group to file
        design_formula_fn = f'{out_dir_curr}/design_formula{dea_diffexpr_out_fn_tag}_all_folds.txt'
        with open(design_formula_fn,'w') as f:
            f.write(f'{formula}\n')
        
        # 12e. Write the list of powered conditions for this group to file
        powered_cond_fn = f'{out_dir_curr}/powered_conditions{dea_diffexpr_out_fn_tag}_all_folds.txt'
        with open(powered_cond_fn,'w') as f:
            for pc in powered_condition_levels_across_folds:
                f.write(f'{pc}\n')

# 13. Save the current state of the random
# number generator for further processing
numpy_random_state_dict = numpy_random_state.get_state(
        legacy=False
        )
numpy_random_state_full_output_fn = f'{output_rand_state_dir}/{rand_state_output_fn}'
with open(numpy_random_state_full_output_fn,'wb') as outfile:
    pickle.dump(numpy_random_state_dict,
            outfile)

# 14. If all possible group values were processed with this script (i.e. if
# this default behavior was not overridden by a single grouping value),
# save the list of validated groups to file (for use by a parallelization
# shell script)
if not grouping_value_override:
    dea_group_list_out_full_fn = dea_diffexpr_out_dir + f'/diffexpr_ready_group_list{dea_diffexpr_out_fn_tag}.txt'
    with open(dea_group_list_out_full_fn,'w') as f:
        for drg in diffexpr_ready_groups:
            f.write(f'{drg}\n')

# 15. Save the number of folds used in DEA setup to file (for use by a
# parallelization shell script)
# 15a. Write the default number of folds used
nfolds_full_fn = f'{nfolds_dir}/{nfolds_out_fn}'
with open(nfolds_full_fn,'w') as f:
    f.write(f'{N_FOLDS_KFCV}\n')
# 15b. Write files for any exceptions used
for group,info in kfcv_exceptions_dict.items():
    nfolds_fn_pieces = nfolds_out_fn.split('.')
    nfolds_fn_exception_curr = f'{nfolds_fn_pieces[0]}__exception_{group}.{nfolds_fn_pieces[1]}'
    nfolds_full_fn_exception_curr = f'{nfolds_dir}/{nfolds_fn_exception_curr}'
    with open(nfolds_full_fn_exception_curr,'w') as f:
        n_folds_to_use_curr = info['n_folds_total']
        n_partitions_to_use_curr = info['n_partitions']
        n_folds_removed_per_partition = info['n_folds_removed_per_partition']
        f.write(f'Total number of folds generated: {n_folds_to_use_curr}\n')
        f.write(f'Number of partitions created from folds: {n_partitions_to_use_curr}\n')
        f.write(f'Number of folds removed per partition: {n_folds_removed_per_partition}\n')

sys.exit()


