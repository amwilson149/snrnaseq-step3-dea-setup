
# 01. Set up environment
# 01a. Activate conda environment
ml anaconda3/2020.11
module  load 'anaconda3/2020.11'
Shell debugging temporarily silenced: export LMOD_SH_DBG_ON=1 for Lmod's output

The following have been reloaded with a version change:
  1) gcc/14.2.0 => gcc/8.3.0

Shell debugging restarted
ml -python
module  load '-python'
Shell debugging temporarily silenced: export LMOD_SH_DBG_ON=1 for Lmod's output
Shell debugging restarted
source /hpc/packages/minerva-centos7/anaconda3/2020.11/etc/profile.d/conda.sh
export CONDA_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/conda'
export _CE_M=''
export _CE_CONDA=''
export CONDA_PYTHON_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/python'

# Copyright (C) 2012 Anaconda, Inc
# SPDX-License-Identifier: BSD-3-Clause

__add_sys_prefix_to_path() {
    # In dev-mode CONDA_EXE is python.exe and on Windows
    # it is in a different relative location to condabin.
    if [ -n "${_CE_CONDA}" ] && [ -n "${WINDIR+x}" ]; then
        SYSP=$(\dirname "${CONDA_EXE}")
    else
        SYSP=$(\dirname "${CONDA_EXE}")
        SYSP=$(\dirname "${SYSP}")
    fi

    if [ -n "${WINDIR+x}" ]; then
        PATH="${SYSP}/bin:${PATH}"
        PATH="${SYSP}/Scripts:${PATH}"
        PATH="${SYSP}/Library/bin:${PATH}"
        PATH="${SYSP}/Library/usr/bin:${PATH}"
        PATH="${SYSP}/Library/mingw-w64/bin:${PATH}"
        PATH="${SYSP}:${PATH}"
    else
        PATH="${SYSP}/bin:${PATH}"
    fi
    \export PATH
}

__conda_hashr() {
    if [ -n "${ZSH_VERSION:+x}" ]; then
        \rehash
    elif [ -n "${POSH_VERSION:+x}" ]; then
        :  # pass
    else
        \hash -r
    fi
}

__conda_activate() {
    if [ -n "${CONDA_PS1_BACKUP:+x}" ]; then
        # Handle transition from shell activated with conda <= 4.3 to a subsequent activation
        # after conda updated to >= 4.4. See issue #6173.
        PS1="$CONDA_PS1_BACKUP"
        \unset CONDA_PS1_BACKUP
    fi

    \local cmd="$1"
    shift
    \local ask_conda
    CONDA_INTERNAL_OLDPATH="${PATH}"
    __add_sys_prefix_to_path
    ask_conda="$(PS1="$PS1" "$CONDA_EXE" $_CE_M $_CE_CONDA shell.posix "$cmd" "$@")" || \return $?
    rc=$?
    PATH="${CONDA_INTERNAL_OLDPATH}"
    \eval "$ask_conda"
    if [ $rc != 0 ]; then
        \export PATH
    fi
    __conda_hashr
}

__conda_reactivate() {
    \local ask_conda
    CONDA_INTERNAL_OLDPATH="${PATH}"
    __add_sys_prefix_to_path
    ask_conda="$(PS1="$PS1" "$CONDA_EXE" $_CE_M $_CE_CONDA shell.posix reactivate)" || \return $?
    PATH="${CONDA_INTERNAL_OLDPATH}"
    \eval "$ask_conda"
    __conda_hashr
}

conda() {
    if [ "$#" -lt 1 ]; then
        "$CONDA_EXE" $_CE_M $_CE_CONDA
    else
        \local cmd="$1"
        shift
        case "$cmd" in
            activate|deactivate)
                __conda_activate "$cmd" "$@"
                ;;
            install|update|upgrade|remove|uninstall)
                CONDA_INTERNAL_OLDPATH="${PATH}"
                __add_sys_prefix_to_path
                "$CONDA_EXE" $_CE_M $_CE_CONDA "$cmd" "$@"
                \local t1=$?
                PATH="${CONDA_INTERNAL_OLDPATH}"
                if [ $t1 = 0 ]; then
                    __conda_reactivate
                else
                    return $t1
                fi
                ;;
            *)
                CONDA_INTERNAL_OLDPATH="${PATH}"
                __add_sys_prefix_to_path
                "$CONDA_EXE" $_CE_M $_CE_CONDA "$cmd" "$@"
                \local t1=$?
                PATH="${CONDA_INTERNAL_OLDPATH}"
                return $t1
                ;;
        esac
    fi
}

if [ -z "${CONDA_SHLVL+x}" ]; then
    \export CONDA_SHLVL=0
    # In dev-mode CONDA_EXE is python.exe and on Windows
    # it is in a different relative location to condabin.
    if [ -n "${_CE_CONDA+x}" ] && [ -n "${WINDIR+x}" ]; then
        PATH="$(\dirname "$CONDA_EXE")/condabin${PATH:+":${PATH}"}"
    else
        PATH="$(\dirname "$(\dirname "$CONDA_EXE")")/condabin${PATH:+":${PATH}"}"
    fi
    \export PATH

    # We're not allowing PS1 to be unbound. It must at least be set.
    # However, we're not exporting it, which can cause problems when starting a second shell
    # via a first shell (i.e. starting zsh from bash).
    if [ -z "${PS1+x}" ]; then
        PS1=
    fi
fi
conda activate CO-deg-setup-env
PS1='(CO-deg-setup-env) '
export PATH='/sc/arion/work/wilsoa28/.conda/envs/CO-deg-setup-env/bin:/hpc/packages/minerva-centos7/anaconda3/2020.11/bin:/hpc/packages/minerva-centos7/gcc/8.3.0_32b/bin:/hpc/packages/minerva-centos7/anaconda3/2020.11/condabin:/hpc/users/wilsoa28/google-cloud-sdk/bin:/hpc/users/wilsoa28/git-filter-repo:/hpc/packages/minerva-rocky9/git/2.46.0/bin:/hpc/packages/minerva-common/vim/8.0/bin:/hpc/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/hpc/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/bin:/usr/bin:/usr/mbin:/local/bin:/usr/local:/usr/ucb:/usr/local/sbin:/usr/sbin:/usr/lpp/mmfs/bin:/hpc/users/wilsoa28/.local/bin:/hpc/users/wilsoa28/bin'
export CONDA_PREFIX='/sc/arion/work/wilsoa28/.conda/envs/CO-deg-setup-env'
export CONDA_SHLVL='1'
export CONDA_DEFAULT_ENV='CO-deg-setup-env'
export CONDA_PROMPT_MODIFIER='(CO-deg-setup-env) '
export CONDA_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/conda'
export _CE_M=''
export _CE_CONDA=''
export CONDA_PYTHON_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/python'

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



