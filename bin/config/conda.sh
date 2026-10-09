#!/usr/bin/env bash

# Sourced by 16S_main.sh and by every step script.
# Only defines functions; environments are created once by setup_conda_envs (called from 16S_main.sh).

# CONDA_EXE is exported by conda's shell hook, so step scripts started as subprocesses can still find conda
conda_base="$("${CONDA_EXE:-conda}" info --base)"
source "${conda_base}/etc/profile.d/conda.sh"

conda_config_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

setup_conda_envs() {
    echo "Creating conda environments."

    local env_names=("16Sdada" "16Scutadapt" "16Skraken")

    conda config --set channel_priority flexible

    for env_name in "${env_names[@]}"; do
        if conda env list | grep -q "$env_name"; then
            echo "Conda environment '$env_name' already exists. Skipping creation."
        else
            echo "Creating conda environment '$env_name'..."
            conda env create -f "${conda_config_dir}/${env_name}_env.yml"
        fi
    done
}

activate_dada() {
    set +u
    conda activate 16Sdada
    set -u
}

activate_cutadapt() {
    set +u
    conda activate 16Scutadapt
    set -u
}

activate_kraken() {
    set +u
    conda activate 16Skraken
    set -u
}
