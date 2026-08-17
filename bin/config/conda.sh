#!/usr/bin/env bash

source "$(conda info --base)/etc/profile.d/conda.sh"

echo "Creating conda environments."

env_names=("16Sdada" "16Scutadapt" "16Skraken")

conda config --set channel_priority flexible

for env_name in "${env_names[@]}"; do
    if conda env list | grep -q "$env_name"; then
        echo "Conda environment '$env_name' already exists. Skipping creation."
    else
        echo "Creating conda environment '$env_name'..."
        case "$env_name" in
            "16Sdada")
                conda env create -f "${script_dir}/config/${env_name}_env.yml"
                ;;
            "16Scutadapt")
                conda env create -f "${script_dir}/config/${env_name}_env.yml"
                ;;
            "16Skraken")
                conda env create -f "${script_dir}/config/${env_name}_env.yml"
                ;;
        esac
    fi
done



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
