#!/usr/bin/env bash

source "$(conda info --base)/etc/profile.d/conda.sh"


activate_dada() {
    set +u
    conda activate 16S_dada
    set -u
}


activate_cutadapt() {
    set +u
    conda activate cutadapt
    set -u
}


activate_kraken() {
    set +u
    conda activate kraken
    set -u
}

# Used like:
# source config/conda.sh

# activate_kraken
