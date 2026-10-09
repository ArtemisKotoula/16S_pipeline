#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
    echo "Usage: bash dada_pipeline.sh <input_directory> <output_directory> <threads> <right_len> <left_len> <dada_db> [sample_sheet]"
    exit 1
fi

# Input arguments
input_dir="$1"
dada2_outDir="$2"
threads="$3"
right_len="$4"
left_len="$5"
dada_db="$6"
sample_sheet="${7:-}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/conda.sh"

activate_dada

echo "DADA pipeline initialized."
echo "Input directory: ${input_dir}"

mkdir -p "${dada2_outDir}"

# Set the number of threads that Rscript can use.
last_cpu=$((threads - 1))

taskset -c 0-${last_cpu} Rscript "${script_dir}/dada2.R" "${input_dir}" "${dada2_outDir}" "${right_len}" "${left_len}" "${dada_db}" "${sample_sheet}"
echo "DADA analysis completed."
