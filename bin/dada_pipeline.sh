#!/user/bin/env bash

set -euo pipefail

# Input argument
input_dir="$1"

if [[ -z "$input_dir" ]]; then
    echo "Usage: bash dada2_pipeline.sh <input_directory>"
    exit 1
fi

activate_dada

echo "DADA pipeline initialized."
echo "Input directory: ${input_dir}"

mkdir -p "${dada2_outDir}"

# Set the number of threads that Rscript can use.
last_cpu=$((threads - 1))

taskset -c 0-${last_cpu} Rscript "${script_dir}/dada2.R" "${input_dir}" "${dada2_outDir}" "${right_len}" "${left_len}" "${dada_db}"
echo "DADA analysis completed."
