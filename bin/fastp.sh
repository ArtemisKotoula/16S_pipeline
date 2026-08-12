#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 <input_dir> <output_dir> <threads> <min_length> <quality_threshold>"
    exit 1
fi

input_dir="$1"
output_dir="$2"
threads="$3"
min_length="$4"
quality_threshold="$5"

activate_dada

mkdir -p "${output_dir}"

echo "Starting fastp quality filtering..."

for sample_dir in "${input_dir}"/*; do

    [[ -d "${sample_dir}" ]] || continue

    sample_name=$(basename "${sample_dir}")

    echo "Processing ${sample_name}"

    mkdir -p "${output_dir}/${sample_name}"

    R1="${sample_dir}/${sample_name}_R1_trimmed.fastq"
    R2="${sample_dir}/${sample_name}_R2_trimmed.fastq"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "ERROR: Missing trimmed FASTQ files for ${sample_name}"
        exit 1
    fi

    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${output_dir}/${sample_name}/${sample_name}_R1_filtered.fastq" \
        -O "${output_dir}/${sample_name}/${sample_name}_R2_filtered.fastq" \
        -l "${min_length}" \
        -r \
        --cut_right_mean_quality "${quality_threshold}" \
        --thread "${threads}" \
        --html "${output_dir}/${sample_name}/${sample_name}_fastp.html" \
        --json "${output_dir}/${sample_name}/${sample_name}_fastp.json"

done

echo "fastp filtering completed."