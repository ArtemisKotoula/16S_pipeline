#!/usr/bin/env bash

# Trim primers from raw reads using cutadapt

set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <input_dir> <output_dir> <forward_primer> <reverse_primer> <threads> <min_length>"
    exit 1
fi

raw_data="$1"
cutadapt_outDir="$2"
r1_primer="$3"
r2_primer="$4"
threads="$5"
min_length="$6"

activate_cutadapt

mkdir -p "${cutadapt_outDir}"

for sample_dir in "${raw_data}"/*; do

    [[ -d "${sample_dir}" ]] || continue

    sample_name=$(basename "${sample_dir}")

    echo "Processing ${sample_name}"

    mkdir -p "${cutadapt_outDir}/${sample_name}"

    R1=$(find "${sample_dir}" -maxdepth 1 -type f \( -name "*R1*.fastq.gz" -o -name "*R1*.fastq" \))
    R2=$(find "${sample_dir}" -maxdepth 1 -type f \( -name "*R2*.fastq.gz" -o -name "*R2*.fastq" \))

    if [[ -z "${R1}" || -z "${R2}" ]]; then
        echo "Missing FASTQ files for ${sample_name}"
        continue
    fi

    cutadapt \
        -j "${threads}" \
        -g "${r1_primer}" \
        -G "${r2_primer}" \
        -m "${min_length}" \
        -o "${cutadapt_outDir}/${sample_name}/${sample_name}_R1_trimmed.fastq" \
        -p "${cutadapt_outDir}/${sample_name}/${sample_name}_R2_trimmed.fastq" \
        "${R1}" "${R2}"

done
