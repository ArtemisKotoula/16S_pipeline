#!/usr/bin/env bash

# Trim primers from raw reads using cutadapt

set -euo pipefail

if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <input_dir> <output_dir> <forward_primer> <reverse_primer> <threads> <min_len>"
    exit 1
fi

raw_data="$1"
cutadapt_outDir="$2"
r1_primer="$3"
r2_primer="$4"
threads="$5"
min_len="$6"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/conda.sh"

activate_cutadapt

mkdir -p "${cutadapt_outDir}"

for sample_dir in "${raw_data}"/*; do

    [[ -d "${sample_dir}" ]] || continue

    sample_name=$(basename "${sample_dir}")

    echo "Processing ${sample_name}"

    # R1/R2 must be a separate token in the file name (e.g. S1_R1.fastq.gz, S1_S1_L001_R1_001.fastq.gz),
    # so that sample names containing "R1"/"R2" (e.g. CTR1) are not matched
    R1_files=()
    R2_files=()
    for fq in "${sample_dir}"/*.fastq "${sample_dir}"/*.fastq.gz; do
        [[ -f "${fq}" ]] || continue
        fq_name=$(basename "${fq}")
        if [[ "${fq_name}" =~ [._-]R1([._-].*)?\.fastq(\.gz)?$ ]]; then
            R1_files+=("${fq}")
        elif [[ "${fq_name}" =~ [._-]R2([._-].*)?\.fastq(\.gz)?$ ]]; then
            R2_files+=("${fq}")
        fi
    done

    if [[ ${#R1_files[@]} -ne 1 || ${#R2_files[@]} -ne 1 ]]; then
        echo "WARNING: Expected exactly one R1 and one R2 FASTQ file for ${sample_name}, found ${#R1_files[@]} R1 and ${#R2_files[@]} R2. Skipping sample."
        continue
    fi

    R1="${R1_files[0]}"
    R2="${R2_files[0]}"

    mkdir -p "${cutadapt_outDir}/${sample_name}"

    cutadapt \
        -j "${threads}" \
        -g "${r1_primer}" \
        -G "${r2_primer}" \
        -m "${min_len}" \
        -o "${cutadapt_outDir}/${sample_name}/${sample_name}_R1_trimmed.fastq" \
        -p "${cutadapt_outDir}/${sample_name}/${sample_name}_R2_trimmed.fastq" \
        "${R1}" "${R2}"

done
