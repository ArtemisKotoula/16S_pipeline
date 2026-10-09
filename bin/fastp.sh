#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "Usage: $0 -mode dada -in input_dir -out output_dir -t threads -min_len min_length -lenR right_len -lenL left_len"
    echo "       $0 -mode kraken -in input_dir -out output_dir -t threads -min_len min_length -qt quality_threshold"
    exit 1
}

mode=""
input_dir=""
output_dir=""
threads=""
min_length=""
right_len=""
left_len=""
quality_threshold=""

echo "Checking arguments and starting fastp quality filtering..."

while [[ $# -gt 0 ]]; do
    case "$1" in
        -mode) mode="$2"; shift 2 ;;
        -in) input_dir="$2"; shift 2 ;;
        -out) output_dir="$2"; shift 2 ;;
        -t) threads="$2"; shift 2 ;;
        -min_len) min_length="$2"; shift 2 ;;
        -lenR) right_len="$2"; shift 2 ;;
        -lenL) left_len="$2"; shift 2 ;;
        -qt) quality_threshold="$2"; shift 2 ;;
        *) usage ;;
    esac
done

if [[ -z "${input_dir}" || -z "${output_dir}" || -z "${threads}" || -z "${min_length}" ]]; then
    usage
fi

# Mode-specific fastp options
case "${mode}" in
    "dada")
        [[ -n "${right_len}" && -n "${left_len}" ]] || usage
        # Trimming at fixed lengths for DADA2
        fastp_opts=(--max_len1 "${right_len}" --max_len2 "${left_len}")
        ;;
    "kraken")
        [[ -n "${quality_threshold}" ]] || usage
        # Quality trimming for Kraken
        fastp_opts=(-r --cut_right_mean_quality "${quality_threshold}")
        ;;
    *) usage ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/conda.sh"

activate_dada

mkdir -p "${output_dir}"

for sample_dir in "${input_dir}"/*; do
    [[ -d "${sample_dir}" ]] || continue

    sample_name=$(basename "${sample_dir}")

    echo "Processing ${sample_name}"

    mkdir -p "${output_dir}/${sample_name}"

    R1="${sample_dir}/${sample_name}_R1_trimmed.fastq"
    R2="${sample_dir}/${sample_name}_R2_trimmed.fastq"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "ERROR: Missing trimmed FASTQ files for ${sample_name}"; exit 1
    fi

    fastp -i "${R1}" -I "${R2}" --thread "${threads}" \
        -o "${output_dir}/${sample_name}/${sample_name}_R1_filtered.fastq" \
        -O "${output_dir}/${sample_name}/${sample_name}_R2_filtered.fastq" \
        -l "${min_length}" "${fastp_opts[@]}" \
        --html "${output_dir}/${sample_name}/${sample_name}_fastp.html" \
        --json "${output_dir}/${sample_name}/${sample_name}_fastp.json"
done

echo "fastp filtering completed."
