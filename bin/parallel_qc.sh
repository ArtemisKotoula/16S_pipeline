#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <input_dir> <FastQC_output_dir> <MultiQC_output_dir> <threads>"
    exit 1
fi

input_dir="$1"
FoutDir="$2"
MoutDir="$3"
threads="$4"

# Activate conda environment
activate_dada

mkdir -p "${FoutDir}" "${MoutDir}"

echo "Running FastQC in parallel..."

# Number of samples/jobs to run simultaneously
jobs=$(( threads / 2 ))
if [[ "${jobs}" -lt 1 ]]; then
    jobs=1
fi

# Run FastQC in parallel
find "${input_dir}" -maxdepth 2 -type f \
    \( -name "*.fastq.gz" -o -name "*.fastq" \) \
    -print0 |
    xargs -0 -n 1 -P "${jobs}" \
    fastqc -t 2 -o "${FoutDir}"

echo "FastQC completed."

echo "Running MultiQC..."

multiqc \
    "${FoutDir}" \
    -o "${MoutDir}"

echo "QC completed."

