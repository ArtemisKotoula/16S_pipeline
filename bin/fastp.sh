#!/usr/bin/env bash

set -euo pipefail

activate_dada
echo "Checking arguments and starting fastp quality filtering..."

case "$mode" in
    "dada")
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -in) input_dir="$2"; shift 2 ;;
                -out) output_dir="$2"; shift 2 ;;
                -t) threads="$2"; shift 2 ;;
                -min_len) min_length="$2"; shift 2 ;;
                -lenR) right_len="$2"; shift 2;;
                -lenL) left_len="$2"; shift 2 ;;
                *)
                    echo "Usage: $0 -in input_dir -out output_dir -t threads -min_len min_length -lenR right_len -lenL left_len"
                    exit 1
                    ;;
            esac
        done
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
                -l "${min_length}" --max_len1 "${right_len}" --max_len2 "${left_len}" \
                --html "${output_dir}/${sample_name}/${sample_name}_fastp.html" \
                --json "${output_dir}/${sample_name}/${sample_name}_fastp.json"
        done

    ;;
    "kraken")

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -in) input_dir="$2"; shift 2 ;;
            -out) output_dir="$2"; shift 2 ;;
            -t) threads="$2"; shift 2 ;;
            -min_len) min_length="$2"; shift 2 ;;
            -qt) quality_threshold="$2"; shift 2 ;;
            *)
                echo "Usage: $0 -in input_dir -out output_dir -t threads -min_len min_length -qt quality_threshold"
                exit 1
                ;;
        esac
    done

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
            -l "${min_length}" -r --cut_right_mean_quality "${quality_threshold}" \
            --html "${output_dir}/${sample_name}/${sample_name}_fastp.html" \
            --json "${output_dir}/${sample_name}/${sample_name}_fastp.json"
    done
    ;;
    *) exit 1 ;;
esac

echo "fastp filtering completed."