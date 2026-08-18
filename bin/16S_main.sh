#!/usr/bin/env bash

set -euo pipefail

mode=${1:-}
rerun=${2:-}

if [[ -z "$mode" ]]; then
    echo "Usage: bash 16S_main.sh [dada|kraken] [--rerun]"
    exit 1
fi

if [[ "$rerun" == "--rerun" ]]; then
    echo "Rerun enabled."
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${script_dir}/config/config.sh"
source "${script_dir}/config/conda.sh"

# LOG FUNCTION
log(){
    echo "[$(date '+%F %T')] --- $1"
}

case "$rerun" in
    "--rerun")
        log "Rerun enabled. Searching for existing output directories."
        prev_dir=$(find "${results_dir}" -maxdepth 1 -type d | sort | tail -1)
        echo "Previous output directory: ${prev_dir}"
        if find "${prev_dir}/02_01_cutadapt" -mindepth 1 -print -quit | grep -q .; then
            echo "Found existing cutadapt output. Will skip preprocessing."
            cutadapt_outDir="${prev_dir}/02_01_cutadapt"
            skip_pre=true
        else
            echo "No existing cutadapt output found. Will start from raw data."
            skip_pre=false
        fi
        ;;
    *)

esac


log "Starting 16S pipeline in $mode mode"

# Create Results Direcotry and Logfile
mkdir -p "${out_dir}"
logfile="${out_dir}/pipeline.log"

exec > >(tee -a "$logfile") 2>&1

# Common Preprocessing Steps

case "$skip_pre" in
    true)
        log "Skipping preprocessing steps. Starting with fastp from previous cutadapt results in ${cutadapt_outDir}."
        ;;
    false)
        log "Starting preprocessing steps."

        log "Starting QC on raw reads"
        source "${script_dir}/parallel_qc.sh" "${raw_data}" "${fastqc_outDir}/raw" "${multiqc_outDir}/raw" "${threads}"
        log "Qc for raw reads completed. Results are in ${fastqc_outDir}/raw and ${multiqc_outDir}/raw"

        log "Starting trimming of primers from raw reads"
        source "${script_dir}/trim.sh" "${raw_data}" "${cutadapt_outDir}" "${r1_primer}" "${r2_primer}" "${threads}" 100
        log "Trimming completed. Results are in ${cutadapt_outDir}"

        log "Starting QC on trimmed reads"
        source "${script_dir}/parallel_qc.sh" "${cutadapt_outDir}" "${fastqc_outDir}/trimmed" "${multiqc_outDir}/trimmed" "${threads}"
        log "Qc for trimmed reads completed. Results are in ${fastqc_outDir}/trimmed and ${multiqc_outDir}/trimmed" 
        ;;
esac

# Choose analysis pipeline
case "$mode" in

    dada)
        log "Starting DADA2 analysis pipeline"

        source "${script_dir}/dada_pipeline.sh" "${cutadapt_outDir}"

        log "DADA2 analysis completed. Results are in ${dada2_outDir}"

        ;;

    kraken)
        # Last preprocess strep for kraken: quality filtering with fastp
        log "Starting quality filtering of trimmed reads for Kraken"
        source "${script_dir}/fastp.sh" "${cutadapt_outDir}" "${fastp_outDir}" "${threads}" "${min_length}" "${quality_threshold}"   
        log "Quality filtering completed. Results are in ${fastp_outDir}"

        log "Generating fastp summary report"
        source "${script_dir}/report_fastp.sh" "${fastp_outDir}" "${report_fastp_outDir}"
        log "fastp summary report generated. Results are in ${report_fastp_outDir}. Samples that do not meet the cutoff criteria have been moved to ${calc_cutoff_outDir}"

        log "Starting QC on filtered reads"
        source "${script_dir}/parallel_qc.sh" "${fastp_outDir}" "${fastqc_outDir}/filtered" "${multiqc_outDir}/filtered" "${threads}"
        log "Qc for filtered reads completed. Results are in ${fastqc_outDir}/filtered and ${multiqc_outDir}/filtered"

        log "Starting Kraken-Krona-Bracken analysis"
        source "${script_dir}/kraken_pipeline.sh" "${fastp_outDir}"
        log "Kraken-Krona-Bracken analysis and visualization completed. Results are in ${kraken_outDir} and ${phyloseq_outDir}"
        ;;

    *)
        echo "Unknown pipeline: ${mode}. Exiting."
        exit 1
        ;;

esac


log "Pipeline completed successfully."
