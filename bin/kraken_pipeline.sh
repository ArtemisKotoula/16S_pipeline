#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 8 || $# -gt 9 ]]; then
    echo "Usage: bash kraken_pipeline.sh <input_directory> <kraken_outDir> <krona_outDir> <bracken_outDir> <phyloseq_outDir> <kraken_db> <threads> <read_length> [sample_sheet]"
    exit 1
fi

# Input arguments
input_dir="$1"
kraken_outDir="$2"
krona_outDir="$3"
bracken_outDir="$4"
phyloseq_outDir="$5"
kraken_db="$6"
threads="$7"
avg_mean_l="$8"
sample_sheet="${9:-}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/conda.sh"

activate_kraken

echo "Kraken pipeline initialized."
echo "Input directory: ${input_dir}"

kraken_res="${kraken_outDir}/kraken_results"

mkdir -p "${kraken_res}"

for sample_dir in "${input_dir}"/*; do
    sample_name=$(basename "${sample_dir}")

    echo "Processing sample: ${sample_name}"

    R1="${input_dir}/${sample_name}/${sample_name}_R1_filtered.fastq"
    R2="${input_dir}/${sample_name}/${sample_name}_R2_filtered.fastq"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "Kraken input files for sample ${sample_name} not found. Instead R1 file is set as: ${R1}. Skipping."
        continue
    fi

    kraken2 --db "${kraken_db}" --paired --threads "${threads}" \
        --report "${kraken_res}/${sample_name}_kraken_report.txt" \
        --output "${kraken_res}/${sample_name}_kraken_output.txt" \
        "${R1}" "${R2}"
done

echo "Kraken analysis completed. Results are in ${kraken_res}"


echo "Starting Krona visualization of Kraken results..."

mkdir -p "${krona_outDir}"

echo "Updating krona taxonomy database..."
# Find Krona's updateTaxonomy.sh
krona_update=$(find "$(dirname "$(which ktImportTaxonomy)")"/.. -name updateTaxonomy.sh -type f -print -quit)

"$krona_update"

# One "file,sample_label" argument per sample.
# The label is the part of the sample name before the first "_" (as in the phyloseq plots),
# unless that would give two samples the same label
kr_files=()
kr_labels=()
for f in "${kraken_res}/"*_kraken_output.txt; do
    sample=$(basename "${f}" "_kraken_output.txt")
    kr_files+=("${f}")
    kr_labels+=("${sample%%_*}")
done

if [[ -n "$(printf '%s\n' "${kr_labels[@]}" | sort | uniq -d)" ]]; then
    echo "Shortened sample names are not unique. Using the full sample names in Krona."
    for i in "${!kr_files[@]}"; do
        kr_labels[i]=$(basename "${kr_files[i]}" "_kraken_output.txt")
    done
fi

kr_report_files=()
for i in "${!kr_files[@]}"; do
    kr_report_files+=("${kr_files[i]},${kr_labels[i]}")
done

ktImportTaxonomy -q 2 -t 3 "${kr_report_files[@]}" -o "${krona_outDir}/krona_all_samples.html"
# -q 2 = read ID is in column 2
#-t 3 = taxonomy ID is in column 3

echo "Krona visualization completed. Results are in ${krona_outDir}/krona_all_samples.html"

###########################################################s

echo "Buiding Bracken database for Kraken results..."
bracken-build -d ${kraken_db} -t ${threads} -l ${avg_mean_l}
echo "Bracken database built in ${kraken_db}"

echo "Starting Bracken analysis on KRAKEN results"

mkdir -p "${bracken_outDir}"

for report_file in "${kraken_res}/"*_kraken_report.txt; do
    sample_name=$(basename "${report_file}" "_kraken_report.txt")
    bracken_output="${bracken_outDir}/${sample_name}.bracken"
    bracken -d "${kraken_db}" -i "${report_file}" -o "${bracken_output}" -r "${avg_mean_l}" -l G
done

# get the sample names from the bracken output files
sample_names=()
for bracken_file in "${bracken_outDir}"/*.bracken; do
    sample_name=$(basename "${bracken_file}" ".bracken")
    sample_names+=("${sample_name}")
done

echo "${sample_names[@]}"
# Join array elements with commas
sample_names=$(IFS=,; echo "${sample_names[*]}")

echo "${sample_names}"

#combine all bracken results and generate genus abundance
combine_bracken_outputs.py --files "${bracken_outDir}"/*.bracken \
--names "${sample_names}" \
--output "${bracken_outDir}/genus_abundance.tsv"        
#modified script, adding "name = f"{name}-{taxid}" after line 106 to avoid duplicate names in the combined output file
# (the stock script exits when the same name has different taxonomy IDs, e.g. "uncultured" in SILVA)

echo "Bracken analysis completed. Results are in ${bracken_outDir}"

############################################################
echo "Starting visualization of Bracken results with phyloseq"

activate_dada

Rscript "${script_dir}/phyloseq.R" "${bracken_outDir}/genus_abundance.tsv" "${phyloseq_outDir}" "${sample_sheet}"

echo "Phyloseq visualization completed. Results are in ${phyloseq_outDir}"

