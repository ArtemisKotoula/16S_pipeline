#!/usr/bin/env bash

####################################3
# fastp json processing script
# create a summary report for all samples
# also calculate the average mean length
#####################################

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi


fastIn_dir="$1"
reportOut_dir="$2"

activate_dada

mkdir -p "${reportOut_dir}"

reportOut_file="${reportOut_dir}/fastp_reads_report.csv"
summaryOut_file="${reportOut_dir}/fastp_summary_report.txt"

###################################3

echo "Generating fastp summary report..."

echo "Sample_Name,Total_Reads,Filtered_Reads,R1_Mean_Length,R2_Mean_Length,Filtered_R1_Mean_Length,Filtered_R2_Mean_Length" > "${reportOut_file}"

for sample_dir in "${fastIn_dir}"/*; do
    [[ -d "${sample_dir}" ]] || continue

    sample_name=$(basename "${sample_dir}")

    echo "Processing ${sample_name}"

    fastp_json="${sample_dir}/${sample_name}_fastp.json"

    echo "fastp json file: ${fastp_json}"

    if [[ ! -f "${fastp_json}" ]]; then
        echo "ERROR: Missing fastp JSON file for ${sample_name}. Exiting."
        exit 1
    fi

    # Extract values from the JSON file using jq
    total_reads=$(jq '.summary.before_filtering.total_reads' "${fastp_json}")
    filtered_reads=$(jq '.summary.after_filtering.total_reads' "${fastp_json}")
    r1_mean_length=$(jq '.summary.before_filtering.read1_mean_length' "${fastp_json}")
    r2_mean_length=$(jq '.summary.before_filtering.read2_mean_length' "${fastp_json}")
    filtered_r1_mean_length=$(jq '.summary.after_filtering.read1_mean_length' "${fastp_json}")
    filtered_r2_mean_length=$(jq '.summary.after_filtering.read2_mean_length' "${fastp_json}")

    # Append to the summary report
    echo "${sample_name},${total_reads},${filtered_reads},${r1_mean_length},${r2_mean_length},${filtered_r1_mean_length},${filtered_r2_mean_length}" >> "${reportOut_file}"
done

# Calculate average mean lengths across all samples
avg_r1_mean_length=$(awk -F',' 'NR>1 {sum+=$4; count++} END {if (count > 0) print sum/count; else print 0}' "${reportOut_file}")
avg_r2_mean_length=$(awk -F',' 'NR>1 {sum+=$5; count++} END {if (count > 0) print sum/count; else print 0}' "${reportOut_file}")
avg_filtered_r1_mean_length=$(awk -F',' 'NR>1 {sum+=$6; count++} END {if (count > 0) print sum/count; else print 0}' "${reportOut_file}")
avg_filtered_r2_mean_length=$(awk -F',' 'NR>1 {sum+=$7; count++} END {if (count > 0) print sum/count; else print 0}' "${reportOut_file}")

# Calculate average total and filtered reads across both R1 and R2
avg_total_reads=$(awk -F',' 'NR>1 {sum+=$2; count++} END {if (count > 0) print sum/count; else print 0}' "${reportOut_file}")
avg_filtered_reads=$(awk -F',' 'NR>1 {sum+=$3; count++} END {if (count > 0) print sum/count; else print 0}' "${reportOut_file}")

avg_mean_l=$(awk -v r1="$avg_filtered_r1_mean_length" -v r2="$avg_filtered_r2_mean_length" 'BEGIN {printf "%.0f", (r1 + r2) / 2}')

# Print the mean lengths to the summary report
{
echo "Average Mean Lengths Across All Samples:"
echo "R1 Mean Length: ${avg_r1_mean_length}"
echo "R2 Mean Length: ${avg_r2_mean_length}"
echo "Filtered R1 Mean Length: ${avg_filtered_r1_mean_length}"
echo "Filtered R2 Mean Length: ${avg_filtered_r2_mean_length}"
echo "Average Total Reads: ${avg_total_reads}"
echo "Average Filtered Reads: ${avg_filtered_reads}"
echo -e "Average Filtered Mean Length R1+R2: ${avg_mean_l}"
} >> "${summaryOut_file}"

echo "fastp summary reports generated at ${reportOut_file} and ${summaryOut_file}"

echo "Calculating cutoff values for filtering"
mkdir -p "${calc_cutoff_outDir}"

# Parameters have been set in config.sh
min_reads_required=$(python "${script_dir}/calc_cutoff.py"\
    --frequency "${frequency}"\
    --confidence "${confidence}"\
    --min_taxon_reads "${min_reads}")

echo "Minimum reads required for a taxon to be considered present: ${min_reads_required}"
echo "Removing samples that do not meet the cutoff criteria and moving them to ${calc_cutoff_outDir}"

# Read the summary report and filter samples based on the cutoff
awk -F',' -v min_reads_rec="${min_reads_required}" -v In_dir="${fastIn_dir}" -v out_dir="${calc_cutoff_outDir}" \
  'NR>1 {
      if ($3 < min_reads_rec) {
          print $1;
          system("mv \"" In_dir "/" $1 "\" \"" out_dir "\"")
      }
  }' "${reportOut_file}"