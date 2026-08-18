#!/usr/bin/env bash

threads=4

# Input raw data directory
raw_data="./data"

# Output directory
results_dir="./16S_results"
out_dir="${results_dir}/16S_$(date '+%Y%m%d_%H%M%S')"

# Primers
r1_primer="CCTACGGGNGGCWGCAG"
r2_primer="GACTACHVGGGTATCTAATCC"

#Read cutoffs for sample filtering 
frequency=0.01 
confidence=0.99
min_reads=10 #min numbers of reads that a taxa must appear in (Default for bracken)

# Kraken database
kraken_db="./dbs/kraken2_silva"

# Results directories for each step
fastqc_outDir="${out_dir}/01_01_fastqc"
multiqc_outDir="${out_dir}/01_02_multiqc"
cutadapt_outDir="${out_dir}/02_01_cutadapt"

#DADA PIPELINE
dada2_outDir="${out_dir}/03_dada2"
# Legths to truncate reads to for DADA2 filtering. 
right_len=260
left_len=220
# DADA database
dada_db="./dbs/dada2_silva_db/silva_nr_v138_train_set.fa.gz"

fastp_outDir="${out_dir}/02_02_fastp"
quality_threshold=20
min_length=100
report_fastp_outDir="${out_dir}/02_03_report_fastp"

calc_cutoff_outDir="${out_dir}/02_04_failed_cutoff_samples"

kraken_outDir="${out_dir}/03_kraken"
krona_outDir="${kraken_outDir}/03_01_krona"
bracken_outDir="${kraken_outDir}/03_02_bracken"

phyloseq_outDir="${out_dir}/04_phyloseq"

# a variable to control whether to skip preprocessing steps (cutadapt and fastp) if previous outputs exist
# false by default, but will be set to true if previous outputs are found and if rerun is specified
skip_pre=false
