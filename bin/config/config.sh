#!/usr/bin/env bash

threads=4

# Input raw data directory
raw_data="./data"

# Output directory 
out_dir="./16S_results/16S_$(date '+%Y%m%d_%H%M%S')"

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

dada2_outDir="${out_dir}/03_dada2"

fastp_outDir="${out_dir}/02_02_fastp"
quality_threshold=20
min_length=100
report_fastp_outDir="${out_dir}/02_03_report_fastp"

calc_cutoff_outDir="${out_dir}/02_04_failed_cutoff_samples"

kraken_outDir="${out_dir}/03_kraken"
krona_outDir="${kraken_outDir}/03_01_krona"
bracken_outDir="${kraken_outDir}/03_02_bracken"


phyloseq_outDir="${out_dir}/04_phyloseq"
