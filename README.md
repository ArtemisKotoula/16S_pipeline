# 16S_pipeline
A bioinformatics pipeline for processing paired-end 16S rRNA amplicon sequencing data — from raw FASTQ reads through quality control, primer trimming, quality filtering, taxonomic classification, and downstream statistical visualization.
 
The pipeline will support three alternative analysis branches after primer trimming:
 
- **Kraken2 / Bracken** branch (fully implemented in this repo) — read-based taxonomic classification, abundance re-estimation, and `phyloseq`-based visualization.
- **DADA2** branch (not fully implemented, viaualization pending) — ASV-based analysis.
- **QIIME2** branch (not implemented) — ASV-based analysis.

## Pipeline Overview

The Overview of the current, complete pipeline

### 1. Preprocessing

1.1 FastQC and MultiQC on raw reads (parallel_qc.sh)

1.2 Primer trimming using cutadapt (trim.sh)

1.3 FastQC and MultiQC on trimmed reads (parallel_qc.sh)

### 2.1 Kraken Branch

2.1.1 Quality FIltering using fastp (fastp.sh)

2.1.2 Summarize fastp reports and remove low depth samples (report_fastp.sh and calc_cutoff.py)

2.1.3 FastQC and MultiQC on filtered reads (parallel_qc.sh)

2.1.4 Kraken2 classification, Krona visualization and Bracken re-estimation (kraken_pipeline.sh)

2.1.5 Statistics and plots visualization with phyloseq (phyloseq.R)

### 2.2 Dada2 Branch

2.2.1 Create quality plots for reads (dada2.R)

2.2.2 Filter and Trim right and left reads to fixed length (dada2.R)

2.2.3 DADA2 classification (dada2.R)


## Repository Structure

| File | Purpose |
|---|---|
| `16S_main.sh` | Main entry point. Runs the shared QC + trimming steps, then dispatches to the pipeline branch. |
| `config/config.sh` | Central configuration: paths, primers, thresholds, and all output directories. |
| `config/conda.sh` | Defines helper functions (`activate_dada`, `activate_cutadapt`, `activate_kraken`) that activate the required conda environments. |
| `config/<CONDA_ENV>_env.yml` | yml files for creating necessary conda enviromntets |
| `parallel_qc.sh` | Runs FastQC in parallel across samples, then aggregates results with MultiQC. Used at three separate stages (raw, trimmed, filtered reads). |
| `trim.sh` | Removes primer sequences from raw reads with `cutadapt`. |
| `fastp.sh` | Quality/length filtering of primer-trimmed reads with `fastp` (used ahead of the Kraken branch). |
| `report_fastp.sh` | Parses all per-sample `fastp` JSON reports into a combined CSV, computes summary statistics, calculates a minimum-read-depth cutoff, and moves under-depth samples out of the analysis set. |
| `calc_cutoff.py` | Computes the minimum number of reads needed to detect a taxon at a given frequency and confidence level (binomial survival function), used by `report_fastp.sh`. |
| `kraken_pipeline.sh` | Runs Kraken2 classification, builds a Krona plot, builds/runs Bracken, combines Bracken output across samples, and calls `phyloseq.R`. |
| `phyloseq.R` | Builds a `phyloseq` object from the combined Bracken table and produces genus barplots, heatmaps, PCoA ordination with PERMANOVA/`betadisper`, and alpha-diversity plots. |
| `dada_pipeline.sh` | Runs the DADA2 classification through calling `dada2.R` |
| `dada2.R` | The R script that filters the reads, learns the error rates, runs the dada algorithm, and assigns taxonomy. |


## Repository / directory layout expected by the scripts
 
`16S_main.sh` sources its configuration from a `config/` subdirectory relative to its own location:
 
```
bin/
|── config/
    ├── 16Scutadapt_env.yml
    ├── 16Sdada_env.yml
    ├── 16Skraken_env.yml
    ├── config.sh
    └── conda.sh
├── 16S_main.sh
├── calc_cutoff.py
├── dada_pipeline.sh
├── dada2.R
├── fastp.sh
├── kraken_pipeline.sh
├── parallel_qc.sh
├── phyloseq.R
├── report_fastp.sh
└── trim.sh
```

 ## Configuration (`config/config.sh`)
 
All paths, primers, and thresholds live in one file:
 
| Variable | Description | Example |
|---|---|---|
| `threads` | Number of threads used across steps | `12` |
| `raw_data` | Path to raw FASTQ input directory | - |
| `out_dir` | Timestamped root output directory for the run | `results/16S_<timestamp>` |
| `r1_primer` / `r2_primer` | Forward/reverse primer sequences trimmed by cutadapt | - |
| `frequency` | Expected minimum taxon frequency for the depth cutoff calculation | `0.01` |
| `confidence` | Required detection probability for the depth cutoff calculation | `0.99` |
| `min_reads` | Minimum reads a taxon must have to be considered present | `10` |
| `kraken_db`, `dada_db` | Paths to pre-built Kraken2 and DADA2 databases| — |
| `quality_threshold` | Mean-quality cutoff passed to fastp (`--cut_right_mean_quality`) | `20` |
| `min_length` | Minimum read length passed to fastp / cutadapt | `100` |
| `right_len`, `left_len` | Lengths at which the right and left reads will be truncated at for DADA2 filtering | 260, 220 |
| `fastqc_outDir`, `multiqc_outDir`, `cutadapt_outDir`, `fastp_outDir`, `report_fastp_outDir`, `calc_cutoff_outDir`, `kraken_outDir`, `krona_outDir`, `bracken_outDir`, `phyloseq_outDir`, `dada2_outDir` | Per-step output subdirectories, all nested under `out_dir` | — |
| `skip_pre` | A boolean variable that is used when run with the --rerun option to skip the preprocessing steps. Is by default "false" and changes automatically. | `false` |
 
Edit these values (in particular `raw_data`, `out_dir`, `kraken_db` / `dada_db`, and the primer sequences) before running the pipeline.

 
## Expected input structure
 
Raw data should be organized as one subdirectory per sample, each containing a forward and reverse read file with `R1`/`R2` in the filename (`.fastq` or `.fastq.gz`):
 
```
raw_data/
├── Sample1/
│   ├── Sample1_R1.fastq.gz
│   └── Sample1_R2.fastq.gz
├── Sample2/
│   ├── Sample2_R1.fastq.gz
│   └── Sample2_R2.fastq.gz
└── ...
```
 
## Usage

```bash
bash 16S_main.sh "kraken|dada" [--rerun]
```
 
The pipeline logs all output to `<out_dir>/pipeline.log` (via `tee`) in addition to the terminal.

There is a `--rerun` option for the pipeline to skip the common preprocessing steps, where the previous run directory is located and if results are located within the cutadapt directory, the pipeline continues without the preprocessing of the raw files.

## Notes and known limitations
 

- **`combine_bracken_outputs.py` requires a local modification.** Per the comment in `kraken_pipeline.sh`, the stock script must be patched to append the taxon ID to the name (`name = f"{name}-{taxid}"`) to avoid duplicate row names when combining outputs — otherwise `phyloseq.R` will fail to build unique taxa names.
- **Sample grouping in `phyloseq.R`** is inferred from sample names via the regex `^[0-9]*([A-Z]+).*` (leading digits stripped, then leading uppercase letters taken as the group). Rename samples accordingly, or adjust this regex if your sample naming convention differs.
- **Primers** in the default config should be updated according to amplicon region.
- The databases for both Kraken and DADA must be built/downloaded separately and their path set in `config.sh` before running.

- **How the depth cutoff is calculated.** `calc_cutoff.py` finds the minimum total read count `N` needed so that a taxon at relative abundance `frequency` has at least a `confidence` probability of getting `min_reads` reads, using a Binomial(N, `frequency`) model. `report_fastp.sh` then moves any sample below that cutoff, from the fastp output directory into a "failed cutoff samples" directory, so that it is not included in the downstream analysis.
- `report_fastp.sh` computes `avg_mean_l` — the average of the post-filtering R1 and R2 mean read lengths across all samples — from the fastp summary statistics it just aggregated. Later, it is used as the read-length parameter for both the Bracken database build (`bracken-build -l ${avg_mean_l}`) and every per-sample Bracken run (`bracken -r ${avg_mean_l}`), so Bracken's abundance re-distribution is matched to the actual (filtered) read length of the dataset rather than a hardcoded value.

- In order to limit DADA to a specified number of threads, the `taskset` command is used. The dada processes are then limited to the first N threads of the system, specified by the `threads` parameter in `config.sh`.

## Future Major Releases

- Complete the visualization part of the dada brach of the pipeline
- Start implemetation of the qiime branch
