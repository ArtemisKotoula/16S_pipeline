#############################################################
# Packages

packages <- c("dada2", "ggplot2", "phyloseq", "Biostrings")
for (package_name in packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name)
  }
  library(package_name, character.only = TRUE); packageVersion(package_name)
}

###############################################################
# Arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("dada2.R needs to be run with an input path, an output directory a right and a left filter length and the dada database path.")
}

path <- args[1]
out_path <- args[2]
right_length <- as.numeric(args[3])
left_length <- as.numeric(args[4])
dada_db <- args[5]

cat("Reading from:", path, "\n")
cat("Writing to:", out_path, "\n")
cat("Filter lengths: right =", right_length, ", left =", left_length, "\n")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
# They were created by cutadapt
# path to all forward and reverse fastq files
fnFs <- sort(list.files(path, pattern="_R1_trimmed.fastq", full.names = TRUE, recursive = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_trimmed.fastq", full.names = TRUE, recursive = TRUE))

head(fnFs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#View the first few sample names
head(sample.names)

QplotF <- plotQualityProfile(fnFs[1:2])

ggsave(
    filename = file.path(out_path, "quality_profileF.png"),
    create.dir = TRUE, plot = QplotF, width = 10, height = 6, dpi = 300,
)
)

QplotR <- plotQualityProfile(fnRs[1:2]) 

ggsave(
    filename = file.path(out_path, "quality_profileR.png"),
    create.dir = TRUE,
)

filtFs <- file.path(out_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(out_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

##running filtering with truncation lengths of 260 for forward reads and 220 for reverse reads

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(right_length, left_length),
              maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 

head(out)

# Learn error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plot error rates for forward reads
png(file.path(out_path, "errF_plot.png"), width = 1200, height = 900, res = 300)
plotErrors(errF, nominalQ = TRUE)
dev.off()


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers [[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# # Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


taxa <- assignTaxonomy(seqtab.nochim, dada_db, multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.csv(
  taxa.print,
  file = file.path(out_path, "taxa_Genus.csv"),
  row.names = FALSE
)
