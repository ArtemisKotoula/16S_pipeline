#############################################################
# Packages
packages <- c("dada2", "ggplot2", "phyloseq", "Biostrings", "vegan")
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
  stop("dada2.R needs to be run with an input path, an output directory a right and a left filter length and the dada database path (and optionally a sample sheet).")
}

path <- args[1]
out_path <- args[2]
right_length <- as.numeric(args[3])
left_length <- as.numeric(args[4])
dada_db <- args[5]
sample_sheet <- if (length(args) >= 6) args[6] else NA

# Fixed seed so that NMDS and the PERMANOVA (adonis2) permutation p-values are reproducible
set.seed(1234)

cat("Reading from:", path, "\n")
cat("Writing to:", out_path, "\n")
cat("Filter lengths: right =", right_length, ", left =", left_length, "\n")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
# They were created by cutadapt
# path to all forward and reverse fastq files
fnFs <- sort(list.files(path, pattern="_R1_filtered.fastq", full.names = TRUE, recursive = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_filtered.fastq", full.names = TRUE, recursive = TRUE))

# Extract sample names, filenames have format: SAMPLENAME_R1_filtered.fastq
# (sample names may contain "_", so only the suffix is removed)
sample.names <- sub("_R1_filtered\\.fastq$", "", basename(fnFs))

# Sample groups
# Groups are read from the sample sheet (tab-separated, columns "sample" and "group") if it exists,
# otherwise they are inferred from the sample names. The regex pattern captures:
# - optional numbers at the beginning
# - followed by uppercase letters (lowercase letters are ignored)
get_groups <- function(samples, sample_sheet) {
  if (!is.na(sample_sheet) && nzchar(sample_sheet) && file.exists(sample_sheet)) {
    cat("Reading sample groups from:", sample_sheet, "\n")
    sheet <- read.delim(sample_sheet, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("sample", "group") %in% colnames(sheet))) {
      stop("Sample sheet must contain the columns 'sample' and 'group': ", sample_sheet)
    }
    missing <- setdiff(samples, sheet$sample)
    if (length(missing) > 0) {
      stop("Samples missing from the sample sheet: ", paste(missing, collapse = ", "))
    }
    return(as.character(sheet$group[match(samples, sheet$sample)]))
  }
  cat("No sample sheet found. Inferring groups from sample names.\n")
  sub("^[0-9]*([A-Z]+).*", "\\1", samples)
}

# Read the groups before running DADA2, so that sample sheet errors are reported immediately
sample.groups <- get_groups(sample.names, sample_sheet)

QplotF <- plotQualityProfile(fnFs[1:2])

ggsave(
  filename = file.path(out_path, "quality_profileF.png"),
  create.dir = TRUE, plot = QplotF, width = 10, height = 6, dpi = 300
)

QplotR <- plotQualityProfile(fnRs[1:2]) 

ggsave(
  filename = file.path(out_path, "quality_profileR.png"),
  create.dir = TRUE, plot = QplotR, width = 10, height = 6, dpi = 300
)

filtFs <- file.path(out_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(out_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(0,0),
              maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 

# Learn error rates for forward and reverse reads
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plot error rates for forward reads
png(file.path(out_path, "errF_plot.png"), width = 1200, height = 900, res = 300)
plotErrors(errF, nominalQ = TRUE)
dev.off()

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Number of reads kept at each step, per sample
write.csv(track, file = file.path(out_path, "read_tracking.csv"), row.names = TRUE)

taxa <- assignTaxonomy(seqtab.nochim, dada_db, multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

write.csv(taxa.print, file = file.path(out_path, "taxa_Genus.csv"), row.names = FALSE)

###############################################################
# ASV outputs, linked by ASV ID
asv_seqs <- colnames(seqtab.nochim)
asv_ids <- sprintf("ASV%d", seq_along(asv_seqs))

# ASV sequences
writeXStringSet(DNAStringSet(setNames(asv_seqs, asv_ids)), file.path(out_path, "ASVs.fasta"))

# ASV counts (rows: ASVs, columns: samples)
asv_counts <- t(seqtab.nochim)
rownames(asv_counts) <- asv_ids
write.table(
  data.frame(ASV_ID = asv_ids, asv_counts, check.names = FALSE),
  file = file.path(out_path, "ASV_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

# ASV taxonomy
asv_taxa <- taxa[match(asv_seqs, rownames(taxa)), , drop = FALSE]
write.table(
  data.frame(ASV_ID = asv_ids, Sequence = asv_seqs, asv_taxa, check.names = FALSE),
  file = file.path(out_path, "ASV_taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

print("Starting phyloseq visualization")

samdf <- data.frame(
  SampleID = sample.names, Group = sample.groups, row.names = sample.names
)

OTU <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
TAX <- tax_table(taxa)
SAM <- sample_data(samdf)

ps <- phyloseq(OTU, TAX, SAM)

ps <- prune_samples(sample_sums(ps) > 0, ps)
# ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Save the phyloseq object for further analysis
saveRDS(ps, file = file.path(out_path, "phyloseq_object.rds"))

# Order samples by Group
sample_order <- order(sample_data(ps)$Group)

ps <- prune_samples(sample_names(ps)[sample_order], ps)

#ALpha DIversity
rich <- plot_richness(
  ps, x = "Group", measures = c("Shannon", "Simpson"), color = "Group"
) +
  geom_boxplot(
    aes(group = Group), alpha = 0.2, outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15, size = 2, alpha = 0.8
  ) +
  facet_wrap(
    ~variable, scales = "free_y"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none", strip.background = element_rect(fill = "white")
  ) +
  labs(
    x = "Group", y = "Diversity"
  )

ggsave(
  filename = file.path(out_path, "alpha_diversity.png"),
  plot = rich, create.dir = TRUE, width = 10, height = 6, dpi = 300
)

# Transform to relative abundance
ps.prop <- transform_sample_counts(ps, function(otu) otu / sum(otu))
ps.genus <- tax_glom(ps.prop, taxrank = "Genus", NArm = FALSE)

# Top 20 taxa
genus_abund <- taxa_sums(ps.genus)
top_genera <- names(sort(genus_abund, decreasing = TRUE) )[1:min(20, length(genus_abund))]
ps.top_genus <- prune_taxa(top_genera,  ps.genus)

# BAR PLOT - TOP 20 GENERA
bar_genus <- plot_bar(
  ps.top_genus, x = "SampleID", fill = "Genus"
) +
  facet_wrap(
    ~Group, scales = "free_x"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.title.x = element_blank(), legend.position = "right"
  ) +
  labs(
    title = "Relative abundance of top 20 bacterial genera", y = "Relative abundance"
  )

ggsave(
  filename = file.path(out_path, "relative_abundance_Genus_top20.pdf"),
  plot = bar_genus, width = 16, height = 9, units = "in", create.dir = TRUE
)
###############################################################
# BRAY-CURTIS NMDS
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")

ord <- plot_ordination(
  ps.prop, ord.nmds.bray, color = "Group", label = "SampleID"
) +
  geom_point(
    size = 4
  ) +
  theme_bw() +
  labs(
    title = "Bray-Curtis NMDS", color = "Group"
  ) +
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path( out_path, "nmds_bray_curtis.pdf"),
  plot = ord, width = 12, height = 8, units = "in", create.dir = TRUE
)
###############################################################
#Distance matrix
bray_dist <- phyloseq::distance(ps.genus, method = "bray")

# PERMANOVA
metadata <- as(sample_data(ps.genus), "data.frame")
adonis_stats <- adonis2(bray_dist ~ Group, data = metadata)

# Extract PERMANOVA results
permanova_R2 <- adonis_stats$R2[1]
permanova_p <- adonis_stats$`Pr(>F)`[1]

# PERMDISP / beta dispersion
disp <- betadisper(bray_dist, metadata$Group)
disp_anova <- anova(disp)

# Extract p-value
permdisp_p <- disp_anova$`Pr(>F)`[1]

# Print results to terminal/log
cat("\n===== PERMANOVA =====\n")
cat("R2 =", permanova_R2, "\n")
cat("p =", permanova_p, "\n")

cat("\n===== PERMDISP =====\n")
cat("p =", permdisp_p, "\n")

#############################################
pcoa <- plot_ordination(
  ps.genus, ordinate(ps.genus, method = "PCoA", distance = bray_dist), color = "Group"
) +
  geom_point(
      size = 3, alpha = 0.8
  ) +
  stat_ellipse(
      aes(group = Group, linetype = "90%"), level = 0.90
  ) +
  stat_ellipse(
      aes(group = Group, linetype = "95%"), level = 0.95
  ) +
  scale_linetype_manual(
      name = "Ellipse", values = c("90%" = "dotted", "95%" = "dashed")
  ) +
  theme_bw()+
  labs(
      title = "PCoA of microbial community composition",
      subtitle = paste0(
          "PERMANOVA: R² = ", round(permanova_R2, 3),
          ", p = ", signif(permanova_p, 3),
          "  |  ",
          "PERMDISP: p = ", signif(permdisp_p, 3)
      ),
      x = "PCoA1",
      y = "PCoA2",
      color = "Group"
  ) +
  theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5)
  )

ggsave(
  filename = file.path(out_path, "pcoa_genus_bray_curtis.pdf"),
  plot = pcoa, width = 12, height = 8, units = "in", create.dir = TRUE
)
###############################################################
# Genus abundance table
genus_table <- psmelt(ps.genus)
write.csv(genus_table,file = file.path(out_path, "genus_relative_abundance.csv"),row.names = FALSE)
###############################################################
# SAMPLE METADATA
write.csv(samdf,file = file.path(out_path, "sample_metadata.csv"),row.names = TRUE)

