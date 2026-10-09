#!/usr/bin/env Rscript
##################################################################
# Visualize the results of the 16S pipeline using phyloseq
# Will read tabels from bracken and dada pipelines and create plots and tables for visualization
##################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Phyloseq needs to be run with both an input path and an output directory (and optionally a sample sheet).")
}

input_path <- args[1]
output_dir <- args[2]
sample_sheet <- if (length(args) >= 3) args[3] else NA

# Fixed seed so that NMDS and the PERMANOVA (adonis2) permutation p-values are reproducible
set.seed(1234)

cat("Reading from:", input_path, "\n", "Writing to:", output_dir, "\n")

# ===============================
# Load packages
# ===============================
library(phyloseq)
library(ggplot2)
library(vegan)
library(grid)
# ===============================
# Read Bracken combined output
# ===============================

bracken <- read.delim(input_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# ===============================
# Build OTU table
# ===============================

# Keep only count columns (*_num)
count_cols <- grep("_num$", colnames(bracken), value = TRUE)

otu <- bracken[, count_cols]

# Remove "_num" from sample names
colnames(otu) <- sub("_num$", "", colnames(otu))

# Genus names, without the "-taxid" suffix added by the modified combine_bracken_outputs.py.
# Only the trailing taxid is removed, so genus names containing hyphens (e.g. Escherichia-Shigella) are kept intact.
# Names without the suffix (stock combine_bracken_outputs.py) are kept as they are
genus <- if ("taxonomy_id" %in% colnames(bracken)) {
  taxid_suffix <- paste0("-", bracken$taxonomy_id)
  ifelse(
    endsWith(bracken$name, taxid_suffix),
    substr(bracken$name, 1, nchar(bracken$name) - nchar(taxid_suffix)),
    bracken$name
  )
} else {
  sub("-[0-9]+$", "", bracken$name)
}

# Unique taxa names: genus name + taxonomy id
taxa_ids <- if ("taxonomy_id" %in% colnames(bracken)) paste(genus, bracken$taxonomy_id, sep = "-") else make.unique(bracken$name)

rownames(otu) <- taxa_ids

OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)

# ===============================
# Build taxonomy table
# ===============================

tax <- data.frame(Kingdom = "Bacteria", Genus = genus, row.names = taxa_ids, stringsAsFactors = FALSE)

TAX <- tax_table(as.matrix(tax))

# ===============================
# Build sample metadata
# ===============================

meta <- data.frame(SampleID = colnames(otu), row.names = colnames(otu), stringsAsFactors = FALSE)

SAM <- sample_data(meta)

# ===============================
# Sample groups
# ===============================
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

# ===============================
# Build phyloseq object
# ===============================

ps <- phyloseq(OTU, TAX, SAM)

# Check object
#print(ps)

# Check taxonomy ranks
#print(rank_names(ps))

# Get sample names
samples <- sample_names(ps)

group <- get_groups(samples, sample_sheet)

# Add group information to phyloseq metadata
sample_data(ps)$Group <- group

# Check
#sample_data(ps)
# ===============================
# Relative abundance
# ===============================

ps.rel <- transform_sample_counts(ps,function(x) x / sum(x))

# ===============================
# Keep top 30 genera + Other
# ===============================

# Calculate total abundance of each genus across all samples
taxa_abundance <- taxa_sums(ps.rel)

# Get the top 30 most abundant genera (or all of them, if there are fewer than 30)
top <- names(sort(taxa_abundance, decreasing = TRUE))[seq_len(min(30, length(taxa_abundance)))]

# Keep the top 30
ps.top <- prune_taxa(top, ps.rel)

# Identify everything outside the top 30
other <- setdiff(taxa_names(ps.rel), top)

# Combine all non-top-30 genera into "Other"
if (length(other) > 0) {
  # Extract abundance matrix for the non-top-30 genera
  other_matrix <- as.matrix(otu_table(prune_taxa(other, ps.rel)))

  # Make sure taxa are rows
  if (!taxa_are_rows(ps.rel)) {
  other_matrix <- t(other_matrix)
  }

  # Sum all non-top-30 genera for each sample
  other_counts <- colSums(other_matrix)

  # Create OTU table for Other
  OTHER_OTU <- otu_table(matrix(other_counts,nrow = 1,dimnames = list("Other", names(other_counts))), taxa_are_rows = TRUE)

  # Create taxonomy table
  OTHER_TAX <- tax_table(matrix("Other",nrow = 1,dimnames = list("Other", "Genus")))

  # Create phyloseq object
  ps.other <- phyloseq(OTHER_OTU,OTHER_TAX,sample_data(ps.rel))

  # Combine top 30 + Other
  ps.top <- merge_phyloseq(ps.top,ps.other)

}

# Make sure Genus is stored as character, NOT factor
tax_table(ps.top)[, "Genus"] <- as.character(tax_table(ps.top)[, "Genus"])

# ===============================
# Plot
# ===============================
# Extract data
plot_data <- psmelt(ps.top)

#Get genus names
genus_names <- unique(as.character(plot_data$Genus))

#  Set Other as last level in Genus factor
genus_names <- c("Other",setdiff(genus_names, "Other"))

plot_data$Genus <- factor(as.character(plot_data$Genus),levels = genus_names)

#Create colors for all genera
genus_colors <- setNames(
grDevices::hcl.colors(
length(genus_names) - 1,
palette = "Dark 3"
),genus_names[genus_names != "Other"])

# Add grey for Other
genus_colors["Other"] <- "grey70"

#Plot Barplot
p <- ggplot(
plot_data, aes(
x = Sample, y = Abundance, fill = Genus)
) +
geom_bar(
stat = "identity", position = "stack"
) +
scale_fill_manual(
values = genus_colors, breaks = genus_names, drop = FALSE
) +
theme_bw() +
labs(
x = "Sample", y = "Relative abundance", fill = "Genus"
) +
guides(
fill = guide_legend(ncol = 1)
) +
theme(
axis.text.x = element_text( angle = 90, hjust = 1, size = 10
),
axis.text.y = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 12)
)

# Save as PDF
ggsave(
  filename = file.path(output_dir, "genus_barplot.pdf"),
  plot = p, width = 25, height = 11, units = "in", create.dir = TRUE
)

hetmap_top <- plot_heatmap(ps.top, method = "NMDS", distance = "bray")

ggsave(
  filename = file.path(output_dir, "heatmap.pdf"),
  width = 22, plot = hetmap_top, height = 11, units = "in", create.dir = TRUE
)

bray_dist <- phyloseq::distance(ps.rel, method = "bray")

ordination <- ordinate(ps.rel, method = "PCoA", distance = bray_dist)

############################################################
# PERMANOVA
metadata <- as(sample_data(ps.rel), "data.frame")

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

pcoa_plot <- plot_ordination(
    ps.rel, ordination, color = "Group"
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
        x = "PCoA1", y = "PCoA2", color = "Group"
    ) +
    theme(
        plot.title = element_text(
            size = 16, face = "bold", hjust = 0.5
        ),
        plot.subtitle = element_text(
            size = 11, hjust = 0.5
        )
    )


ggsave(
    filename = file.path(output_dir, "pcoa_plot.pdf"),
    plot = pcoa_plot, width = 16, height = 9, units = "in", create.dir = TRUE
)

##################################################
#heatmap with not only top
ps.heat <- filter_taxa(ps.rel, function(x) mean(x) > 0.001, prune = TRUE)

heat <- plot_heatmap(ps.heat, method = "NMDS", distance = "bray", sample.label = "Group")

ggsave(
  filename = file.path(output_dir, "heatmap_all.pdf"),
  plot = heat, width = 16, height = 9, units = "in", create.dir = TRUE
)

##################################################3
# mupltiple alpha diversity plots
alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

alpha$Group <- sample_data(ps)$Group

# Shannon
p1 <- ggplot(alpha, aes(Group, Shannon, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1) +
  theme_bw() +
  ggtitle("Shannon diversity")

# Simpson
p2 <- ggplot(alpha, aes(Group, Simpson, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1) +
  theme_bw() +
  ggtitle("Simpson diversity")

pdf(file.path(output_dir, "alpha_diversity_all.pdf"), width = 12, height = 6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
print(p1, vp = viewport(layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.col = 2))
dev.off()

cat("Phyloseq analysis completed. Plots saved to:", output_dir, "\n")
