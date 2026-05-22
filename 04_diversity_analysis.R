# ============================================================
# VIRAL DIVERSITY ANALYSIS
# Alpha and Beta diversity metrics for virome data
#
# Input files (copy from pipeline output to INPUT_DIR):
#   *_all_viruses.txt   : Layer 1 known virus abundance
#   *_diamond.tsv       : Layer 2 divergent virus hits
#   *_summary.txt       : Read counts per layer per sample
#
# Output (saved to OUTPUT_DIR):
#   layer1_alpha_diversity.csv       : Alpha metrics for Layer 1
#   layer2_alpha_diversity.csv       : Alpha metrics for Layer 2
#   *_layer1_top20.pdf/png           : Top 20 viruses per sample (Layer 1)
#   *_layer2_top20.pdf/png           : Top 20 viruses per sample (Layer 2)
#   virome_composition_3layers.pdf   : Proportional bar chart of all layers
#   beta_bray_curtis.csv             : Bray-Curtis dissimilarity matrix
#   beta_jaccard.csv                 : Jaccard distance matrix
#   beta_pcoa_bray_curtis.pdf/png    : PCoA ordination plot
#   beta_heatmap_bray_curtis.pdf     : Heatmap of Bray-Curtis distances
#
# NOTE: Beta diversity requires >= 2 samples.
# ============================================================

# === PACKAGE INSTALLATION ===
if (!requireNamespace("vegan",   quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr",   quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("fossil",  quietly = TRUE)) install.packages("fossil")

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(fossil)

# === PATHS ===
INPUT_DIR  <- "."                  # directory containing pipeline output files
OUTPUT_DIR <- "./diversity_results"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ============================================================
# FUNCTION: calculate alpha diversity metrics
# ============================================================
calc_alpha <- function(counts, sample_name) {
  counts   <- counts[counts > 0]
  richness <- length(counts)
  shannon  <- diversity(counts, index = "shannon")
  simpson  <- diversity(counts, index = "simpson")
  chao1    <- chao1(counts)
  evenness <- shannon / log(richness)  # Pielou evenness

  data.frame(
    sample   = sample_name,
    richness = richness,
    shannon  = round(shannon, 4),
    simpson  = round(simpson, 4),
    chao1    = round(chao1,   4),
    evenness = round(evenness, 4)
  )
}

# ============================================================
# LAYER 1: known viruses (_all_viruses.txt)
# Columns: accession, length, reads
# ============================================================
cat("=== Processing Layer 1: known viruses ===\n")

layer1_files   <- list.files(INPUT_DIR, pattern = "_all_viruses.txt", full.names = TRUE)
layer1_metrics <- data.frame()
layer1_data    <- list()

for (f in layer1_files) {
  sample         <- gsub("_all_viruses.txt", "", basename(f))
  df             <- read.table(f, header = FALSE, col.names = c("accession", "length", "reads"))
  layer1_data[[sample]] <- df
  layer1_metrics <- rbind(layer1_metrics, calc_alpha(df$reads, sample))
}

cat("Layer 1 alpha diversity metrics:\n")
print(layer1_metrics)
write.csv(layer1_metrics, file.path(OUTPUT_DIR, "layer1_alpha_diversity.csv"), row.names = FALSE)

# ============================================================
# LAYER 2: divergent viruses (_diamond.tsv)
# Columns: qseqid sseqid pident length evalue bitscore stitle
# ============================================================
cat("\n=== Processing Layer 2: divergent viruses ===\n")

layer2_files   <- list.files(INPUT_DIR, pattern = "_diamond.tsv", full.names = TRUE)
layer2_metrics <- data.frame()
layer2_data    <- list()

for (f in layer2_files) {
  sample <- gsub("_diamond.tsv", "", basename(f))
  df     <- read.table(f, header = FALSE, sep = "\t",
                       col.names = c("qseqid", "sseqid", "pident", "length",
                                     "evalue", "bitscore", "stitle"))

  # Extract species name from brackets if present
  df$species <- gsub(".*\\[(.+)\\].*", "\\1", df$stitle)
  df$species[!grepl("\\[", df$stitle)] <- df$stitle[!grepl("\\[", df$stitle)]

  counts_df             <- df %>% count(species, name = "reads")
  layer2_data[[sample]] <- counts_df
  layer2_metrics        <- rbind(layer2_metrics, calc_alpha(counts_df$reads, sample))
}

cat("Layer 2 alpha diversity metrics:\n")
print(layer2_metrics)
write.csv(layer2_metrics, file.path(OUTPUT_DIR, "layer2_alpha_diversity.csv"), row.names = FALSE)

# ============================================================
# PLOTS — Layer 1 alpha diversity
# ============================================================
cat("\n=== Generating Layer 1 plots ===\n")

layer1_long <- layer1_metrics %>%
  pivot_longer(cols = c(richness, shannon, simpson, chao1, evenness),
               names_to = "metric", values_to = "value")

p1 <- ggplot(layer1_long, aes(x = sample, y = value, fill = sample)) +
  geom_bar(stat = "identity") +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Alpha Diversity — Known viruses (Layer 1)", x = "Sample", y = "Value")

ggsave(file.path(OUTPUT_DIR, "layer1_alpha_diversity.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(OUTPUT_DIR, "layer1_alpha_diversity.png"), plot = p1, width = 10, height = 8, dpi = 300)

# Top 20 viruses per sample (relative abundance)
for (sample in names(layer1_data)) {
  df <- layer1_data[[sample]] %>%
    arrange(desc(reads)) %>%
    head(20) %>%
    mutate(rel_abundance = reads / sum(layer1_data[[sample]]$reads) * 100)

  p2 <- ggplot(df, aes(x = reorder(accession, rel_abundance), y = rel_abundance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_bw() +
    labs(title = paste("Top 20 viruses — Layer 1 —", sample),
         x = "Virus", y = "Relative abundance (%)")

  ggsave(file.path(OUTPUT_DIR, paste0(sample, "_layer1_top20.pdf")), plot = p2, width = 10, height = 8)
  ggsave(file.path(OUTPUT_DIR, paste0(sample, "_layer1_top20.png")), plot = p2, width = 10, height = 8, dpi = 300)
}

# ============================================================
# PLOTS — Layer 2 top 20
# ============================================================
cat("=== Generating Layer 2 plots ===\n")

for (sample in names(layer2_data)) {
  df <- layer2_data[[sample]] %>%
    arrange(desc(reads)) %>%
    head(20) %>%
    mutate(rel_abundance = reads / sum(layer2_data[[sample]]$reads) * 100)

  p3 <- ggplot(df, aes(x = reorder(species, rel_abundance), y = rel_abundance)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    coord_flip() +
    theme_bw() +
    labs(title = paste("Top 20 divergent viruses — Layer 2 —", sample),
         x = "Species", y = "Relative abundance (%)")

  ggsave(file.path(OUTPUT_DIR, paste0(sample, "_layer2_top20.pdf")), plot = p3, width = 10, height = 8)
  ggsave(file.path(OUTPUT_DIR, paste0(sample, "_layer2_top20.png")), plot = p3, width = 10, height = 8, dpi = 300)
}

# ============================================================
# PLOTS — 3-layer composition summary
# ============================================================
cat("=== Generating 3-layer composition plot ===\n")

summary_files <- list.files(INPUT_DIR, pattern = "_summary.txt", full.names = TRUE)
summary_data  <- data.frame()

for (f in summary_files) {
  df           <- read.table(f, header = TRUE, sep = "\t")
  summary_data <- rbind(summary_data, df)
}

summary_long <- summary_data %>%
  pivot_longer(cols = c(layer1_known, layer2_divergent, layer3_unclassified),
               names_to = "layer", values_to = "reads") %>%
  group_by(sample) %>%
  mutate(percentage = reads / sum(reads) * 100)

summary_long$layer <- factor(summary_long$layer,
  levels = c("layer1_known", "layer2_divergent", "layer3_unclassified"),
  labels = c("Known viruses", "Divergent viruses", "Unclassified"))

p4 <- ggplot(summary_long, aes(x = sample, y = percentage, fill = layer)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("steelblue", "darkorange", "gray60")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Virome composition by classification layer",
       x = "Sample", y = "Percentage of reads (%)", fill = "Layer")

ggsave(file.path(OUTPUT_DIR, "virome_composition_3layers.pdf"), plot = p4, width = 10, height = 6)
ggsave(file.path(OUTPUT_DIR, "virome_composition_3layers.png"), plot = p4, width = 10, height = 6, dpi = 300)

# ============================================================
# BETA DIVERSITY (requires >= 2 samples)
# ============================================================
if (length(layer1_data) >= 2) {
  cat("\n=== Calculating beta diversity ===\n")

  if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
  library(pheatmap)

  all_species <- unique(unlist(lapply(layer1_data, function(df) df$accession)))

  comm_matrix <- do.call(rbind, lapply(names(layer1_data), function(s) {
    df    <- layer1_data[[s]]
    reads <- df$reads[match(all_species, df$accession)]
    reads[is.na(reads)] <- 0
    reads
  }))
  rownames(comm_matrix) <- names(layer1_data)
  colnames(comm_matrix) <- all_species

  # Bray-Curtis dissimilarity
  bray_dist <- vegdist(comm_matrix, method = "bray")
  cat("Bray-Curtis dissimilarity:\n")
  print(as.matrix(bray_dist))
  write.csv(as.matrix(bray_dist), file.path(OUTPUT_DIR, "beta_bray_curtis.csv"))

  # Jaccard distance
  jacc_dist <- vegdist(comm_matrix, method = "jaccard", binary = TRUE)
  cat("Jaccard distance:\n")
  print(as.matrix(jacc_dist))
  write.csv(as.matrix(jacc_dist), file.path(OUTPUT_DIR, "beta_jaccard.csv"))

  # PCoA with Bray-Curtis
  pcoa    <- cmdscale(bray_dist, eig = TRUE, k = 2)
  pcoa_df <- data.frame(sample = rownames(pcoa$points),
                        PC1 = pcoa$points[, 1],
                        PC2 = pcoa$points[, 2])
  var_exp <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

  p5 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, label = sample)) +
    geom_point(size = 4, color = "steelblue") +
    geom_text(vjust = -0.8, size = 3.5) +
    theme_bw() +
    labs(title = "PCoA — Bray-Curtis (Layer 1)",
         x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"))

  ggsave(file.path(OUTPUT_DIR, "beta_pcoa_bray_curtis.pdf"), plot = p5, width = 8, height = 6)
  ggsave(file.path(OUTPUT_DIR, "beta_pcoa_bray_curtis.png"), plot = p5, width = 8, height = 6, dpi = 300)

  # Heatmap of Bray-Curtis distances
  pdf(file.path(OUTPUT_DIR, "beta_heatmap_bray_curtis.pdf"), width = 8, height = 6)
  pheatmap(as.matrix(bray_dist),
           main = "Bray-Curtis dissimilarity",
           color = colorRampPalette(c("white", "steelblue"))(50),
           display_numbers = TRUE,
           number_format = "%.2f")
  dev.off()

  cat("  Beta diversity results saved to:", OUTPUT_DIR, "\n")

} else {
  cat("\n  Only 1 sample available — beta diversity skipped.\n")
  cat("  Add more samples to calculate Bray-Curtis, Jaccard, and PCoA.\n")
}

cat("\n=== Analysis complete. Results in:", OUTPUT_DIR, "===\n")
