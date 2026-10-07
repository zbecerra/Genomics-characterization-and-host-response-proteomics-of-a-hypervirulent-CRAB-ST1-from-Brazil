# Your original R script. Changes: (1) hard-coded personal path replaced by a relative one (BASE_DIR); (2) column labels corrected
# (ATCC_A-C = ATCC 19606, FAGO_D-F = BR5; the original plot labels were reversed); (3) regulation direction computed from the
# group means instead of the 'Difference ATCC_FAGO' column. Regenerate Figure 5 with this version and re-check the figure and text.
# Input: proteome/reguladas_pvalue.txt (table of differentially abundant proteins with LFQ intensities, 3 replicates per strain; produced outside this repository).
# ======================================================================
# Proteome heatmap - A. baumannii BR5 vs ATCC19606
# Matching reference figure style exactly
# ======================================================================

library(pheatmap)
library(RColorBrewer)
library(grid)

BASE_DIR <- "proteome"   # folder with reguladas_pvalue.txt (edit as needed)

# ── 1. Load data ───────────────────────────────────────────────────────────
df <- read.delim(file.path(BASE_DIR, "reguladas_pvalue.txt"),
                 sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

cat("Proteins loaded:", nrow(df), "\n")

# ── 2. Build matrix ────────────────────────────────────────────────────────
intensity_cols <- c("ATCC_A","ATCC_B","ATCC_C","FAGO_D","FAGO_E","FAGO_F")

# Use protein ID as row names
df$protein_id <- sub(";.*", "", df$`T: Majority protein IDs`)

mat <- as.matrix(df[, intensity_cols])
rownames(mat) <- df$protein_id

# ── 3. Z-score normalize per row (each protein) ───────────────────────────
mat_scaled <- t(scale(t(mat)))
cat("Z-score range:", round(min(mat_scaled, na.rm=TRUE), 2),
    "to", round(max(mat_scaled, na.rm=TRUE), 2), "\n")

# ── 4. Column annotation (sample groups) ──────────────────────────────────
ann_col <- data.frame(
  Strain = c("ATCC19606","ATCC19606","ATCC19606","BR5","BR5","BR5"),
  row.names = intensity_cols
)

# Rename columns to match reference figure labels
colnames(mat_scaled) <- c("ATCC_A","ATCC_B","ATCC_C","FAGO_D","FAGO_E","FAGO_F")

# Custom x-axis labels matching reference
# CORRECTED: ATCC_* = ATCC 19606 (non-virulent), FAGO_* = BR5 (virulent)
col_labels <- c(
  "ATCC_A" = "A. baumannii ATCC19606",
  "ATCC_B" = "A. baumannii ATCC19606",
  "ATCC_C" = "A. baumannii ATCC19606",
  "FAGO_D" = "A. baumannii BR5",
  "FAGO_E" = "A. baumannii BR5",
  "FAGO_F" = "A. baumannii BR5"
)

# ── 5. Colors ──────────────────────────────────────────────────────────────
# Matching reference: blue → white → red
col_breaks <- seq(-2, 2, length.out = 101)
col_colors <- colorRampPalette(c("#0000ff","#673dff","#9265ff","#b38bff",
                                  "#f7f7f7",
                                  "#ff9e81","#ff7b5a","#ff5332","#ff0500"))(100)

ann_colors <- list(
  Strain = c("ATCC19606" = "#ffffff", "BR5" = "#ffffff")
)

# ── 6. Draw heatmap using pheatmap ────────────────────────────────────────
pdf(file.path(BASE_DIR, "proteome_heatmap.pdf"), width = 7, height = 10)

pheatmap(
  mat_scaled,
  color            = col_colors,
  breaks           = col_breaks,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  annotation_col   = ann_col,
  annotation_colors = ann_colors,
  labels_col       = c("A. baumannii\nATCC19606",
                        "A. baumannii\nATCC19606",
                        "A. baumannii\nATCC19606",
                        "A. baumannii\nBR5",
                        "A. baumannii\nBR5",
                        "A. baumannii\nBR5"),   # CORRECTED order: columns ATCC_A-C, FAGO_D-F
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  fontsize_row     = 10,
  fontsize_col     = 10,
  fontsize         = 11,
  border_color     = "white",
  treeheight_row   = 60,
  treeheight_col   = 40,
  legend_breaks    = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
  legend_labels    = c("-1.5","-1","-0.5","0","0.5","1","1.5"),
  main             = ""
)

dev.off()

# Also save PNG
png(file.path(BASE_DIR, "proteome_heatmap.png"),
    width = 7, height = 7, units = "in", res = 900)

pheatmap(
  mat_scaled,
  color            = col_colors,
  breaks           = col_breaks,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  annotation_col   = ann_col,
  annotation_colors = ann_colors,
  labels_col       = c("A. baumannii\nATCC19606",
                        "A. baumannii\nATCC19606",
                        "A. baumannii\nATCC19606",
                        "A. baumannii\nBR5",
                        "A. baumannii\nBR5",
                        "A. baumannii\nBR5"),   # CORRECTED order: columns ATCC_A-C, FAGO_D-F
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  fontsize_row     = 10,
  fontsize_col     = 10,
  fontsize         = 11,
  border_color     = "NA",
  treeheight_row   = 60,
  treeheight_col   = 40,
  legend_breaks    = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
  legend_labels    = c("-1.5","-1","-0.5","0","0.5","1","1.5"),
  main             = ""
)

dev.off()

cat("Heatmap saved!\n")

# ── 7. Print summary ───────────────────────────────────────────────────────
# CORRECTED: direction from the intensities themselves (no dependence on a difference-column name or sign convention)
df$direction <- ifelse(rowMeans(df[, c("ATCC_A","ATCC_B","ATCC_C")]) > rowMeans(df[, c("FAGO_D","FAGO_E","FAGO_F")]),
                       "Up in ATCC19606", "Up in BR5")
cat("\nRegulation summary:\n")
print(table(df$direction))
