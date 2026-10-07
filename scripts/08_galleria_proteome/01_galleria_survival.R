# Verified script: your original R script; only the hard-coded personal path was replaced by a relative one (see setwd below).
# Input: gm_a.baumannii_BR5.xlsx (sheet 'A. baumannii': time + number of live larvae per replicate, columns <strain>_R1..R3; the workbook is not included in this repository yet).
# Output: survival curves, pairwise Fleming-Harrington(0,1) FDR matrix, supplementary pairwise p-value workbook.
# ======================================================================
# Galleria mellonella survival analysis – Acinetobacter baumannii BR5
# ======================================================================

# ----------------------------
# 1. Load packages and define alfa
# ----------------------------
library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(coin)
library(viridis)
library(scales)
library(dendextend)
library(ggplot2)
library(reshape2)
library(cluster)
library(ggdendro)   # needed for dendro_data()

# Alfa value
alpha <- 0.05

# ----------------------------
# 2. Read and prepare data
# ----------------------------
# setwd("/path/Users\Name\Documents/")
setwd("gm_a.baumannii_BR5")   # folder containing gm_a.baumannii_BR5.xlsx (edit as needed)

library(readxl)

excel_sheets("gm_a.baumannii_BR5.xlsx")

data_counts <- read_excel("gm_a.baumannii_BR5.xlsx", sheet = "A. baumannii")

names(data_counts)[1] <- "time"
data_counts$time <- as.numeric(data_counts$time)
data_counts <- data_counts %>% filter(!is.na(time)) %>% arrange(time)

if (!any(data_counts$time == 0)) {
  zero_row <- data.frame(time = 0)
  for (col in names(data_counts)[-1]) zero_row[[col]] <- 10
  data_counts <- bind_rows(zero_row, data_counts) %>% arrange(time)
}

cat("Data dimensions:", dim(data_counts), "\n")
cat("Time points:", unique(data_counts$time), "\n")

# ----------------------------
# 3. Convert counts to individual survival data
# ----------------------------
replicate_to_survival <- function(time_vec, count_vec, start_n = 10) {
  time_full  <- c(0, time_vec)
  count_full <- c(start_n, count_vec)
  
  # NEW: catch data entry errors before they silently vanish
  if (any(diff(count_full) > 0)) {
    warning("Non-monotonic counts detected — some intervals will be skipped.")
  }
  
  res <- data.frame(time = numeric(), event = integer())
  for (i in 2:length(time_full)) {
    n_deaths <- count_full[i-1] - count_full[i]
    if (n_deaths > 0)
      res <- rbind(res, data.frame(time = rep(time_full[i], n_deaths), event = 1))
  }
  final_alive <- tail(count_full, 1)
  if (final_alive > 0)
    res <- rbind(res, data.frame(time = rep(tail(time_full, 1), final_alive), event = 0))
  return(res)
}

all_individuals <- list()
time_vals <- data_counts$time

for (col_name in names(data_counts)[-1]) {
  strain   <- sub("_R[0-9]+$", "", col_name)
  counts   <- data_counts[[col_name]]
  surv_data <- replicate_to_survival(time_vals, counts, start_n = 10)
  surv_data$group     <- strain
  surv_data$replicate <- col_name
  all_individuals[[col_name]] <- surv_data
}

final_data       <- bind_rows(all_individuals)
if (nrow(final_data) == 0) stop("No survival data generated. Check count columns.")
final_data$group <- factor(final_data$group)

cat("\nTotal individuals:", nrow(final_data), "\n")
print(table(final_data$group, final_data$event))

# ----------------------------
# 4. Survival analysis: KM, log-rank, FH(0,1)
# ----------------------------
fit    <- survfit(Surv(time, event) ~ group, data = final_data)

logrank_all <- survdiff(Surv(time, event) ~ group, data = final_data)
p_logrank   <- 1 - pchisq(logrank_all$chisq, length(logrank_all$n) - 1)
cat("\nOverall log-rank p-value:", p_logrank, "\n")

global_fh <- logrank_test(Surv(time, event) ~ group, data = final_data,
                          type = "Fleming-Harrington", rho = 0, gamma = 1)
stat_fh <- as.numeric(global_fh@statistic@teststatistic)
df_fh   <- length(levels(final_data$group)) - 1
p_fh    <- pchisq(stat_fh, df = df_fh, lower.tail = FALSE)
cat("Global Fleming-Harrington(0,1) p-value:", p_fh, "\n")

# ----------------------------
# 5. Pairwise comparisons (Fleming-Harrington + FDR correction)
# ----------------------------
groups <- levels(final_data$group)
n_grp  <- length(groups)

raw_p <- matrix(NA, n_grp, n_grp, dimnames = list(groups, groups))
fdr_p <- raw_p

# --- FIXED safe_fh_test: removes early variance check ---
safe_fh_test <- function(data_pair) {
  g <- unique(data_pair$group)
  if (length(g) != 2) return(NA_real_)
  
  # Try Fleming-Harrington test first
  tryCatch({
    test <- logrank_test(
      Surv(time, event) ~ group,
      data = data_pair,
      type = "Fleming-Harrington",
      rho = 0,
      gamma = 1
    )
    stat_p <- as.numeric(test@statistic@teststatistic)
    # If test statistic is 0 or NA, p = 1 (identical curves)
    if (is.na(stat_p) || stat_p == 0) return(1.0)
    return(pchisq(stat_p^2, df = 1, lower.tail = FALSE))
  }, error = function(e) {
    # Fallback to standard log-rank (handles zero-variance groups correctly)
    sd <- survdiff(Surv(time, event) ~ group, data = data_pair)
    return(1 - pchisq(sd$chisq, df = 1))
  })
}

for (i in 1:(n_grp - 1)) {
  for (j in (i + 1):n_grp) {
    pair_data <- final_data %>% filter(group %in% c(groups[i], groups[j]))
    p_val     <- safe_fh_test(pair_data)
    raw_p[i, j] <- raw_p[j, i] <- p_val
  }
}
diag(raw_p) <- 1

upper_vals <- raw_p[upper.tri(raw_p)]
non_na_idx <- which(!is.na(upper_vals))
fdr_vals_temp <- rep(NA_real_, length(upper_vals))
if (length(non_na_idx) > 0) {
  fdr_vals_temp[non_na_idx] <- p.adjust(upper_vals[non_na_idx], method = "fdr")
}
fdr_p[upper.tri(fdr_p)] <- fdr_vals_temp
fdr_p[lower.tri(fdr_p)] <- t(fdr_p)[lower.tri(fdr_p)]
diag(fdr_p) <- 1

cat("\nPairwise FDR matrix created.\n")

sig_pairs <- which(fdr_p < alpha & upper.tri(fdr_p), arr.ind = TRUE)
if (nrow(sig_pairs) > 0) {
  cat("\nSignificantly different strain pairs (FDR <", alpha, "):\n")
  for (i in 1:nrow(sig_pairs)) {
    s1 <- rownames(fdr_p)[sig_pairs[i, 1]]
    s2 <- colnames(fdr_p)[sig_pairs[i, 2]]
    pv <- fdr_p[sig_pairs[i, 1], sig_pairs[i, 2]]
    cat(s1, "vs", s2, "-> p =", format(pv, scientific = TRUE), "\n")
  }
} else {
  cat("No significant differences at FDR <", alpha, "\n")
}
cat("\nNumber of significantly different pairs:", sum(fdr_p[upper.tri(fdr_p)] < alpha, na.rm = TRUE), "\n")

# ----------------------------
# 6. Heatmap of pairwise FDR-adjusted p-values
# ----------------------------
fdr_p[is.na(fdr_p)] <- 1
p_matrix_long <- melt(fdr_p, varnames = c("Strain1", "Strain2"), value.name = "p_value")
p_matrix_long <- p_matrix_long[p_matrix_long$Strain1 != p_matrix_long$Strain2, ]

p_heat <- ggplot(p_matrix_long, aes(x = Strain1, y = Strain2, fill = p_value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "log10",
    name = "FDR p-value",
    direction = -1
  ) +
  labs(title = "Pairwise strain comparisons", subtitle = "(Fleming-Harrington test)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.x = element_blank(),   # removes "Strain1"
    axis.title.y = element_blank(),   # removes "Strain2"
    panel.grid = element_blank()
  ) +
  coord_fixed()

print(p_heat)
ggsave("pairwise_heatmap_a.baumannii_BR5.png", plot = p_heat, width = 8, height = 6, dpi = 600)


###
# Replace NA with non-significant
fdr_p[is.na(fdr_p)] <- 1

# Convert to long format
p_matrix_long <- melt(
  fdr_p,
  varnames = c("Strain1", "Strain2"),
  value.name = "p_value"
)

# Optional: remove diagonal
p_matrix_long <- subset(p_matrix_long, Strain1 != Strain2)

# Categorize p-values
p_matrix_long$significance <- cut(
  p_matrix_long$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("p < 0.001", "p < 0.01", "p < 0.05", "NS")
)

# Set legend order
p_matrix_long$significance <- factor(
  p_matrix_long$significance,
  levels = c("p < 0.001", "p < 0.01", "p < 0.05", "NS")
)

# Plot
p_heat <- ggplot(
  p_matrix_long,
  aes(x = Strain1, y = Strain2, fill = significance)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "p < 0.001" = "#5A009D",  # purple
      "p < 0.01"  = "#E72968",  # pink
      "p < 0.05"  = "#F4A12F",  # orange
      "NS"        = "#E6E6E6"   # light gray
    ),
    name = "FDR p-value"
  ) +
  coord_fixed() +
  labs(
    title = "Pairwise strain comparisons",
    subtitle = "(Fleming-Harrington test, FDR corrected)",
    x = NULL,
    y = NULL
  ) +
theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      size = 11
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(0.6, "cm"),
    legend.key.height = unit(0.6, "cm")
  )

print(p_heat)

ggsave(
  "pairwise_heatmap_a.baumannii_BR5.png",
  p_heat,
  width = 9,
  height = 7,
  dpi = 700
)

# Color options
# magma
# plasma 
# cividis
# mako 
# inferno 
# viridis
# rocket
# turbo

# ----------------------------
# 7. Survival plots (export)
# ----------------------------

# Optional - selected strains
library(survival)
library(survminer)
library(dplyr)

# Strains to plot
selected <- c("PBS", "BR5","ATCC_19606")

# Keep only strains present in the dataset
selected <- intersect(selected, levels(final_data$group))

if (length(selected) >= 2) {
  
  # Manual colors
  strain_colors <- c(
    "PBS"        = "#6E7B8B",
    "BR5"        = "#CD0000",
    "ATCC_19606" = "#00B2EE"
  )
  
  # Subset data
  sub_data <- final_data %>%
    filter(group %in% selected) %>%
    mutate(group = factor(group, levels = selected))
  
  # Survival fit
  sub_fit <- survfit(
    Surv(time, event) ~ group,
    data = sub_data
  )
  
  # Plot with x‑axis limited to 36 hours
  # Plot with x‑axis limited to 36 hours
  p_sub <- ggsurvplot(
    fit = sub_fit,
    data = sub_data,
    pval = FALSE,
    pval.method = FALSE,
    conf.int = FALSE,
    risk.table = FALSE,
    break.time.by = 6,           # Ticks at 0, 6, 12, 18, 24, 30, 36
    xlim = c(0, 36),             # <-- ADD THIS LINE
    xlab = "Time (hours)",
    ylab = "Survival (%)",
    title = "gm_a.baumannii",
    legend.title = "Strain",
    legend = "right",
    palette = unname(strain_colors[selected]),
    size = 1,
    surv.scale = "percent"
  )
  
  print(p_sub)
  
  ggsave(
    filename = "gm_a.baumannii_BR5.png",
    plot = p_sub$plot,
    width = 7,
    height = 5,
    dpi = 700
  )
}

# ----------------------------
# 8. Export statistical results for supplementary material
# ----------------------------

library(writexl)

# Export the FDR-adjusted p-value matrix
# Replace NA with 1 (non-significant) for clarity
fdr_p_export <- fdr_p
fdr_p_export[is.na(fdr_p_export)] <- 1

# Convert to data frame (keeping row/col names)
df_export <- as.data.frame(fdr_p_export)

# Add a column for strain names so it's easier to read in Excel
df_export <- cbind(Strain = rownames(df_export), df_export)

# Write to Excel
write_xlsx(df_export, "Supplementary_Table_S8_Pairwise_pvalues_a.baumannii_BR5.xlsx")

cat("\n✅ Supplementary Table S8 exported: Supplementary_Table_S8_Pairwise_pvalues_a.baumannii_BR5.xlsx\n")