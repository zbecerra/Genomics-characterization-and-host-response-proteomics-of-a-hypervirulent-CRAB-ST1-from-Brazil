# =============================================================================
# Recombination hotspot figure - A. baumannii ST1 core-gene alignment
# Requires: recomb_density_by_window_v2.tsv (download from your cluster first)
# Packages: ggplot2, ggrepel (for non-overlapping labels)
# =============================================================================

# install.packages(c("ggplot2", "ggrepel"))  # run once if not already installed
library(ggplot2)
library(ggrepel)

# -----------------------------------------------------------------------------
# 1. Load density data
#    Download first: scp user@server:path/to/recomb_density_by_window_v2.tsv .
#    Then set your working directory / path below
# -----------------------------------------------------------------------------
density <- read.delim("recomb_density_by_window_v2.tsv", sep = "\t")
density$Midpoint <- (density$Window_Start + density$Window_End) / 2

# -----------------------------------------------------------------------------
# 2. Hotspot annotations (from today's gene-lookup analysis)
#    category: "AMR", "Virulence", or "Other"
# -----------------------------------------------------------------------------
hotspots <- data.frame(
  Midpoint = c(830000, 810000, 790000, 1190000, 570000,
               1650000, 590000, 1670000, 2050000, 2070000,
               2010000, 1930000, 1830000),
  Gene     = c("pilY1", "gspD", "pilC", "lpxH", "lpxP/lpxL",
               "adeG", "pgaA/pgaB", "sul1", "lpxK", "lpxO / abaF",
               "lpxC", "barB", "gsp"),
  Events   = c(296, 119, 68, 67, 61,
               40, 38, 38, 36, 32,
               33, 31, 28),
  Category = c("Virulence", "Virulence", "Virulence", "Virulence", "Virulence",
               "AMR", "Virulence", "AMR", "Virulence", "AMR",
               "Virulence", "Virulence", "Virulence")
)

# -----------------------------------------------------------------------------
# 3. Color scheme (matches your existing figure palette style)
# -----------------------------------------------------------------------------
cat_colors <- c("AMR" = "#D6604D", "Virulence" = "#4393C3", "Other" = "#999999")

# -----------------------------------------------------------------------------
# 4. Build the figure
#    All gene labels sit in one fixed row above the plot, connected to their
#    actual peak by a straight vertical dotted line. Labels only shift
#    horizontally (never vertically) to avoid overlapping each other.
# -----------------------------------------------------------------------------
label_y <- max(density$Event_Count) * 1.18   # fixed height for every label

p <- ggplot(density, aes(x = Midpoint, y = Event_Count)) +
  # density area across the whole alignment
  geom_area(fill = "grey85", color = "grey50", linewidth = 0.3) +

  # straight vertical connector from each real peak up to the label row
  geom_segment(data = hotspots,
               aes(x = Midpoint, xend = Midpoint, y = Events, yend = label_y,
                   color = Category),
               linewidth = 0.4, linetype = "dotted", show.legend = FALSE) +

  # marker at the real peak
  geom_point(data = hotspots, aes(x = Midpoint, y = Events, color = Category),
             size = 2.2) +

  # gene labels, all at the same fixed height, nudged only horizontally
  geom_text_repel(data = hotspots,
                   aes(x = Midpoint, y = label_y, label = Gene, color = Category),
                   size = 3.2, fontface = "italic",
                   direction = "x",           # horizontal movement only
                   segment.color = NA,        # suppress ggrepel's own leader line
                   box.padding = 0.3,
                   max.overlaps = 20,
                   show.legend = FALSE) +

  scale_color_manual(values = cat_colors, name = "Category") +
  scale_x_continuous(labels = scales::comma, expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, label_y * 1.08)) +   # room for the label row
  coord_cartesian(xlim = c(400000, max(density$Window_End))) +  # crop empty left margin

  labs(
    title = "Recombination hotspots across the A. baumannii ST1 core-gene alignment",
    subtitle = "Position reflects Panaroo core-gene concatenation order, not physical chromosomal position",
    x = "Position in core-gene alignment (bp)",
    y = "Recombination events (Gubbins/ClonalFrameML, n = 897 genomes)"
  ) +

  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(size = 9, color = "grey40", face = "italic"),
    legend.position = "top"
  )

print(p)

# -----------------------------------------------------------------------------
# 5. Save
# -----------------------------------------------------------------------------
ggsave("Figure_recombination_hotspots.png", p, width = 10, height = 5, dpi = 600)
ggsave("Figure_recombination_hotspots.pdf", p, width = 10, height = 5)
