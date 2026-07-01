library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

pdf(NULL)  # suppress Rplots.pdf side-effect

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(script_dir)

summary_tsv <- "qc_metrics_summary.tsv"
file_paths  <- list(
  PacBio = "CCS15-PB_qc.txt",
  ONT    = "CCS15-ONT_qc.txt"
)

# ── Summary metrics data ─────────────────────────────────────────────────────

summary_data <- read.table(summary_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
summary_data$Platform <- ifelse(grepl("PB", summary_data$sample_name), "PacBio", "ONT")
summary_data$Platform <- factor(summary_data$Platform, levels = c("ONT", "PacBio"))

# ── FiberSeq QC data ─────────────────────────────────────────────────────────

all_data <- lapply(names(file_paths), function(sample) {
  read.table(file_paths[[sample]], header = TRUE, sep = "\t") %>%
    mutate(sample = sample)
}) %>% bind_rows() %>%
  filter(value != "UNK")

statistics_to_plot <- c("fiber_length", "cpg_count", "m6a_count", "nuc_count", "nuc_length")

# Summary stats (pre-cap) for median lines
summary_stats_individual <- all_data %>%
  filter(statistic %in% statistics_to_plot) %>%
  mutate(value_numeric = as.numeric(value)) %>%
  group_by(statistic, sample) %>%
  summarise(
    total_count = sum(count),
    mean        = weighted.mean(value_numeric, count),
    median      = median(rep(value_numeric, count)),
    std_dev     = sqrt(sum(count * (value_numeric - weighted.mean(value_numeric, count))^2) / sum(count)),
    min         = min(value_numeric),
    max         = max(value_numeric),
    .groups     = "drop"
  )

# Apply display cutoffs and aggregate
cutoffs <- list(
  fiber_length = 60000,
  cpg_count    = 1000,
  m6a_count    = 4000,
  nuc_count    = 500,
  nuc_length   = 500
)

capped_data <- all_data %>%
  filter(statistic %in% statistics_to_plot) %>%
  mutate(
    value_numeric = as.numeric(value),
    value_capped  = case_when(
      statistic == "fiber_length" & value_numeric > cutoffs$fiber_length ~ cutoffs$fiber_length,
      statistic == "cpg_count"    & value_numeric > cutoffs$cpg_count    ~ cutoffs$cpg_count,
      statistic == "m6a_count"    & value_numeric > cutoffs$m6a_count    ~ cutoffs$m6a_count,
      statistic == "nuc_count"    & value_numeric > cutoffs$nuc_count    ~ cutoffs$nuc_count,
      statistic == "nuc_length"   & value_numeric > cutoffs$nuc_length   ~ cutoffs$nuc_length,
      TRUE ~ value_numeric
    )
  ) %>%
  group_by(statistic, sample, value_capped) %>%
  summarise(count = sum(count), .groups = "drop")

# ── Shared theme helpers ──────────────────────────────────────────────────────

box_theme <- theme_minimal() +
  theme(
    text          = element_text(family = "sans", size = 14),
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text     = element_text(size = 10, color = "black"),
    axis.ticks.y  = element_line(color = "black"),
    panel.grid    = element_blank(),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  )

custom_sample_colors <- c("PacBio" = "#CC6677", "ONT" = "#0077BB")

# ── Summary boxplots (top row) ────────────────────────────────────────────────

yield_plot <- ggplot(summary_data, aes(x = Platform, y = Yield, fill = Platform)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "black") +
  scale_fill_manual(values = c(ONT = "#332288", PacBio = "#332288")) +
  geom_jitter(shape = 21, fill = "#332288", color = "black",
              size = 2, alpha = 0.9, width = 0.4, height = 0) +
  labs(title = "Yield (Gb)", y = NULL, x = "") +
  box_theme

coverage_plot <- ggplot(summary_data, aes(x = Platform, y = Mean.Coverage, fill = Platform)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "black") +
  scale_fill_manual(values = c(ONT = "#117733", PacBio = "#117733")) +
  geom_jitter(shape = 21, fill = "#117733", color = "black",
              size = 2, alpha = 0.9, width = 0.4, height = 0) +
  labs(title = "Mean Coverage", y = NULL, x = "") +
  box_theme

n50_plot <- ggplot(summary_data, aes(x = Platform, y = N50, fill = Platform)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "black") +
  scale_fill_manual(values = c(ONT = "#882255", PacBio = "#882255")) +
  geom_jitter(shape = 21, fill = "#882255", color = "black",
              size = 2, alpha = 0.9, width = 0.4, height = 0) +
  labs(title = "N50 (bp)", y = NULL, x = "") +
  box_theme

# ── Density plots (bottom row) ────────────────────────────────────────────────

facet_titles <- c(
  "fiber_length" = "Fiber Length",
  "cpg_count"    = "5mC Count",
  "m6a_count"    = "6mA Count",
  "nuc_count"    = "Nucleosome Count",
  "nuc_length"   = "Nucleosome Length"
)

# ── Median label adjustment variables ─────────────────────────────────────────
# Tweak these to reposition / resize the stacked top-right median labels

median_label_size    <- 3.2  # font size (ggplot units ≈ pt × 0.35)
median_label_vjust_1 <- 1.3  # distance down from top for label 1 (ONT)  — larger = lower
median_label_vjust_2 <- 2.8  # distance down from top for label 2 (PacBio) — larger = lower
median_label_hjust   <- 1.05 # horizontal alignment: 1 = flush right edge, >1 = slight inset

# 147 bp reference line label (nuc_length panel only) — stacked with median labels
nuc147_label_size  <- 3.2  # font size for "147 bp" label
nuc147_label_vjust <- 4.3  # distance down from top of panel (below the two median labels)

# ── Build median label data ───────────────────────────────────────────────────

stat_units <- c(
  fiber_length = "bp",
  cpg_count    = "bases",
  m6a_count    = "bases",
  nuc_count    = "nucleosomes",
  nuc_length   = "bp"
)

# ONT on top, PacBio below; vjust_val controls vertical stacking
median_label_data <- summary_stats_individual %>%
  filter(statistic %in% statistics_to_plot, sample %in% c("PacBio", "ONT")) %>%
  mutate(
    unit       = stat_units[statistic],
    label_text = paste0(formatC(round(median), format = "d", big.mark = ","), " ", unit),
    vjust_val  = ifelse(sample == "ONT", median_label_vjust_1, median_label_vjust_2)
  )

# 147 bp reference line + label (filtered to nuc_length facet)
nuc147_vline <- data.frame(statistic = "nuc_length", xintercept = 147)
nuc147_label <- data.frame(statistic = "nuc_length", x = 147, label = "147 bp")

# ── Density plot ──────────────────────────────────────────────────────────────

density_plot <- ggplot(capped_data, aes(x = value_capped, fill = sample)) +
  geom_density(aes(weight = count), alpha = 0.6, position = "identity") +
  # Median dashed lines
  geom_vline(
    data     = summary_stats_individual %>% filter(sample %in% c("PacBio", "ONT")),
    aes(xintercept = median, color = sample),
    linetype = "dashed", linewidth = 0.5
  ) +
  # Median value labels — stacked in top-right corner of each panel
  geom_text(
    data         = median_label_data,
    aes(x = Inf, y = Inf, label = label_text, color = sample, vjust = vjust_val),
    hjust        = median_label_hjust,
    size         = median_label_size,
    fontface     = "bold",
    inherit.aes  = FALSE,
    show.legend  = FALSE
  ) +
  scale_fill_manual(values  = custom_sample_colors) +
  scale_color_manual(values = custom_sample_colors) +
  facet_wrap(
    ~ factor(statistic, levels = statistics_to_plot),
    labeller = as_labeller(facet_titles),
    scales   = "free",
    ncol     = 3
  ) +
  labs(x = "Value", y = "Density", fill = "Sample", color = "Sample") +
  theme(
    text             = element_text(family = "sans", size = 14),
    axis.text.x      = element_text(size = 10, color = "black"),
    axis.text.y      = element_text(size = 10, color = "black"),
    axis.title       = element_text(face = "bold", size = 12),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey97", color = "black"),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid       = element_blank(),
    legend.position  = "bottom"
  )

# ── Combine with patchwork ────────────────────────────────────────────────────

top_row <- (
  (yield_plot + labs(tag = "A")) | coverage_plot | n50_plot
) + plot_layout(widths = c(1, 1, 1))

combined_plot <- top_row / (density_plot + labs(tag = "B")) +
  plot_layout(heights = c(1, 1.6)) &
  theme(plot.tag = element_text(face = "bold", size = 14))

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave("combined_qc_summary.svg",
       plot = combined_plot, width = 210, height = 250, units = "mm")

ggsave("combined_qc_summary.png",
       plot = combined_plot, width = 210, height = 250, units = "mm")

# ── Report ────────────────────────────────────────────────────────────────────

fmt_int  <- function(x) formatC(round(x), format = "d", big.mark = ",")
fmt_dec2 <- function(x) sprintf("%.2f", x)

yield_med   <- median(summary_data$Yield)
yield_range <- range(summary_data$Yield)
cov_med     <- median(summary_data$Mean.Coverage)
cov_range   <- range(summary_data$Mean.Coverage)
n50_med     <- median(summary_data$N50)
n50_range   <- range(summary_data$N50)

get_med <- function(stat, samp) {
  summary_stats_individual$median[
    summary_stats_individual$statistic == stat &
    summary_stats_individual$sample    == samp
  ]
}

pb_len  <- get_med("fiber_length", "PacBio")
pb_m6a  <- get_med("m6a_count",    "PacBio")
pb_cpg  <- get_med("cpg_count",    "PacBio")
pb_nuc  <- get_med("nuc_count",    "PacBio")
ont_len <- get_med("fiber_length", "ONT")
ont_m6a <- get_med("m6a_count",    "ONT")
ont_cpg <- get_med("cpg_count",    "ONT")
ont_nuc <- get_med("nuc_count",    "ONT")

writeLines(c(
  "== Per-sample aligned QC (Fig. 2A; from qc_metrics_summary.tsv) ==",
  sprintf("N samples: %d", nrow(summary_data)),
  sprintf("Median aligned data yield: %s Gb (range: [%s; %s])",
          fmt_dec2(yield_med), fmt_dec2(yield_range[1]), fmt_dec2(yield_range[2])),
  sprintf("Median mean coverage: %sX (range: [%s; %s])",
          fmt_dec2(cov_med), fmt_dec2(cov_range[1]), fmt_dec2(cov_range[2])),
  sprintf("Median N50 read length: %s bp (range: [%s; %s])",
          fmt_int(n50_med), fmt_int(n50_range[1]), fmt_int(n50_range[2])),
  "",
  "== Per-read fiber-seq medians (Fig. 2B; CCS15 only) ==",
  sprintf("CCS15 PacBio: median read length %s bp; per read median %s 6mA and %s 5mC calls; %s nucleosomes",
          fmt_int(pb_len),  fmt_int(pb_m6a),  fmt_int(pb_cpg),  fmt_int(pb_nuc)),
  sprintf("CCS15 ONT:    median read length %s bp; per read median %s 6mA and %s 5mC calls; %s nucleosomes",
          fmt_int(ont_len), fmt_int(ont_m6a), fmt_int(ont_cpg), fmt_int(ont_nuc)),
  ""
), "report.txt")
