library(ggplot2)
library(dplyr)
library(tidyr)

# Working directory and data file paths
base_dir <- file.path(getwd(), "analyses", "fiber_seq_qc")
file_paths <- list(
  PacBio = file.path(base_dir, "CCS15-PB_qc.txt"),
  ONT = file.path(base_dir, "CCS15-ONT_qc.txt")
)

# Combine data with a sample column
all_data <- lapply(names(file_paths), function(sample) {
  read.table(file_paths[[sample]], header = TRUE, sep = "\t") %>% 
    mutate(sample = sample)
}) %>% bind_rows()

# Remove rows where `value` column contains "UNK"
all_data <- all_data %>%
  filter(value != "UNK")

# Filter for nucleosome length statistic only
nuc_length_data_raw <- all_data %>%
  filter(statistic == "nuc_length") %>%
  mutate(value_numeric = as.numeric(value))

# Calculate summary statistics BEFORE capping for individual samples
summary_stats_individual <- nuc_length_data_raw %>%
  group_by(sample) %>%
  summarise(
    total_count = sum(count),
    mean = weighted.mean(value_numeric, count),
    median = median(rep(value_numeric, count)),
    std_dev = sqrt(sum(count * (value_numeric - weighted.mean(value_numeric, count))^2) / sum(count)),
    min = min(value_numeric),
    max = max(value_numeric),
    .groups = 'drop'
  )

# Calculate summary statistics BEFORE capping for merged samples
summary_stats_merged <- nuc_length_data_raw %>%
  summarise(
    sample = "Merged",
    total_count = sum(count),
    mean = weighted.mean(value_numeric, count),
    median = median(rep(value_numeric, count)),
    std_dev = sqrt(sum(count * (value_numeric - weighted.mean(value_numeric, count))^2) / sum(count)),
    min = min(value_numeric),
    max = max(value_numeric)
  )

# Combine individual and merged statistics
summary_stats_all <- bind_rows(summary_stats_individual, summary_stats_merged)

# Print summary statistics
print("Summary Statistics for Nucleosome Length BEFORE Capping:")
print(summary_stats_all, n = Inf)

# Extract median values for ONT and PacBio
median_pacbio <- summary_stats_individual %>% filter(sample == "PacBio") %>% pull(median)
median_ont <- summary_stats_individual %>% filter(sample == "ONT") %>% pull(median)

# Apply 500 bp cutoff and aggregate counts
nuc_length_data <- nuc_length_data_raw %>%
  mutate(value_capped = ifelse(value_numeric > 500, 500, value_numeric)) %>%
  # Group by sample and capped value to sum counts
  group_by(sample, value_capped) %>%
  summarise(count = sum(count), .groups = 'drop')

# Custom sample colors
custom_sample_colors <- c("PacBio" = "#C34043", "ONT" = "#658594")

# Create density plot
density_plot <- ggplot(nuc_length_data, aes(x = value_capped, weight = count, fill = sample)) +
  geom_density(alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 147, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = median_pacbio, linetype = "dashed", color = "#C34043", linewidth = 0.5) +
  geom_vline(xintercept = median_ont, linetype = "dashed", color = "#658594", linewidth = 0.5) +
  scale_fill_manual(values = custom_sample_colors) +
  labs(x = "Nucleosome Length (bp)", y = "Density", fill = "Sample", 
       title = "Nucleosome Length") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(
    text = element_text(family = "sans", size = 14),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave("figures/nucleosome_length_density.svg", plot = density_plot, width = 210, height = 100, units = "mm")
print(density_plot)

