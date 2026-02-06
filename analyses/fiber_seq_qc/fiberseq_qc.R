# Load necessary libraries
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

# Define statistics to include in analysis
statistics_to_plot <- c("fiber_length", "cpg_count", "m6a_count", "nuc_count")

# Calculate summary statistics BEFORE capping for individual samples
summary_stats_individual <- all_data %>%
  filter(statistic %in% statistics_to_plot) %>%
  mutate(value_numeric = as.numeric(value)) %>%
  group_by(statistic, sample) %>%
  summarise(
    total_count = sum(count),
    mean = weighted.mean(value_numeric, count),
    median = median(rep(value_numeric, count)),
    std_dev = sqrt(sum(count * (value_numeric - weighted.mean(value_numeric, count))^2) / sum(count)),
    min = min(value_numeric),
    max = max(value_numeric),
    .groups = 'drop'
  )

# Print median values to the console
cat("Median Values (BEFORE Capping):\n")
print(summary_stats_individual %>%
        select(statistic, sample, median))

# Define cutoffs for each statistic
cutoffs <- list(
  fiber_length = 60000,  # 60 kb
  cpg_count = 1000,      # 1k
  m6a_count = 4000,      # 4k
  nuc_count = 500        # 500
)

# Apply cutoffs and aggregate counts
capped_data <- all_data %>%
  filter(statistic %in% statistics_to_plot) %>%
  mutate(
    value_numeric = as.numeric(value),
    value_capped = case_when(
      statistic == "fiber_length" & value_numeric > cutoffs$fiber_length ~ cutoffs$fiber_length,
      statistic == "cpg_count" & value_numeric > cutoffs$cpg_count ~ cutoffs$cpg_count,
      statistic == "m6a_count" & value_numeric > cutoffs$m6a_count ~ cutoffs$m6a_count,
      statistic == "nuc_count" & value_numeric > cutoffs$nuc_count ~ cutoffs$nuc_count,
      TRUE ~ value_numeric
    )
  ) %>%
  # Group by statistic, sample, and capped value to sum counts
  group_by(statistic, sample, value_capped) %>%
  summarise(count = sum(count), .groups = 'drop')

# Plot settings
facet_titles <- c(
  "fiber_length" = " Length",
  "cpg_count" = "5mC Count",
  "m6a_count" = "6mA Count",
  "nuc_count" = "Nucleosome Count"
)

custom_sample_colors <- c("PacBio" = "#C34043", "ONT" = "#658594", "Merged" = "#7A6FCA")

# Generate the density plot
density_plot <- ggplot(capped_data, aes(x = value_capped, fill = sample)) +
  # Density plots
  geom_density(
    aes(weight = count),
    alpha = 0.6, position = "identity"
  ) +
  # Add dashed vertical lines for median values
  geom_vline(
    data = summary_stats_individual %>% filter(sample %in% c("PacBio", "ONT")),
    aes(xintercept = median, color = sample),
    linetype = "dashed", size = 0.5
  ) +
  scale_fill_manual(values = custom_sample_colors) +
  scale_color_manual(values = custom_sample_colors) +  # Match dashed lines to sample colors
  facet_wrap(
    ~ factor(statistic, levels = statistics_to_plot), 
    labeller = as_labeller(facet_titles), 
    scales = "free"
  ) +
  labs(x = "Value", y = "Density", fill = "Sample", color = "Sample") +
  theme(
    text = element_text(family = "sans", size = 14),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey97", color = "black"),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

# Save the plot
ggsave("figures/fiber_seq_qc_with_medians.svg", plot = density_plot, width = 210, height = 100, units = "mm")
print(density_plot)