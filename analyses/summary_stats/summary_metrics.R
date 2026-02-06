library(ggplot2)
library(patchwork)

# Read data
data <- read.table(file.path(getwd(), "analyses/summary_stats/qc_metrics_summary.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Add a 'Platform' column based on the 'sample_name'
data$Platform <- ifelse(grepl("PB", data$sample_name), "PacBio", "ONT")

# Ensure 'Platform' is treated as a factor
data$Platform <- factor(data$Platform, levels = c("ONT", "PacBio"))

# Print summary statistics for each metric
print(summary(data$Yield))
print(summary(data$Mean.Coverage))
print(summary(data$N50))

# Create individual boxplots with filled jitter points, bold titles, and platform on the x-axis
yield_plot <- ggplot(data, aes(x = Platform, y = Yield, fill = Platform)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "black") +
  scale_fill_manual(values = c("#2D4F67", "#1A2F3E")) +
  geom_jitter(shape = 21, fill = "#1A2F3E", color = "black", size = 2, alpha = 0.9, width = 0.4,height = 0) +
  labs(title = "Yield [Gb]", y = NULL, x = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

coverage_plot <- ggplot(data, aes(x = Platform, y = Mean.Coverage, fill = Platform)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "black") +
  scale_fill_manual(values = c("#76946A", "#4E703E")) +
  geom_jitter(shape = 21, fill = "#4E703E", color = "black", size = 2, alpha = 0.9, width = 0.4, height = 0) +
  labs(title = "Mean Coverage", y = NULL, x = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

n50_plot <- ggplot(data, aes(x = Platform, y = N50, fill = Platform)) +
  geom_boxplot(alpha = 0.7, outlier.colour = "black") +
  scale_fill_manual(values = c("#C34043", "#7A2121")) +
  geom_jitter(shape = 21, fill = "#7A2121", color = "black", size = 2, alpha = 0.9, width = 0.4, height = 0) +
  labs(title = "N50", y = NULL, x = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.ticks.y = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Combine the plots with patchwork
final_plot <- (yield_plot | coverage_plot | n50_plot) +
  plot_layout(guides = "collect")

# Save plot
ggsave("figures/summary_boxplots_grouped.svg", plot = final_plot, width = 210, height = 100, units = "mm")

# Display the combined plot (optional)
print(final_plot)
