# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read data for COLO829 Structural Variants into a dataframe
data_COLO829 <- read.table(header = TRUE, text = "
Caller Precision Recall F1_Score
'ONT Severus Paired' 0.68 0.8225806 0.7445255
'PacBio Severus Paired' 0.6575342 0.7741935 0.7111111
'ONT Severus TO' 0.117506 0.7903226 0.2045929
'PacBio Severus TO' 0.1753731 0.7580645 0.2848485
'ONT Severus Benchmark' 0.720 0.868 0.787
'PacBio Severus Benchmark' 0.787 0.868 0.825
")

# Set factor levels for Caller to ensure the desired order
data_COLO829$Caller <- factor(
  data_COLO829$Caller,
  levels = c("ONT Severus Benchmark","ONT Severus Paired", "PacBio Severus Benchmark", "PacBio Severus Paired", "ONT Severus TO", 
             "PacBio Severus TO")
)

# Reshape the data
data_COLO829_long <- data_COLO829 %>%
  pivot_longer(
    cols = c(Precision, Recall, F1_Score),
    names_to = "Metric",
    values_to = "Value"
  )

# Create plot for COLO829
p2 <- ggplot(data_COLO829_long, aes(x = Metric, y = Value, fill = Caller)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", linewidth = 0.25, width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", Value)),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3.5
  ) +
  labs(
    title = "COLO829 Structural Variants",
    x = "Metric",
    y = "Value"
  ) +
  scale_fill_manual(
    values = c(
      "ONT Severus Paired" = "#9CABCA",
      "PacBio Severus Paired" = "#FFA066",
      "ONT Severus TO" = "#76946A",
      "PacBio Severus TO" = "#D36E70",
      "ONT Severus Benchmark" = "#DCA561",
      "PacBio Severus Benchmark" = "#2D4F67"
    ),
    name = "Caller"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 1, vjust = 1),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "top"
  )

# Read data for HG008 Structural Variants into a dataframe
data_HG008 <- read.table(header = TRUE, text = "
Caller Precision Recall F1_Score
'ONT Severus Paired' 0.6296296296296297 0.8947368421052632 0.7391304347826088
'PacBio Severus Paired' 0.9104477611940298 0.9172932330827067 0.9138576779026217
'ONT Severus TO' 0.3034300791556728 0.8646616541353384 0.44921874999999994
'PacBio Severus TO' 0.4057971014492754 0.8421052631578947 0.5476772616136919
'ONT Severus Benchmark' 0.72 0.91 0.8043478
'PacBio Severus Benchmark' 0.93 0.94 0.9349593
")

# Set factor levels for Caller to ensure the desired order
data_HG008$Caller <- factor(
  data_HG008$Caller,
  levels = c("ONT Severus Benchmark","ONT Severus Paired", "PacBio Severus Benchmark", "PacBio Severus Paired", "ONT Severus TO", 
             "PacBio Severus TO")
)

# Reshape the data
data_HG008_long <- data_HG008 %>%
  pivot_longer(
    cols = c(Precision, Recall, F1_Score),
    names_to = "Metric",
    values_to = "Value"
  )

# Create plot for HG008
p1 <- ggplot(data_HG008_long, aes(x = Metric, y = Value, fill = Caller)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", linewidth = 0.25, width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", Value)),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3.5
  ) +
  labs(
    title = "HG008 Structural Variants",
    x = "Metric",
    y = "Value"
  ) +
  scale_fill_manual(
    values = c(
      "ONT Severus Paired" = "#9CABCA",
      "PacBio Severus Paired" = "#FFA066",
      "ONT Severus TO" = "#76946A",
      "PacBio Severus TO" = "#D36E70",
      "ONT Severus Benchmark" = "#DCA561",
      "PacBio Severus Benchmark" = "#2D4F67"
    ),
    name = "Caller"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(hjust = 1, vjust = 1),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "top"
  )

ggsave("HG008_SVs.svg", plot = p1, width = 210, height = 100, units = "mm")
ggsave("COLO829_SVs.svg", plot = p2, width = 210, height = 100, units = "mm")

