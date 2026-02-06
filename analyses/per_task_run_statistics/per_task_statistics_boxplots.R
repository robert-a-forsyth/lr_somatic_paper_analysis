library(ggplot2)
library(dplyr)
library(patchwork) 
setwd("/Users/u0155044/Documents/lr_somatic_paper_analysis")
data <- read.table(file.path(getwd(), "improved_report.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# get task names without sample names in parentheses
data$Task <- gsub(".*:", "", gsub(" \\(.*\\)$", "", data$name))
data$Task <- factor(data$Task, levels = unique(data$Task))

# Function to convert duration to seconds
convert_to_seconds <- function(duration) {
  h <- as.numeric(sub(".*?(\\d+(\\.\\d+)?)h.*", "\\1", duration))
  m <- as.numeric(sub(".*?(\\d+(\\.\\d+)?)m.*", "\\1", duration))
  s <- as.numeric(sub(".*?(\\d+(\\.\\d+)?)s.*", "\\1", duration))
  
  h[is.na(h)] <- 0
  m[is.na(m)] <- 0
  s[is.na(s)] <- 0
  
  return(h * 3600 + m * 60 + s)
}

data$DurationSeconds <- convert_to_seconds(data$duration)

data$PercentCPU <- as.numeric(sub("%", "", data$X.cpu))

# Convert peak_vmem to GB
convert_to_gb <- function(vmem) {
  mb <- as.numeric(sub(" MB", "", vmem))
  gb <- as.numeric(sub(" GB", "", vmem))
  mb[!is.na(mb)] <- mb[!is.na(mb)] / 1024  # Convert MB values to GB
  gb[is.na(mb)] <- gb[is.na(mb)]           # Keep GB values as-is
  return(ifelse(is.na(mb), gb, mb))        # Combine non-NA values
}

data$PeakVmemGB <- convert_to_gb(data$peak_vmem)

custom_theme <- theme(
  text = element_text(family = "sans", size = 14),   
  axis.text.x = element_blank(),                    
  axis.text.y = element_text(family = "sans", size = 12),
  axis.title = element_text(family = "sans", size = 12,face="bold"), 
  plot.title = element_text(hjust = 0.5, family = "sans", size = 12), 
  legend.text = element_text(size = 12),             
  legend.title = element_text(size = 12),            
  panel.background = element_rect(fill = NA, color = "black"), 
  panel.grid = element_blank()                       
)

data$DurationHours <- data$DurationSeconds / 3600  

data$PercentCPU <- as.numeric(data$PercentCPU)
data$cpus <- as.numeric(data$cpus)

data$PercentCPU[is.na(data$PercentCPU)] <- 0
data$cpus[is.na(data$cpus)] <- 1  

# Calculate the percentage of requested CPU used
data$PercentRequestedCPU <- (data$PercentCPU / 100) / data$cpus

# 1. Boxplot for DurationHours
p1 <- ggplot(data, aes(x = Task, y = DurationHours)) +
  geom_boxplot(fill = "#7FB4CA", color = "black") + # Black boxplot outlines
  labs(x = NULL, y = "Duration (Hours)") +          # Remove x-axis title
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 1)) + # Adjust limits based on your data
  custom_theme

# 2. Boxplot for PercentCPU
p2 <- ggplot(data, aes(x = Task, y = PercentCPU)) +
  geom_boxplot(fill = "#98BB6C", color = "black") +
  labs(x = NULL, y = "CPU Usage (%)") +             # Remove x-axis title
  scale_y_continuous(limits = c(0, 3000), breaks = seq(0, 3000, 1000)) + 
  custom_theme

# 3. Boxplot for PeakVmemGB
p3 <- ggplot(data, aes(x = Task, y = PeakVmemGB)) +
  geom_boxplot(fill = "#D27E99", color = "black") +
  labs(x = NULL, y = "Peak VMem (GB)") +             # Remove x-axis title
  custom_theme

# 4. Boxplot for PercentRequestedCPU
p4 <- ggplot(data, aes(x = Task, y = PercentRequestedCPU)) +
  geom_boxplot(fill = "#E69F00", color = "black") +
  geom_hline(yintercept = 1, color = "red", linetype = "dotted") + 
  labs(x = "Task", y = "Fraction of Requested CPU") +
  scale_y_continuous(limits = c(0, max(data$PercentRequestedCPU, na.rm = TRUE) * 1.1)) + 
  custom_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

combined_plot <- p1 / p2 / p3 / p4 +
  plot_layout(heights = c(1, 1, 1, 1)) 

ggsave("figures/per_task_statistics_boxplots_raw.svg", combined_plot, width = 210, height = 260, units = "mm")

# Calculate overall statistics for PercentRequestedCPU
data$TotalRealTimeSeconds <- convert_to_seconds(data$realtime) # Assuming 'realtime' column exists and is structured like 'duration'
data$CPUHours <- (data$TotalRealTimeSeconds / 3600) * (data$PercentCPU / 100)
total_cpu_hours <- sum(data$CPUHours, na.rm = TRUE)
cat("Total CPU Hours:", total_cpu_hours, "\n")
overall_stats <- data %>%
  summarise(
    MedianPercentRequestedCPU = median(PercentRequestedCPU, na.rm = TRUE),
    MeanPercentRequestedCPU = mean(PercentRequestedCPU, na.rm = TRUE),
    SDPercentRequestedCPU = sd(PercentRequestedCPU, na.rm = TRUE)
  )

cat("Overall Median PercentRequestedCPU:", overall_stats$MedianPercentRequestedCPU, "\n")
cat("Overall Mean PercentRequestedCPU:", overall_stats$MeanPercentRequestedCPU, "\n")
cat("Overall Standard Deviation PercentRequestedCPU:", overall_stats$SDPercentRequestedCPU, "\n")

stats_by_task <- data %>%
  group_by(Task) %>%
  summarise(
    MedianPercentRequestedCPU = median(PercentRequestedCPU, na.rm = TRUE),
    MeanPercentRequestedCPU = mean(PercentRequestedCPU, na.rm = TRUE),
    SDPercentRequestedCPU = sd(PercentRequestedCPU, na.rm = TRUE)
  )

print(stats_by_task)
