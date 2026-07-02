library(ggplot2)
library(dplyr)
library(patchwork)
library(rsvg)
library(grid)
library(egg)

pdf(NULL)  # prevent Rplots.pdf side-effect

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(script_dir)

data <- read.table(
  "improved_report.txt",
  header = FALSE, sep = "\t", stringsAsFactors = FALSE,
  col.names = c("name", "status", "duration", "pct_cpu", "peak_vmem", "cpus", "full_name")
)
data <- data[data$status %in% c("CACHED", "COMPLETED"), ]

full_path <- gsub(" \\(.*\\)$", "", data$name)
leaf_task <- sub(".*:", "", full_path)

data$Task <- dplyr::case_when(
  grepl(":PREPARE_REFERENCE_FILES:", full_path)         ~ "PREPARE_REFERENCE_FILES",
  leaf_task == "METAEXTRACT"                            ~ "PREPARE_REFERENCE_FILES",
  grepl(":DEEPVARIANT:", full_path)                     ~ "DEEPVARIANT",
  grepl(":DEEPSOMATIC:", full_path)                     ~ "DEEPSOMATIC",
  grepl("^FIBERTOOLSRS_", leaf_task)                    ~ "FIBERTOOLS",
  grepl("^LONGPHASE_", leaf_task) & grepl("SOMATIC",  leaf_task) ~ "LONGPHASE_SOMATIC",
  grepl("^LONGPHASE_", leaf_task) & grepl("GERMLINE", leaf_task) ~ "LONGPHASE_GERMLINE",
  leaf_task == "LONGPHASE_HAPLOTAG"                     ~ "LONGPHASE_HAPLOTAG",
  leaf_task %in% c("CRAMINO_PRE", "CRAMINO_POST")       ~ "CRAMINO",
  leaf_task %in% c("NANOPLOT_PRE", "NANOPLOT_POST")     ~ "NANOPLOT",
  grepl("^SAMTOOLS_", leaf_task)                        ~ "SAMTOOLS",
  grepl("^BCFTOOLS_", leaf_task)                        ~ "BCFTOOLS",
  leaf_task == "STANDARDIZE_AF"                         ~ "BCFTOOLS",
  TRUE                                                  ~ leaf_task
)

task_group_map <- c(
  CLAIR3 = "clair", CLAIRS = "clair", CLAIRSTO = "clair",
  SV_VEP = "vep", SOMATIC_VEP = "vep", GERMLINE_VEP = "vep",
  DEEPSOMATIC = "deep", DEEPVARIANT = "deep",
  FIBERTOOLS = "fiber",
  LONGPHASE_HAPLOTAG = "longphase", LONGPHASE_SOMATIC = "longphase", LONGPHASE_GERMLINE = "longphase"
)
task_df <- data.frame(Task = unique(data$Task), stringsAsFactors = FALSE)
task_df$group <- ifelse(task_df$Task %in% names(task_group_map),
                        task_group_map[task_df$Task], task_df$Task)
task_df$len <- nchar(task_df$Task)
task_df$group_min <- ave(task_df$len, task_df$group, FUN = min)
task_df <- task_df[order(task_df$group_min, task_df$group, task_df$len), ]
data$Task <- factor(data$Task, levels = task_df$Task)

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
data$DurationHours   <- data$DurationSeconds / 3600
data$PercentCPU      <- as.numeric(sub("%", "", data$pct_cpu))

convert_to_gb <- function(vmem) {
  mb <- as.numeric(sub(" MB", "", vmem))
  gb <- as.numeric(sub(" GB", "", vmem))
  mb[!is.na(mb)] <- mb[!is.na(mb)] / 1024
  gb[is.na(mb)]  <- gb[is.na(mb)]
  return(ifelse(is.na(mb), gb, mb))
}

data$PeakVmemGB <- convert_to_gb(data$peak_vmem)

data$PercentCPU[is.na(data$PercentCPU)] <- 0
data$cpus <- as.numeric(data$cpus)
data$cpus[is.na(data$cpus)] <- 1

data$PercentRequestedCPU <- (data$PercentCPU / 100) / data$cpus

# Extract sample tag from full_name column, e.g. "... (CCS15-PB)" → "CCS15-PB"
data$Sample <- sub(".*\\(([^)]+)\\)$", "\\1", data$full_name)

axis_label_size <- 7

custom_theme <- theme(
  text = element_text(family = "sans", size = 14),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.y = element_text(family = "sans", size = axis_label_size),
  axis.title = element_text(family = "sans", size = axis_label_size, face = "bold"),
  plot.title = element_text(hjust = 0.5, family = "sans", size = 12),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  panel.background = element_rect(fill = NA, color = "black"),
  panel.grid = element_blank(),
  plot.margin = margin(1, 2, 1, 2, "pt")
)

box_lwd    <- 0.3
outlier_sz <- 0.6

p1 <- ggplot(data, aes(x = Task, y = DurationHours)) +
  geom_boxplot(fill = "#7FB4CA", color = "black",
               linewidth = box_lwd, outlier.size = outlier_sz) +
  labs(x = NULL, y = "Duration (Hours)") +
  custom_theme

p2 <- ggplot(data, aes(x = Task, y = PercentCPU)) +
  geom_boxplot(fill = "#98BB6C", color = "black",
               linewidth = box_lwd, outlier.size = outlier_sz) +
  labs(x = NULL, y = "CPU Usage (%)") +
  custom_theme

p3 <- ggplot(data, aes(x = Task, y = PeakVmemGB)) +
  geom_boxplot(fill = "#D27E99", color = "black",
               linewidth = box_lwd, outlier.size = outlier_sz) +
  labs(x = NULL, y = "Peak VMem (GB)") +
  custom_theme

p4 <- ggplot(data, aes(x = Task, y = PercentRequestedCPU)) +
  geom_boxplot(fill = "#E69F00", color = "black",
               linewidth = box_lwd, outlier.size = outlier_sz) +
  geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  labs(x = "Task", y = "Requested CPU (Frac.)") +
  scale_y_continuous(limits = c(0, max(data$PercentRequestedCPU, na.rm = TRUE) * 1.1)) +
  custom_theme +
  theme(axis.text.x = element_text(size = axis_label_size, angle = 45, hjust = 1))

# Combined A4 figure: lrsomatic logo (A) on top, boxplots (B) on bottom
img_grob <- rasterGrob(
  rsvg_nativeraster("lrsomatic_2.0.svg", width = 2400),
  interpolate = TRUE
)
logo_panel <- wrap_elements(full = img_grob)

small_text <- theme(
  text         = element_text(size = 7.5),
  axis.text.y  = element_text(size = axis_label_size),
  axis.title   = element_text(size = axis_label_size, face = "bold"),
  plot.title   = element_text(size = 7.5),
  legend.text  = element_text(size = 7.5),
  legend.title = element_text(size = 7.5),
  plot.margin  = margin(0, 2.3, 0, 1.7, "pt")
)
boxplot_panel <- egg::ggarrange(
  p1 + small_text,
  p2 + small_text,
  p3 + small_text,
  p4  + small_text,
  ncol = 1
)

combined_a4 <- logo_panel / wrap_elements(boxplot_panel) +
  plot_layout(heights = c(0.7, 1.3)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("per_task_statistics_combined_a4.svg",
       combined_a4, width = 210, height = 297, units = "mm")
ggsave("per_task_statistics_combined_a4.png",
       combined_a4, width = 6.5, height = 9, units = "in", dpi = 300)

# Report: resource statistics
fmt_mem <- function(gb) {
  if (gb < 1) sprintf("%.0f MB", gb * 1024) else sprintf("%.2f GB", gb)
}

summarise_resources <- function(df, label) {
  med_cpus     <- median(df$cpus, na.rm = TRUE)
  min_cpus     <- min(df$cpus, na.rm = TRUE)
  max_cpus     <- max(df$cpus, na.rm = TRUE)
  med_vmem     <- median(df$PeakVmemGB, na.rm = TRUE)
  min_vmem     <- min(df$PeakVmemGB, na.rm = TRUE)
  max_vmem     <- max(df$PeakVmemGB, na.rm = TRUE)
  med_req_frac <- median(df$PercentRequestedCPU, na.rm = TRUE)
  c(
    sprintf("== %s ==", label),
    sprintf("N tasks: %d", nrow(df)),
    sprintf("Median CPUs requested: %g (range %g - %g)", med_cpus, min_cpus, max_cpus),
    sprintf("Median peak VMem: %s (range %s - %s)",
            fmt_mem(med_vmem), fmt_mem(min_vmem), fmt_mem(max_vmem)),
    sprintf("Median fraction of requested CPU used: %.0f%%", 100 * med_req_frac),
    ""
  )
}

excluded_samples <- c("CCS15-PB", "CCS15-PB-TO")
data_no_ccs15pb  <- data[!data$Sample %in% excluded_samples, ]

writeLines(
  c(
    summarise_resources(data,            "All tasks"),
    summarise_resources(data_no_ccs15pb, "Excluding PacBio CCS15 (CCS15-PB, CCS15-PB-TO)")
  ),
  "summary.txt"
)
