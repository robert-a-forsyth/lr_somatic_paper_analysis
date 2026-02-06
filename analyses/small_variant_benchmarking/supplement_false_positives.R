library(ggplot2)
library(dplyr)
library(vcfR)
library(tidyr)

# Set working directory
setwd("/Users/u0155044/Documents/lr_somatic_paper_analysis/analyses/small_variant_benchmarking")

# Path to intermediate files
intermediate_dir <- "/Users/u0155044/Documents/lr_somatic_paper_analysis/analyses/small_variant_benchmarking/intermediate_files"

# Define samples
samples <- c("COLO829-PB", "COLO829-ONT", "COLO829-PB-TO", "COLO829-ONT-TO")

# Function to extract QUAL scores from VCF
extract_qual_from_vcf <- function(vcf_path) {
  tryCatch({
    vcf <- read.vcfR(vcf_path, verbose = FALSE)
    
    qual_scores <- as.numeric(vcf@fix[, "QUAL"])
    
    # Determine variant type
    ref <- vcf@fix[, "REF"]
    alt <- vcf@fix[, "ALT"]
    var_type <- ifelse(nchar(ref) == 1 & nchar(alt) == 1, "SNV", "Indel")
    
    data.frame(
      QUAL = qual_scores,
      VariantType = var_type,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(paste("Error reading", vcf_path, ":", e$message))
    return(NULL)
  })
}

# Process each sample
all_data <- list()

for (sample in samples) {
  message(paste("Processing", sample, "..."))
  
  # Paths to isec results
  isec_indels <- file.path(intermediate_dir, paste0("isec_indels_", sample))
  isec_snvs <- file.path(intermediate_dir, paste0("isec_snvs_", sample))
  
  # Read TP and FP for indels (0003.vcf = TP, 0001.vcf = FP)
  tp_indel_file <- paste0(isec_indels, "/0003.vcf")
  if (file.exists(tp_indel_file)) {
    tp_indels <- extract_qual_from_vcf(tp_indel_file)
    if (!is.null(tp_indels) && nrow(tp_indels) > 0) {
      tp_indels$Classification <- "True Positive"
      tp_indels$Sample <- sample
      all_data[[length(all_data) + 1]] <- tp_indels
    }
  }
  
  fp_indel_file <- paste0(isec_indels, "/0001.vcf")
  if (file.exists(fp_indel_file)) {
    fp_indels <- extract_qual_from_vcf(fp_indel_file)
    if (!is.null(fp_indels) && nrow(fp_indels) > 0) {
      fp_indels$Classification <- "False Positive"
      fp_indels$Sample <- sample
      all_data[[length(all_data) + 1]] <- fp_indels
    }
  }
  
  # Read TP and FP for SNVs
  tp_snv_file <- paste0(isec_snvs, "/0003.vcf")
  if (file.exists(tp_snv_file)) {
    tp_snvs <- extract_qual_from_vcf(tp_snv_file)
    if (!is.null(tp_snvs) && nrow(tp_snvs) > 0) {
      tp_snvs$Classification <- "True Positive"
      tp_snvs$Sample <- sample
      all_data[[length(all_data) + 1]] <- tp_snvs
    }
  }
  
  fp_snv_file <- paste0(isec_snvs, "/0001.vcf")
  if (file.exists(fp_snv_file)) {
    fp_snvs <- extract_qual_from_vcf(fp_snv_file)
    if (!is.null(fp_snvs) && nrow(fp_snvs) > 0) {
      fp_snvs$Classification <- "False Positive"
      fp_snvs$Sample <- sample
      all_data[[length(all_data) + 1]] <- fp_snvs
    }
  }
}

# Combine all data
qual_data <- bind_rows(all_data) %>%
  filter(!is.na(QUAL))

# Create output directory
dir.create("quality_score_plots", showWarnings = FALSE)

# QUAL distribution faceted by variant type, showing TP vs FP (vertical violins)
p <- ggplot(qual_data, aes(x = Sample, y = QUAL, fill = Classification)) +
  geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(width = 0.9)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5, 
               position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 14, linetype = "dashed", color = "red", linewidth = 0.5,alpha = 0.5) +
  facet_wrap(~VariantType, ncol = 2) +
  scale_fill_manual(values = c("True Positive" = "#2E7D32", "False Positive" = "#C62828")) +
  labs(
       x = "Sample",
       y = "QUAL Score",
       fill = "Classification") +
  theme(
    text = element_text(family = "sans", size = 14),
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey97", color = "black"),
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )
ggsave("/Users/u0155044/Documents/lr_somatic_paper_analysis/figures/supplementary_true_false_positives.svg", plot = p, width = 210, height = 100, units = "mm")


