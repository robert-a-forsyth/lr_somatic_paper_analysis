# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read data from the file
# Replace the path with the actual path to your input file
data <- read.table("/Users/u0155044/Documents/lr_somatic_figures/CCS15-ONT.segments.txt", header = TRUE, sep = "\t")

# Define adjustable parameters
allele_width <- 0.4

# Remove 'chr' prefix from the chromosome column, and round major/minor allele copy numbers
data <- data %>%
  mutate(
    chr = gsub("^chr", "", chr),  # Remove 'chr' prefix if present
    rounded_nMajor = round(nMajor),  # Round major allele values
    rounded_nMinor = round(nMinor),  # Round minor allele values
    chr = factor(chr, levels = c(as.character(1:22), "X"), ordered = TRUE)  # Ensure proper chromosome ordering
  )

# Compute cumulative offsets for chromosomes
chromosome_offsets <- data %>%
  group_by(chr) %>%
  summarise(
    chr_length = max(endpos) - min(startpos) + 1,  # Length of the chromosome
    chr_max_pos = max(endpos)                     # Max position within the chromosome
  ) %>%
  mutate(
    offset = lag(cumsum(chr_length), default = 0),       # Cumulative offset for each chromosome
    boundary_position = offset + chr_max_pos            # Boundary position for chromosome
  )

# Add cumulative positions for plotting
plot_data <- data %>%
  left_join(chromosome_offsets, by = "chr") %>%
  mutate(
    cum_start = startpos + offset,      # Adjusted start position with offset
    cum_end = endpos + offset,          # Adjusted end position with offset
    major_ymin = rounded_nMajor,        # Bottom of the major allele rectangle
    major_ymax = rounded_nMajor + allele_width,  # Top of the major allele rectangle
    minor_ymin = rounded_nMinor - allele_width,  # Bottom of the minor allele rectangle
    minor_ymax = rounded_nMinor                # Top of the minor allele rectangle
  )

# Prepare chromosome labels and positions for the x-axis
chromosome_labels <- chromosome_offsets %>%
  mutate(label_position = boundary_position)  # Position of the label at the chromosome boundary

# Prepare boundary lines for chromosomes
chromosome_lines <- chromosome_offsets %>%
  select(boundary_position)

# Custom theme with options for title and subtitle sizes
custom_theme <- theme(
  text = element_text(family = "sans", size = 14),     # Base font size
  axis.text.x = element_text(hjust = 1, size = 10, color = "black"),  # X-axis text
  axis.text.y = element_text(size = 10, color = "black"),             # Y-axis text
  axis.title = element_text(face = "bold", size = 12),                # Axis titles
  plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),   # Main title
  plot.subtitle = element_text(hjust = 0.5, size = 10, face = "bold"), # Subtitle
  strip.text = element_text(size = 12, face = "bold"),                # Bold facet (subplot) titles
  strip.background = element_rect(fill = "grey97", color = "black"),  # Facet title background
  panel.background = element_rect(fill = NA, color = "black"),        # Black box around the panel
  panel.grid.major.y = element_line(color = "grey", size = 0.5, alpha = 0.7),      # Horizontal gridlines for whole numbers
  panel.grid.minor.y = element_blank(),                               # Remove minor gridlines
  panel.grid.major.x = element_blank(),                               # Remove vertical gridlines
  legend.position = "none"                                            # Remove legend from the plot
)

# Generate the main plot (no legend included directly)
final_plot <- ggplot() +
  # Add rectangles for major alleles
  geom_rect(
    data = plot_data,
    aes(xmin = cum_start, xmax = cum_end, ymin = major_ymin, ymax = major_ymax, fill = "Major Allele")
  ) +
  # Add rectangles for minor alleles
  geom_rect(
    data = plot_data,
    aes(xmin = cum_start, xmax = cum_end, ymin = minor_ymin, ymax = minor_ymax, fill = "Minor Allele")
  ) +
  # Add dotted lines for chromosome boundaries
  geom_vline(
    data = chromosome_lines,
    aes(xintercept = boundary_position),
    linetype = "dotted",
    color = "grey",
    alpha = 0.7
  ) +
  # Add chromosome labels at the boundaries
  scale_x_continuous(
    breaks = chromosome_labels$label_position,
    labels = chromosome_labels$chr
  ) +
  # Set y-axis ticks at every whole number with gridlines at each tick
  scale_y_continuous(
    breaks = seq(0, 10, 1)  # Ticks at every whole number from 0 to 10
  ) +
  # Set colors for Major and Minor Alleles
  scale_fill_manual(
    values = c("Major Allele" = "#C34043", "Minor Allele" = "#466cb8"),
    name = "Allele Type"
  ) +
  # Add axis labels, title, and subtitle
  labs(
    x = "Chromosome",
    y = "Allele Value",
    title = "ASCAT ONT Copy Number Profile",           # Add main title
    subtitle = "Ploidy: 65%     Purity: 5.09"            # Add subtitle
  ) +
  # Set y-axis limits from 0 to 10
  coord_cartesian(ylim = c(0, 10)) +
  # Apply the custom theme
  custom_theme

# Save the main plot as an SVG
svg(filename = "CCS15_ONT_ASCAT.svg", width = 8.5, height = 3)
print(final_plot)
dev.off()


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Read data from the file
# Replace the path with the actual path to your input file
data <- read.table("/Users/u0155044/Documents/lr_somatic_figures/CCS15-PB.segments.txt", header = TRUE, sep = "\t")

# Define adjustable parameters
allele_width <- 0.4

# Remove 'chr' prefix from the chromosome column, and round major/minor allele copy numbers
data <- data %>%
  mutate(
    chr = gsub("^chr", "", chr),  # Remove 'chr' prefix if present
    rounded_nMajor = round(nMajor),  # Round major allele values
    rounded_nMinor = round(nMinor),  # Round minor allele values
    chr = factor(chr, levels = c(as.character(1:22), "X"), ordered = TRUE)  # Ensure proper chromosome ordering
  )

# Compute cumulative offsets for chromosomes
chromosome_offsets <- data %>%
  group_by(chr) %>%
  summarise(
    chr_length = max(endpos) - min(startpos) + 1,  # Length of the chromosome
    chr_max_pos = max(endpos)                     # Max position within the chromosome
  ) %>%
  mutate(
    offset = lag(cumsum(chr_length), default = 0),       # Cumulative offset for each chromosome
    boundary_position = offset + chr_max_pos            # Boundary position for chromosome
  )

# Add cumulative positions for plotting
plot_data <- data %>%
  left_join(chromosome_offsets, by = "chr") %>%
  mutate(
    cum_start = startpos + offset,      # Adjusted start position with offset
    cum_end = endpos + offset,          # Adjusted end position with offset
    major_ymin = rounded_nMajor,        # Bottom of the major allele rectangle
    major_ymax = rounded_nMajor + allele_width,  # Top of the major allele rectangle
    minor_ymin = rounded_nMinor - allele_width,  # Bottom of the minor allele rectangle
    minor_ymax = rounded_nMinor                # Top of the minor allele rectangle
  )

# Prepare chromosome labels and positions for the x-axis
chromosome_labels <- chromosome_offsets %>%
  mutate(label_position = boundary_position)  # Position of the label at the chromosome boundary

# Prepare boundary lines for chromosomes
chromosome_lines <- chromosome_offsets %>%
  select(boundary_position)

# Custom theme with options for title and subtitle sizes
custom_theme <- theme(
  text = element_text(family = "sans", size = 14),     # Base font size
  axis.text.x = element_text(hjust = 1, size = 10, color = "black"),  # X-axis text
  axis.text.y = element_text(size = 10, color = "black"),             # Y-axis text
  axis.title = element_text(face = "bold", size = 12),                # Axis titles
  plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),   # Main title
  plot.subtitle = element_text(hjust = 0.5, size = 10, face = "bold"), # Subtitle
  strip.text = element_text(size = 12, face = "bold"),                # Bold facet (subplot) titles
  strip.background = element_rect(fill = "grey97", color = "black"),  # Facet title background
  panel.background = element_rect(fill = NA, color = "black"),        # Black box around the panel
  panel.grid.major.y = element_line(color = "grey", size = 0.5, alpha = 0.7),      # Horizontal gridlines for whole numbers
  panel.grid.minor.y = element_blank(),                               # Remove minor gridlines
  panel.grid.major.x = element_blank(),                               # Remove vertical gridlines
  legend.position = "none"                                            # Remove legend from the plot
)

# Generate the main plot (no legend included directly)
final_plot <- ggplot() +
  # Add rectangles for major alleles
  geom_rect(
    data = plot_data,
    aes(xmin = cum_start, xmax = cum_end, ymin = major_ymin, ymax = major_ymax, fill = "Major Allele")
  ) +
  # Add rectangles for minor alleles
  geom_rect(
    data = plot_data,
    aes(xmin = cum_start, xmax = cum_end, ymin = minor_ymin, ymax = minor_ymax, fill = "Minor Allele")
  ) +
  # Add dotted lines for chromosome boundaries
  geom_vline(
    data = chromosome_lines,
    aes(xintercept = boundary_position),
    linetype = "dotted",
    color = "grey",
    alpha = 0.7
  ) +
  # Add chromosome labels at the boundaries
  scale_x_continuous(
    breaks = chromosome_labels$label_position,
    labels = chromosome_labels$chr
  ) +
  # Set y-axis ticks at every whole number with gridlines at each tick
  scale_y_continuous(
    breaks = seq(0, 10, 1)  # Ticks at every whole number from 0 to 10
  ) +
  # Set colors for Major and Minor Alleles
  scale_fill_manual(
    values = c("Major Allele" = "#C34043", "Minor Allele" = "#466cb8"),
    name = "Allele Type"
  ) +
  # Add axis labels, title, and subtitle
  labs(
    x = "Chromosome",
    y = "Allele Value",
    title = "ASCAT PacBio Copy Number Profile",           # Add main title
    subtitle = "Ploidy: 81%     Purity: 3.44"            # Add subtitle
  ) +
  # Set y-axis limits from 0 to 10
  coord_cartesian(ylim = c(0, 10)) +
  # Apply the custom theme
  custom_theme

# Save the main plot as an SVG
svg(filename = "CCS15_PB_ASCAT.svg", width = 8.5, height = 3)
print(final_plot)
dev.off()



