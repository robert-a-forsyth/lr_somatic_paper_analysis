# make_figure4.R
# Faceted A4 comparison of ASCAT, Wakhan, and PURPLE copy number profiles
# for CCS15 (ONT and PacBio). Output: PDF + SVG in Fig4/.
#
# Run from inside Fig4/:  Rscript make_figure4.R
# Override data root:     DATA_DIR=../Run/out Rscript make_figure4.R
# Override PURPLE file:   PURPLE_SEG=/path/to/CCS15M.purple.segment.tsv Rscript make_figure4.R
# By default PURPLE_SEG is read from Fig4/CCS15M.purple.segment.tsv (bundled in repo)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ── Data root ─────────────────────────────────────────────────────────────────
DATA_DIR <- Sys.getenv("DATA_DIR", "../Run/out")

# ── Path helpers ──────────────────────────────────────────────────────────────
ascat_path <- function(s) {
  file.path(DATA_DIR, s, "ascat", paste0(s, ".segments.txt"))
}

wakhan_path <- function(s, hp) {
  pat  <- file.path(DATA_DIR, s, "wakhan", "solution_1", "bed_output",
                    sprintf("%s_*_copynumbers_segments_HP_%d.bed", s, hp))
  hits <- Sys.glob(pat)
  if (length(hits) != 1)
    stop("Expected 1 Wakhan file for ", s, " HP", hp,
         "; got ", length(hits), " matching: ", pat)
  hits
}

# ── Input files ───────────────────────────────────────────────────────────────
ASCAT_ONT      <- ascat_path("CCS15-ONT")
ASCAT_ONT_TO   <- ascat_path("CCS15-ONT-TO")
ASCAT_PB       <- ascat_path("CCS15-PB")
ASCAT_PB_TO    <- ascat_path("CCS15-PB-TO")

WAKHAN_ONT1    <- wakhan_path("CCS15-ONT",    1); WAKHAN_ONT2    <- wakhan_path("CCS15-ONT",    2)
WAKHAN_ONT_TO1 <- wakhan_path("CCS15-ONT-TO", 1); WAKHAN_ONT_TO2 <- wakhan_path("CCS15-ONT-TO", 2)
WAKHAN_PB1     <- wakhan_path("CCS15-PB",     1); WAKHAN_PB2     <- wakhan_path("CCS15-PB",     2)
WAKHAN_PB_TO1  <- wakhan_path("CCS15-PB-TO",  1); WAKHAN_PB_TO2  <- wakhan_path("CCS15-PB-TO",  2)

PURPLE_SEG     <- Sys.getenv("PURPLE_SEG", "CCS15M.purple.segment.tsv")

# ── Output files ──────────────────────────────────────────────────────────────
OUT_PDF    <- "combined_cnv_plots.pdf"
OUT_SVG    <- "combined_cnv_plots.svg"
OUT_TO_PDF <- "combined_cnv_plots_to.pdf"
OUT_TO_SVG <- "combined_cnv_plots_to.svg"

# ── Panel labels ──────────────────────────────────────────────────────────────
PANEL_META <- list(
  purple        = list(title    = "Oncoanalyser (PURPLE) - Illumina",
                       subtitle = "Ploidy: 3.45 | Purity: 63%"),
  ascat_pb      = list(title    = "ASCAT - PacBio",
                       subtitle = "Ploidy: 3.46 | Purity: 80% | GoF: 97.6%"),
  wakhan_pb     = list(title    = "Wakhan - PacBio",
                       subtitle = "Ploidy: 3.42 | Purity: 76% | Confidence: 0.91"),
  ascat_ont     = list(title    = "ASCAT - ONT",
                       subtitle = "Ploidy: 1.85 | Purity: 85% | GoF: 95.8%"),
  wakhan_ont    = list(title    = "Wakhan - ONT",
                       subtitle = "Ploidy: 3.44 | Purity: 65% | Confidence: 0.96"),
  ascat_pb_to   = list(title    = "ASCAT - PacBio (Tumor-Only)",
                       subtitle = "Ploidy: 3.50 | Purity: 79% | GoF: 97.7%"),
  wakhan_pb_to  = list(title    = "Wakhan - PacBio (Tumor-Only)",
                       subtitle = "Ploidy: 1.64 | Purity: 33% | Confidence: 0.88"),
  ascat_ont_to  = list(title    = "ASCAT - ONT (Tumor-Only)",
                       subtitle = "Ploidy: 5.14 | Purity: 58% | GoF: 97.1%"),
  wakhan_ont_to = list(title    = "Wakhan - ONT (Tumor-Only)",
                       subtitle = "Ploidy: 3.28 | Purity: 44% | Confidence: 0.94")
)

# ── Constants ─────────────────────────────────────────────────────────────────
CHR_LEVELS   <- c(as.character(1:22), "X")
MAJOR_COLOR  <- "#C34043"
MINOR_COLOR  <- "#466cb8"
ALLELE_WIDTH <- 0.4
Y_MAX        <- 6

# ── Helper: normalise chromosome column ───────────────────────────────────────
normalize_chr <- function(x) {
  factor(gsub("^chr", "", as.character(x)), levels = CHR_LEVELS)
}

# ── Helper: canonical chromosome offsets ─────────────────────────────────────
compute_offsets <- function(list_of_dfs) {
  bind_rows(list_of_dfs) %>%
    filter(!is.na(chr)) %>%
    group_by(chr) %>%
    summarise(chr_span = max(end), .groups = "drop") %>%
    arrange(chr) %>%
    mutate(
      offset            = lag(cumsum(as.numeric(chr_span)), default = 0),
      boundary_position = offset + as.numeric(chr_span),
      mid_position      = offset + as.numeric(chr_span) / 2
    )
}

# ── Helper: shared panel theme ────────────────────────────────────────────────
panel_theme <- theme(
  text               = element_text(family = "sans", size = 13),
  axis.text.y        = element_text(size = 10, color = "black"),
  axis.title.y       = element_text(face = "bold", size = 11),
  axis.title.x       = element_text(face = "bold", size = 11),
  plot.title         = element_text(hjust = 0.5, size = 13, face = "bold"),
  plot.subtitle      = element_text(hjust = 0.5, size = 10),
  panel.background   = element_rect(fill = NA, color = "black"),
  panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_blank(),
  legend.position    = "bottom",
  legend.key.size    = unit(0.45, "cm"),
  legend.text        = element_text(size = 10),
  legend.title       = element_text(size = 11, face = "bold"),
  plot.margin        = margin(10, 4, 2, 14)
)

no_x_axis <- theme(
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x = element_blank()
)

# ── Helper: build one CN panel ────────────────────────────────────────────────
make_cnv_panel <- function(df, offsets, title, subtitle, show_x_axis = TRUE) {
  plot_data <- df %>%
    filter(!is.na(chr)) %>%
    left_join(offsets, by = "chr") %>%
    mutate(
      cum_start    = start + offset,
      cum_end      = end   + offset,
      major_capped = pmin(major, Y_MAX),
      minor_capped = pmin(minor, Y_MAX),
      major_ymin   = major_capped,
      major_ymax   = major_capped + ALLELE_WIDTH,
      minor_ymin   = minor_capped - ALLELE_WIDTH,
      minor_ymax   = minor_capped
    )

  x_scale <- scale_x_continuous(
    breaks = offsets$mid_position,
    labels = offsets$chr,
    expand = c(0, 0)
  )

  p <- ggplot() +
    geom_rect(data = plot_data,
              aes(xmin = cum_start, xmax = cum_end,
                  ymin = major_ymin, ymax = major_ymax,
                  fill = "Major allele")) +
    geom_rect(data = plot_data,
              aes(xmin = cum_start, xmax = cum_end,
                  ymin = minor_ymin, ymax = minor_ymax,
                  fill = "Minor allele")) +
    geom_vline(data = offsets,
               aes(xintercept = boundary_position),
               linetype = "dotted", color = "grey60", alpha = 0.6) +
    scale_fill_manual(
      values = c("Major allele" = MAJOR_COLOR, "Minor allele" = MINOR_COLOR),
      name   = "Allele"
    ) +
    x_scale +
    scale_y_continuous(breaks = seq(0, Y_MAX, 1)) +
    coord_cartesian(ylim = c(0, Y_MAX)) +
    labs(title = title, subtitle = subtitle,
         x = "Chromosome", y = "Copy number") +
    panel_theme

  if (!show_x_axis) p <- p + no_x_axis
  p
}

# ── Data loaders ──────────────────────────────────────────────────────────────
load_ascat <- function(path) {
  read.table(path, header = TRUE, sep = "\t") %>%
    transmute(
      chr   = normalize_chr(chr),
      start = startpos,
      end   = endpos,
      major = as.integer(round(nMajor)),
      minor = as.integer(round(nMinor))
    ) %>%
    filter(!is.na(chr))
}

load_purple <- function(path) {
  read.table(path, header = TRUE, sep = "\t") %>%
    transmute(
      chr   = normalize_chr(chromosome),
      start = start,
      end   = end,
      major = as.integer(round(majorAlleleCopyNumber)),
      minor = as.integer(round(minorAlleleCopyNumber))
    ) %>%
    filter(!is.na(chr))
}

load_wakhan <- function(hp1_path, hp2_path) {
  col_names <- c("chr", "start", "end", "coverage",
                 "copynumber_state", "confidence", "svs_breakpoints_ids")
  read_bed  <- function(p) {
    read.table(p, header = FALSE, sep = "\t",
               comment.char = "#", col.names = col_names)
  }
  hp1 <- read_bed(hp1_path)
  hp2 <- read_bed(hp2_path)

  joined <- inner_join(
    hp1 %>% select(chr, start, end, cn_hp1 = copynumber_state),
    hp2 %>% select(chr, start, end, cn_hp2 = copynumber_state),
    by = c("chr", "start", "end")
  )

  dropped <- min(nrow(hp1), nrow(hp2)) - nrow(joined)
  if (dropped > 0)
    message("  Note: ", dropped, " rows dropped (interval mismatch between HP1/HP2).")

  joined %>%
    transmute(
      chr   = normalize_chr(chr),
      start = start,
      end   = end,
      major = as.integer(round(pmax(cn_hp1, cn_hp2))),
      minor = as.integer(round(pmin(cn_hp1, cn_hp2)))
    ) %>%
    filter(!is.na(chr))
}

# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading ASCAT segments...")
dat_ascat_pb     <- load_ascat(ASCAT_PB)
dat_ascat_ont    <- load_ascat(ASCAT_ONT)
dat_ascat_pb_to  <- load_ascat(ASCAT_PB_TO)
dat_ascat_ont_to <- load_ascat(ASCAT_ONT_TO)

message("Loading Wakhan haplotype segments...")
dat_wakhan_pb     <- load_wakhan(WAKHAN_PB1,     WAKHAN_PB2)
dat_wakhan_ont    <- load_wakhan(WAKHAN_ONT1,    WAKHAN_ONT2)
dat_wakhan_pb_to  <- load_wakhan(WAKHAN_PB_TO1,  WAKHAN_PB_TO2)
dat_wakhan_ont_to <- load_wakhan(WAKHAN_ONT_TO1, WAKHAN_ONT_TO2)

message("Loading PURPLE segments...")
dat_purple <- load_purple(PURPLE_SEG)

# ── Canonical chromosome offsets ──────────────────────────────────────────────
message("Computing chromosome offsets...")
offsets <- compute_offsets(list(
  dat_purple,
  dat_ascat_pb,     dat_wakhan_pb,
  dat_ascat_ont,    dat_wakhan_ont,
  dat_ascat_pb_to,  dat_wakhan_pb_to,
  dat_ascat_ont_to, dat_wakhan_ont_to
))

# ── Assemble 5-panel figure ───────────────────────────────────────────────────
make_figure <- function(p_purple, p_ascat_pb, p_wakhan_pb,
                        p_ascat_ont, p_wakhan_ont) {
  p_purple    /
  p_ascat_pb  /
  p_ascat_ont /
  p_wakhan_pb /
  p_wakhan_ont +
    plot_layout(heights = rep(1, 5), guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16,
                                  family = "sans", margin = margin(b = 6)),
        plot.tag   = element_text(size = 14, family = "sans")
      )
    ) &
    theme(
      legend.position   = "bottom",
      plot.tag.position = c(0, 0.93),
      plot.tag          = element_text(face = "bold", size = 14, family = "sans")
    )
}

# ── Figure 1: Tumor-Normal ────────────────────────────────────────────────────
message("Building tumor-normal figure...")
figure_tn <- make_figure(
  make_cnv_panel(dat_purple,     offsets,
    PANEL_META$purple$title,     PANEL_META$purple$subtitle,     show_x_axis = FALSE),
  make_cnv_panel(dat_ascat_pb,   offsets,
    PANEL_META$ascat_pb$title,   PANEL_META$ascat_pb$subtitle,   show_x_axis = FALSE),
  make_cnv_panel(dat_wakhan_pb,  offsets,
    PANEL_META$wakhan_pb$title,  PANEL_META$wakhan_pb$subtitle,  show_x_axis = FALSE),
  make_cnv_panel(dat_ascat_ont,  offsets,
    PANEL_META$ascat_ont$title,  PANEL_META$ascat_ont$subtitle,  show_x_axis = FALSE),
  make_cnv_panel(dat_wakhan_ont, offsets,
    PANEL_META$wakhan_ont$title, PANEL_META$wakhan_ont$subtitle, show_x_axis = TRUE)
)

# ── Figure 2: Tumor-Only ──────────────────────────────────────────────────────
message("Building tumor-only figure...")
figure_to <- make_figure(
  make_cnv_panel(dat_purple,       offsets,
    PANEL_META$purple$title,       PANEL_META$purple$subtitle,       show_x_axis = FALSE),
  make_cnv_panel(dat_ascat_pb_to,  offsets,
    PANEL_META$ascat_pb_to$title,  PANEL_META$ascat_pb_to$subtitle,  show_x_axis = FALSE),
  make_cnv_panel(dat_wakhan_pb_to, offsets,
    PANEL_META$wakhan_pb_to$title, PANEL_META$wakhan_pb_to$subtitle, show_x_axis = FALSE),
  make_cnv_panel(dat_ascat_ont_to, offsets,
    PANEL_META$ascat_ont_to$title, PANEL_META$ascat_ont_to$subtitle, show_x_axis = FALSE),
  make_cnv_panel(dat_wakhan_ont_to, offsets,
    PANEL_META$wakhan_ont_to$title, PANEL_META$wakhan_ont_to$subtitle, show_x_axis = TRUE)
)

# ── Save ──────────────────────────────────────────────────────────────────────
message("Saving figures...")
ggsave(OUT_PDF,    figure_tn, width = 8.27, height = 9.5, units = "in")
ggsave(OUT_SVG,    figure_tn, width = 8.27, height = 9.5, units = "in")
ggsave(OUT_TO_PDF, figure_to, width = 8.27, height = 9.5, units = "in")
ggsave(OUT_TO_SVG, figure_to, width = 8.27, height = 9.5, units = "in")
message("Done.\n  ", OUT_PDF, "\n  ", OUT_SVG,
        "\n  ", OUT_TO_PDF, "\n  ", OUT_TO_SVG)
