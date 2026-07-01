library(tidyverse)
library(patchwork)

tol <- list(
  indigo = "#332288",
  cyan   = "#88CCEE",
  teal   = "#44AA99",
  green  = "#117733",
  olive  = "#999933",
  sand   = "#DDCC77",
  rose   = "#CC6677",
  wine   = "#882255",
  purple = "#AA4499"
)

base_sz <- 8

theme_kana <- function(bs = base_sz) {
  theme_bw(base_size = bs, base_family = "sans") +
    theme(
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "#3B3B4F", linewidth = 0.7, fill = NA),
      panel.background   = element_rect(fill = "white"),
      strip.background   = element_rect(fill = "#EBEBEB", color = "#AAAAAA", linewidth = 0.5),
      strip.text         = element_text(size = bs, face = "bold", color = "#444444"),
      strip.placement    = "outside",
      axis.text          = element_text(size = bs, color = "#2A2A37"),
      axis.title         = element_text(size = bs, face = "bold", color = "#1F1F28"),
      axis.ticks         = element_line(color = "#54546D", linewidth = 0.4),
      axis.line          = element_blank(),
      legend.position    = "bottom",
      legend.key.size    = unit(0.45, "cm"),
      legend.key         = element_rect(fill = NA, color = NA),
      legend.text        = element_text(size = bs, color = "#2A2A37"),
      legend.title       = element_text(size = bs, face = "bold", color = "#1F1F28"),
      legend.background  = element_rect(fill = "white", color = NA),
      legend.margin      = margin(3, 3, 3, 3),
      plot.title         = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      plot.margin        = margin(6, 8, 4, 6)
    )
}

# ── PLOT A: Benchmarking PON ─────────────────────────────────────────────────

df <- read_tsv("benchmarking_results.tsv")

keep_callers <- c("clairs", "clairsto", "deepsomatic", "consensus", "summed")

caller_pal <- c(
  "Clair-S"     = tol$indigo,
  "ClairS-TO"   = tol$cyan,
  "DeepSomatic" = tol$teal,
  "Consensus"   = tol$rose,
  "Union"       = tol$purple
)

df_filt <- df %>%
  filter(Caller %in% keep_callers) %>%
  mutate(
    CallerLabel = case_when(
      Caller == "clairs"      ~ "Clair-S",
      Caller == "clairsto"    ~ "ClairS-TO",
      Caller == "deepsomatic" ~ "DeepSomatic",
      Caller == "consensus"   ~ "Consensus",
      Caller == "summed"      ~ "Union"
    ),
    CallerLabel = factor(CallerLabel, levels = names(caller_pal)),
    TechLabel = case_when(
      Technology == "PacBio" & Mode == "paired"     ~ "PB",
      Technology == "PacBio" & Mode == "tumor_only" ~ "PB-TO",
      Technology == "ONT"    & Mode == "paired"     ~ "ONT",
      Technology == "ONT"    & Mode == "tumor_only" ~ "ONT-TO"
    ),
    TechLabel = factor(TechLabel, levels = c("PB", "ONT", "PB-TO", "ONT-TO"))
  )

df_long <- df_filt %>%
  pivot_longer(cols = c(Precision, Recall, F1), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric      = factor(Metric, levels = c("Precision", "Recall", "F1")),
    VariantType = recode(VariantType, "indel" = "Indels", "all" = "All"),
    VariantType = factor(VariantType, levels = c("SNV", "Indels", "All"))
  )

pA <- ggplot(df_long, aes(x = Metric, y = Value, fill = CallerLabel)) +
  geom_col(
    position = position_dodge(width = 0.88), width = 0.82,
    color = "#1F1F28", linewidth = 0.18
  ) +
  geom_text(
    aes(label = sprintf("%.2f", Value)),
    position = position_dodge(width = 0.88),
    angle = 90, hjust = -0.15, vjust = 0.5,
    size = 2.2, fontface = "bold", color = "black",
    nudge_y = 0.03
  ) +
  facet_grid(TechLabel ~ VariantType, switch = "y") +
  scale_y_continuous(
    limits = c(0, 1.3), breaks = seq(0, 1, 0.2),
    expand = c(0, 0), position = "right"
  ) +
  scale_fill_manual(values = caller_pal) +
  labs(x = NULL, y = "Score", fill = "Caller") +
  theme_kana() +
  guides(fill = guide_legend(nrow = 1, override.aes = list(linewidth = 0.4)))

# ── SV shared palette ────────────────────────────────────────────────────────

caller_order_sv <- c(
  "ONT Severus Benchmark",   "ONT Severus Paired",
  "PacBio Severus Benchmark","PacBio Severus Paired",
  "ONT Severus TO",          "PacBio Severus TO"
)

sv_pal <- c(
  "ONT Severus Benchmark"    = tol$indigo,
  "ONT Severus Paired"       = tol$cyan,
  "PacBio Severus Benchmark" = tol$wine,
  "PacBio Severus Paired"    = tol$sand,
  "ONT Severus TO"           = tol$teal,
  "PacBio Severus TO"        = tol$rose
)

# ── PLOT B: HG008 Structural Variants ───────────────────────────────────────

data_HG008 <- read.table(header = TRUE, text = "
Caller Precision Recall F1_Score
'ONT Severus Paired' 0.6296296296296297 0.8947368421052632 0.7391304347826088
'PacBio Severus Paired' 0.9104477611940298 0.9172932330827067 0.9138576779026217
'ONT Severus TO' 0.3034300791556728 0.8646616541353384 0.44921874999999994
'PacBio Severus TO' 0.4057971014492754 0.8421052631578947 0.5476772616136919
'ONT Severus Benchmark' 0.72 0.91 0.8043478
'PacBio Severus Benchmark' 0.93 0.94 0.9349593
") %>%
  mutate(Caller = factor(Caller, levels = caller_order_sv)) %>%
  pivot_longer(cols = c(Precision, Recall, F1_Score), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric, "F1_Score" = "F1"),
         Metric = factor(Metric, levels = c("Precision", "Recall", "F1")))

pB <- ggplot(data_HG008, aes(x = Metric, y = Value, fill = Caller)) +
  geom_col(
    position = position_dodge(width = 0.88), width = 0.82,
    color = "#1F1F28", linewidth = 0.18
  ) +
  geom_text(
    aes(label = sprintf("%.2f", Value)),
    position = position_dodge(width = 0.88),
    angle = 90, hjust = -0.15, vjust = 0.5,
    size = 2.2, fontface = "bold", color = "black",
    nudge_y = 0.03
  ) +
  scale_fill_manual(values = sv_pal, name = "Caller") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2),
                     breaks = seq(0, 1, 0.2)) +
  labs(x = "Metric", y = "Score", title = "HG008 (SV)") +
  theme_kana() +
  theme(plot.title = element_text(size = base_sz, face = "bold", hjust = 0.5,
                                  color = "#1F1F28", margin = margin(b = 4))) +
  guides(fill = guide_legend(nrow = 3, override.aes = list(linewidth = 0.4)))

# ── PLOT C: COLO829 Structural Variants ─────────────────────────────────────

data_COLO829 <- read.table(header = TRUE, text = "
Caller Precision Recall F1_Score
'ONT Severus Paired' 0.68 0.8225806 0.7445255
'PacBio Severus Paired' 0.6575342 0.7741935 0.7111111
'ONT Severus TO' 0.117506 0.7903226 0.2045929
'PacBio Severus TO' 0.1753731 0.7580645 0.2848485
'ONT Severus Benchmark' 0.720 0.868 0.787
'PacBio Severus Benchmark' 0.787 0.868 0.825
") %>%
  mutate(Caller = factor(Caller, levels = caller_order_sv)) %>%
  pivot_longer(cols = c(Precision, Recall, F1_Score), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric, "F1_Score" = "F1"),
         Metric = factor(Metric, levels = c("Precision", "Recall", "F1")))

pC <- ggplot(data_COLO829, aes(x = Metric, y = Value, fill = Caller)) +
  geom_col(
    position = position_dodge(width = 0.88), width = 0.82,
    color = "#1F1F28", linewidth = 0.18
  ) +
  geom_text(
    aes(label = sprintf("%.2f", Value)),
    position = position_dodge(width = 0.88),
    angle = 90, hjust = -0.15, vjust = 0.5,
    size = 2.2, fontface = "bold", color = "black",
    nudge_y = 0.03
  ) +
  scale_fill_manual(values = sv_pal, name = "Caller") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2),
                     breaks = seq(0, 1, 0.2)) +
  labs(x = "Metric", y = "Score", title = "COLO829 (SV)") +
  theme_kana() +
  theme(plot.title = element_text(size = base_sz, face = "bold", hjust = 0.5,
                                  color = "#1F1F28", margin = margin(b = 4))) +
  guides(fill = guide_legend(nrow = 3, override.aes = list(linewidth = 0.4)))

# ── Assemble with patchwork ──────────────────────────────────────────────────

combined <- (pA /
  ((pB | pC) +
     plot_layout(guides = "collect") &
     theme(legend.position = "bottom")) +
  plot_layout(heights = c(58, 42)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16,
                                family = "sans", margin = margin(b = 6))
    )
  )) &
  theme(plot.tag = element_text(size = 14, family = "sans", face = "bold"))

ggsave("combined_figure.pdf", combined,
       width = 6.5, height = 9, units = "in", device = cairo_pdf)
