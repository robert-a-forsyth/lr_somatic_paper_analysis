# make_figure5.R
# Combines four phased methylation locus SVGs into a labelled vertical figure
# by assembling SVG XML directly with xml2. No rasterisation at any step.
# Output: SVG + PDF (via rsvg).
#
# Run from inside Fig5/ after make_figure5_methylation.sh has produced the SVGs:
#   Rscript make_figure5.R

suppressPackageStartupMessages({
  library(xml2)
  library(rsvg)
})

# ── Input SVGs (produced by make_figure5_methylation.sh) ─────────────────────
SVG_FILES <- c(
  "CCS15-PB_5mC_CG_locus_phased.svg",   # A
  "CCS15-PB_6mA_A_locus_phased.svg",    # B
  "CCS15-ONT_5mC_CG_locus_phased.svg",  # C
  "CCS15-ONT_6mA_A_locus_phased.svg"    # D
)

TITLES <- c(
  "CCS15 PacBio 5mC Methylation",
  "CCS15 PacBio 6mA Methylation",
  "CCS15 ONT 5mC Methylation",
  "CCS15 ONT 6mA Methylation"
)

LABELS <- c("A", "B", "C", "D")

# ── Layout constants (all in mm, matching A4 = 210 × 297 mm) ─────────────────
PAGE_W   <- 210
MARGIN   <-   5
GAP      <-   3
LABEL_H  <-   7
FONT_PT  <-   9
LABEL_PT <-  10

# ── Output paths ──────────────────────────────────────────────────────────────
OUT_SVG <- "methylation_locus_panels.svg"
OUT_PDF <- "methylation_locus_panels.pdf"

# ── Helper: extract viewBox from an SVG document ──────────────────────────────
get_viewbox <- function(doc) {
  vb <- xml_attr(doc, "viewBox")
  if (!is.na(vb)) {
    v <- as.numeric(strsplit(trimws(vb), "[[:space:],]+")[[1]])
    return(list(x = v[1], y = v[2], w = v[3], h = v[4]))
  }
  w <- as.numeric(gsub("[^0-9.]", "", xml_attr(doc, "width")))
  h <- as.numeric(gsub("[^0-9.]", "", xml_attr(doc, "height")))
  list(x = 0, y = 0, w = w, h = h)
}

# ── Helper: prefix all IDs and their references inside a cloned document ──────
prefix_ids <- function(doc, prefix) {
  all_nodes <- xml_find_all(doc, ".//*")

  ids <- vapply(all_nodes, function(n) {
    id <- xml_attr(n, "id")
    if (is.na(id)) "" else id
  }, character(1))
  ids <- ids[nchar(ids) > 0]

  if (length(ids) == 0) return(invisible(NULL))

  for (n in all_nodes) {
    id <- xml_attr(n, "id")
    if (!is.na(id) && id %in% ids)
      xml_set_attr(n, "id", paste0(prefix, id))
  }

  for (n in all_nodes) {
    attrs <- xml_attrs(n)
    for (attr_name in names(attrs)) {
      val     <- attrs[[attr_name]]
      new_val <- val
      for (old_id in ids) {
        new_val <- gsub(paste0("url\\(#", old_id, "\\)"),
                        paste0("url(#", prefix, old_id, ")"), new_val)
        new_val <- gsub(paste0("^#", old_id, "$"),
                        paste0("#", prefix, old_id),          new_val)
      }
      if (!identical(new_val, val))
        xml_set_attr(n, attr_name, new_val)
    }
  }
}

# ── Read all SVGs ─────────────────────────────────────────────────────────────
message("Reading SVGs...")
docs <- lapply(SVG_FILES, read_xml)
vbs  <- lapply(docs, get_viewbox)

# ── Compute layout ────────────────────────────────────────────────────────────
panel_w  <- PAGE_W - 2 * MARGIN
panel_hs <- sapply(vbs, function(v) panel_w * (v$h / v$w))
total_h  <- MARGIN +
  sum(LABEL_H + panel_hs) +
  (length(SVG_FILES) - 1) * GAP +
  MARGIN

# ── Build root SVG ────────────────────────────────────────────────────────────
message("Assembling combined SVG...")
root <- xml_new_root("svg",
  xmlns             = "http://www.w3.org/2000/svg",
  `xmlns:xlink`     = "http://www.w3.org/1999/xlink",
  version           = "1.1",
  width             = sprintf("%.4fmm", PAGE_W),
  height            = sprintf("%.4fmm", total_h),
  viewBox           = sprintf("0 0 %.4f %.4f", PAGE_W, total_h)
)

xml_add_child(root, "rect",
  x = "0", y = "0",
  width  = as.character(PAGE_W),
  height = as.character(total_h),
  fill   = "white"
)

y <- MARGIN

for (i in seq_along(SVG_FILES)) {
  vb <- vbs[[i]]
  ph <- panel_hs[[i]]

  xml_add_child(root, "text",
    LABELS[i],
    x             = sprintf("%.4f", MARGIN),
    y             = sprintf("%.4f", y + LABEL_H * 0.75),
    `font-family` = "Helvetica, Arial, sans-serif",
    `font-size`   = sprintf("%dpt", LABEL_PT),
    `font-weight` = "bold",
    `text-anchor` = "start",
    fill          = "black"
  )

  xml_add_child(root, "text",
    TITLES[i],
    x             = sprintf("%.4f", MARGIN + LABEL_PT * 0.85),
    y             = sprintf("%.4f", y + LABEL_H * 0.75),
    `font-family` = "Helvetica, Arial, sans-serif",
    `font-size`   = sprintf("%dpt", FONT_PT),
    `font-weight` = "bold",
    `text-anchor` = "start",
    fill          = "black"
  )

  y <- y + LABEL_H

  nested <- xml_add_child(root, "svg",
    x                   = sprintf("%.4f", MARGIN),
    y                   = sprintf("%.4f", y),
    width               = sprintf("%.4f", panel_w),
    height              = sprintf("%.4f", ph),
    viewBox             = sprintf("%.4f %.4f %.4f %.4f", vb$x, vb$y, vb$w, vb$h),
    preserveAspectRatio = "xMidYMid meet"
  )

  doc_clone <- read_xml(as.character(docs[[i]]))
  prefix_ids(doc_clone, sprintf("p%d_", i))

  for (child in xml_children(doc_clone)) {
    xml_add_child(nested, child)
  }

  y <- y + ph + GAP
}

# ── Write outputs ─────────────────────────────────────────────────────────────
message("Writing SVG...")
write_xml(root, OUT_SVG)

message("Converting to PDF...")
rsvg_pdf(OUT_SVG, OUT_PDF)

message("Done.\n  ", OUT_SVG, "\n  ", OUT_PDF)
