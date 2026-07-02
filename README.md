# LRSomatic Paper Analysis

Scripts and data to reproduce all figures in the LRSomatic paper.

## Repository layout

```
lr_somatic_paper_analysis/
├── Run/                      Pipeline execution (samplesheet, Nextflow command)
│   ├── run.sh                Nextflow command to regenerate pipeline outputs
│   ├── data.csv              Samplesheet (12 samples)
│   └── download.sh           Download public BAMs (HG008, COLO829)
├── Fig1/                     Resource usage boxplots — self-contained
├── Fig2/                     QC summary — self-contained
├── Fig3/                     Variant benchmarking (SNV/indel + SV)
├── Fig4/                     Copy number variation profiles
├── Fig5/                     Methylation locus plots
├── Table1/                   Disk usage report (pre-computed; work dirs deleted)
├── environment.yml           Full dependency list
└── README.md                 This file
```

## Prerequisites

### R packages

Requires R ≥ 4.4. On VSC/HPC systems, load with `module load R/4.4.2-gfbf-2024a` (or newer).

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "patchwork",
  "rsvg", "xml2", "egg", "tidyverse"
))
```

### Conda environments

See `environment.yml` for full setup instructions. Required environments:

| Environment | Used by | Purpose |
|---|---|---|
| `methylartist` | Fig5 | Methylation locus plots |
| `minda` | Fig3 SV benchmarking | SV truth-set comparison |
| `vcf_tsv` | Fig3 SV benchmarking | VCF to TSV conversion |

### System tools

- `bcftools` and `tabix` (available via `module load BCFtools` on VSC/HPC)

### User-provided files

These large files are **not bundled** and must be provided:

| File | Used by | Where to place |
|---|---|---|
| `Homo_sapiens_assembly38_masked_noALT.fasta` | Fig3, Fig5 | `Run/data/` or set `REF=` env var |
| `gencode.v46.basic.annotation.sorted.gtf.gz` | Fig5 | `Run/data/` or set `GTF=` env var |

See `environment.yml` for download sources.

## Pipeline outputs

Figure scripts expect pipeline outputs at `../Run/out/` relative to each figure directory (i.e. `Run/out/` within this repo). If outputs are elsewhere, override with `DATA_DIR`:

```bash
DATA_DIR=/path/to/pipeline/out Rscript make_figureX.R
# or for shell scripts:
DATA_DIR=/path/to/pipeline/out bash script.sh
```

To regenerate pipeline outputs from scratch, run `Run/run.sh` with Nextflow (requires BAM files from `Run/download.sh` for public samples, and access to the CCS15 cell-line BAMs).

## Regenerating figures

All figure scripts should be run from inside their own directory.

### Figure 1 — Pipeline resource usage

```bash
cd Fig1
Rscript per_task_statistics_boxplots.R
# Outputs: per_task_statistics_combined_a4.{svg,png}
```

Inputs are bundled (`improved_report.txt`, `lrsomatic_2.0.svg`).

### Figure 2 — QC summary

```bash
cd Fig2
Rscript combined_summary_metrics.R
# Outputs: combined_qc_summary.{svg,png}
```

Inputs are bundled (`qc_metrics_summary.tsv`, `CCS15-PB_qc.txt`, `CCS15-ONT_qc.txt`).

To regenerate `qc_metrics_summary.tsv` from pipeline outputs:

```bash
cd Fig2
DATA_DIR=../Run/out bash summary_stats.sh
```

### Figure 3 — Variant benchmarking

**SNV/indel figure** (uses pre-computed `benchmarking_results.tsv`):

```bash
cd Fig3
Rscript make_figure3.R
# Output: combined_figure.pdf
```

**To re-run SNV/indel benchmarking** (requires `REF` and `bcftools`):

```bash
cd Fig3
REF=/path/to/reference.fasta bash benchmark_colo829.sh
```

**To re-run SV benchmarking** (requires `minda`, `vcf_tsv` conda envs):

```bash
cd Fig3/sv_benchmarking/COLO829
DATA_DIR=../../../Run/out \
  VAF_SEVERUS_COMP=/path/to/severus_pb.vcf.gz \
  bash benchmarking.sh

cd Fig3/sv_benchmarking/HG008
DATA_DIR=../../../Run/out bash benchmarking.sh
```

Note: `VAF_SEVERUS_COMP` (COLO829 script) can be extracted from `sv_benchmarking/variant_calls_and_benchmarks.tar.gz`.

### Figure 4 — Copy number variation profiles

```bash
cd Fig4
Rscript make_figure4.R
# Outputs: combined_cnv_plots.{pdf,svg}, combined_cnv_plots_to.{pdf,svg}
```

Reads from `DATA_DIR` (default `../Run/out`) for ASCAT and Wakhan segments. PURPLE file (`CCS15M.purple.segment.tsv`) is bundled in `Fig4/`.

### Figure 5 — Methylation locus plots

Step 1: Generate per-locus SVGs (requires `methylartist` conda env, BAM files, REF, and GTF):

```bash
cd Fig5
DATA_DIR=../Run/out \
  REF=/path/to/Homo_sapiens_assembly38_masked_noALT.fasta \
  GTF=/path/to/gencode.v46.basic.annotation.sorted.gtf.gz \
  bash make_figure5_methylation.sh
```

Step 2: Combine SVGs into final figure:

```bash
cd Fig5
Rscript make_figure5.R
# Outputs: methylation_locus_panels.{svg,pdf}
```

### Table 1 — Disk usage

Pre-computed outputs (`disk_usage_table.tsv`, `final_sample_disk_usage.tsv`) are present. The Nextflow work directories required to re-run `disk_usage_report.py` are no longer available.
