#!/usr/bin/env bash
# Run from inside Fig5/
# Requires: conda env 'methylartist' with methylartist installed
#
# Override defaults:
#   DATA_DIR=../Run/out
#   REF=../Run/data/Homo_sapiens_assembly38_masked_noALT.fasta
#   GTF=../Run/data/gencode.v46.basic.annotation.sorted.gtf.gz

set -euo pipefail

: "${DATA_DIR:=../Run/out}"
: "${REF:=../Run/data/Homo_sapiens_assembly38_masked_noALT.fasta}"
: "${GTF:=../Run/data/gencode.v46.basic.annotation.sorted.gtf.gz}"

ONT_BAM="${DATA_DIR}/CCS15-ONT/bamfiles/CCS15-ONT_tumor.bam"
PB_BAM="${DATA_DIR}/CCS15-PB/bamfiles/CCS15-PB_tumor.bam"
LOCUS="chr14:100825400-100826941"

# 5mC ONT phased SVG
conda run -n methylartist python methylartist locus \
  -b "$ONT_BAM" \
  -i "$LOCUS" \
  -r "$REF" \
  -g "$GTF" \
  -n CG --motifsize 2 -m m --phased --svg \
  --width 8.27 --height 2.92 \
  -o CCS15-ONT_5mC_CG_locus_phased.svg

# 5mC PB phased SVG
conda run -n methylartist python methylartist locus \
  -b "$PB_BAM" \
  -i "$LOCUS" \
  -r "$REF" \
  -g "$GTF" \
  -n CG --motifsize 2 -m m --phased --svg \
  --width 8.27 --height 2.92 \
  -o CCS15-PB_5mC_CG_locus_phased.svg

# 6mA ONT phased SVG
conda run -n methylartist python methylartist locus \
  -b "$ONT_BAM" \
  -i "$LOCUS" \
  -r "$REF" \
  -g "$GTF" \
  -n A --motifsize 1 -m a --infer_implicit --force_implicit --phased --svg \
  --width 8.27 --height 2.92 \
  -o CCS15-ONT_6mA_A_locus_phased.svg

# 6mA PB phased SVG
conda run -n methylartist python methylartist locus \
  -b "$PB_BAM" \
  -i "$LOCUS" \
  -r "$REF" \
  -g "$GTF" \
  -n A --motifsize 1 -m a --infer_implicit --force_implicit --phased --svg \
  --width 8.27 --height 2.92 \
  -o CCS15-PB_6mA_A_locus_phased.svg
