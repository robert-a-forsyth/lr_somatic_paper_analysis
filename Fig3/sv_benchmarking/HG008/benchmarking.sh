#!/usr/bin/env bash
# Run from inside Fig3/sv_benchmarking/HG008/
# Requires conda environments: minda, vcf_tsv
# Requires: bcftools, tabix, python3
#
# Configurable via environment variables:
#   CONDA_SH     - path to conda.sh (default: ~/miniconda3/etc/profile.d/conda.sh)
#   MINDA_PATH   - path to minda.py (default: searches PATH for minda.py)
#   DATA_DIR     - pipeline output root (default: ../../../Run/out)
#   SV_BENCH_DIR - sv_benchmarking directory containing truth sets (default: ..)
#
# Usage:
#   cd Fig3/sv_benchmarking/HG008
#   DATA_DIR=/path/to/pipeline/out bash benchmarking.sh

set -euo pipefail

: "${CONDA_SH:=${HOME}/miniconda3/etc/profile.d/conda.sh}"
: "${MINDA_PATH:=minda.py}"
: "${DATA_DIR:=../../../Run/out}"
: "${SV_BENCH_DIR:=..}"

source "$CONDA_SH"
conda activate minda

RAW_HG008_ONT="${DATA_DIR}/HG008-ONT/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
RAW_HG008_ONT_TO="${DATA_DIR}/HG008-ONT-TO/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
RAW_HG008_PB="${DATA_DIR}/HG008-PB/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
RAW_HG008_PB_TO="${DATA_DIR}/HG008-PB-TO/variants/severus/somatic_SVs/severus_somatic.vcf.gz"

TRUTH_SET_RAW="${SV_BENCH_DIR}/GRCh38_HG008-T-V0.4_somatic-stvar_PASS.draftbenchmark.vcf.gz"
BED="${SV_BENCH_DIR}/GRCh38_HG008-T-V0.4_somatic-stvar.draftbenchmark.bed"

VAF_HG008_ONT=ONT_vaf.vcf.gz
VAF_HG008_ONT_TO=ONT_TO_vaf.vcf.gz
VAF_HG008_PB=PB_vaf.vcf.gz
VAF_HG008_PB_TO=PB_TO_vaf.vcf.gz

bcftools filter -i "VAF>.1" "$RAW_HG008_ONT_TO" -Oz -o "$VAF_HG008_ONT_TO"
bcftools filter -i "VAF>.1" "$RAW_HG008_ONT"    -Oz -o "$VAF_HG008_ONT"
bcftools filter -i "VAF>.1" "$RAW_HG008_PB"     -Oz -o "$VAF_HG008_PB"
bcftools filter -i "VAF>.1" "$RAW_HG008_PB_TO"  -Oz -o "$VAF_HG008_PB_TO"

tabix -p vcf "$VAF_HG008_ONT_TO"
tabix -p vcf "$VAF_HG008_ONT"
tabix -p vcf "$VAF_HG008_PB"
tabix -p vcf "$VAF_HG008_PB_TO"

REGION_HG008_ONT=ONT_vaf_region.vcf.gz
REGION_HG008_ONT_TO=ONT_TO_vaf_region.vcf.gz
REGION_HG008_PB=PB_vaf_region.vcf.gz
REGION_HG008_PB_TO=PB_TO_vaf_region.vcf.gz

bcftools view -R "$BED" "$VAF_HG008_ONT_TO" -Oz -o "$REGION_HG008_ONT_TO"
bcftools view -R "$BED" "$VAF_HG008_ONT"    -Oz -o "$REGION_HG008_ONT"
bcftools view -R "$BED" "$VAF_HG008_PB"     -Oz -o "$REGION_HG008_PB"
bcftools view -R "$BED" "$VAF_HG008_PB_TO"  -Oz -o "$REGION_HG008_PB_TO"

TRUTH_SET=filtered_truth.vcf.gz
bcftools view -R "$BED" "$TRUTH_SET_RAW" -Oz -o "$TRUTH_SET"

### RUN BENCHMARK
"$MINDA_PATH" truthset --base "$TRUTH_SET_RAW" \
    --vcfs "$VAF_HG008_ONT" "$VAF_HG008_ONT_TO" "$VAF_HG008_PB" "$VAF_HG008_PB_TO" \
    --min_size 50 --tolerance 500 --out_dir minda_out

"$MINDA_PATH" truthset --base "$TRUTH_SET" \
    --vcfs "$REGION_HG008_ONT" "$REGION_HG008_ONT_TO" "$REGION_HG008_PB" "$REGION_HG008_PB_TO" \
    --min_size 50 --tolerance 500 --out_dir minda_out_2

### GET TSV FILES
conda activate vcf_tsv

vcf2tsvpy --input_vcf "$VAF_HG008_ONT"    --out_tsv HG008_ONT_vaf.tsv;    sed -i '1d' HG008_ONT_vaf.tsv
vcf2tsvpy --input_vcf "$VAF_HG008_ONT_TO" --out_tsv HG008_ONT_TO_vaf.tsv; sed -i '1d' HG008_ONT_TO_vaf.tsv
vcf2tsvpy --input_vcf "$VAF_HG008_PB"     --out_tsv HG008_PB_vaf.tsv;     sed -i '1d' HG008_PB_vaf.tsv
vcf2tsvpy --input_vcf "$VAF_HG008_PB_TO"  --out_tsv HG008_PB_TO_vaf.tsv;  sed -i '1d' HG008_PB_TO_vaf.tsv
