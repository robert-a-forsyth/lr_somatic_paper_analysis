#!/usr/bin/env bash
# Run from inside Fig3/sv_benchmarking/COLO829/
# Requires conda environments: minda, vcf_tsv
# Requires: bcftools, python3
#
# Configurable via environment variables:
#   CONDA_SH      - path to conda.sh (default: ~/miniconda3/etc/profile.d/conda.sh)
#   MINDA_PATH    - path to minda.py (default: searches PATH for minda.py)
#   DATA_DIR      - pipeline output root (default: ../../../Run/out)
#   SV_BENCH_DIR  - sv_benchmarking directory containing truth sets and add_vaf.py (default: ..)
#   VAF_SEVERUS_COMP - path to a comparison Severus PB VCF (required for minda benchmarking)
#                      Extract from sv_benchmarking/variant_calls_and_benchmarks.tar.gz if needed
#
# Usage:
#   cd Fig3/sv_benchmarking/COLO829
#   DATA_DIR=/path/to/pipeline/out VAF_SEVERUS_COMP=/path/to/severus_pb.vcf.gz bash benchmarking.sh

set -euo pipefail

: "${CONDA_SH:=${HOME}/miniconda3/etc/profile.d/conda.sh}"
: "${MINDA_PATH:=minda.py}"
: "${DATA_DIR:=../../../Run/out}"
: "${SV_BENCH_DIR:=..}"

source "$CONDA_SH"
conda activate minda

RAW_COLO_ONT="${DATA_DIR}/COLO829-ONT/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
RAW_COLO_ONT_TO="${DATA_DIR}/COLO829-ONT-TO/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
RAW_COLO_PB="${DATA_DIR}/COLO829-PB/variants/severus/somatic_SVs/severus_somatic.vcf.gz"
RAW_COLO_PB_TO="${DATA_DIR}/COLO829-PB-TO/variants/severus/somatic_SVs/severus_somatic.vcf.gz"

TRUTH_SET="${SV_BENCH_DIR}/truthset_somaticSVs_COLO829_hg38lifted_chr_alt_vaf.vcf.gz"
VAF_SCRIPT_PATH="${SV_BENCH_DIR}/add_vaf.py"

### ADD VAF TO RAW FILES
VAF_COLO_ONT=ONT_vaf.vcf.gz
VAF_COLO_ONT_TO=ONT_TO_vaf.vcf.gz
VAF_COLO_PB=PB_vaf.vcf.gz
VAF_COLO_PB_TO=PB_TO_vaf.vcf.gz

python3 "$VAF_SCRIPT_PATH" "$RAW_COLO_ONT"    "$VAF_COLO_ONT"
python3 "$VAF_SCRIPT_PATH" "$RAW_COLO_ONT_TO" "$VAF_COLO_ONT_TO"
python3 "$VAF_SCRIPT_PATH" "$RAW_COLO_PB"     "$VAF_COLO_PB"
python3 "$VAF_SCRIPT_PATH" "$RAW_COLO_PB_TO"  "$VAF_COLO_PB_TO"

### RUN BENCHMARK
# VAF_SEVERUS_COMP is a comparison Severus PB VCF; extract from
# ../variant_calls_and_benchmarks.tar.gz if not available elsewhere
: "${VAF_SEVERUS_COMP:?VAF_SEVERUS_COMP must be set (path to severus_pb comparison vcf.gz)}"

"$MINDA_PATH" truthset --base "$TRUTH_SET" \
    --vcfs "$VAF_SEVERUS_COMP" "$VAF_COLO_ONT" "$VAF_COLO_ONT_TO" "$VAF_COLO_PB" "$VAF_COLO_PB_TO" \
    --min_size 50 --tolerance 500 --out_dir minda_out --vaf 0.1

### GET TSV FILES
conda activate vcf_tsv

vcf2tsvpy --input_vcf "$VAF_COLO_ONT"    --out_tsv COLO_ONT_vaf.tsv;    sed -i '1d' COLO_ONT_vaf.tsv
vcf2tsvpy --input_vcf "$VAF_COLO_ONT_TO" --out_tsv COLO_ONT_TO_vaf.tsv; sed -i '1d' COLO_ONT_TO_vaf.tsv
vcf2tsvpy --input_vcf "$VAF_COLO_PB"     --out_tsv COLO_PB_vaf.tsv;     sed -i '1d' COLO_PB_vaf.tsv
vcf2tsvpy --input_vcf "$VAF_COLO_PB_TO"  --out_tsv COLO_PB_TO_vaf.tsv;  sed -i '1d' COLO_PB_TO_vaf.tsv

### FILTER VAF < 0.9 AND RE-BENCHMARK
VAF_COLO_ONT2=ONT_vaf2.vcf.gz
VAF_COLO_ONT_TO2=ONT_TO_vaf2.vcf.gz
VAF_COLO_PB2=PB_vaf2.vcf.gz
VAF_COLO_PB_TO2=PB_TO_vaf2.vcf.gz

bcftools filter -i "INFO/VAF<.9" "$VAF_COLO_ONT"    -Oz -o "$VAF_COLO_ONT2"
bcftools filter -i "INFO/VAF<.9" "$VAF_COLO_ONT_TO" -Oz -o "$VAF_COLO_ONT_TO2"
bcftools filter -i "INFO/VAF<.9" "$VAF_COLO_PB"     -Oz -o "$VAF_COLO_PB2"
bcftools filter -i "INFO/VAF<.9" "$VAF_COLO_PB_TO"  -Oz -o "$VAF_COLO_PB_TO2"

conda activate minda

"$MINDA_PATH" truthset --base "$TRUTH_SET" \
    --vcfs "$VAF_COLO_ONT2" "$VAF_COLO_ONT_TO2" "$VAF_COLO_PB2" "$VAF_COLO_PB_TO2" \
    --min_size 50 --tolerance 500 --out_dir minda_out2 --vaf 0.1
