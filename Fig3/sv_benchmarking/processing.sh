#!/usr/bin/env bash
# Run from inside Fig3/sv_benchmarking/
# Filters COLO829 PacBio paired Severus SV calls (PASS filter, VAF >= 0.1)
#
# Configurable via environment variables:
#   DATA_DIR - pipeline output root (default: ../../Run/out)
#
# Usage:
#   cd Fig3/sv_benchmarking
#   DATA_DIR=/path/to/pipeline/out bash processing.sh

: "${DATA_DIR:=../../Run/out}"

COLO829_PB_PAIR_SEVERUS="${DATA_DIR}/COLO829-PB/variants/severus/somatic_SVs/severus_somatic.vcf.gz"

bcftools filter -i 'FILTER="PASS"' "$COLO829_PB_PAIR_SEVERUS" \
| bcftools filter -i 'VAF>=0.1' | bcftools view -s COLO829-PB_tumor -o COLO829-PB_severus.filtered.vcf
