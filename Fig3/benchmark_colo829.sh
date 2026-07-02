#!/usr/bin/env bash
# ==============================================================
# COLO829 Somatic Variant Benchmarking
#
# Benchmarks ClairS/ClairSTO and DeepSomatic against a curated
# truth set across PacBio and ONT, in paired and tumour-only mode.
#
# Output: benchmarking_results.tsv
#   Columns: Sample, Technology, Mode, Caller, VariantType,
#            TP, FP, FN, Precision, Recall, F1
#
# Usage: bash benchmark_colo829.sh
# ==============================================================
set -euo pipefail

module load BCFtools

# ---- Paths ---------------------------------------------------
# DATA_DIR: root of lr_somatic pipeline output (default: ../Run/out)
# REF: path to reference genome FASTA (must be provided by user; not bundled)
: "${DATA_DIR:=../Run/out}"
: "${REF:=../Run/data/Homo_sapiens_assembly38_masked_noALT.fasta}"

TRUTH_RAW="COLO-829--COLO-829BL.snv.indel.final.v6.annotated.final.vcf.gz"

WORKDIR="benchmarking"
OUTTSV="benchmarking_results.tsv"

# ---- ClairS / ClairSTO QUAL thresholds -----------------------
# Applied to the QUAL column of ClairS/ClairSTO calls only.
# Set to 0 to disable filtering for a given category.

# PacBio Paired Mode (ClairS)
PACBIO_PAIRED_INDEL_QUAL=8
PACBIO_PAIRED_SNV_QUAL=8

# PacBio Tumor-Only Mode (ClairS-TO)
PACBIO_TO_INDEL_QUAL=8
PACBIO_TO_SNV_QUAL=8

# ONT Paired Mode (ClairS)
ONT_PAIRED_INDEL_QUAL=8
ONT_PAIRED_SNV_QUAL=8

# ONT Tumor-Only Mode (ClairS-TO)
ONT_TO_INDEL_QUAL=8
ONT_TO_SNV_QUAL=8

mkdir -p "$WORKDIR"

# ---- TSV header ----------------------------------------------
printf "Sample\tTechnology\tMode\tCaller\tVariantType\tTP\tFP\tFN\tPrecision\tRecall\tF1\n" > "$OUTTSV"


# ==============================================================
# Helper functions
# ==============================================================

# Filter to PASS/., left-align and split multi-allelic, index
prepare_vcf() {
    local raw="$1"
    local out="$2"
    bcftools view -f 'PASS,.' "$raw" \
        | bcftools norm -m -any -c w -f "$REF" -Oz -o "$out"
    bcftools index -t "$out"
}

# Split a VCF into snvs / indels / all subtypes, index each.
# Optional $3/$4: minimum QUAL threshold for SNVs / indels (0 = no filter).
# "all" is rebuilt as the concat of the (optionally filtered) snvs + indels.
split_by_type() {
    local vcf="$1"
    local prefix="$2"
    local snv_qual="${3:-0}"
    local indel_qual="${4:-0}"

    if [[ "$snv_qual" -gt 0 ]]; then
        bcftools view --type snps   -i "QUAL>=${snv_qual}"   "$vcf" -Oz -o "${prefix}_snvs.vcf.gz"
    else
        bcftools view --type snps   "$vcf" -Oz -o "${prefix}_snvs.vcf.gz"
    fi
    bcftools index -t "${prefix}_snvs.vcf.gz"

    if [[ "$indel_qual" -gt 0 ]]; then
        bcftools view --type indels -i "QUAL>=${indel_qual}" "$vcf" -Oz -o "${prefix}_indels.vcf.gz"
    else
        bcftools view --type indels "$vcf" -Oz -o "${prefix}_indels.vcf.gz"
    fi
    bcftools index -t "${prefix}_indels.vcf.gz"

    # "all" = quality-filtered SNVs + quality-filtered indels merged
    bcftools concat --allow-overlaps "${prefix}_snvs.vcf.gz" "${prefix}_indels.vcf.gz" \
        | bcftools sort -Oz -o "${prefix}_all.vcf.gz"
    bcftools index -t "${prefix}_all.vcf.gz"
}

# Count non-header variant records in a VCF
vcf_count() {
    local vcf="$1"
    if [[ -f "$vcf" ]]; then
        bcftools view -H "$vcf" | wc -l
    else
        echo 0
    fi
}

# Check that both VCFs share the same chromosome prefix (chr vs no-chr)
# Emits a WARNING to stderr if they differ — isec results would be empty
check_chr_prefix() {
    local vcf1="$1" vcf2="$2" label="$3"
    local chr1 chr2
    chr1=$(bcftools view -H "$vcf1" | head -1 | cut -f1 || true)
    chr2=$(bcftools view -H "$vcf2" | head -1 | cut -f1 || true)
    if [[ -z "$chr1" || -z "$chr2" ]]; then return; fi
    local p1 p2
    p1=$(echo "$chr1" | grep -c "^chr" || true)
    p2=$(echo "$chr2" | grep -c "^chr" || true)
    if [[ "$p1" != "$p2" ]]; then
        echo "WARNING [$label]: chr-prefix mismatch — '$chr1' vs '$chr2'. isec hits will be zero." >&2
    fi
}

# Run bcftools isec between a caller VCF and a truth VCF.
# Sets globals BENCH_TP, BENCH_FP, BENCH_FN.
run_isec() {
    local caller_vcf="$1"
    local truth_vcf="$2"
    local isec_dir="$3"

    mkdir -p "$isec_dir"
    bcftools isec -p "$isec_dir" -Oz "$caller_vcf" "$truth_vcf" 2>/dev/null

    # 0000 = caller only (FP), 0001 = truth only (FN), 0002 = both (TP from caller POV)
    BENCH_TP=$(vcf_count "${isec_dir}/0002.vcf.gz")
    BENCH_FP=$(vcf_count "${isec_dir}/0000.vcf.gz")
    BENCH_FN=$(vcf_count "${isec_dir}/0001.vcf.gz")
}

# Compute metrics and append one TSV row
write_row() {
    local sample="$1" tech="$2" mode="$3" caller="$4" vtype="$5"
    local TP="$6" FP="$7" FN="$8"
    python3 - "$sample" "$tech" "$mode" "$caller" "$vtype" "$TP" "$FP" "$FN" <<'PYEOF'
import sys
sample, tech, mode, caller, vtype = sys.argv[1:6]
TP, FP, FN = int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
prec = TP / (TP + FP) if (TP + FP) > 0 else 0.0
rec  = TP / (TP + FN) if (TP + FN) > 0 else 0.0
f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
print(f"{sample}\t{tech}\t{mode}\t{caller}\t{vtype}\t{TP}\t{FP}\t{FN}\t{prec:.4f}\t{rec:.4f}\t{f1:.4f}")
PYEOF
}

# Benchmark one caller across snvs/indels/all, append rows to TSV
benchmark_caller() {
    local sample="$1" tech="$2" mode="$3" caller_label="$4"
    local caller_prefix="$5"   # path prefix; files are ${prefix}_snvs/indels/all.vcf.gz
    local isec_base="$6"       # directory prefix for isec temp outputs

    for vtype in snvs indels all; do
        local caller_vcf="${caller_prefix}_${vtype}.vcf.gz"
        local truth_vcf="$WORKDIR/truth_${vtype}.vcf.gz"
        local isec_dir="${isec_base}_${vtype}"

        run_isec "$caller_vcf" "$truth_vcf" "$isec_dir"

        local display_type
        case "$vtype" in
            snvs)   display_type="SNV"   ;;
            indels) display_type="indel" ;;
            all)    display_type="all"   ;;
        esac

        write_row "$sample" "$tech" "$mode" "$caller_label" \
            "$display_type" "$BENCH_TP" "$BENCH_FP" "$BENCH_FN" >> "$OUTTSV"
    done
}


# ==============================================================
# Step 1: Filter and normalise truth VCF
# ==============================================================
echo "==> Preparing truth VCF..."

TRUTH_FILTERED="$WORKDIR/truth_filtered.vcf.gz"
bcftools view \
    -i 'INFO/HighConfidence=1 && INFO/num_callers >= 2 && FORMAT/AF[1:0]>=0.05 && FORMAT/DP[0]>4 && FORMAT/DP[1]>4' \
    "$TRUTH_RAW" \
    | bcftools norm -m -any -c w -f "$REF" -Oz -o "$TRUTH_FILTERED"
bcftools index -t "$TRUTH_FILTERED"

split_by_type "$TRUTH_FILTERED" "$WORKDIR/truth"

echo "    SNV:   $(vcf_count "$WORKDIR/truth_snvs.vcf.gz")"
echo "    indel: $(vcf_count "$WORKDIR/truth_indels.vcf.gz")"
echo "    all:   $(vcf_count "$WORKDIR/truth_all.vcf.gz")"


# ==============================================================
# process_sample: prepare, intersect/union, split, benchmark
#
# Args:
#   $1  sample name (e.g. COLO829-PB)
#   $2  technology  (PacBio | ONT)
#   $3  mode        (paired | tumor_only)
#   $4  clairs_combined   — somatic.vcf.gz for TO, or "" for paired
#   $5  clairs_snvs_raw   — snvs.vcf.gz for paired, or "" for TO
#   $6  clairs_indels_raw — indel.vcf.gz for paired, or "" for TO
#   $7  deepsomatic_raw
#   $8  clairs_snv_qual   — QUAL threshold for ClairS SNVs (0 = off)
#   $9  clairs_indel_qual — QUAL threshold for ClairS indels (0 = off)
# ==============================================================
process_sample() {
    local sample="$1" tech="$2" mode="$3"
    local clairs_combined="$4"
    local clairs_snvs_raw="$5"
    local clairs_indels_raw="$6"
    local ds_raw="$7"
    local clairs_snv_qual="${8:-0}"
    local clairs_indel_qual="${9:-0}"

    local sdir="$WORKDIR/$sample"
    mkdir -p "$sdir"

    echo ""
    echo "==> [$sample | $tech | $mode] Preparing caller VCFs..."

    # ---- ClairS / ClairSTO normalisation ----
    local clairs_norm="$sdir/clairs_norm.vcf.gz"
    if [[ -n "$clairs_snvs_raw" ]]; then
        # Paired mode: concatenate the separate snvs + indel files first
        local clairs_merged_raw="$sdir/clairs_merged_raw.vcf.gz"
        bcftools concat --allow-overlaps "$clairs_snvs_raw" "$clairs_indels_raw" \
            | bcftools sort -Oz -o "$clairs_merged_raw"
        bcftools index -t "$clairs_merged_raw"
        prepare_vcf "$clairs_merged_raw" "$clairs_norm"
    else
        # Tumour-only: somatic.vcf.gz already contains both SNV and indel
        prepare_vcf "$clairs_combined" "$clairs_norm"
    fi

    # ---- DeepSomatic normalisation ----
    local ds_norm="$sdir/deepsomatic_norm.vcf.gz"
    prepare_vcf "$ds_raw" "$ds_norm"

    # Sanity check: warn if chr prefix differs between callers and truth
    check_chr_prefix "$clairs_norm" "$WORKDIR/truth_all.vcf.gz" "${sample}:clairs-vs-truth"
    check_chr_prefix "$ds_norm"     "$WORKDIR/truth_all.vcf.gz" "${sample}:deepsomatic-vs-truth"

    # ---- Summed / union (all variants from either caller, deduplicated) ----
    echo "    Building summed (union)..."
    local summed_norm="$sdir/summed_norm.vcf.gz"
    bcftools concat --allow-overlaps "$clairs_norm" "$ds_norm" \
        | bcftools sort \
        | bcftools norm -d any -c w -f "$REF" -Oz -o "$summed_norm"
    bcftools index -t "$summed_norm"

    # ---- Split all four call sets by variant type ----
    echo "    Splitting by variant type..."
    # ClairS/ClairSTO: apply per-technology/mode QUAL thresholds
    split_by_type "$sdir/clairs_norm.vcf.gz"     "$sdir/clairs"     "$clairs_snv_qual" "$clairs_indel_qual"
    # DeepSomatic, summed: no additional QUAL threshold
    split_by_type "$sdir/deepsomatic_norm.vcf.gz" "$sdir/deepsomatic"
    split_by_type "$sdir/summed_norm.vcf.gz"      "$sdir/summed"

    # ---- Consensus (intersection of QUAL-filtered callers) ----
    echo "    Building consensus (intersection of QUAL-filtered)..."
    local caller_isec_dir="$sdir/caller_isec"
    mkdir -p "$caller_isec_dir"
    bcftools isec -p "$caller_isec_dir" -Oz \
        "$sdir/clairs_all.vcf.gz" \
        "$sdir/deepsomatic_all.vcf.gz" 2>/dev/null
    # 0002.vcf.gz = clairs_all records also present in deepsomatic_all
    local consensus_norm="$sdir/consensus_norm.vcf.gz"
    cp "${caller_isec_dir}/0002.vcf.gz" "$consensus_norm"
    bcftools index -t "$consensus_norm"
    split_by_type "$sdir/consensus_norm.vcf.gz" "$sdir/consensus"

    # ---- Benchmark ----
    echo "    Benchmarking..."
    local clairs_display
    [[ "$mode" == "paired" ]] && clairs_display="clairs" || clairs_display="clairsto"

    benchmark_caller "$sample" "$tech" "$mode" "$clairs_display"  \
        "$sdir/clairs"      "$sdir/bench_clairs"
    benchmark_caller "$sample" "$tech" "$mode" "deepsomatic"      \
        "$sdir/deepsomatic" "$sdir/bench_deepsomatic"
    benchmark_caller "$sample" "$tech" "$mode" "consensus"        \
        "$sdir/consensus"   "$sdir/bench_consensus"
    benchmark_caller "$sample" "$tech" "$mode" "summed"           \
        "$sdir/summed"      "$sdir/bench_summed"

    echo "    Done."
}


# ==============================================================
# Run all four COLO829 samples
# ==============================================================

# COLO829-PB  (PacBio, paired)
process_sample \
    "COLO829-PB" "PacBio" "paired" \
    "" \
    "$DATA_DIR/COLO829-PB/variants/clairs/snvs.vcf.gz" \
    "$DATA_DIR/COLO829-PB/variants/clairs/indel.vcf.gz" \
    "$DATA_DIR/COLO829-PB/variants/deepsomatic/COLO829-PB.vcf.gz" \
    "$PACBIO_PAIRED_SNV_QUAL" "$PACBIO_PAIRED_INDEL_QUAL"

# COLO829-PB-TO  (PacBio, tumour-only)
process_sample \
    "COLO829-PB-TO" "PacBio" "tumor_only" \
    "$DATA_DIR/COLO829-PB-TO/variants/clairsto/somatic.vcf.gz" \
    "" "" \
    "$DATA_DIR/COLO829-PB-TO/variants/deepsomatic/COLO829-PB-TO.vcf.gz" \
    "$PACBIO_TO_SNV_QUAL" "$PACBIO_TO_INDEL_QUAL"

# COLO829-ONT  (ONT, paired)
process_sample \
    "COLO829-ONT" "ONT" "paired" \
    "" \
    "$DATA_DIR/COLO829-ONT/variants/clairs/snvs.vcf.gz" \
    "$DATA_DIR/COLO829-ONT/variants/clairs/indel.vcf.gz" \
    "$DATA_DIR/COLO829-ONT/variants/deepsomatic/COLO829-ONT.vcf.gz" \
    "$ONT_PAIRED_SNV_QUAL" "$ONT_PAIRED_INDEL_QUAL"

# COLO829-ONT-TO  (ONT, tumour-only)
process_sample \
    "COLO829-ONT-TO" "ONT" "tumor_only" \
    "$DATA_DIR/COLO829-ONT-TO/variants/clairsto/somatic.vcf.gz" \
    "" "" \
    "$DATA_DIR/COLO829-ONT-TO/variants/deepsomatic/COLO829-ONT-TO.vcf.gz" \
    "$ONT_TO_SNV_QUAL" "$ONT_TO_INDEL_QUAL"


# ==============================================================
echo ""
echo "==> Benchmarking complete."
echo "    Results: $OUTTSV"
NROWS=$(( $(wc -l < "$OUTTSV") - 1 ))
echo "    Rows: $NROWS  (expected 48)"
if [[ "$NROWS" -ne 48 ]]; then
    echo "    WARNING: unexpected row count — check stderr for errors." >&2
fi
