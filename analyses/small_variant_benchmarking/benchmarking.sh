#!/bin/bash

set -euo pipefail  # Ensure the script exits on error and unset variables are not used

# Working directory
WORKDIR="/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/small_variant_benchmarking"
cd "$WORKDIR"

# Directories for Input files
INPUT_DIRS=(
    "/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-PB"
    "/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-ONT-TO"
    "/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-ONT"
    "/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-PB-TO"
)

RAW_TRUTH_VCF="COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf.gz"

# Functions
index_vcf_if_needed() {
    local vcf_file=$1
        echo "Indexing $vcf_file..."
        bcftools index "$vcf_file" -f
}

calculate_metrics() {
    local ISEC_DIR=$1
    if [[ -f "${ISEC_DIR}/0002.vcf" ]]; then
        TP=$(grep -vc '^#' "${ISEC_DIR}/0002.vcf") # Shared variants (TP)
    else
        TP=0
    fi
    if [[ -f "${ISEC_DIR}/0001.vcf" ]]; then
        FP=$(grep -vc '^#' "${ISEC_DIR}/0001.vcf") # Variants unique to file B (FP)
    else
        FP=0
    fi
    if [[ -f "${ISEC_DIR}/0000.vcf" ]]; then
        FN=$(grep -vc '^#' "${ISEC_DIR}/0000.vcf") # Variants unique to file A (FN)
    else
        FN=0
    fi

    local PRECISION RECALL F1
    if [[ $((TP + FP)) -eq 0 ]]; then
        PRECISION=0
    else
        PRECISION=$(awk "BEGIN {printf \"%.6f\", $TP / ($TP + $FP)}")
    fi
    if [[ $((TP + FN)) -eq 0 ]]; then
        RECALL=0
    else
        RECALL=$(awk "BEGIN {printf \"%.6f\", $TP / ($TP + $FN)}")
    fi
    if [[ $(awk "BEGIN {print ($PRECISION + $RECALL)}") == 0 ]]; then
        F1=0
    else
        F1=$(awk "BEGIN {printf \"%.6f\", 2 * $PRECISION * $RECALL / ($PRECISION + $RECALL)}")
    fi

    echo -e "$TP\t$FP\t$FN\t$PRECISION\t$RECALL\t$F1"
}

# Ensure indexes exist for truth VCF
index_vcf_if_needed "$RAW_TRUTH_VCF"

# Preprocess truth VCF
FILTERED_TRUTH="filtered_truth.vcf.gz"
TRUTH_INDELS="truth_indels.vcf.gz"
TRUTH_SNVS="truth_snvs.vcf.gz"
bcftools view -i 'INFO/HighConfidence=1 && INFO/num_callers >= 2 && FORMAT/AF[1:0]>=0.05 && FORMAT/DP[0]>4 && FORMAT/DP[1]>4' "$RAW_TRUTH_VCF" -Oz -o "$FILTERED_TRUTH"
index_vcf_if_needed "$FILTERED_TRUTH"
bcftools view -v indels "$FILTERED_TRUTH" -Oz -o "$TRUTH_INDELS"
index_vcf_if_needed "$TRUTH_INDELS"
bcftools view -v snps "$FILTERED_TRUTH" -Oz -o "$TRUTH_SNVS"
index_vcf_if_needed "$TRUTH_SNVS"

# Process each input directory
for DIR in "${INPUT_DIRS[@]}"; do
    echo "Processing directory: $DIR"

    OUTPUT_PREFIX="$(basename "$DIR")"
    ISEC_COMBINED="isec_combined_${OUTPUT_PREFIX}"
    ISEC_INDELS="isec_indels_${OUTPUT_PREFIX}"
    ISEC_SNVS="isec_snvs_${OUTPUT_PREFIX}"
    OUTPUT_TSV="${OUTPUT_PREFIX}_variant_analysis_results.tsv"

    if [[ "$DIR" == *"-TO" ]]; then
        RAW_SOMATIC="$DIR/variants/clairsto/somatic.vcf.gz"
        FILTERED_SOMATIC="filtered_${OUTPUT_PREFIX}_somatic.vcf.gz"
        FILTERED_INDELS="filtered_${OUTPUT_PREFIX}_indels.vcf.gz"
        FILTERED_SNVS="filtered_${OUTPUT_PREFIX}_snvs.vcf.gz"
        
        # Ensure index for somatic file
        index_vcf_if_needed "$RAW_SOMATIC"
        
        # Filter PASS variants and split into SNVs and Indels (add QUAL >= 13 for PacBio)
        if [[ "$DIR" == *"-PB-TO" ]]; then
            bcftools filter -i 'FILTER="PASS" && FORMAT/AF>=0.05 && QUAL>=13' "$RAW_SOMATIC" -Oz -o "$FILTERED_SOMATIC"
        else
            bcftools filter -i 'FILTER="PASS" && FORMAT/AF>=0.05' "$RAW_SOMATIC" -Oz -o "$FILTERED_SOMATIC"
        fi
        index_vcf_if_needed "$FILTERED_SOMATIC"
        bcftools view -v indels "$FILTERED_SOMATIC" -Oz -o "$FILTERED_INDELS"
        index_vcf_if_needed "$FILTERED_INDELS"
        bcftools view -v snps "$FILTERED_SOMATIC" -Oz -o "$FILTERED_SNVS"
        index_vcf_if_needed "$FILTERED_SNVS"

    else
        RAW_INDELS="$DIR/variants/clairs/indel.vcf.gz"
        RAW_SNVS="$DIR/variants/clairs/snvs.vcf.gz"
        FILTERED_INDELS="filtered_${OUTPUT_PREFIX}_indels.vcf.gz"
        FILTERED_SNVS="filtered_${OUTPUT_PREFIX}_snvs.vcf.gz"

        # Ensure indexes for input VCF files
        index_vcf_if_needed "$RAW_INDELS"
        index_vcf_if_needed "$RAW_SNVS"

        # Filter PASS variants (add QUAL >= 13 for PacBio data)
        if [[ "$DIR" == *"-PB" ]]; then
            bcftools filter -i 'FILTER="PASS" && FORMAT/AF>=0.05 && QUAL>=13' "$RAW_INDELS" -Oz -o "$FILTERED_INDELS"
            index_vcf_if_needed "$FILTERED_INDELS"
            bcftools filter -i 'FILTER="PASS" && FORMAT/AF>=0.05 && QUAL>=13' "$RAW_SNVS" -Oz -o "$FILTERED_SNVS"
            index_vcf_if_needed "$FILTERED_SNVS"
        else
            bcftools filter -i 'FILTER="PASS" && FORMAT/AF>=0.05' "$RAW_INDELS" -Oz -o "$FILTERED_INDELS"
            index_vcf_if_needed "$FILTERED_INDELS"
            bcftools filter -i 'FILTER="PASS" && FORMAT/AF>=0.05' "$RAW_SNVS" -Oz -o "$FILTERED_SNVS"
            index_vcf_if_needed "$FILTERED_SNVS"
        fi
    fi

    # Combine SNVs and indels
    COMBINED_FILTERED="combined_${OUTPUT_PREFIX}_variants.vcf.gz"
    bcftools concat -o "$COMBINED_FILTERED" -O z "$FILTERED_SNVS" "$FILTERED_INDELS" -a
    index_vcf_if_needed "$COMBINED_FILTERED"

    # Intersections
    mkdir -p "$ISEC_COMBINED" "$ISEC_INDELS" "$ISEC_SNVS"
    bcftools isec "$FILTERED_TRUTH" "$COMBINED_FILTERED" -p "$ISEC_COMBINED"
    bcftools isec "$TRUTH_INDELS" "$FILTERED_INDELS" -p "$ISEC_INDELS"
    bcftools isec "$TRUTH_SNVS" "$FILTERED_SNVS" -p "$ISEC_SNVS"

    # Output metrics
    {
        echo -e "Category\tTP\tFP\tFN\tPrecision\tRecall\tF1"
        echo -e "Combined\t$(calculate_metrics "$ISEC_COMBINED")"
        echo -e "Indels\t$(calculate_metrics "$ISEC_INDELS")"
        echo -e "SNVs\t$(calculate_metrics "$ISEC_SNVS")"
    } > "$OUTPUT_TSV"

    echo "Results saved to $OUTPUT_TSV"
done

echo "All analyses completed."