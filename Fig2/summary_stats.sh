#!/bin/bash

# Define the output directory and results file
OUTPUT_DIR="gr38_CCS15/data/lr_somatic/out/"
RESULTS_FILE="qc_metrics_summary.tsv"

# Write the header row in the results file
echo -e "sample_name\tYield\tN50\tMean Coverage" > $RESULTS_FILE

# Loop through each sample directory in the output directory
find "$OUTPUT_DIR" -type f -path "*/qc/tumor/cramino_aln/*_cramino.txt" | while read -r FILE; do
    # Extract the sample name from the path
    SAMPLE_NAME=$(basename "$FILE" | sed 's/_cramino\.txt//')

    # Extract the values from the file
    YIELD=$(grep -w "Yield \[Gb\]" "$FILE" | awk 'NR==1 {print $NF}')
    N50=$(grep -w "N50" "$FILE" | awk 'NR==1 {print $NF}')
    MEAN_COVERAGE=$(grep "Mean coverage" "$FILE" | awk 'NR==1 {print $NF}')

    # Append the results to the TSV file
    echo -e "${SAMPLE_NAME}\t${YIELD}\t${N50}\t${MEAN_COVERAGE}" >> "$RESULTS_FILE"
done

echo "Extracted metrics saved to $RESULTS_FILE."