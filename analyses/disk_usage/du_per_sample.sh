#!/bin/bash

# Input file path
input_file="/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/pipeline_info/execution_trace_2025-12-04_16-47-42.txt"

# Base directory to calculate disk usage
base_dir="/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/work"

# Log file path
log_file="disk_usage_debug.log"
> "$log_file" # Clear the log file before starting

# Initialize dictionaries to track disk usage and hashes for each sample
declare -A sample_hash_map
declare -A sample_disk_usage

# Read input file and build a sample-to-hash mapping
while IFS=$'\t' read -r task_id hash _ name _; do
    if [[ $task_id != "task_id" ]]; then
        # Extract sample name from parentheses in the "name" column
        sample=$(echo "$name" | grep -oP '\(\K[^\)]+') # Extract text inside parentheses
        sample_hash_map["$sample"]+="$hash "
    fi
done < "$input_file"

# Iterate over the sample-hash mapping to calculate disk usage
for sample in "${!sample_hash_map[@]}"; do
    echo "Processing sample: $sample" | tee -a "$log_file"
    echo "---------------------------------" | tee -a "$log_file"
    total_usage=0  # Running total for this sample

    for hash_prefix in ${sample_hash_map[$sample]}; do
        # Use find to locate the **full hash directory**
        echo "Searching for directories starting with hash prefix: $hash_prefix" | tee -a "$log_file"
        
        # Dynamic directory search for the correct full hash directory
        dir=$(find "$base_dir" -type d -path "*/$hash_prefix*" -mindepth 2 -maxdepth 2 -print -quit 2>/dev/null)
        
        if [[ -n "$dir" ]]; then
            echo "Matched directory: $dir to hash prefix $hash_prefix (sample: $sample)" | tee -a "$log_file"

            # Calculate disk usage for the parent directory only
            usage=$(du -sb "$dir" 2>/dev/null | awk '{ print $1 }')
            total_usage=$((total_usage + usage))
            human_readable_usage=$(numfmt --to=iec "$usage")
            running_total_hr=$(numfmt --to=iec "$total_usage")

            # Output detailed log
            echo "Directory: $dir" | tee -a "$log_file"
            echo "Sample: $sample" | tee -a "$log_file"
            echo "Usage: $human_readable_usage" | tee -a "$log_file"
            echo "Running Total (Sample): $running_total_hr" | tee -a "$log_file"
            echo "---------------------------------" | tee -a "$log_file"
        else
            # Log no-directory-found errors
            echo "No directory found starting with hash prefix: $hash_prefix" | tee -a "$log_file"
            echo "[DEBUG] Checking manually in base_dir: $base_dir"
            find "$base_dir" -type d -name "${hash_prefix}*" -mindepth 1 -maxdepth 3 | tee -a "$log_file"
        fi
    done
    # Store the total usage for the sample
    sample_disk_usage["$sample"]=$total_usage
done

# Display the final disk usage summary for each sample and the total across all samples
echo "Final Summary:" | tee -a "$log_file"
grand_total=0
for sample in "${!sample_disk_usage[@]}"; do
    sample_total=${sample_disk_usage[$sample]}
    grand_total=$((grand_total + sample_total))
    total_human_readable=$(numfmt --to=iec "$sample_total")
    echo "Sample: $sample, Total Disk Usage: $total_human_readable" | tee -a "$log_file"
done

# Print the grand total
grand_total_human_readable=$(numfmt --to=iec "$grand_total")
echo "---------------------------------" | tee -a "$log_file"
echo "Grand Total (Overall Disk Usage): $grand_total_human_readable" | tee -a "$log_file"