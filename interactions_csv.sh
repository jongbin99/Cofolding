#!/bin/bash

# Define the base directory where working directories are located
BASE_WORKING_DIR="/mnt/nfs/home/jkim/work/ifp/work_dir"

# Define the output file
OUTPUT_FILE="combined_interactions.csv"

# Find all .interactions.csv files inside subdirectories of working_####
csv_files=($(find "$BASE_WORKING_DIR"/working_*/* -type f -name "*.interactions.csv"))

# Check if there are any CSV files found
if [[ ${#csv_files[@]} -eq 0 ]]; then
    exit 1
fi

# Get the header from the first CSV file
head -n 1 "${csv_files[0]}" > "$OUTPUT_FILE"

# Append all other CSV contents (excluding the header) to the output file
for file in "${csv_files[@]}"; do
    tail -n +2 "$file" >> "$OUTPUT_FILE"
done

