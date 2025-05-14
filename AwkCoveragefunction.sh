#!/bin/bash

INPUT_DIR="/emine190-workingdir/output/coverage/Cov_genes"
OUTPUT_FILE="/emine190-workingdir/output/coverage/combined_coverage_genesover0.12.txt"


bedfiles="/emine190-workingdir/output/Bedfiles"
good_coverage_sample="/emine190-workingdir/output/files_withgoodCov"

# Output file can be modified to show what cutoff for total coverage is there. 
# Clear the output file if it exists
> "$OUTPUT_FILE"

# Loop through all .bed files in the specified directory
for file in "$INPUT_DIR"/*.bed; do
  if [[ -f "$file" ]]; then
    filename=$(basename "$file") #get filename
    id="${filename%_coverage_genes.bed}"
    coverage=$(awk '{total += $8; count += ($3 - $2)} END {print total / count}' "${file}")
    if [[ $(echo "$coverage > 0.12" | bc -l) -eq 1 ]]; then #check if coverage > 0.1
      echo "$filename $coverage" >> "$OUTPUT_FILE"
      cp "$bedfiles/$id".bed "$good_coverage_sample"
    fi
  fi
done
