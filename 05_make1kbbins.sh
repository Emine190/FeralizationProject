#!/bin/bash

#This script is used to make the 1kbp bins with the reads for each 1kbp combined from each original read. i.e 1 1 345 8 followed by 1 346 1000 14 now looks like 1 1 1000 22 that the lst column now is added together from all the previous to make each row now 1kbp.

INPUT_DIR="/emine190-workingdir/output/files_withgoodCov"
GENOME_FILE="/emine190-workingdir/original-files/genome_bins.bed"
OUTPUT_DIR="/emine190-workingdir/output/1000bpbins"
OUTPUTmap="/emine190-workingdir/output/map1000bpbins"
module load samtools
module load bedtools

# Loop through all .bed files in the directory
for file in "$INPUT_DIR"/*.bed; do
    if [[ -f "$file" ]]; then
        # Get the filename without the path and extension
        filename=$(basename "$file")
        filename_no_ext="${filename%_*}" # Remove .bed extension
        output_file="${OUTPUT_DIR}/${filename_no_ext}_binned.bed"
        output_map="${OUTPUTmap}/${filename_no_ext}_map.bed"
        echo "Processing $file..."

        # Use bedtools intersect to map data onto genome bins
        bedtools intersect -a "$GENOME_FILE" -b "$file" -wa -wb > "$output_file"
        bedtools map -a "$GENOME_FILE" -b "$file" -c 5 -o mean > "$output_map"
        echo "Created binned output for $filename → $output_file"
    fi
done

echo "All files processed!"
