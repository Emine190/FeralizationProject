#!/bin/bash

# --- Configuration ---
MAIN_FOLDER="/feral_methylation/MeDiP/" # Path to the main folder containing sample subfolders
OUTPUT_FOLDER="/emine190-workingdir/output/fastp_output" # Path to the output folder. Create it if it doesn't exist
OUTPUT_FOLDERjsonhtml="/emine190-workingdir/output/fastp_output/jsonandhtml"

THREADS=10             # Number of threads to use



# --- load fastp ---
module load fastp
module load samtools

# --- Create output folder if it doesn't exist
#mkdir -p "$OUTPUT_FOLDER"

# --- Loop through sample subfolders ---
find "$MAIN_FOLDER" -mindepth 1 -maxdepth 1 -type d | while read SAMPLE_FOLDER; do
  SAMPLE_NAME=$(basename "$SAMPLE_FOLDER")

  # Extract the "true" sample name by removing "sample_" prefix (if present)
  true_sample_name="${SAMPLE_NAME#Sample_}" # Removes "sample_" from the beginning

  # Find FASTQ files using the "true" sample name
  READ1=$(find "$SAMPLE_FOLDER" -name "${true_sample_name}*L001_R1_001.fastq.gz" -o -name "${true_sample_name}*_R1_001.fastq.gz" | head -n 1)
  READ2=$(find "$SAMPLE_FOLDER" -name "${true_sample_name}*L001_R2_001.fastq.gz" -o -name "${true_sample_name}*_R2_001.fastq.gz" | head -n 1)


  if [ -z "$READ1" ] || [ -z "$READ2" ]; then
    echo "Error: FASTQ files not found for sample $SAMPLE_NAME (true name: $true_sample_name) in $SAMPLE_FOLDER"
    continue
  fi
  fastp -i "$READ1" -I "$READ2" -o "$OUTPUT_FOLDER/${true_sample_name}.R1.fastq.gz" -O "$OUTPUT_FOLDER/${true_sample_name}.R2.fastq.gz" -j "$OUTPUT_FOLDERjsonhtml/${true_sample_name}.json" -h "$OUTPUT_FOLDERjsonhtml/${true_sample_name}.html" --dont_overwrite --average_qual 25 --trim_poly_g --trim_poly_x --qualified_quality_phred 20 -t "$THREADS"

done

echo "fastp processing complete for all samples."
