#!/bin/bash


# --- Configuration ---
REF_GENOME="/referencegenomes/Gallus_gallus"  # Path to your reference genome
SAMPLE_DIR="/emine190-workingdir/output/fastp_output"      # Main directory containing ALL FASTQ files
OUTPUT_DIR="/emine190-workingdir/output/BwaAligned"          # Directory to store the output BAM files

THREADS=10                                          # Number of threads for BWA

# load in the modules 

module load bwa
module load samtools
module load bcftools

#mkdir -p "$OUTPUT_DIR"

# Find all R1 files and loop through them (since R2 should match)
find "$SAMPLE_DIR" -maxdepth 1 -name "*.R1.fastq.gz" | while read READ1; do

  # Extract SAMPLE_NAME from R1 filename (improved)
  SAMPLE_NAME=$(basename "$READ1" | sed 's/.R1.fastq.gz//')  # Remove suffix

  # Find matching R2 file (more robust)
  READ2=$(find "$SAMPLE_DIR" -maxdepth 1 -name "${SAMPLE_NAME}.R2.fastq.gz" | head -n 1)

  if [ -z "$READ2" ]; then
    echo "Error: R2 file not found for sample $SAMPLE_NAME"
    continue
  fi

  # Check if "ERC" is in the sample name (if needed)
  if [[ "$SAMPLE_NAME" == *"ERC"* ]]; then
    echo "Processing sample $SAMPLE_NAME..."

    OUTPUT_BAM="$OUTPUT_DIR/${SAMPLE_NAME}.bam"

    echo "Aligning and sorting sample $SAMPLE_NAME..."

    bwa mem -t "$THREADS" -M "$REF_GENOME" "$READ1" "$READ2" | \
      samtools view -Sbh - | \
      samtools sort -T "$OUTPUT_DIR/${SAMPLE_NAME}.tmp" -o "$OUTPUT_BAM"

    samtools index "$OUTPUT_BAM"

    echo "Finished processing sample $SAMPLE_NAME"
  else
    echo "Skipping sample $SAMPLE_NAME (does not contain ERC)."
  fi
done

echo "Alignment and sorting complete."
