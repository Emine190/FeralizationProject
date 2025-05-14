#!/bin/bash

# --- Configuration ---
INPUT="/emine190-workingdir/output/withreadgroup"  # Directory containing the input BAM files
OUTPUT="/emine190-workingdir/output/Dedup_bam"     # Directory for output BAM and metrics files

#####load modules 

module load samtools
module load picard

# --- Loop through each BAM file in the input directory ---
find "$INPUT" -maxdepth 1 -name "*.bam" -print0 | while IFS= read -r -d $'\0' f; do #Handles filenames with spaces

  # Extract sample ID (filename without "_sorted.bam")
  name="${f%.bam}"
  id="${name##*/}"

 if [[ "$name" == *"BDA"* ]]; then
    echo "Processing sample $name..."

deduplicate_bam="$OUTPUT/$id"_dedup.bam
dedup_metrics="$OUTPUT/$id"_dedup_metrics.txt


  # --- Picard MarkDuplicates --- Run Markduplicates straigth away as thwe module is loaded. 
  picard MarkDuplicates \
    I="$f" \
    O="$deduplicate_bam" \
    M="$dedup_metrics" \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT #Often good practice to add this

  echo "Finished MarkDuplicates for sample $id"
samtools index "$deduplicate_bam" #this is done to index the file after removing duplicates.

  else
    echo "Skipping sample $name (does not contain BDA)."
  fi
done

echo "MarkDuplicates complete."
