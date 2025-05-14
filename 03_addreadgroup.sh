#!/bin/bash
# This scripts is a solution if the trouble with read groups persists.
# --- Pointers to all the correct folders ---

INPUT="/emine190-workingdir/output/BwaAligned"  # Directory containing the input BAM files
OUTPUT="/emine190-workingdir/output/withreadgroup"     # Directory for output BAM and metrics files




#--- Modules ---
module load samtools
module load bcftools
module load picard

# A find to identify the files that fill the criteria that it is a bam file and can be modified to just be specific files. 
find "$INPUT" -maxdepth 1 -name "*.bam" -print0 | while IFS= read -r -d $'\0' f; do #Handles filenames with spaces

  # Extract sample ID
  name="${f%.bam}"
  id="${name##*/}"
withreadgroup_bam="$OUTPUT/$id".bam

picard AddOrReplaceReadGroups \
       I="$f" \
       O="$withreadgroup_bam"\
       RGID=4 \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=20


done
