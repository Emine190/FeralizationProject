#!/bin/bash


#-- All pointers for static files or directories --
INPUT="/emine190-workingdir/output/Dedup_bam"

REF_transcripts="/referencegenomes/Gallus_gallustranscripts.bed"

REF_genes="" 

REF_exons=""

OUTPUT_Cov_Exons="emine190-workingdir/output/coverage/Cov_exons" #path to output
OUTPUT_Cov_Transcript="/emine190-workingdir/output/coverage/Cov_transcript" #path to output
OUTPUT_Cov_Genes="/emine190-workingdir/output/coverage/Cov_genes" #path to output

OUTPUT_stats="/emine190-workingdir/output/Bam_stats"
OUTPUT="/emine190-workingdir/output/Bedfiles"
#-- load modules --

module load bcftools
module load samtools
module load bedtools


find "$INPUT" -maxdepth 1 -name "*.bam" -print0 | while IFS= read -r -d $'\0' f; do #Handles filenames with spaces

name="${f%.bam}"
  id="${name##*/}"

output_bed="$OUTPUT/$id".bed

#samtools stats $f > $OUTPUT_stats/$id"_bam_stats.txt"

echo $id
#bedtools bamtobed -i "$f" > "$output_bed"

#bedtools coverage -a $REF_transcripts -b $output_bed > "$OUTPUT_Cov_Transcript/$id"_coverage_transcripts.bed

bedtools coverage -a $REF_genes -b $output_bed > "$OUTPUT_Cov_Genes/$id"_coverage_genes.bed

#bedtools coverage -a $REF_exons -b $output_bed > "$OUTPUT_Cov_Exons/$id"_coverage_exons.bed

echo "done with $id starting Next sample"

done
