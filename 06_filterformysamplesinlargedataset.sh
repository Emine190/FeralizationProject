#!/bin/bash

# Input data file
DATA_FILE="/emine190-workingdir/original-files/cc5k.autosomes.selectFERAL.renamesamples.biallelic.phase.filtermafmiss.genotypematrixNUMformat.txt"

# --- Configuration ---
# *** USER: SET THESE VARIABLES ***

# Path to the main large data file (e.g., genotypes)
INPUT_DATA_FILE="/emine190-workingdir/original-files/cc5k.autosomes.selectFERAL.renamesamples.biallelic.phase.filtermafmiss.genotypematrixNUMformat.txt" 

# Path to the text file containing the list of Sample IDs (headers) to keep
# (One ID per line)
ID_LIST_FILE="/emine190-workingdir/output/map1000bpbins/samples_to_keep.txt" 

# Path for the output file with the filtered columns
OUTPUT_FILE="/emine190-workingdir/original-files/cc5k.autosomes.selectFERAL.renamesamples.biallelic.phase.filtermafmiss.genotypematrixNUMformat_emine190version.txt"

# *** END USER SETTINGS ***

# Check if input files exist
if [ ! -f "$INPUT_DATA_FILE" ]; then
  echo "Error: Input data file '$INPUT_DATA_FILE' not found."
  exit 1
fi

if [ ! -f "$ID_LIST_FILE" ]; then
  echo "Error: ID list file '$ID_LIST_FILE' not found."
  exit 1
fi

echo "Input Data File: $INPUT_DATA_FILE"
echo "ID List File:    $ID_LIST_FILE"
echo "Output File:     $OUTPUT_FILE"
echo "Processing..."

# Use awk to perform the filtering
# Explanation of the awk script:
#   -F'\t'             : Set input field separator to Tab
#   BEGIN{OFS="\t"}    : Set output field separator to Tab
#   'FILENAME==ARGV[1]': Action block executed only for the first file (ID list)
#     { wanted_ids[$1]=1; next } : Store each ID from the first file in the 'wanted_ids' array.
#                                   'next' skips to the next line/file.
#   'FNR==1'           : Action block executed only for the first line of the second file (header)
#     { ... }            : Print the first header ($1). Then loop through the rest ($2 to $NF).
#                          If a header is in 'wanted_ids', print it and store its column index 'i'
#                          in the 'keep_idx' array. Finally, print a newline and 'next'.
#   'FNR>1'            : Action block executed for all data lines (after header) of the second file
#     { ... }            : Print the first field ($1). Then loop through the stored indices in 'keep_idx'.
#                          For each stored index 'j', retrieve the actual column number 'col_index'.
#                          Print the data from that specific column '$col_index'. Finally, print a newline.

awk -F'\t' '
BEGIN { OFS="\t" }
FILENAME==ARGV[1] {
    wanted_ids[$1]=1
    next
}
FNR==1 {
    printf "%s", $1 # Print first column header (e.g., "SNP")
    num_wanted = 0
    for (i=2; i<=NF; i++) { # Loop through sample headers (column 2 onwards)
        if ($i in wanted_ids) {
            printf "%s%s", OFS, $i # Print wanted sample header
            keep_idx[++num_wanted] = i # Store the column index to keep
        }
    }
    printf "\n" # End the header line
    next # Move to the next line of the input data file
}
FNR>1 {
    printf "%s", $1 # Print first column data (e.g., SNP ID)
    for (j=1; j<=num_wanted; j++) {
        col_index = keep_idx[j]
        printf "%s%s", OFS, $(col_index) # Print data from the wanted column
    }
    printf "\n" # End the data line
}
' "$ID_LIST_FILE" "$INPUT_DATA_FILE" > "$OUTPUT_FILE"
#  ^----------------^------------------^------------- Pass ID list first, then data file

# Check awk's exit status
if [ $? -eq 0 ]; then
  echo "Processing complete. Filtered data saved to '$OUTPUT_FILE'."
else
  echo "Error during awk processing. Check input files and formats."
  # Optional: rm -f "$OUTPUT_FILE" # Remove potentially incomplete output
  exit 1
fi

exit 0
