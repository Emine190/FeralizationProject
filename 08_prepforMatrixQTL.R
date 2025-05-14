# ============================================================================
# Script 1: Preprocess Individual BED Files to MatrixEQTL Input Format
# (Version 3: Assumes Identical Bins, Specific Filename Structure, Filters for Autosomes)
# ============================================================================
#
# Description:
# Reads multiple BED files (chr, start, end, score).
# **Assumption**: All input BED files contain the exact same set of bins.
# Filters to keep only bins on **numerical chromosomes (autosomes)**.
# Extracts sample IDs (assuming simple .bed suffix removal based on user input).
# Compiles data into two files:
#   1. Autosomal methylation score matrix (bins x samples).
#   2. Autosomal bin location file (bin_id, chr, start, end).
#
# Requirements: R, data.table package.
# Input Format: BED files (tab-sep, no header, 4 cols: chr, start, end, score).
#
# ============================================================================

# --- 1. Load Libraries ---

library(data.table)

# --- 2. Define Parameters ---
# *** SET THESE VARIABLES ***

# Directory containing the input BED files
bed_file_directory <- "/cfs/klemming/projects/supr/sllstore2017078/emine190-workingdir/output/map1000bpbins/"

# Set working directory for outputs
setwd("/cfs/klemming/projects/supr/sllstore2017078/emine190-workingdir/output/QTLRanalysis/")
message("Output files will be saved to: ", getwd())

# Pattern to identify relevant BED files within the directory.
# Use "\\.bed$" to match any file ending in .bed
bed_file_pattern <- "\\.bed$"

# Base name for the output files
output_basename <- "methylation"

# *** END USER SETTINGS ***

# Construct full output filenames reflecting the filtering
output_matrix_file <- paste0(output_basename, "_dataNOsexchromo.txt")
output_location_file <- paste0(output_basename, "_locationsNosexchromo.txt")

# --- 3. List Sample Files and Extract Sample IDs ---
message("\nSearching for BED files matching pattern '", bed_file_pattern, "'...")
bed_files <- list.files(path = bed_file_directory,
                        pattern = bed_file_pattern,
                        full.names = TRUE,
                        ignore.case = TRUE) # Allow case-insensitivity for pattern

if (length(bed_files) == 0) {
  stop(paste("FATAL ERROR: No files found matching pattern '", bed_file_pattern,
             "' in directory '", bed_file_directory,
             "'. Please check the path and pattern.", sep=""))
}

message("Found ", length(bed_files), " files potentially representing samples.")

# Extract sample IDs by removing the .bed suffix (as per user's provided script)
# If a more complex pattern is needed (like the VE-3292- prefix/suffix), revise this gsub.
sample_ids <- gsub("\\.bed$", "", basename(bed_files), ignore.case = TRUE)

# Check for duplicate sample IDs
if (any(duplicated(sample_ids))) {
    warning("Duplicate sample IDs detected after extraction! Check file naming consistency.")
    print(sample_ids[duplicated(sample_ids)])
}
names(bed_files) <- sample_ids # Name the file list with sample IDs

message("\nExtracted Sample IDs and corresponding files:")
print(data.frame(SampleID = sample_ids, FileName = basename(bed_files)))
message("-> IMPORTANT: Verify these Sample IDs match those in your genotype data header!")


# --- 4. Define Bins based on the First File AND Filter for Autosomes ---
message("\nDefining bins based on the first file (assuming all files have identical bins)...")
first_file_path <- bed_files[1]
message("Reading bins from: ", basename(first_file_path))

tryCatch({
    unique_bins <- fread(first_file_path, sep = "\t", header = FALSE,
                         col.names = c("chr", "start", "end", "score"))

    # Basic validation
    if (nrow(unique_bins) == 0) stop("First file is empty.")
    if (ncol(unique_bins) != 4) stop("First file does not have 4 columns.")

    # Ensure numeric types for coordinates
    if(!is.numeric(unique_bins$start)) unique_bins[, start := as.numeric(as.character(start))]
    if(!is.numeric(unique_bins$end)) unique_bins[, end := as.numeric(as.character(end))]
    if (any(is.na(unique_bins$start)) || any(is.na(unique_bins$end))) {
        stop("FATAL ERROR: NA values found in start/end coordinates of the first file. Cannot define bins.")
    }

    # Keep only location columns for defining bins
    unique_bins <- unique_bins[, .(chr, start, end)]

    # ** FILTERING STEP: Keep only numerical chromosomes (autosomes) **
    message("Filtering bins to keep only autosomes (numerical chromosomes)...")
    # Create a standardized chromosome representation (remove "chr" prefix if present)
    unique_bins[, chr_standardized := gsub("chr", "", chr, ignore.case = TRUE)]
    # Check if the standardized chromosome name can be converted to an integer
    unique_bins[, is_numeric_chr := !is.na(suppressWarnings(as.integer(chr_standardized)))]

    original_bin_count <- nrow(unique_bins)
    # Subset the data table to keep only rows where is_numeric_chr is TRUE
    unique_bins <- unique_bins[is_numeric_chr == TRUE]
    filtered_bin_count <- nrow(unique_bins)

    if (filtered_bin_count == 0) {
        stop("FATAL ERROR: No bins remaining after filtering for numerical chromosomes. Check chromosome naming in BED files (e.g., '1', '2' or 'chr1', 'chr2').")
    }
    message("Filtering complete: Kept ", filtered_bin_count, " autosomal bins out of ", original_bin_count, " total defined bins.")

    # Remove temporary columns used for filtering
    unique_bins[, chr_standardized := NULL]
    unique_bins[, is_numeric_chr := NULL]

    # Create the unique bin ID for the remaining autosomal bins
    unique_bins[, bin_id := paste(chr, as.integer(start), as.integer(end), sep = "_")]

    # Sort the remaining autosomal bins
    # Convert chromosome (now known to be numeric) to integer for proper sorting
    unique_bins[, chr_num := as.integer(gsub("chr", "", chr, ignore.case = TRUE))] # Re-extract numeric part for sorting
    setorder(unique_bins, chr_num, start)
    unique_bins[, chr_num := NULL] # Remove temporary sort column
    setkey(unique_bins, bin_id) # Key for efficient joining later

    message("Successfully defined and filtered ", nrow(unique_bins), " autosomal bins.")

}, error = function(e) {
    stop("FATAL ERROR processing the first file ('", basename(first_file_path), "'): ", e$message)
})


# --- 5. Write Autosomal Methylation Location File ---
message("\nWriting filtered (autosomal) methylation location file...")
# Select final columns for the location file (already sorted)
location_output <- unique_bins[, .(bin_id, chr, start, end)]

# Write the file
fwrite(location_output, file = output_location_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
message("Autosomal location file saved to: ", output_location_file)


# --- 6. Read Scores for All Samples (for Autosomal Bins) ---
sample_data_list <- list() # To store scores per sample
processed_samples_count <- 0
skipped_files <- c()

message("\nReading scores from all sample files...")
# Get the set of AUTOSOMAL bin IDs we want to keep scores for
autosomal_bin_ids <- unique_bins$bin_id

for (i in seq_along(bed_files)) {
  sample_id <- names(bed_files)[i]
  file_path <- bed_files[i]
  message("Processing scores: ", sample_id, " (", i, "/", length(bed_files), ") - ", basename(file_path))

  tryCatch({
    # Read the current sample file
    dt <- fread(file_path, sep = "\t", header = FALSE,
                col.names = c("chr", "start", "end", "score"),
                check.names = FALSE, showProgress = FALSE)

    # Basic validation and cleaning (as before)
    if (nrow(dt) == 0 || ncol(dt) != 4) {
        message("  WARNING: File empty or has incorrect columns. Skipping.")
        skipped_files <- c(skipped_files, basename(file_path))
        sample_data_list[[sample_id]] <- data.table(bin_id=character(0), score=numeric(0))
        next
    }
    if(!is.numeric(dt$score)) dt[, score := as.numeric(as.character(score))]
    dt <- dt[!is.na(score)]
    if(!is.numeric(dt$start)) dt[, start := as.numeric(as.character(start))]
    if(!is.numeric(dt$end)) dt[, end := as.numeric(as.character(end))]
    dt <- dt[!is.na(start) & !is.na(end)]
     if(nrow(dt) == 0) {
        message("  WARNING: No valid rows remain after cleaning NAs. Skipping.")
        skipped_files <- c(skipped_files, basename(file_path))
        sample_data_list[[sample_id]] <- data.table(bin_id=character(0), score=numeric(0))
        next
     }

    # Create bin IDs consistently
    dt[, bin_id := paste(chr, as.integer(start), as.integer(end), sep = "_")]

    # ** FILTER SCORES **: Keep only scores for the autosomal bins defined earlier
    dt_filtered <- dt[bin_id %in% autosomal_bin_ids, .(bin_id, score)]

    if (nrow(dt_filtered) == 0) {
         message("  WARNING: No scores found corresponding to the defined autosomal bins in this file.")
         # Still store an empty table to avoid errors later, but note it
         sample_data_list[[sample_id]] <- data.table(bin_id=character(0), score=numeric(0))
    } else {
         setkey(dt_filtered, bin_id)
         sample_data_list[[sample_id]] <- dt_filtered
    }

    processed_samples_count <- processed_samples_count + 1

  }, error = function(e) {
    message("  ERROR processing scores file: ", basename(file_path), " - ", e$message)
    skipped_files <- c(skipped_files, basename(file_path))
    sample_data_list[[sample_id]] <- data.table(bin_id=character(0), score=numeric(0)) # Store empty table on error
  })
} # End loop through files

# Report issues
if(length(skipped_files) > 0) {
    message("\nWarning: Issues encountered with the following files (may result in NAs or exclusion): ", paste(unique(skipped_files), collapse=", "))
}
if (processed_samples_count != length(bed_files)) {
    message("\nWarning: Attempted to process scores for ", processed_samples_count, " files, but found ", length(bed_files), " total files matching pattern.")
}


# --- 7. Create Autosomal Methylation Score Matrix ---
message("\nCreating filtered (autosomal) methylation score matrix...")

# Use the ORDERED unique AUTOSOMAL bins defined earlier
autosomal_bin_ids_ordered <- unique_bins[, .(bin_id)] # Get bin_ids in the final sorted order
setkey(autosomal_bin_ids_ordered, bin_id)

# Start with the ordered autosomal bin IDs as the base
final_matrix_dt <- autosomal_bin_ids_ordered

# Iterate through each sample for which we have score data
processed_sample_ids <- names(sample_data_list)
for (sample_id in processed_sample_ids) {
  message("  Adding scores for sample: ", sample_id)
  sample_scores <- sample_data_list[[sample_id]] # This now contains only scores for autosomal bins

  if (is.null(sample_scores) || nrow(sample_scores) == 0) {
       message("    No autosomal scores found/processed for this sample. Adding NA column.")
       final_matrix_dt[, (sample_id) := NA_real_]
       next
  }

  setkey(sample_scores, bin_id)

  # Perform a left join: keep all AUTOSOMAL bins (in order), add scores
  final_matrix_dt <- merge(final_matrix_dt, sample_scores, by = "bin_id", all.x = TRUE)
  setnames(final_matrix_dt, "score", sample_id)
}

message("\nAutosomal matrix construction complete.")
message("Final matrix dimensions: ", nrow(final_matrix_dt), " autosomal bins x ", (ncol(final_matrix_dt)), " columns (including bin_id).")
message("Number of samples in matrix header: ", length(processed_sample_ids))

# Final check
if (!all(names(bed_files) %in% colnames(final_matrix_dt))) {
    missing_samples <- setdiff(names(bed_files), colnames(final_matrix_dt))
    warning("The following expected samples are missing as columns in the final matrix, likely due to processing errors: ", paste(missing_samples, collapse=", "))
}


# --- 8. Save the Autosomal Methylation Matrix File ---
message("\nWriting filtered (autosomal) methylation score matrix file...")
fwrite(final_matrix_dt, file = output_matrix_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "NA")
message("Autosomal methylation score matrix saved to: ", output_matrix_file)

message("\nPreprocessing finished successfully (autosomal data only).")
message("Ensure sample IDs in '", output_matrix_file, "' header match sample IDs in your genotype file header!")
message("Ensure your genotype data used in the next step also excludes non-autosomal markers for consistency.")
