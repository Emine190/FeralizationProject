# --- 1. Load Libraries ---

#install.packages("MatrixEQTL", repos = "http://cloud.r-project.org/")
#install.packages("qqman", repos = "http://cloud.r-project.org/")

library(MatrixEQTL)

# For plotting

library(qqman) # Loaded later if needed for plotting


# --- 2. Define File Paths and Parameters ---
# *** USER: SET THESE VARIABLES ***

# -- Input Files --
# Output from Script preprocesseQTL (Preprocessing)
methylation_data_file <- "/emine190-workingdir/output/QTLRanalysis/methylation_dataNOsexchromo.txt"
methylation_loc_file <- "/emine190-workingdir/output/QTLRanalysis/methylation_locationsNosexchromo.txt"

# Your Genotype Files
genotype_data_file <- "/emine190-workingdir/original-files/cc5k.autosomes.selectFERAL.renamesamples.biallelic.phase.filtermafmiss.genotypematrixNUMformat_emine190version.txt"
marker_loc_file <- "/emine190-workingdir/original-files/cc5k.autosomes.selectFERAL.renamesamples.biallelic.phase.filtermafmiss.snplocation.txt"

setwd("/emine190-workingdir/output/QTLRanalysis/CisQTL/")

# Optional Covariates File (samples x covariates, tab-separated, with header)
# Set to NULL if no covariates are used.
# Example: covariates_file <- "covariates.txt"
covariates_file <- NULL

# -- Output Files --
# Base name for MatrixEQTL results files
output_results_basename <- "meqtl_resultsV2" # Will create e.g., meqtl_results_cis.txt

# -- Analysis Parameters --
# Model type: modelLINEAR, modelANOVA, or modelLINEAR_CROSS
# modelLINEAR is standard for additive effects on quantitative traits.
useModel <- modelLINEAR

# P-value threshold for reporting results in the output files.
# NOTE: This is NOT the significance threshold after multiple testing correction.
# Set pvOutputThreshold_trans = 0 to perform ONLY cis-analysis (faster).
pvOutputThreshold_cis <- 1e-3  # Report cis-QTLs with p < 0.001
pvOutputThreshold_trans = 0    # DOne in a different script

# Distance threshold for defining cis-QTLs (in base pairs).
# Markers within this distance upstream/downstream of a bin's start/end are considered cis.
cisDist <- 1e6 # 1 Megabase (1,000,000 bp)

# *** END USER SETTINGS ***

# Construct full output filenames
output_file_cis <- paste0(output_results_basename, "_cis.txt")
output_file_trans <- paste0(output_results_basename, "_trans.txt") # Will only be created if pvOutputThreshold_trans > 0
manhattan_plot_file <- paste0(output_results_basename, "_cis_manhattan.png")
qq_plot_file <- paste0(output_results_basename, "_qqplot.pdf")


# --- 3. Load Methylation Data (Phenotypes) ---
message("\nLoading methylation data (phenotypes)...")
meth = SlicedData$new();
meth$fileDelimiter = "\t";      # Tab delimited file
meth$fileOmitCharacters = "NA"; # Define missing data character
meth$fileSkipRows = 1;          # Skip header row (sample IDs)
meth$fileSkipColumns = 1;       # Skip first column (bin IDs)
meth$fileSliceSize = 2000;      # Read data in chunks (adjust if memory issues occur)
tryCatch({
    meth$LoadFile(methylation_data_file);
}, error = function(e) {
    stop("FATAL ERROR loading methylation data file '", methylation_data_file, "': ", e$message)
})
message("Methylation data loaded: ", meth$nCols(), " samples, ", meth$nRows(), " bins.")

# --- 4. Load Genotype Data ---
message("Loading genotype data...")
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # Tab delimited file
snps$fileOmitCharacters = "NA"; # Define missing data character
snps$fileSkipRows = 1;          # Skip header row (sample IDs)
snps$fileSkipColumns = 1;       # Skip first column (SNP IDs)
snps$fileSliceSize = 2000;      # Read data in chunks
tryCatch({
    snps$LoadFile(genotype_data_file);
}, error = function(e) {
    stop("FATAL ERROR loading genotype data file '", genotype_data_file, "': ", e$message)
})
message("Genotypes loaded: ", snps$nCols(), " samples, ", snps$nRows(), " markers.")

# --- 5. Load Covariates (Optional) ---
cvrt = SlicedData$new() # Create the object even if not used
if (!is.null(covariates_file) && file.exists(covariates_file)) {
    message("Loading covariates data...")
    cvrt$fileDelimiter = "\t";      # Tab delimited file
    cvrt$fileOmitCharacters = "NA"; # Define missing data character
    cvrt$fileSkipRows = 1;          # Skip header row (covariate names)
    cvrt$fileSkipColumns = 1;       # Skip first column (sample IDs)
    tryCatch({
        cvrt$LoadFile(covariates_file);
    }, error = function(e) {
        stop("FATAL ERROR loading covariates file '", covariates_file, "': ", e$message)
    })
    message("Covariates loaded: ", cvrt$nCols(), " covariates for ", cvrt$nRows(), " samples.")
    # Important check: Covariate samples must match phenotype/genotype samples
     if (!identical(cvrt$columnNames, meth$columnNames)) {
          # Try to reorder covariates to match methylation data
          message("Attempting to reorder covariates to match methylation sample order...")
          common_samples <- intersect(cvrt$columnNames, meth$columnNames)
          if (length(common_samples) < length(meth$columnNames)) {
              warning("Not all methylation samples found in covariates file!")
          }
          if (length(common_samples) == 0) {
              stop("FATAL ERROR: No common samples between covariates and methylation data.")
          }
          cvrt$ColumnSubsample(match(common_samples, cvrt$columnNames))
          meth$ColumnSubsample(match(common_samples, meth$columnNames))
          snps$ColumnSubsample(match(common_samples, snps$columnNames))
           if (!identical(cvrt$columnNames, meth$columnNames)) { # Check again after reordering attempt
                 stop("FATAL ERROR: Could not align sample order between covariates and methylation/genotype data after subsetting.")
           }
           message("Successfully subsetted and aligned data to ", length(common_samples), " common samples.")
      }

} else {
    if(!is.null(covariates_file)) message("Covariates file specified but not found: ", covariates_file, ". Proceeding without covariates.")
    # No covariates loaded, cvrt object remains empty
}


# --- 6. Check Sample Alignment ---
message("Checking sample alignment between genotype and methylation data...")
# MatrixEQTL requires samples to be in the *exact* same order.
if (!identical(snps$columnNames, meth$columnNames)) {
    message("Sample names or order mismatch between genotype and methylation files. Attempting to reconcile...")
    common_samples <- intersect(snps$columnNames, meth$columnNames)
    if (length(common_samples) == 0) {
        stop("FATAL ERROR: No common samples found between genotype and methylation data. Check headers!")
    }
    message("Found ", length(common_samples), " common samples. Subsetting data...")
    snps$ColumnSubsample(match(common_samples, snps$columnNames))
    meth$ColumnSubsample(match(common_samples, meth$columnNames))
    # Also subset covariates if they were loaded and matched originally
    if(cvrt$nCols() > 0 && identical(names(cvrt$.slices), common_samples)) {
       cvrt$ColumnSubsample(match(common_samples, cvrt$columnNames))
    } else if (cvrt$nCols() > 0) {
       # This case should ideally be handled during covariate loading alignment check
       warning("Covariate samples might be out of sync after reconciling genotype/methylation samples. Check carefully.")
    }

    if (!identical(snps$columnNames, meth$columnNames)) {
      stop("FATAL ERROR: Failed to reconcile sample order between genotype and methylation data.")
    }
     message("Successfully aligned data to ", length(common_samples), " common samples.")
} else {
    message("Sample names and order match perfectly between genotype and methylation data (", snps$nCols(), " samples).")
}
# Final check for covariates if used
if (cvrt$nCols() > 0 && !identical(cvrt$columnNames, meth$columnNames)) {
   stop("FATAL ERROR: Final check failed - Covariate samples do not match aligned genotype/methylation samples.")
}


# --- 7. Load Location Data ---
message("Loading location files...")
tryCatch({
    snpspos <- read.table(marker_loc_file, header = TRUE, stringsAsFactors = FALSE, sep="\t", check.names = FALSE)
    # Expect columns: snp_id, chr, pos
    if(!all(c("snp_id", "chr", "pos") %in% colnames(snpspos))) {
        stop("Marker location file ('", marker_loc_file, "') must contain columns named 'snp_id', 'chr', and 'pos'.")
    }
}, error = function(e) {
    stop("FATAL ERROR loading marker location file '", marker_loc_file, "': ", e$message)
})

tryCatch({
    methpos <- read.table(methylation_loc_file, header = TRUE, stringsAsFactors = FALSE, sep="\t", check.names = FALSE)
    # Expect columns: bin_id, chr, start, end
     if(!all(c("bin_id", "chr", "start", "end") %in% colnames(methpos))) {
        stop("Methylation location file ('", methylation_loc_file, "') must contain columns named 'bin_id', 'chr', 'start', and 'end'.")
    }
}, error = function(e) {
    stop("FATAL ERROR loading methylation location file '", methylation_loc_file, "': ", e$message)
})

message("Location files loaded.")
# Note: MatrixEQTL internally checks if the order of IDs in location files matches data files.
# If Script 1 (preprocessing) was run correctly, the methylation location order should match.
# You need to ensure your marker location file rows are in the same order as your genotype data file rows.


# --- 8. Run MatrixEQTL Analysis ---
message("\nStarting MatrixEQTL analysis...")
# Define output file for trans results (can be NULL if trans analysis is disabled)
output_file_trans_set <- if(pvOutputThreshold_trans > 0) output_file_trans else NULL

# Error covariance matrix: use numeric() for identity matrix assumption
errorCovariance <- numeric()

me <- Matrix_eQTL_main(
    snps = snps,
    gene = meth,                # Use 'gene' terminology for phenotypes
    cvrt = cvrt,                # Covariates SlicedData object (can be empty)
    output_file_name = output_file_trans_set, # Trans results output file
    pvOutputThreshold = pvOutputThreshold_trans, # Trans p-value threshold
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,             # Print progress messages
    output_file_name.cis = output_file_cis, # Cis results output file
    pvOutputThreshold.cis = pvOutputThreshold_cis, # Cis p-value threshold
    snpspos = snpspos,          # Data frame with SNP locations (snp_id, chr, pos)
    genepos = methpos,          # Data frame with bin locations (bin_id, chr, start, end) -> MatrixEQTL expects 'geneid' not 'bin_id' here, it uses first col.
    cisDist = cisDist,          # Cis distance threshold
    pvalue.hist = "qqplot",     # Generate a Q-Q plot of all p-values tested for cis
    min.pv.by.genesnp = FALSE,  # Set to TRUE to only record the top SNP per gene (bin) - usually FALSE
    noFDRsaveMemory = FALSE     # Calculate FDR (requires more memory, set TRUE only if needed)
)

# Rename the generated qq-plot file if it exists
default_qq_name <- if(pvOutputThreshold_cis > 0 && file.exists("qqplot.pdf")) "qqplot.pdf" else NULL # Default name used by MatrixEQTL if cis is run
if (!is.null(default_qq_name) && file.exists(default_qq_name)) {
    file.rename(default_qq_name, qq_plot_file)
    message("Q-Q plot saved to: ", qq_plot_file)
} else if (pvOutputThreshold_cis > 0) {
     message("Q-Q plot was requested but 'qqplot.pdf' not found. Check MatrixEQTL output.")
}


message("\nMatrixEQTL analysis complete.")

# --- 9. Results Summary ---
cat('Analysis done in: ', sprintf("%.2f", me$time.in.sec), ' seconds\n')
cat('Found ', me$cis$neqtls, ' cis-meQTLs passing the reporting threshold (p < ', pvOutputThreshold_cis, ')\n', sep='')
cat('cis-meQTL results saved to: ', output_file_cis, '\n')
if (pvOutputThreshold_trans > 0) {
    cat('Found ', me$trans$neqtls, ' trans-meQTLs passing the reporting threshold (p < ', pvOutputThreshold_trans, ')\n', sep='')
    cat('trans-meQTL results saved to: ', output_file_trans, '\n')
} else {
    cat('Trans-meQTL analysis was not performed (pvOutputThreshold_trans = 0).\n')
}

# Reminder about significance
message("\nREMINDER: The results files contain associations passing the reporting p-value threshold.")
message("For significance, filter these results based on the 'FDR' column (e.g., FDR < 0.05 or 0.10).")
