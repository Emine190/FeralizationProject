# Install and load necessary packages if you haven't already
library(MatrixEQTL)
library(qqman)
library(ggplot2)

# --- 1. Define Input File Names and Parameters ---
methylation_file = "/emine190-workingdir/output/QTLRanalysis/methylation_dataNOsexchromo.txt"
genotype_file = "/emine190-workingdir/original-files/cc5k.autosomes.selectFERAL.renamesamples.biallelic.phase.filtermafmiss.genotypematrixNUMformat_emine190version.txt"
covariate_file = NULL
output_file_name = "trans_mqtl_results.txt"

# Define cis-window (set to NULL for trans-QTL only)
cis_window = NULL

# Set working directory
setwd("/emine190-workingdir/output/QTLRanalysis/TransQTL")

# --- 2. Load Data ---
me_df = read.table(methylation_file, header = TRUE, row.names = 1)
ge_df = read.table(genotype_file, header = TRUE, row.names = 1)
cv_df = if (!is.null(covariate_file)) read.table(covariate_file, header = TRUE, row.names = 1) else NULL

# --- 3. Ensure Sample Order Consistency ---
common_samples = intersect(colnames(me_df), colnames(ge_df))
if (!is.null(cv_df)) {
  common_samples = intersect(common_samples, colnames(cv_df))
}

me_df = me_df[, common_samples]
ge_df = ge_df[, common_samples]
if (!is.null(cv_df)) {
  cv_df = cv_df[, common_samples]
}

# --- 4. Prepare Data for MatrixEQTL with Slicing ---
# Convert to matrix
methylation_matrix = as.matrix(me_df)
genotype_matrix = as.matrix(ge_df)
covariate_matrix = if (!is.null(cv_df)) as.matrix(cv_df) else matrix()

# Create SlicedData objects
methylation_data = SlicedData$new()
methylation_data$CreateFromMatrix(methylation_matrix)
methylation_data$ResliceCombined(sliceSize = 3000)  # Adjust sliceSize as needed

genotype_data = SlicedData$new()
genotype_data$CreateFromMatrix(genotype_matrix)
genotype_data$ResliceCombined(sliceSize = 8000)  # Adjust sliceSize as needed

covariate_data = SlicedData$new()
covariate_data$CreateFromMatrix(covariate_matrix)

# --- 5. Run MatrixEQTL ---
model = modelLINEAR
#output_file = file.path(getwd(), output_file_name)
errorCovariance = numeric()

mqtl_results = Matrix_eQTL_main(
  snps = genotype_data,
  gene = methylation_data,
  output_file_name = output_file_name,
  pvOutputThreshold = 1e-5,
  useModel = model,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsave = FALSE
)

cat('Trans-QTL analysis done.\n')
cat(paste('Results saved to:', output_file, '\n'))

# --- 6. Optional Summary ---
if (!is.null(mqtl_results)) {
  cat('\n--- Summary of significant trans-QTLs ---\n')
  print(mqtl_results$trans$eqtls)
}
