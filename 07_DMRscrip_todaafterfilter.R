# Load necessary libraries
#if (!require("data.table")) install.packages("data.table")
#if (!require("dplyr")) install.packages("dplyr")
#if (!require("ggplot2")) install.packages("ggplot2")
#if (!require("matrixStats")) install.packages("matrixStats")
#if (!require("tidyr")) install.packages("tidyr")

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(tidyr)

#------------------------------
# Define file path and load all files
#------------------------------
# Path to BED files
file_path <- "/emine190-workingdir/output/map1000bpbins/"  
files <- list.files(path = file_path, pattern = "*.bed", full.names = TRUE)
setwd("/emine190-workingdir/output/DMRfromR/")
# Read all files and combine into one dataframe
data_list <- lapply(files, function(f) {
  dt <- fread(f, col.names = c("chrom", "start", "end", "score"))
  sample_id <- gsub(".bed", "", basename(f))
  dt[, sample := sample_id]
  return(dt)
})

# Merge all samples based on genomic coordinates
merged_data <- Reduce(function(x, y) merge(x, y, by = c("chrom", "start", "end"), all = TRUE), data_list)

#If data is skewed test to log transform it
#merged_data$mean_methylation <- log2(merged_data$mean_methylation + 1)

# Convert sample columns into a matrix of methylation scores
methylation_matrix <- as.matrix(merged_data[, 4:ncol(merged_data), with = FALSE])

#------------------------------
# Calculate Mean and Variability
#------------------------------
# Compute mean and standard deviation across all samples
merged_data$mean_methylation <- rowMeans(methylation_matrix, na.rm = TRUE)
merged_data$sd_methylation <- rowSds(methylation_matrix, na.rm = TRUE)

#------------------------------
# Identify DMRs Based on Outliers
#------------------------------
# Define thresholds for outliers:
# - Top/bottom 5% of mean methylation
# - Top 5% of variability

mean_cutoff_high <- quantile(merged_data$mean_methylation, 0.95, na.rm = TRUE)
mean_cutoff_low <- quantile(merged_data$mean_methylation, 0.05, na.rm = TRUE)
sd_cutoff <- quantile(merged_data$sd_methylation, 0.95, na.rm = TRUE)

# Filter out DMRs based on outliers
dmr <- merged_data %>%
  filter(mean_methylation > mean_cutoff_high | 
         mean_methylation < mean_cutoff_low | 
         sd_methylation > sd_cutoff)

#------------------------------
# Save Results
#------------------------------
# Save DMRs to BED file
fwrite(dmr[, c("chrom", "start", "end", "mean_methylation", "sd_methylation")], 
       "dmr.bed", sep = "\t", col.names = FALSE)

#------------------------------
# Visualization (Volcano-like plot)
#------------------------------
ggplot(merged_data, aes(x = mean_methylation, y = sd_methylation)) +
  geom_point(alpha = 0.5) +
  geom_point(data = dmr, aes(x = mean_methylation, y = sd_methylation), color = "red", size = 2) +
  geom_vline(xintercept = c(mean_cutoff_low, mean_cutoff_high), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = sd_cutoff, linetype = "dashed", color = "green") +
  labs(title = "Methylation Outliers",
       x = "Mean Methylation",
       y = "Methylation Variability (SD)") +
  theme_minimal()
ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6, units = "in")
#------------------------------
# DONE!
#------------------------------
cat("\nNumber of DMRs identified: ", nrow(dmr), "\n")
