# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# --- 0. Load Gene Annotations ---
gene_annotation_file <- "/emine190-workingdir/original-files/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.113.gtftobedgenes-job8657516.bed" # Replace with your file name
gene_annotations <- read.table(gene_annotation_file, header = FALSE, stringsAsFactors = FALSE)
colnames(gene_annotations) <- c("chromosome", "start", "end", "gene_id", "V5", "V6")
gene_annotations <- gene_annotations %>% select(gene_id, chromosome, start, end)

# --- 1. Define Input File Names and Parameters ---
cis_qtl_results_file <- "/emine190-workingdir/output/QTLRanalysis/CisQTL/significant_meqtl_results.txt" # Replace with your cis-mQTL results file
snp_positions_file <- "/emine190-workingdir/original-files/cc5k.autosomes.selectFERAL.renamesamples.biallelic.phase.filtermafmiss.snplocation.txt"
methylation_locations_file <- "/emine190-workingdir/output/QTLRanalysis/methylation_locationsNosexchromo.txt"
cis_window <- 1e6 # 1 Mb for filtering cis-QTLs

setwd("/emine190-workingdir/output/QTLRanalysis/CisQTL")

# --- 2. Load Data ---
cis_mqtls <- read.table(cis_qtl_results_file, header = TRUE, stringsAsFactors = FALSE)
snp_positions <- read.table(snp_positions_file, header = TRUE, stringsAsFactors = FALSE)
methylation_locations <- read.table(methylation_locations_file, header = TRUE, stringsAsFactors = FALSE)

# --- 3. Prepare Methylation Locations ---
methylation_locations_processed <- methylation_locations %>%
    mutate(gene_id = bin_id, gene_chr = chr, gene_start = start, gene_end = end, gene_midpoint = (start + end) / 2) %>%
    select(gene_id, gene_chr, gene_start, gene_end, gene_midpoint)

# --- 4. Merge Data ---
cis_mqtls_merged <- cis_mqtls %>%
    inner_join(snp_positions, by = c("SNP" = "snp_id")) %>%
    rename(snp_id = SNP, snp_chr = chr, snp_position = pos) %>%
    inner_join(methylation_locations_processed, by = c("gene" = "gene_id")) %>%
    filter(FDR < 0.05) # Assuming your file is already filtered for FDR < 0.05

# --- 5. Manhattan Plot (for cis-QTLs) ---
manhattan_plot_cis <- ggplot(cis_mqtls_merged, aes(x = snp_position / 1e6, y = -log10(pvalue), color = snp_chr)) +
    geom_point(size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(data = gene_annotations, aes(xintercept = start / 1e6), color = "grey70", linetype = "dotted", size = 0.2) +
    labs(title = "Cis-QTL Manhattan Plot with Gene Locations", x = "SNP Position (Mb)", y = "-log10(p-value)", color = "SNP Chromosome") +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 6), panel.grid.minor = element_blank())

print(manhattan_plot_cis)
ggsave("manhattan_plot_cis.pdf", plot = manhattan_plot_cis, width = 10, height = 6, units = "in", dpi = 300)

# --- 6. QQ Plot (for cis-QTL p-values) ---
qqplot_data_cis <- data.frame(observed = sort(-log10(cis_mqtls_merged$pvalue)), expected = sort(-log10(ppoints(nrow(cis_mqtls_merged)))))
qqplot_cis <- ggplot(qqplot_data_cis, aes(x = expected, y = observed)) +
    geom_point(size = 1.5) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = "QQ Plot of Cis-QTL p-values", x = "Expected -log10(p-value)", y = "Observed -log10(p-value)") +
    theme_bw()
print(qqplot_cis)
ggsave("qqplot_cis.pdf", plot = qqplot_cis, width = 6, height = 6, units = "in", dpi = 300)

# --- 8. LocusZoom-like Plots (for most significant gene per chromosome) ---
top_gene_per_chr <- cis_mqtls_merged %>%
    group_by(snp_chr) %>%
    slice_min(FDR, n = 1) %>%
    ungroup()

cat("\nGenerating LocusZoom plots for the most significant gene per chromosome...\n")

for (chr in unique(top_gene_per_chr$snp_chr)) {
    top_qtl_chr <- top_gene_per_chr %>% filter(snp_chr == chr)

    if (nrow(top_qtl_chr) > 0) {
        representative_gene <- top_qtl_chr$gene[1] # Take the first if multiple with same min FDR

        gene_region_mqtl <- methylation_locations_processed %>% filter(gene_id == representative_gene)

        if (nrow(gene_region_mqtl) > 0) {
            plot_region_start <- gene_region_mqtl$gene_midpoint - cis_window
            plot_region_end <- gene_region_mqtl$gene_midpoint + cis_window

            locus_data <- cis_mqtls_merged %>%
                filter(gene == representative_gene, snp_chr == gene_region_mqtl$gene_chr,
                       snp_position >= plot_region_start, snp_position <= plot_region_end) %>%
                mutate(distance_to_gene = snp_position - gene_region_mqtl$gene_midpoint)

            relevant_genes <- gene_annotations %>%
                filter(chromosome == gene_region_mqtl$gene_chr, start <= plot_region_end, end >= plot_region_start)

            if (nrow(locus_data) > 0) {
                locuszoom_plot_cis <- ggplot(locus_data, aes(x = distance_to_gene / 1e3, y = -log10(pvalue))) +
                    geom_point(aes(color = -log10(FDR)), size = 2) +
                    scale_color_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
                    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
                    geom_segment(data = relevant_genes, aes(x = (start - gene_region_mqtl$gene_midpoint) / 1e3, xend = (end - gene_region_mqtl$gene_midpoint) / 1e3, y = -2, yend = -2), linewidth = 1, color = "black") +
                    geom_text_repel(data = relevant_genes, aes(x = ((start + end) / 2 - gene_region_mqtl$gene_midpoint) / 1e3, y = -3, label = gene_id), size = 2, segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0.1) +
                    labs(title = paste("Cis-QTL Locus Zoom (Chr:", chr, ", Gene:", representative_gene, ")"), x = "Distance from Gene Midpoint (kb)", y = "-log10(p-value)") +
                    coord_cartesian(ylim = c(0, max(-log10(locus_data$pvalue)) + 1), clip = "off") +
                    theme_bw() +
                    theme(plot.margin = unit(c(1, 1, 2, 1), "lines"))
                print(locuszoom_plot_cis)
                filename <- paste0("locuszoom_cis_chr", chr, "_top_gene_", representative_gene, ".pdf")
                ggsave(filename, plot = locuszoom_plot_cis, width = 8, height = 6, units = "in", dpi = 300)
            }
        }
    }
}

# --- 9. Summary Table of Significant cis-QTLs (FDR < 0.05 already filtered) ---
significant_cis_qtls_table <- cis_mqtls_merged %>%
    arrange(FDR)

print("Summary Table of Significant Cis-QTLs (FDR < 0.05):")
#print(significant_cis_qtls_table)

# Save the significant QTLs table to a txt file
#write.table(significant_cis_qtls_table, "significant_cis_qtls_table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Optional: Using kableExtra for better table formatting
if (requireNamespace("kableExtra", quietly = TRUE)) {
    library(kableExtra)
    formatted_table_cis <- significant_cis_qtls_table %>%
    mutate(distance = snp_position - gene_midpoint) %>%
    dplyr::select(gene, snp_id, snp_chr, snp_position, gene_chr, gene_start, gene_end,
                  distance, beta, pvalue, FDR) %>%
    kable() %>%
    kable_styling(bootstrap_options = c("striped", "hover"))
    print(formatted_table_cis)

    # Save the kableExtra formatted table as an HTML file
    save_kable(formatted_table_cis, file = "formatted_significant_cis_qtls_table.html")

    # Save the kableExtra formatted table as a TXT file of the kable output
    capture.output(formatted_table_cis, file = "formatted_significant_cis_qtls_table_kable.txt")

} else {
    cat("\nInstall 'kableExtra' to save the formatted table as HTML or TXT.\n")
}
