# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))
# Need to run bulk RNA-seq analysis first
# save.image(paste0(onedrive_dir, "Influenza Analysis/bulk_RNA_analysis.RData"))
load(paste0(onedrive_dir, "Influenza Analysis/bulk_RNA_analysis.RData"))

# Data frame that captures which sc pseudobulk genes pass different logFC thresholds for different bulk differential expression analyses
# The basic idea is, for each day (D2, D5, D8, D28), we look at whether the gene passes the threshold for:
# log2FC 0.1, log2FC 0.3, log2FC 0.585, log2FC 1, logFC 2
# log2FC -0.1, log2FC -0.3, log2FC -0.585, log2FC -1, logFC -2
high_pos_pseudobulk_sc_genes_bulk_passing_df <- data.frame(gene = character(), cell_types = character(), D2_0.2 = logical(), D2_0.3 = logical(), D2_0.585 = logical(), D2_1 = logical(), D2_2 = logical(),
                                                  D5_0.2 = logical(), D5_0.3 = logical(), D5_0.585 = logical(), D5_1 = logical(), D5_2 = logical(),
                                                  D8_0.2 = logical(), D8_0.3 = logical(), D8_0.585 = logical(), D8_1 = logical(), D8_2 = logical(),
                                                  D28_0.2 = logical(), D28_0.3 = logical(), D28_0.585 = logical(), D28_1 = logical(), D28_2 = logical(),
                                                  D2_negative_0.2 = logical(), D2_negative_0.3 = logical(), D2_negative_0.585 = logical(), D2_negative_1 = logical(), D2_negative_2 = logical(),
                                                  D5_negative_0.2 = logical(), D5_negative_0.3 = logical(), D5_negative_0.585 = logical(), D5_negative_1 = logical(), D5_negative_2 = logical(),
                                                  D8_negative_0.2 = logical(), D8_negative_0.3 = logical(), D8_negative_0.585 = logical(), D8_negative_1 = logical(), D8_negative_2 = logical(),
                                                  D28_negative_0.2 = logical(), D28_negative_0.3 = logical(), D28_negative_0.585 = logical(), D28_negative_1 = logical(), D28_negative_2 = logical())

high_neg_pseudobulk_sc_genes_bulk_passing_df <- data.frame(gene = character(), cell_types = character(), D2_0.2 = logical(), D2_0.3 = logical(), D2_0.585 = logical(), D2_1 = logical(), D2_2 = logical(),
                                                      D5_0.2 = logical(), D5_0.3 = logical(), D5_0.585 = logical(), D5_1 = logical(), D5_2 = logical(),
                                                      D8_0.2 = logical(), D8_0.3 = logical(), D8_0.585 = logical(), D8_1 = logical(), D8_2 = logical(),
                                                      D28_0.2 = logical(), D28_0.3 = logical(), D28_0.585 = logical(), D28_1 = logical(), D28_2 = logical(),
                                                      D2_negative_0.2 = logical(), D2_negative_0.3 = logical(), D2_negative_0.585 = logical(), D2_negative_1 = logical(), D2_negative_2 = logical(),
                                                      D5_negative_0.2 = logical(), D5_negative_0.3 = logical(), D5_negative_0.585 = logical(), D5_negative_1 = logical(), D5_negative_2 = logical(),
                                                      D8_negative_0.2 = logical(), D8_negative_0.3 = logical(), D8_negative_0.585 = logical(), D8_negative_1 = logical(), D8_negative_2 = logical(),
                                                      D28_negative_0.2 = logical(), D28_negative_0.3 = logical(), D28_negative_0.585 = logical(), D28_negative_1 = logical(), D28_negative_2 = logical())

low_pos_pseudobulk_sc_genes_bulk_passing_df <- data.frame(gene = character(), cell_types = character(), D2_0.2 = logical(), D2_0.3 = logical(), D2_0.585 = logical(), D2_1 = logical(), D2_2 = logical(),
                                                           D5_0.2 = logical(), D5_0.3 = logical(), D5_0.585 = logical(), D5_1 = logical(), D5_2 = logical(),
                                                           D8_0.2 = logical(), D8_0.3 = logical(), D8_0.585 = logical(), D8_1 = logical(), D8_2 = logical(),
                                                           D28_0.2 = logical(), D28_0.3 = logical(), D28_0.585 = logical(), D28_1 = logical(), D28_2 = logical(),
                                                           D2_negative_0.2 = logical(), D2_negative_0.3 = logical(), D2_negative_0.585 = logical(), D2_negative_1 = logical(), D2_negative_2 = logical(),
                                                           D5_negative_0.2 = logical(), D5_negative_0.3 = logical(), D5_negative_0.585 = logical(), D5_negative_1 = logical(), D5_negative_2 = logical(),
                                                           D8_negative_0.2 = logical(), D8_negative_0.3 = logical(), D8_negative_0.585 = logical(), D8_negative_1 = logical(), D8_negative_2 = logical(),
                                                           D28_negative_0.2 = logical(), D28_negative_0.3 = logical(), D28_negative_0.585 = logical(), D28_negative_1 = logical(), D28_negative_2 = logical())

low_neg_pseudobulk_sc_genes_bulk_passing_df <- data.frame(gene = character(), cell_types = character(), D2_0.2 = logical(), D2_0.3 = logical(), D2_0.585 = logical(), D2_1 = logical(), D2_2 = logical(),
                                                           D5_0.2 = logical(), D5_0.3 = logical(), D5_0.585 = logical(), D5_1 = logical(), D5_2 = logical(),
                                                           D8_0.2 = logical(), D8_0.3 = logical(), D8_0.585 = logical(), D8_1 = logical(), D8_2 = logical(),
                                                           D28_0.2 = logical(), D28_0.3 = logical(), D28_0.585 = logical(), D28_1 = logical(), D28_2 = logical(),
                                                           D2_negative_0.2 = logical(), D2_negative_0.3 = logical(), D2_negative_0.585 = logical(), D2_negative_1 = logical(), D2_negative_2 = logical(),
                                                           D5_negative_0.2 = logical(), D5_negative_0.3 = logical(), D5_negative_0.585 = logical(), D5_negative_1 = logical(), D5_negative_2 = logical(),
                                                           D8_negative_0.2 = logical(), D8_negative_0.3 = logical(), D8_negative_0.585 = logical(), D8_negative_1 = logical(), D8_negative_2 = logical(),
                                                           D28_negative_0.2 = logical(), D28_negative_0.3 = logical(), D28_negative_0.585 = logical(), D28_negative_1 = logical(), D28_negative_2 = logical())

high_pos_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_high_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                  raw_high_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                  raw_high_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                  raw_high_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                  high_pos_pseudobulk_sc_genes_bulk_passing_df, sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$sc_log2FC > 0,])

high_neg_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_high_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                      raw_high_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                      raw_high_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                      raw_high_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                      high_neg_pseudobulk_sc_genes_bulk_passing_df, sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$sc_log2FC < 0,])

# Validation threshold: genes that pass 0.2 FC (or -0.2 FC) for D28 bulk
# 28 genes
high_passing_pos_gene_df <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]
high_passing_pos_gene_df <- fill_in_special_notes_pos(high_passing_pos_gene_df)
write.table(high_passing_pos_gene_df, file = paste0(onedrive_dir, "Influenza Analysis/high_passing_pos_df.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
high_passing_pos_genes <- high_passing_pos_gene_df$gene
# 86 genes
high_passing_neg_gene_df <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]
high_passing_neg_gene_df <- fill_in_special_notes_neg(high_passing_neg_gene_df)
write.table(high_passing_neg_gene_df, file = paste0(onedrive_dir, "Influenza Analysis/high_passing_neg_df.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
high_passing_neg_genes <- high_passing_neg_gene_df$gene

# Check these validated genes on low viral load individuals
low_pos_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_low_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                          low_pos_pseudobulk_sc_genes_bulk_passing_df, sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% high_passing_pos_genes,])

low_neg_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_low_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                          low_neg_pseudobulk_sc_genes_bulk_passing_df, sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% high_passing_neg_genes,])

# Genes that pass 0.2 FC (or -0.2 FC) for D28 bulk (low)
# 6 genes
low_pos_pseudobulk_sc_genes_bulk_passing_df <- fill_in_special_notes_pos(low_pos_pseudobulk_sc_genes_bulk_passing_df, viral_load = "LVL")
write.table(low_pos_pseudobulk_sc_genes_bulk_passing_df, file = paste0(onedrive_dir, "Influenza Analysis/low_passing_pos_df.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
low_passing_pos_genes <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
# 16 genes
low_neg_pseudobulk_sc_genes_bulk_passing_df <- fill_in_special_notes_neg(low_neg_pseudobulk_sc_genes_bulk_passing_df, viral_load = "LVL")
write.table(low_neg_pseudobulk_sc_genes_bulk_passing_df, file = paste0(onedrive_dir, "Influenza Analysis/low_passing_neg_df.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
low_passing_neg_genes <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene



# FLAGGED GENES - those genes which have negative FC in SC data but positive FC in bulk data (or vice versa)
flagged_high_pos_genes <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
flagged_high_pos_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_high_pos_genes,]

flagged_high_neg_genes <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
flagged_high_neg_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_high_neg_genes,]

flagged_low_pos_genes <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
flagged_low_pos_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_low_pos_genes,]

flagged_low_neg_genes <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
flagged_low_neg_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_low_neg_genes,]




