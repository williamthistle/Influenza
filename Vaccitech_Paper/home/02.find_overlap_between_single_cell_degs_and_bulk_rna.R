# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))
# Need to run bulk RNA-seq analysis first
load(paste0(onedrive_dir, "Influenza Analysis/bulk_RNA_analysis.RData"))

# Data frame that captures which sc pseudobulk genes pass different logFC thresholds for different bulk differential expression analyses
# The basic idea is, for each day (D2, D5, D8, D28), we look at whether the gene passes the threshold for:
# log2FC 0.1, log2FC 0.3, log2FC 0.585, log2FC 1, logFC 2
# log2FC -0.1, log2FC -0.3, log2FC -0.585, log2FC -1, logFC -2
high_pos_pseudobulk_sc_genes_bulk_passing_df <- data.frame(gene = character(), D2_0.2 = logical(), D2_0.3 = logical(), D2_0.585 = logical(), D2_1 = logical(), D2_2 = logical(),
                                                  D5_0.2 = logical(), D5_0.3 = logical(), D5_0.585 = logical(), D5_1 = logical(), D5_2 = logical(),
                                                  D8_0.2 = logical(), D8_0.3 = logical(), D8_0.585 = logical(), D8_1 = logical(), D8_2 = logical(),
                                                  D28_0.2 = logical(), D28_0.3 = logical(), D28_0.585 = logical(), D28_1 = logical(), D28_2 = logical(),
                                                  D2_negative_0.2 = logical(), D2_negative_0.3 = logical(), D2_negative_0.585 = logical(), D2_negative_1 = logical(), D2_negative_2 = logical(),
                                                  D5_negative_0.2 = logical(), D5_negative_0.3 = logical(), D5_negative_0.585 = logical(), D5_negative_1 = logical(), D5_negative_2 = logical(),
                                                  D8_negative_0.2 = logical(), D8_negative_0.3 = logical(), D8_negative_0.585 = logical(), D8_negative_1 = logical(), D8_negative_2 = logical(),
                                                  D28_negative_0.2 = logical(), D28_negative_0.3 = logical(), D28_negative_0.585 = logical(), D28_negative_1 = logical(), D28_negative_2 = logical())

high_neg_pseudobulk_sc_genes_bulk_passing_df <- data.frame(gene = character(), D2_0.2 = logical(), D2_0.3 = logical(), D2_0.585 = logical(), D2_1 = logical(), D2_2 = logical(),
                                                      D5_0.2 = logical(), D5_0.3 = logical(), D5_0.585 = logical(), D5_1 = logical(), D5_2 = logical(),
                                                      D8_0.2 = logical(), D8_0.3 = logical(), D8_0.585 = logical(), D8_1 = logical(), D8_2 = logical(),
                                                      D28_0.2 = logical(), D28_0.3 = logical(), D28_0.585 = logical(), D28_1 = logical(), D28_2 = logical(),
                                                      D2_negative_0.2 = logical(), D2_negative_0.3 = logical(), D2_negative_0.585 = logical(), D2_negative_1 = logical(), D2_negative_2 = logical(),
                                                      D5_negative_0.2 = logical(), D5_negative_0.3 = logical(), D5_negative_0.585 = logical(), D5_negative_1 = logical(), D5_negative_2 = logical(),
                                                      D8_negative_0.2 = logical(), D8_negative_0.3 = logical(), D8_negative_0.585 = logical(), D8_negative_1 = logical(), D8_negative_2 = logical(),
                                                      D28_negative_0.2 = logical(), D28_negative_0.3 = logical(), D28_negative_0.585 = logical(), D28_negative_1 = logical(), D28_negative_2 = logical())

low_pos_pseudobulk_sc_genes_bulk_passing_df <- data.frame(gene = character(), D2_0.2 = logical(), D2_0.3 = logical(), D2_0.585 = logical(), D2_1 = logical(), D2_2 = logical(),
                                                           D5_0.2 = logical(), D5_0.3 = logical(), D5_0.585 = logical(), D5_1 = logical(), D5_2 = logical(),
                                                           D8_0.2 = logical(), D8_0.3 = logical(), D8_0.585 = logical(), D8_1 = logical(), D8_2 = logical(),
                                                           D28_0.2 = logical(), D28_0.3 = logical(), D28_0.585 = logical(), D28_1 = logical(), D28_2 = logical(),
                                                           D2_negative_0.2 = logical(), D2_negative_0.3 = logical(), D2_negative_0.585 = logical(), D2_negative_1 = logical(), D2_negative_2 = logical(),
                                                           D5_negative_0.2 = logical(), D5_negative_0.3 = logical(), D5_negative_0.585 = logical(), D5_negative_1 = logical(), D5_negative_2 = logical(),
                                                           D8_negative_0.2 = logical(), D8_negative_0.3 = logical(), D8_negative_0.585 = logical(), D8_negative_1 = logical(), D8_negative_2 = logical(),
                                                           D28_negative_0.2 = logical(), D28_negative_0.3 = logical(), D28_negative_0.585 = logical(), D28_negative_1 = logical(), D28_negative_2 = logical())

low_neg_pseudobulk_sc_genes_bulk_passing_df <- data.frame(gene = character(), D2_0.2 = logical(), D2_0.3 = logical(), D2_0.585 = logical(), D2_1 = logical(), D2_2 = logical(),
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
                                                                                  high_pos_pseudobulk_sc_genes_bulk_passing_df, pos_sc_pseudobulk_genes)

high_neg_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_high_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                      raw_high_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                      raw_high_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                      raw_high_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                      high_neg_pseudobulk_sc_genes_bulk_passing_df, neg_sc_pseudobulk_genes)

# Validation threshold: genes that pass 0.2 FC (or -0.2 FC) for D28 bulk
# 28 genes
high_passing_pos_genes <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
# 86 genes
high_passing_neg_genes <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene

# Check these validated genes on low viral load individuals
low_pos_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_low_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                          low_pos_pseudobulk_sc_genes_bulk_passing_df, high_passing_pos_genes)

low_neg_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_low_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                          raw_low_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                          low_neg_pseudobulk_sc_genes_bulk_passing_df, high_passing_neg_genes)

# Validation threshold: genes that pass 0.2 FC (or -0.2 FC) for D28 bulk
# 6 genes
low_passing_pos_genes <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
# 16 genes
low_passing_neg_genes <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene

# This is a pretty low threshold, so let's go up to 0.3 FC / 0.585 FC
# 22 genes
high_passing_pos_genes_0.3 <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.3 == TRUE,]$gene
# 65 genes
high_passing_neg_genes_0.3 <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.3 == TRUE,]$gene
# 5 genes
low_passing_pos_genes_0.3 <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.3 == TRUE,]$gene
# 16 genes
low_passing_neg_genes_0.3 <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.3 == TRUE,]$gene

# 3 genes
high_passing_pos_genes_0.585 <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.585 == TRUE,]$gene
# 2 genes
high_passing_neg_genes_0.585 <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.585 == TRUE,]$gene
# 2 genes
low_passing_pos_genes_0.585 <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.585 == TRUE,]$gene
# 1 gene
low_passing_neg_genes_0.585 <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.585 == TRUE,]$gene

# Subset dataframes to passing genes
high_passing_pos_gene_df <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$gene %in% high_passing_pos_genes,]
high_passing_neg_gene_df <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$gene %in% high_passing_neg_genes,]
low_passing_pos_gene_df <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$gene %in% low_passing_pos_genes,]
low_passing_neg_gene_df <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$gene %in% low_passing_neg_genes,]

### HIGH DOWNSTREAM ANALYSIS ###
# Positive Genes that pass D5/D8/D28 (no genes pass all 4 time points)
d5_d8_d28_high_passing_pos_gene_df <- high_passing_pos_gene_df[high_passing_pos_gene_df$D5_0.2 == TRUE & high_passing_pos_gene_df$D8_0.2 == TRUE & high_passing_pos_gene_df$D28_0.2 == TRUE,]
# Positive genes that pass D8/D28 (still just IL10RA)
d8_d28_high_passing_pos_gene_df <- high_passing_pos_gene_df[high_passing_pos_gene_df$D8_0.2 == TRUE & high_passing_pos_gene_df$D28_0.2 == TRUE,]
# Positive genes that pass D5/D28 (prevalent during acute phase and then upregulated again a month later?)
d5_d28_high_passing_pos_gene_df <- high_passing_pos_gene_df[high_passing_pos_gene_df$D5_0.2 == TRUE & high_passing_pos_gene_df$D28_0.2 == TRUE,]

# Negative Genes that pass D2/D5/D8/D28
# 23 genes
d2_d5_d8_d28_high_passing_neg_gene_df <- high_passing_neg_gene_df[high_passing_neg_gene_df$D2_negative_0.2 == TRUE & high_passing_neg_gene_df$D5_negative_0.2 == TRUE & high_passing_neg_gene_df$D8_negative_0.2 == TRUE & high_passing_neg_gene_df$D28_negative_0.2 == TRUE,]
# Negative Genes that pass D5/D8/D28
# 30 genes
d5_d8_d28_high_passing_neg_gene_df <- high_passing_neg_gene_df[high_passing_neg_gene_df$D5_negative_0.2 == TRUE & high_passing_neg_gene_df$D8_negative_0.2 == TRUE & high_passing_neg_gene_df$D28_negative_0.2 == TRUE,]
# Positive genes that pass D8/D28
# 35 genes
d8_d28_high_passing_neg_gene_df <- high_passing_neg_gene_df[high_passing_neg_gene_df$D8_negative_0.2 == TRUE & high_passing_neg_gene_df$D28_negative_0.2 == TRUE,]
# Positive genes that pass D5/D28 (prevalent during acute phase and then downregulated again a month later?)
# 37 genes
d5_d28_high_passing_neg_gene_df <- high_passing_neg_gene_df[high_passing_neg_gene_df$D5_negative_0.2 == TRUE & high_passing_neg_gene_df$D28_negative_0.2 == TRUE,]



### LOW DOWNSTREAM ANALYSIS ###
# No positive genes found for D2/D5/D8, so looking at multiple time points doesn't make sense







# FLAGGED GENES - those genes which have negative FC in SC data but positive FC in bulk data (or vice versa)
flagged_high_pos_genes <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
flagged_high_pos_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_high_pos_genes,]

flagged_high_neg_genes <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
flagged_high_neg_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_high_neg_genes,]

flagged_low_pos_genes <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
flagged_low_pos_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_low_pos_genes,]

flagged_low_neg_genes <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
flagged_low_neg_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_low_neg_genes,]




