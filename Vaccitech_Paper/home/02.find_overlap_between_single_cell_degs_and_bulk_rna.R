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

low_pos_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_low_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                           raw_low_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                           raw_low_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                           raw_low_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                          low_pos_pseudobulk_sc_genes_bulk_passing_df, pos_sc_pseudobulk_genes)

low_neg_pseudobulk_sc_genes_bulk_passing_df <- find_degs_across_time_points_for_gene_list(raw_low_placebo_period_2_D2_vs_D_minus_1_results,
                                                                                           raw_low_placebo_period_2_D5_vs_D_minus_1_results,
                                                                                           raw_low_placebo_period_2_D8_vs_D_minus_1_results,
                                                                                           raw_low_placebo_period_2_D28_vs_D_minus_1_results,
                                                                                          low_neg_pseudobulk_sc_genes_bulk_passing_df, neg_sc_pseudobulk_genes)

# Validation threshold: genes that pass 0.2 FC (or -0.2 FC) for D28 bulk
# 28 genes
high_passing_pos_genes <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
# 87 genes
high_passing_neg_genes <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
# 14 genes
low_passing_pos_genes <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
# 27 genes
low_passing_neg_genes <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene

# This is a pretty low threshold, so let's go up to 0.3 FC / 0.585 FC
# Honestly, we may just use 0.1 FC for IL10RA since it's an interesting gene biologically and has consistent signal across days
# 22 genes
high_passing_pos_genes_0.3 <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.3 == TRUE,]$gene
# 65 genes
high_passing_neg_genes_0.3 <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.3 == TRUE,]$gene
# 13 genes
low_passing_pos_genes_0.3 <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.3 == TRUE,]$gene
# 27 genes
low_passing_neg_genes_0.3 <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.3 == TRUE,]$gene

# 3 genes
high_passing_pos_genes_0.585 <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.585 == TRUE,]$gene
# 2 genes
high_passing_neg_genes_0.585 <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.585 == TRUE,]$gene
# 4 genes
low_passing_pos_genes_0.585 <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_0.585 == TRUE,]$gene
# 2 genes
low_passing_neg_genes_0.585 <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.585 == TRUE,]$gene

# Subset dataframes to passing genes
high_passing_pos_gene_df <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$gene %in% high_passing_pos_genes,]
high_passing_neg_gene_df <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$gene %in% high_passing_neg_genes,]$gene
low_passing_pos_gene_df <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$gene %in% low_passing_pos_genes,]$gene
low_passing_neg_gene_df <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$gene %in% low_passing_neg_genes,]$gene







# FLAGGED GENES - those genes which have negative FC in SC data but positive FC in bulk data (or vice versa)
flagged_high_pos_genes <- high_pos_pseudobulk_sc_genes_bulk_passing_df[high_pos_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
flagged_high_pos_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_high_pos_genes,]

flagged_high_neg_genes <- high_neg_pseudobulk_sc_genes_bulk_passing_df[high_neg_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
flagged_high_neg_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_high_neg_genes,]

flagged_low_pos_genes <- low_pos_pseudobulk_sc_genes_bulk_passing_df[low_pos_pseudobulk_sc_genes_bulk_passing_df$D28_negative_0.2 == TRUE,]$gene
flagged_low_pos_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_low_pos_genes,]

# ANXA1 is both negative and positive in SC data - maybe interesting?
flagged_low_neg_genes <- low_neg_pseudobulk_sc_genes_bulk_passing_df[low_neg_pseudobulk_sc_genes_bulk_passing_df$D28_0.2 == TRUE,]$gene
flagged_low_neg_gene_df <- sc_pseudobulk_gene_table[sc_pseudobulk_gene_table$Gene_Name %in% flagged_low_neg_genes,]




