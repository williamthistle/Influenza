# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Create MetaIntegrator objects for bulk for all days (HVL and LVL)
high_D2_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D2", "2_D_minus_1")
high_D5_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D5", "2_D_minus_1")
high_D8_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D8", "2_D_minus_1")
high_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D28", "2_D_minus_1")

low_D2_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_placebo_counts, low_placebo_metadata, "2_D2", "2_D_minus_1")
low_D5_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_placebo_counts, low_placebo_metadata, "2_D5", "2_D_minus_1")
low_D8_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_placebo_counts, low_placebo_metadata, "2_D8", "2_D_minus_1")
low_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_placebo_counts, low_placebo_metadata, "2_D28", "2_D_minus_1")

# Because of the way we calculate AUC, in our case an AUC of under 0.3 is equally valuable as AUC of over 0.7
# We test genes one at a time as positive genes in our gene set signature
# If we get an AUC of under 0.3, that means that the same gene would score an AUC of over 0.7 as a negative gene
# in our gene set signature
# Our single cell data gives us single cell granularity on the genes being expressed. This also applies to the ATAC-seq data and our peaks! mintchip can help us verify our choices,
# but the ATAC-seq gives us the cell type!
# TODO: List cell types for genes when reporting
auc_df <- data.frame(Filtering_Assay = character(), Filtering_Method = character(), Discovery_Assay = character(), 
                     Discovery_Dataset = character(), Pos_Genes = integer(), Neg_Genes = integer(), Total_Passing_Genes = integer(), 
                     Total_Genes = integer(), Percentage_of_Passing_Genes = double(), stringsAsFactors = FALSE)
auc_names <- c("Filtering_Assay", "Filtering_Method", "Discovery_Assay", "Discovery_Dataset", "Passing_Pos_Genes", "Passing_Neg_Genes", "Total_Passing_Genes", "Total_Genes", "Percentage_of_Passing_Genes")

# AUCS ON D28 HVL BULK
# Note that our initial gene lists (from pseudobulk filtering / MAGICAL) are further curated using all D28 HVL bulk data

# SC (pseudobulk filtering)
high_sc_bulk_D28_sc_pseudobulk_gene_info <- find_aucs_of_interest(sc_pseudobulk_genes, high_D28_bulk_metaintegrator_obj, "sc_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_sc_bulk_D28_sc_pseudobulk_gene_info[[1]], "gene_auc")

# SC (MAGICAL filtering)
high_sc_bulk_D28_sc_magical_gene_info <- find_aucs_of_interest(sc_magical_genes, high_D28_bulk_metaintegrator_obj, "sc_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_sc_bulk_D28_sc_magical_gene_info[[1]], "gene_auc")

# Multiome (pseudobulk filtering
high_multiome_bulk_D28_multiome_pseudobulk_gene_info <- find_aucs_of_interest(multiome_pseudobulk_genes, high_D28_bulk_metaintegrator_obj, "multiome_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_multiome_bulk_D28_multiome_pseudobulk_gene_info[[1]], "gene_auc")

# Multiome (MAGICAL filtering)
high_multiome_bulk_D28_multiome_magical_gene_info <- find_aucs_of_interest(multiome_magical_genes, high_D28_bulk_metaintegrator_obj, "multiome_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_multiome_bulk_D28_multiome_magical_gene_info[[1]], "gene_auc")

# Most likely, we will next combine SC pseudobulk genes and multiome pseudobulk genes for a combined search across HVL D2/D5/D8 and LVL D28 (not much signal in LVL D2/D5/D8)
# MAGICAL genes may be interesting to look at further down the line
# Are we seeing interesting behavior out of the MAGICAL genes for D2 / D5 / D8 / D28 or for LVL D28?
# Which genes come up in both SC and multiome?
intersecting_genes_between_sc_and_multiome <- intersect(high_sc_bulk_D28_sc_pseudobulk_gene_info[[2]], high_multiome_bulk_D28_multiome_pseudobulk_gene_info[[2]])
combined_final_gene_list <- unique(c(high_sc_bulk_D28_sc_pseudobulk_gene_info[[2]], high_multiome_bulk_D28_multiome_pseudobulk_gene_info[[2]]))

# Find AUCs for combined gene list for HVL D2/D5/D8 and LVL D2/D5/D8/D28
high_bulk_D2_combined_info <- find_aucs_of_interest(combined_final_gene_list, high_D2_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D2 Bulk for HVL Subjects", high_bulk_D2_combined_info[[1]], "gene_auc")

high_bulk_D5_combined_info <- find_aucs_of_interest(combined_final_gene_list, high_D5_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D5 Bulk for HVL Subjects", high_bulk_D5_combined_info[[1]], "gene_auc")

high_bulk_D8_combined_info <- find_aucs_of_interest(combined_final_gene_list, high_D8_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D8 Bulk for HVL Subjects", high_bulk_D8_combined_info[[1]], "gene_auc")

# This one is a little redundant with the individual D28 sc/multiome objects but it may still be convenient for something!
high_bulk_D28_combined_info <- find_aucs_of_interest(combined_final_gene_list, high_D28_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_bulk_D28_combined_info[[1]], "gene_auc")

low_bulk_D2_combined_info <- find_aucs_of_interest(combined_final_gene_list, low_D2_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D2 Bulk for LVL Subjects", low_bulk_D2_combined_info[[1]], "gene_auc")

low_bulk_D5_combined_info <- find_aucs_of_interest(combined_final_gene_list, low_D5_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D5 Bulk for LVL Subjects", low_bulk_D5_combined_info[[1]], "gene_auc")

low_bulk_D8_combined_info <- find_aucs_of_interest(combined_final_gene_list, low_D8_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D8 Bulk for LVL Subjects", low_bulk_D8_combined_info[[1]], "gene_auc")

low_bulk_D28_combined_info <- find_aucs_of_interest(combined_final_gene_list, low_D28_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D28 Bulk for LVL Subjects", low_bulk_D28_combined_info[[1]], "gene_auc")




# Find list of significant genes in bulk data - all (LRT)
placebo_period_2_LRT_metadata <- placebo_metadata[placebo_metadata$time_point == "2_D28" | placebo_metadata$time_point == "2_D8" | 
                                                              placebo_metadata$time_point == "2_D5" | placebo_metadata$time_point == "2_D2" |
                                                              placebo_metadata$time_point == "2_D_minus_1",]
placebo_period_2_LRT_metadata <- placebo_period_2_LRT_metadata[placebo_period_2_LRT_metadata$subject_id 
                                                                         %in% names(table(placebo_period_2_LRT_metadata$subject_id)
                                                                                    [table(placebo_period_2_LRT_metadata$subject_id) == 5]),]
placebo_period_2_LRT_counts <- placebo_counts[rownames(placebo_period_2_LRT_metadata)]
placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = placebo_period_2_LRT_counts,
                                                             colData = placebo_period_2_LRT_metadata,
                                                             design = ~ subject_id + time_point)
placebo_period_2_LRT_analysis <- DESeq(placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id)
placebo_period_2_LRT_analysis_results <- results(placebo_period_2_LRT_analysis, alpha = 0.05)
placebo_period_2_LRT_analysis_results <- placebo_period_2_LRT_analysis_results[order(placebo_period_2_LRT_analysis_results$padj),]
placebo_period_2_LRT_analysis_results <- subset(placebo_period_2_LRT_analysis_results, padj < 0.05)

# Find list of significant genes in bulk data - high (LRT)
high_placebo_period_2_LRT_metadata <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D28" | high_placebo_metadata$time_point == "2_D8" | 
                                                              high_placebo_metadata$time_point == "2_D5" | high_placebo_metadata$time_point == "2_D2" |
                                                              high_placebo_metadata$time_point == "2_D_minus_1",]
high_placebo_period_2_LRT_metadata <- high_placebo_period_2_LRT_metadata[high_placebo_period_2_LRT_metadata$subject_id 
                                                                         %in% names(table(high_placebo_period_2_LRT_metadata$subject_id)
                                                                                    [table(high_placebo_period_2_LRT_metadata$subject_id) == 5]),]
high_placebo_period_2_LRT_counts <- high_placebo_counts[rownames(high_placebo_period_2_LRT_metadata)]
high_placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = high_placebo_period_2_LRT_counts,
                                                             colData = high_placebo_period_2_LRT_metadata,
                                                             design = ~ subject_id + time_point)
high_placebo_period_2_LRT_analysis <- DESeq(high_placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id)
high_placebo_period_2_LRT_analysis_results <- results(high_placebo_period_2_LRT_analysis, alpha = 0.05)
high_placebo_period_2_LRT_analysis_results <- high_placebo_period_2_LRT_analysis_results[order(high_placebo_period_2_LRT_analysis_results$padj),]
high_placebo_period_2_LRT_analysis_results <- subset(high_placebo_period_2_LRT_analysis_results, padj < 0.05)

# Find list of significant genes in bulk data - low (LRT)
low_placebo_period_2_LRT_metadata <- low_placebo_metadata[low_placebo_metadata$time_point == "2_D28" | low_placebo_metadata$time_point == "2_D8" | 
                                                              low_placebo_metadata$time_point == "2_D5" | low_placebo_metadata$time_point == "2_D2" |
                                                              low_placebo_metadata$time_point == "2_D_minus_1",]
low_placebo_period_2_LRT_metadata <- low_placebo_period_2_LRT_metadata[low_placebo_period_2_LRT_metadata$subject_id 
                                                                         %in% names(table(low_placebo_period_2_LRT_metadata$subject_id)
                                                                                    [table(low_placebo_period_2_LRT_metadata$subject_id) == 5]),]
low_placebo_period_2_LRT_counts <- low_placebo_counts[rownames(low_placebo_period_2_LRT_metadata)]
low_placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = low_placebo_period_2_LRT_counts,
                                                             colData = low_placebo_period_2_LRT_metadata,
                                                             design = ~ subject_id + time_point)
low_placebo_period_2_LRT_analysis <- DESeq(low_placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id)
low_placebo_period_2_LRT_analysis_results <- results(low_placebo_period_2_LRT_analysis, alpha = 0.05)
low_placebo_period_2_LRT_analysis_results <- low_placebo_period_2_LRT_analysis_results[order(low_placebo_period_2_LRT_analysis_results$padj),]
low_placebo_period_2_LRT_analysis_results <- subset(low_placebo_period_2_LRT_analysis_results, padj < 0.05)

# Grab list of genes that have AUC > 0.7 and AUC < 0.3 in D28 bulk RNA-seq
sc_pos_genes <- all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$all_bulk_D28_gene_auc > 0.7,]$gene_name
sc_neg_genes <- all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$all_bulk_D28_gene_auc < 0.3,]$gene_name

### ALL BULK PLACEBO DATA ###
# See which of these genes are significant in LRT data
sc_pos_genes_LRT_pass <- c()
for(gene in sc_pos_genes) {
  if(gene %in% rownames(placebo_period_2_LRT_analysis_results)) {
    sc_pos_genes_LRT_pass <- c(sc_pos_genes_LRT_pass, gene)
  }
}

placebo_period_2_LRT_analysis_betas <- coef(placebo_period_2_LRT_analysis)
placebo_period_2_LRT_analysis_betas <- placebo_period_2_LRT_analysis_betas[, -c(1:22)]
placebo_period_2_LRT_analysis_betas <- placebo_period_2_LRT_analysis_betas[rownames(placebo_period_2_LRT_analysis_betas) %in% sc_pos_genes_LRT_pass,]
placebo_period_2_LRT_analysis_thr <- 1
colnames(placebo_period_2_LRT_analysis_betas) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(placebo_period_2_LRT_analysis_betas, breaks=seq(from=-placebo_period_2_LRT_analysis_thr, to=placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14)

# See which of these genes are significant in LRT data - why are only 21/32 found to be significant? Because ALL vs HIGH?
sc_neg_genes_LRT_pass <- c()
for(gene in sc_neg_genes) {
  if(gene %in% rownames(placebo_period_2_LRT_analysis_results)) {
    sc_neg_genes_LRT_pass <- c(sc_neg_genes_LRT_pass, gene)
  }
}

placebo_period_2_LRT_analysis_betas_neg <- coef(placebo_period_2_LRT_analysis)
placebo_period_2_LRT_analysis_betas_neg <- placebo_period_2_LRT_analysis_betas_neg[, -c(1:22)]
placebo_period_2_LRT_analysis_betas_neg <- placebo_period_2_LRT_analysis_betas_neg[rownames(placebo_period_2_LRT_analysis_betas_neg) %in% sc_neg_genes_LRT_pass,]
colnames(placebo_period_2_LRT_analysis_betas_neg) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(placebo_period_2_LRT_analysis_betas_neg, breaks=seq(from=-1.5, to=1.5, length=101),
         cluster_col=FALSE, fontsize_col=14)

### HIGH BULK PLACEBO DATA ###
# See which of these genes are significant in LRT data
high_sc_pos_genes_LRT_pass <- c()
for(gene in sc_pos_genes) {
  if(gene %in% rownames(high_placebo_period_2_LRT_analysis_results)) {
    high_sc_pos_genes_LRT_pass <- c(high_sc_pos_genes_LRT_pass, gene)
  }
}

high_placebo_period_2_LRT_analysis_betas <- coef(high_placebo_period_2_LRT_analysis)
high_placebo_period_2_LRT_analysis_betas <- high_placebo_period_2_LRT_analysis_betas[, -c(1:13)]
high_placebo_period_2_LRT_analysis_betas <- high_placebo_period_2_LRT_analysis_betas[rownames(high_placebo_period_2_LRT_analysis_betas) %in% high_sc_pos_genes_LRT_pass,]
high_placebo_period_2_LRT_analysis_thr <- 1.8
colnames(high_placebo_period_2_LRT_analysis_betas) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(high_placebo_period_2_LRT_analysis_betas, breaks=seq(from=-high_placebo_period_2_LRT_analysis_thr, to=high_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize = 30, width = 16, height = 12, filename = "C:/Users/willi/Desktop/testing_low_pos.png", cex = 0.7)

# See which of these genes are significant in LRT data - why are only 21/32 found to be significant? Because ALL vs HIGH?
high_sc_neg_genes_LRT_pass <- c()
for(gene in sc_neg_genes) {
  if(gene %in% rownames(high_placebo_period_2_LRT_analysis_results)) {
    high_sc_neg_genes_LRT_pass <- c(high_sc_neg_genes_LRT_pass, gene)
  }
}

high_placebo_period_2_LRT_analysis_betas_neg <- coef(high_placebo_period_2_LRT_analysis)
high_placebo_period_2_LRT_analysis_betas_neg <- high_placebo_period_2_LRT_analysis_betas_neg[, -c(1:13)]
high_placebo_period_2_LRT_analysis_betas_neg <- high_placebo_period_2_LRT_analysis_betas_neg[rownames(high_placebo_period_2_LRT_analysis_betas_neg) %in% high_sc_neg_genes_LRT_pass,]
colnames(high_placebo_period_2_LRT_analysis_betas_neg) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(high_placebo_period_2_LRT_analysis_betas_neg, breaks=seq(from=-2, to=2, length=101),
         cluster_col=FALSE, fontsize = 22, width = 16, height = 12, filename = "C:/Users/willi/Desktop/testing_high_neg.png", cex = 0.8)


# Grab DEGs for 2_D28 vs 2_D_minus_1 for bulk (Wald test)
placebo_period_2_D28_vs_D_minus_1_bulk_results <- run_deseq_bulk_analysis_time_series("placebo", placebo_counts,placebo_metadata,
                                                                                      "2_D28", "2_D_minus_1", "~/")



# Look at how data compendium papers use convalescent data

# Find intersecting genes between bulk and pseudobulk
# What about genes that have AUC < 0.3 or AUC > 0.7 for ALL time points? How many of those are there?
# Are those genes interesting or not interesting?
# Should we include subsets of times?
# See prediction in bulk datasets using pos and neg signatures, I guess?
# See how many other genes have same trend as BAG1 (pos AUC in pseudo and fold change is generally negative or vice versa)



# Genes that passed D28 bulk and all bulk RNA-seq (AUC > 0.7) was nothing
#all_pass_sc_pos_genes <- all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$all_bulk_D28_gene_auc > 0.7,]$gene_name
#all_pass_sc_pos_genes <- intersect(all_pass_sc_pos_genes, all_bulk_D2_sc_pseudobulk_gene_aucs[all_bulk_D2_sc_pseudobulk_gene_aucs$all_bulk_D2_gene_auc > 0.7,]$gene_name)
#all_pass_sc_pos_genes <- intersect(all_pass_sc_pos_genes, all_bulk_D5_sc_pseudobulk_gene_aucs[all_bulk_D5_sc_pseudobulk_gene_aucs$all_bulk_D5_gene_auc > 0.7,]$gene_name)
#all_pass_sc_pos_genes <- intersect(all_pass_sc_pos_genes, all_bulk_D8_sc_pseudobulk_gene_aucs[all_bulk_D8_sc_pseudobulk_gene_aucs$all_bulk_D8_gene_auc > 0.7,]$gene_name)

# Genes that passed D28 bulk and all bulk RNA-seq (AUC < 0.3) were "OST4"    "C9orf78" "TMA7"    "UQCR11"  "BAG1" 
#all_pass_sc_neg_genes <- all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$all_bulk_D28_gene_auc < 0.3,]$gene_name
#all_pass_sc_neg_genes <- intersect(all_pass_sc_neg_genes, all_bulk_D2_sc_pseudobulk_gene_aucs[all_bulk_D2_sc_pseudobulk_gene_aucs$all_bulk_D2_gene_auc < 0.3,]$gene_name)
#all_pass_sc_neg_genes <- intersect(all_pass_sc_neg_genes, all_bulk_D5_sc_pseudobulk_gene_aucs[all_bulk_D5_sc_pseudobulk_gene_aucs$all_bulk_D5_gene_auc < 0.3,]$gene_name)
#all_pass_sc_neg_genes <- intersect(all_pass_sc_neg_genes, all_bulk_D8_sc_pseudobulk_gene_aucs[all_bulk_D8_sc_pseudobulk_gene_aucs$all_bulk_D8_gene_auc < 0.3,]$gene_name)

# Genes that passed D28 bulk and all high bulk RNA-seq (AUC > 0.7)
high_sc_pos_genes <- high_non_sc_bulk_D28_sc_pseudobulk_gene_aucs[high_non_sc_bulk_D28_sc_pseudobulk_gene_aucs$high_bulk_D28_minus_training_gene_auc > 0.7,]$gene_name
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D2_sc_pseudobulk_gene_aucs[high_bulk_D2_sc_pseudobulk_gene_aucs$high_bulk_D2_gene_auc > 0.7,]$gene_name)
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]$gene_name)
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D8_sc_pseudobulk_gene_aucs[high_bulk_D8_sc_pseudobulk_gene_aucs$high_bulk_D8_gene_auc > 0.7,]$gene_name)
# Genes that passed D28 bulk and all high bulk RNA-seq (AUC < 0.3)
high_sc_neg_genes <- high_non_sc_bulk_D28_sc_pseudobulk_gene_aucs[high_non_sc_bulk_D28_sc_pseudobulk_gene_aucs$high_bulk_D28_minus_training_gene_auc < 0.3,]$gene_name
high_sc_neg_genes <- intersect(high_sc_neg_genes, high_bulk_D2_sc_pseudobulk_gene_aucs[high_bulk_D2_sc_pseudobulk_gene_aucs$high_bulk_D2_gene_auc < 0.3,]$gene_name)
high_sc_neg_genes <- intersect(high_sc_neg_genes, high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc < 0.3,]$gene_name)
high_sc_neg_genes <- intersect(high_sc_neg_genes, high_bulk_D8_sc_pseudobulk_gene_aucs[high_bulk_D8_sc_pseudobulk_gene_aucs$high_bulk_D8_gene_auc < 0.3,]$gene_name)
# No genes passed pseudobulk and all low bulk RNA-seq (AUC > 0.7)
low_sc_pos_genes <- high_bulk_D28_sc_pseudobulk_gene_aucs[high_bulk_D28_sc_pseudobulk_gene_aucs$high_bulk_D28_gene_auc > 0.7,]$gene_name
low_sc_pos_genes <- intersect(low_sc_pos_genes, low_bulk_D2_sc_pseudobulk_gene_aucs[low_bulk_D2_sc_pseudobulk_gene_aucs$low_bulk_D2_gene_auc > 0.7,]$gene_name)
low_sc_pos_genes <- intersect(low_sc_pos_genes, low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc > 0.7,]$gene_name)
low_sc_pos_genes <- intersect(low_sc_pos_genes, low_bulk_D8_sc_pseudobulk_gene_aucs[low_bulk_D8_sc_pseudobulk_gene_aucs$low_bulk_D8_gene_auc > 0.7,]$gene_name)
low_sc_pos_genes <- intersect(low_sc_pos_genes, low_bulk_D28_sc_pseudobulk_gene_aucs[low_bulk_D28_sc_pseudobulk_gene_aucs$low_bulk_D28_gene_auc > 0.7,]$gene_name)
# No genes passed pseudobulk and all low bulk RNA-seq (AUC < 0.3)
low_sc_neg_genes <- high_bulk_D28_sc_pseudobulk_gene_aucs[high_bulk_D28_sc_pseudobulk_gene_aucs$high_bulk_D28_gene_auc < 0.3,]$gene_name
low_sc_neg_genes <- intersect(low_sc_neg_genes, low_bulk_D2_sc_pseudobulk_gene_aucs[low_bulk_D2_sc_pseudobulk_gene_aucs$low_bulk_D2_gene_auc < 0.3,]$gene_name)
low_sc_neg_genes <- intersect(low_sc_neg_genes, low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc < 0.3,]$gene_name)
low_sc_neg_genes <- intersect(low_sc_neg_genes, low_bulk_D8_sc_pseudobulk_gene_aucs[low_bulk_D8_sc_pseudobulk_gene_aucs$low_bulk_D8_gene_auc < 0.3,]$gene_name)
low_sc_neg_genes <- intersect(low_sc_neg_genes, low_bulk_D28_sc_pseudobulk_gene_aucs[low_bulk_D28_sc_pseudobulk_gene_aucs$low_bulk_D28_gene_auc < 0.3,]$gene_name)
# Genes that passed 

# Test to see if dumb way of doing FC and p-value thresholding results in "better" alignment with what people generally think


# See if those are usually found in MAGICAL

# For SC data
# Sample 1 - 318475c9c36293a5 - HVL
# Sample 2 - fac0a26252aed8d9 - HVL
# Sample 3 - e5c7a7c9f56d67f1 - HVL
# Sample 4 - 2418b55aca2c78d6 - LVL
# Sample 5 - 82302d3a293e261b - LVL
# Sample 6 - 76690d37e925ba09 - HVL

# 80633 HVL
# 35496 LVL

# Humanbase for DEGs in Bulk
# Humanbase for DEGs in cell types (pseudobulk or no? No, right?)

# Testing on bulk D5 - no wonder there isn't perfect overlap, right? Test bulk for D28 and see what you find!
# Interesting idea - compare AUCs in our gene list(s) to other datasets (e.g., Data Compendium discovery or our own bulk RNA-seq)
# Do I need to remove 0 qPCR value person from LVL for BULK as well? Should we just completely remove this individual from consideration?




# Interesting idea - compare AUCs in our gene list(s) to other datasets (e.g., Data Compendium discovery or our own bulk RNA-seq)
# # IDEA - separate into HVL and LVL and see whether it has higher AUC for higher viral load ppl
comparison_aucs <- sc_pseudobulk_aucs
comparison_aucs$discovery_aucs <- sc_discovery_flu_aucs$flu_discovery_gene_auc
comparison_aucs <- comparison_aucs[,c(1,2,4,3)]



curated_sc_pseudobulk_gene_aucs
sc_bulk_D28_sc_pseudobulk_gene_aucs

###### MISC UNORDERED STUFF ######

# DESeq2 analysis comparing 4 SC subjects to the 9 non-SC subjects for each time point demonstrates that there is no significant difference between the groups of samples with respect to gene expression
testing_sc_vs_other_metadata <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D28" | high_placebo_metadata$time_point == "2_D_minus_1",]
testing_sc_vs_other_subject_tag <- testing_sc_vs_other_metadata$subject_id %in% sc_subjects
testing_sc_vs_other_subject_tag <- replace(testing_sc_vs_other_subject_tag, testing_sc_vs_other_subject_tag == FALSE, "OTHER")
testing_sc_vs_other_subject_tag <- replace(testing_sc_vs_other_subject_tag, testing_sc_vs_other_subject_tag == TRUE, "SC")
testing_sc_vs_other_metadata$subject_tag <- testing_sc_vs_other_subject_tag
testing_sc_vs_other_metadata <- testing_sc_vs_other_metadata[testing_sc_vs_other_metadata$subject_id  %in% names(table(testing_sc_vs_other_metadata$subject_id)[table(testing_sc_vs_other_metadata$subject_id) == 2]),]
# D28 testing
testing_sc_vs_other_metadata_d28 <- testing_sc_vs_other_metadata[testing_sc_vs_other_metadata$time_point == "2_D28",]
testing_sc_vs_other_counts_d28 <- high_placebo_counts[rownames(testing_sc_vs_other_metadata_d28)]
testing_sc_vs_other_analysis_d28 <- DESeqDataSetFromMatrix(countData = testing_sc_vs_other_counts_d28, colData = testing_sc_vs_other_metadata_d28, design = ~ subject_tag)
testing_sc_vs_other_analysis_d28 <- DESeq(testing_sc_vs_other_analysis_d28)
testing_sc_vs_other_analysis_results_d28 <- results(testing_sc_vs_other_analysis_d28, contrast = c("subject_tag", "SC", "OTHER"), alpha = 0.05, lfcThreshold = 0.1)
testing_sc_vs_other_analysis_results_pseudobulk_gene_subset_d28 <- testing_sc_vs_other_analysis_results_d28[rownames(testing_sc_vs_other_analysis_results_d28) %in% sc_pseudobulk_genes,]
testing_sc_vs_other_analysis_results_pseudobulk_gene_subset_d28 <- testing_sc_vs_other_analysis_results_pseudobulk_gene_subset_d28[order(testing_sc_vs_other_analysis_results_pseudobulk_gene_subset_d28$padj),]
# D_minus_1 testing
testing_sc_vs_other_metadata_d_minus_1 <- testing_sc_vs_other_metadata[testing_sc_vs_other_metadata$time_point == "2_D_minus_1",]
testing_sc_vs_other_counts_d_minus_1 <- high_placebo_counts[rownames(testing_sc_vs_other_metadata_d_minus_1)]
testing_sc_vs_other_analysis_d_minus_1 <- DESeqDataSetFromMatrix(countData = testing_sc_vs_other_counts_d_minus_1, colData = testing_sc_vs_other_metadata_d_minus_1, design = ~ subject_tag)
testing_sc_vs_other_analysis_d_minus_1 <- DESeq(testing_sc_vs_other_analysis_d_minus_1)
testing_sc_vs_other_analysis_results_d_minus_1 <- results(testing_sc_vs_other_analysis_d_minus_1, contrast = c("subject_tag", "SC", "OTHER"), alpha = 0.05, lfcThreshold = 0.1)
testing_sc_vs_other_analysis_results_pseudobulk_gene_subset_d_minus_1 <- testing_sc_vs_other_analysis_results_d_minus_1[rownames(testing_sc_vs_other_analysis_results_d_minus_1) %in% sc_pseudobulk_genes,]
testing_sc_vs_other_analysis_results_pseudobulk_gene_subset_d_minus_1 <- testing_sc_vs_other_analysis_results_pseudobulk_gene_subset_d_minus_1[order(testing_sc_vs_other_analysis_results_pseudobulk_gene_subset_d_minus_1$padj),]

# TODO: Should do the same for multiome!

# Below, we do analysis using pseudobulk (HVL) and specific subjects analyzed using single-cell / multiome (HVL)
# We don't really use these analyses currently

# Create log transformed pseudobulk count tables (HVL only!)
sc_pseudobulk_counts_log_transformed <- grab_transformed_pseudobulk_counts(sc_pseudobulk_dir, possible_cell_types)
multiome_14_pseudobulk_counts_log_transformed <- grab_transformed_pseudobulk_counts(multiome_pseudobulk_dir, possible_cell_types)

# Create MetaIntegrator objects using pseudobulk count tables (HVL only!)
sc_pseudobulk_metaintegrator_obj <- create_metaintegrator_obj("mine", sc_pseudobulk_counts_log_transformed)
multiome_pseudobulk_metaintegrator_obj <- create_metaintegrator_obj("mine", multiome_14_pseudobulk_counts_log_transformed)

# Create MetaIntegrator objects for the specific subjects that we used for single-cell / multiome (HVL only!)
sc_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_sc_placebo_counts, high_sc_placebo_metadata, "2_D28", "2_D_minus_1")
multiome_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_multiome_placebo_counts, high_multiome_placebo_metadata, "2_D28", "2_D_minus_1")

# Create MetaIntegrator objects for the specific subjects that we didn't use for single-cell / multiome (HVL only!)
# I can use these to see which pseudobulk genes have AUC > 0.7 for the subjects not used in single cell / multiome analysis
# However, because there's no real difference between the subjects chosen for single cell / multiome analysis and the other subjects,
# as demonstrated above, it doesn't really make sense to test only on these subjects
non_sc_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_non_sc_placebo_counts, high_non_sc_placebo_metadata, "2_D28", "2_D_minus_1")
non_multiome_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_non_multiome_placebo_counts, high_non_multiome_placebo_metadata, "2_D28", "2_D_minus_1")

# AUCS ON PSEUDOBULK DATA - I am not really using these data for anything currently

# Single Cell (Pseudobulk)
sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_pseudobulk_genes, sc_pseudobulk_metaintegrator_obj, "sc_paired", "sc_pseudobulk"))
curated_sc_pseudobulk_gene_aucs <- sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc < 0.3 | sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Single Cell", "Total Pseudobulk", sc_pseudobulk_gene_aucs, "sc_pseudobulk_gene_auc")

# Single Cell (MAGICAL)
sc_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_magical_genes, sc_pseudobulk_metaintegrator_obj, "sc_paired", "sc_magical"))
curated_sc_magical_gene_aucs <- sc_magical_gene_aucs[sc_magical_gene_aucs$sc_magical_gene_auc < 0.3 | sc_magical_gene_aucs$sc_magical_gene_auc > 0.7,]
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Single Cell", "Total Pseudobulk", sc_magical_gene_aucs, "sc_magical_gene_auc")

# Multiome (Pseudobulk)
multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(multiome_pseudobulk_genes, multiome_pseudobulk_metaintegrator_obj, "multiome_paired", "multiome_pseudobulk"))
curated_multiome_pseudobulk_gene_aucs <- multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc < 0.3 | multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc > 0.7,]
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Multiome", "Total Pseudobulk", multiome_pseudobulk_gene_aucs, "multiome_pseudobulk_gene_auc")

# Multiome (MAGICAL)
multiome_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(multiome_magical_genes, multiome_pseudobulk_metaintegrator_obj, "multiome_paired", "multiome_magical"))
curated_multiome_magical_gene_aucs <- multiome_magical_gene_aucs[multiome_magical_gene_aucs$multiome_magical_gene_auc < 0.3 | multiome_magical_gene_aucs$multiome_magical_gene_auc > 0.7,]
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Multiome", "Total Pseudobulk", multiome_magical_gene_aucs, "multiome_magical_gene_auc")

# AUCS ON BULK (SAME SUBJECTS AS SINGLE CELL) - I am not really using these data for anything currently

# Single Cell (Pseudobulk)
sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_pseudobulk_genes, sc_D28_bulk_metaintegrator_obj, "sc_paired", "sc_bulk_D28"))
curated_sc_bulk_D28_sc_pseudobulk_gene_aucs <- sc_bulk_D28_sc_pseudobulk_gene_aucs[sc_bulk_D28_sc_pseudobulk_gene_aucs$sc_bulk_D28_gene_auc < 0.3 | sc_bulk_D28_sc_pseudobulk_gene_aucs$sc_bulk_D28_gene_auc > 0.7,]
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Single Cell", "D28 Bulk for Single Cell Subjects", sc_bulk_D28_sc_pseudobulk_gene_aucs, "sc_bulk_D28_gene_auc")

# Single Cell (MAGICAL)
sc_bulk_D28_sc_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_magical_genes, sc_D28_bulk_metaintegrator_obj, "sc_paired", "sc_bulk_D28_magical"))
curated_sc_bulk_D28_sc_magical_gene_aucs <- sc_bulk_D28_sc_magical_gene_aucs[sc_bulk_D28_sc_magical_gene_aucs$sc_bulk_D28_magical_gene_auc < 0.3 | sc_bulk_D28_sc_magical_gene_aucs$sc_bulk_D28_magical_gene_auc > 0.7,]
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Single Cell", "D28 Bulk for Single Cell Subjects", sc_bulk_D28_sc_magical_gene_aucs, "sc_bulk_D28_magical_gene_auc")

# Multiome (Pseudobulk)
multiome_bulk_D28_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(multiome_pseudobulk_genes, multiome_D28_bulk_metaintegrator_obj, "multiome_paired", "multiome_bulk_D28"))
curated_multiome_bulk_D28_multiome_pseudobulk_gene_aucs <- multiome_bulk_D28_multiome_pseudobulk_gene_aucs[multiome_bulk_D28_multiome_pseudobulk_gene_aucs$multiome_bulk_D28_gene_auc < 0.3 | multiome_bulk_D28_multiome_pseudobulk_gene_aucs$multiome_bulk_D28_gene_auc > 0.7,]
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Single Cell", "D28 Bulk for Multiome Subjects", multiome_bulk_D28_multiome_pseudobulk_gene_aucs, "sc_bulk_D28_gene_auc")

# Multiome (MAGICAL)
multiome_bulk_D28_multiome_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_magical_genes, multiome_D28_bulk_metaintegrator_obj, "multiome_paired", "multiome_bulk_D28_magical"))
curated_multiome_bulk_D28_multiome_magical_gene_aucs <- multiome_bulk_D28_multiome_magical_gene_aucs[multiome_bulk_D28_multiome_magical_gene_aucs$multiome_bulk_D28_magical_gene_auc < 0.3 | multiome_bulk_D28_multiome_magical_gene_aucs$multiome_bulk_D28_magical_gene_auc > 0.7,]
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Single Cell", "D28 Bulk for Multiome Subjects", multiome_bulk_D28_multiome_magical_gene_aucs, "multiome_bulk_D28_magical_gene_auc")

### CALCULATE AUCS ON D28 HVL BULK (MINUS TRAINING SAMPLES) - see explanation above on why we don't care about this analysis
high_non_sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_pseudobulk_genes, non_sc_D28_bulk_metaintegrator_obj, "sc_paired", "high_bulk_D28_minus_training"))
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Single Cell", "D28 Bulk for HVL Subjects (Minus Training)", high_non_sc_bulk_D28_sc_pseudobulk_gene_aucs, "high_bulk_D28_minus_training_gene_auc")

### CALCULATE AUCS ON ALL D28 BULK DATA
# We don't actively use this analysis because it's much better to separate HVL and LVL samples
all_bulk_D2_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_pseudobulk_genes, all_D2_bulk_metaintegrator_obj, "sc_paired", "all_bulk_D2"))
all_bulk_D5_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_pseudobulk_genes, all_D5_bulk_metaintegrator_obj, "sc_paired", "all_bulk_D5"))
all_bulk_D8_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_pseudobulk_genes, all_D8_bulk_metaintegrator_obj, "sc_paired", "all_bulk_D8"))
all_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(sc_pseudobulk_genes, all_D28_bulk_metaintegrator_obj, "sc_paired", "all_bulk_D28"))

