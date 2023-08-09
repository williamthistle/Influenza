# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Create MetaIntegrator objects for bulk for all days (HVL and LVL)
high_D2_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D2", "2_D_minus_1")
high_D5_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D5", "2_D_minus_1")
high_D8_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D8", "2_D_minus_1")
high_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D28", "2_D_minus_1")

high_D28_non_sc_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_non_sc_placebo_counts, high_non_sc_placebo_metadata, "2_D28", "2_D_minus_1")
high_D28_non_multiome_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_non_multiome_placebo_counts, high_non_multiome_placebo_metadata, "2_D28", "2_D_minus_1")


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
auc_df <- data.frame(Filtering_Assay = character(), Filtering_Method = character(), Discovery_Assay = character(), 
                     Discovery_Dataset = character(), Pos_Genes = integer(), Neg_Genes = integer(), Total_Passing_Genes = integer(), 
                     Total_Genes = integer(), Percentage_of_Passing_Genes = double(), stringsAsFactors = FALSE)
auc_names <- c("Filtering_Assay", "Filtering_Method", "Discovery_Assay", "Discovery_Dataset", "Passing_Pos_Genes", "Passing_Neg_Genes", "Total_Passing_Genes", "Total_Genes", "Percentage_of_Passing_Genes")

# AUCS ON D28 HVL BULK
# Note that our initial gene lists (from pseudobulk filtering / MAGICAL) are further curated using all D28 HVL bulk data

# SC (pseudobulk filtering)
high_sc_bulk_D28_sc_pseudobulk_gene_info <- find_aucs_of_interest(sc_pseudobulk_gene_table, high_D28_non_sc_bulk_metaintegrator_obj, "sc_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_sc_bulk_D28_sc_pseudobulk_gene_info[[1]], "gene_auc")

# SC (MAGICAL filtering)
#high_sc_bulk_D28_sc_magical_gene_info <- find_aucs_of_interest(sc_magical_gene_table, high_D28_bulk_metaintegrator_obj, "sc_paired")
#auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_sc_bulk_D28_sc_magical_gene_info[[1]], "gene_auc")

# Multiome (pseudobulk filtering)
high_multiome_bulk_D28_multiome_pseudobulk_gene_info <- find_aucs_of_interest(multiome_pseudobulk_gene_table, high_D28_non_multiome_bulk_metaintegrator_obj, "multiome_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_multiome_bulk_D28_multiome_pseudobulk_gene_info[[1]], "gene_auc")

# Multiome (MAGICAL filtering)
#high_multiome_bulk_D28_multiome_magical_gene_info <- find_aucs_of_interest(multiome_magical_genes, high_D28_bulk_metaintegrator_obj, "multiome_paired")
#auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_multiome_bulk_D28_multiome_magical_gene_info[[1]], "gene_auc")

# Next, we combine SC pseudobulk genes and multiome pseudobulk genes for a combined search across HVL D2/D5/D8 and LVL D28 (not much signal in LVL D2/D5/D8)
# MAGICAL genes may be interesting to look at further down the line
# Are we seeing interesting behavior out of the MAGICAL genes for D2 / D5 / D8 / D28 or for LVL D28?
# Which genes come up in both SC and multiome?
intersecting_genes_between_sc_and_multiome <- intersect(high_sc_bulk_D28_sc_pseudobulk_gene_info[[2]], high_multiome_bulk_D28_multiome_pseudobulk_gene_info[[2]])
combined_final_gene_list <- unique(c(high_sc_bulk_D28_sc_pseudobulk_gene_info[[2]], high_multiome_bulk_D28_multiome_pseudobulk_gene_info[[2]]))

# Find AUCs for combined gene list for HVL D2/D5/D8
high_bulk_D2_combined_info <- find_aucs_of_interest(combined_final_gene_list, high_D2_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D2 Bulk for HVL Subjects", high_bulk_D2_combined_info[[1]], "gene_auc")

high_bulk_D5_combined_info <- find_aucs_of_interest(combined_final_gene_list, high_D5_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D5 Bulk for HVL Subjects", high_bulk_D5_combined_info[[1]], "gene_auc")

high_bulk_D8_combined_info <- find_aucs_of_interest(combined_final_gene_list, high_D8_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D8 Bulk for HVL Subjects", high_bulk_D8_combined_info[[1]], "gene_auc")

# Combines info from SC and multiome D28 analyses above
high_bulk_D28_combined_info <- find_aucs_of_interest(combined_final_gene_list, high_D28_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D28 Bulk for HVL Subjects", high_bulk_D28_combined_info[[1]], "gene_auc")

# Find AUCs for combined gene list for LVL D2/D5/D8/D28
low_bulk_D2_combined_info <- find_aucs_of_interest(combined_final_gene_list, low_D2_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D2 Bulk for LVL Subjects", low_bulk_D2_combined_info[[1]], "gene_auc")

low_bulk_D5_combined_info <- find_aucs_of_interest(combined_final_gene_list, low_D5_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D5 Bulk for LVL Subjects", low_bulk_D5_combined_info[[1]], "gene_auc")

low_bulk_D8_combined_info <- find_aucs_of_interest(combined_final_gene_list, low_D8_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D8 Bulk for LVL Subjects", low_bulk_D8_combined_info[[1]], "gene_auc")

low_bulk_D28_combined_info <- find_aucs_of_interest(combined_final_gene_list, low_D28_bulk_metaintegrator_obj, "combined_paired")
auc_df <- add_auc_row(auc_df, auc_names, "Combined", "Cell Type Pseudobulk", "Bulk RNA-Seq", "D28 Bulk for LVL Subjects", low_bulk_D28_combined_info[[1]], "gene_auc")

# Find list of significant genes in bulk data - high (LRT)
high_LRT_analysis_results_info <- run_deseq2_LRT(high_placebo_counts, high_placebo_metadata)

# Find list of significant genes in bulk data - low (LRT)
low_LRT_analysis_results_info <- run_deseq2_LRT(low_placebo_counts, low_placebo_metadata)

# Plot heatmaps for our combined gene list (pos and neg) from SC and multiome (HVL D28 bulk curated) 
# TODO: Decide whether we should look at fold change of genes that aren't significant
# We have already decided these genes are important due to AUC > 0.7, so why do we need them to be significant in LRT analysis?
# This could be even more relevant for LVL since nothing is significant (but maybe we still care about fold change!)
plot_lrt_heatmap(high_bulk_D28_combined_info[[3]], high_LRT_analysis_results_info, 1.8, "C:/Users/willi/Desktop/testing_high_pos.png")
plot_lrt_heatmap(high_bulk_D28_combined_info[[4]], high_LRT_analysis_results_info, 2, "C:/Users/willi/Desktop/testing_high_neg.png")

# Find genes that have AUC > 0.7 for D2 / D5 / D8 (HVL)
high_pos_auc_all_time_points <- high_bulk_D28_combined_info[[3]]
high_pos_auc_all_time_points <- intersect(high_pos_auc_all_time_points, high_bulk_D2_combined_info[[3]])
high_pos_auc_all_time_points <- intersect(high_pos_auc_all_time_points, high_bulk_D5_combined_info[[3]])
high_pos_auc_all_time_points <- intersect(high_pos_auc_all_time_points, high_bulk_D8_combined_info[[3]])

# Find genes that have AUC < 0.3 for D2 / D5 / D8 (HVL)
high_neg_auc_all_time_points <- high_bulk_D28_combined_info[[4]]
high_neg_auc_all_time_points <- intersect(high_neg_auc_all_time_points, high_bulk_D2_combined_info[[4]])
high_neg_auc_all_time_points <- intersect(high_neg_auc_all_time_points, high_bulk_D5_combined_info[[4]])
high_neg_auc_all_time_points <- intersect(high_neg_auc_all_time_points, high_bulk_D8_combined_info[[4]])

# Find genes that have AUC > 0.7 for D28 (LVL)
low_pos_auc_d28 <- high_bulk_D28_combined_info[[3]]
low_pos_auc_d28 <- intersect(low_pos_auc_d28, low_bulk_D28_combined_info[[3]])

# Find genes that have AUC < 0.3 for D28 (LVL)
low_neg_auc_d28 <- high_bulk_D28_combined_info[[4]]
low_neg_auc_d28 <- intersect(low_neg_auc_d28, low_bulk_D28_combined_info[[4]])



# Look at how data compendium papers use convalescent data

# Should we include subsets of times?
# See prediction in bulk datasets using pos and neg signatures, I guess?

# Humanbase for DEGs in Bulk
# Humanbase for DEGs in cell types (pseudobulk or no? No, right?)

# Interesting idea - compare AUCs in our gene list(s) to other datasets (e.g., Data Compendium discovery or our own bulk RNA-seq)

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
# TODO: Need to edit possible_cell_types
#sc_pseudobulk_counts_log_transformed <- grab_transformed_pseudobulk_counts(sc_pseudobulk_dir, possible_cell_types)
#multiome_14_pseudobulk_counts_log_transformed <- grab_transformed_pseudobulk_counts(multiome_pseudobulk_dir, possible_cell_types)

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

