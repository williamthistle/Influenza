library(data.table)
library(DESeq2)
library(MetaIntegrator)

# Set up 
base_dir <- "~/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
source(paste0(base_dir, "pseudobulk_analysis_helper.R"))
source(paste0(base_dir, "Data Compendium/Compendium_Functions.R"))
setup_bulk_analysis()
sample_metadata <- read.table(paste0(base_dir, "sample_metadata.tsv"), sep = "\t", header = TRUE)
cell_types <- c("CD4_Naive", "CD8_Naive", "CD4_Memory", "CD8_Memory", "cDC", "HSPC", "pDC", "Platelet", "Plasmablast", "Proliferating", "NK", "T_Naive", "CD14_Mono", "CD16_Mono", "MAIT")
single_cell_magical_dir <- "C:/Users/willi/OneDrive - Princeton University/Influenza Analysis/Single Cell RNA-Seq/MAGICAL Analyses/Placebo 6 Sample (Run by Aliza)/"
single_cell_pseudobulk_dir <- paste0(single_cell_magical_dir, "scRNA/pseudo_bulk/")
multiome_magical_dir <- "C:/Users/willi/OneDrive - Princeton University/Influenza Analysis/True Multiome/MAGICAL Analyses/14 Placebo Sample (Final)/"
multiome_pseudobulk_dir <- paste0(multiome_magical_dir, "scRNA_pseudobulk/")
set.seed(2000)

# Tables containing results for single cell and multiome RNA-seq processing
# Includes genes that passed pseudobulk filtering and genes that passed MAGICAL filtering
# Note that LR includes latent variable subject in differential expression analysis
single_cell_pseudobulk_gene_table <- read.table(paste0(single_cell_magical_dir, "D28_D1_MAGICAL_12_sample_sc_sc_genes.txt"), sep = "\t", header = TRUE)
single_cell_magical_gene_table <- read.table(paste0(single_cell_magical_dir, "D28_D1_MAGICAL_higher_fc_threshold_results.txt"), sep = "\t", header = TRUE)

multiome_pseudobulk_gene_table <- read.table(paste0(multiome_magical_dir, "D28_D1_MAGICAL_14_sample_multiome_sc_genes.txt"), sep = "\t", header = TRUE)
multiome_magical_gene_table <- read.table(paste0(multiome_magical_dir, "D28_D1_MAGICAL_14_sample_multiome.txt"), sep = "\t", header = TRUE)

multiome_pseudobulk_gene_LR_table <- read.table(paste0(multiome_magical_dir, "D28_D1_MAGICAL_14_sample_multiome_sc_genes_LR.txt"), sep = "\t", header = TRUE)
multiome_magical_gene_LR_table <- read.table(paste0(multiome_magical_dir, "D28_D1_MAGICAL_14_sample_multiome_LR.txt"), sep = "\t", header = TRUE)

# Grab gene lists from result tables and report number of genes
single_cell_pseudobulk_genes <- unique(single_cell_pseudobulk_gene_table$gene)
print(paste0("Number of genes that pass pseudobulk (scRNA): ", length(single_cell_pseudobulk_genes)))
single_cell_magical_genes <- unique(single_cell_magical_gene_table$Gene_symbol)
print(paste0("Number of genes that pass MAGICAL (scRNA): ", length(single_cell_magical_genes)))

multiome_pseudobulk_genes <- unique(multiome_pseudobulk_gene_table$gene)
print(paste0("Number of genes that pass pseudobulk (multiome): ", length(multiome_pseudobulk_genes)))
multiome_magical_genes <- unique(multiome_magical_gene_table$Gene_symbol)
print(paste0("Number of genes that pass MAGICAL (multiome): ", length(multiome_magical_genes)))

multiome_pseudobulk_genes_LR <- unique(multiome_pseudobulk_gene_LR_table$gene)
print(paste0("Number of genes that pass pseudobulk (multiome LR): ", length(multiome_pseudobulk_genes_LR)))
multiome_magical_genes_LR <- unique(multiome_magical_gene_LR_table$Gene_symbol)
print(paste0("Number of genes that pass MAGICAL (multiome LR): ", length(multiome_magical_genes_LR)))

# Create log transformed pseudobulk count tables
single_cell_pseudobulk_counts_log_transformed <- grab_transformed_pseudobulk_counts(single_cell_pseudobulk_dir, cell_types)
multiome_14_pseudobulk_counts_log_transformed <- grab_transformed_pseudobulk_counts(multiome_pseudobulk_dir, cell_types)

# Create MetaIntegrator objects using pseudobulk count tables
single_cell_pseudobulk_metaintegrator_obj <- create_metaintegrator_obj("single cell", single_cell_pseudobulk_counts_log_transformed)
multiome_pseudobulk_metaintegrator_obj <- create_metaintegrator_obj("multiome", multiome_14_pseudobulk_counts_log_transformed)

# Calculate individual AUCs for our gene lists on their respective pseudobulk data
sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(single_cell_pseudobulk_genes, single_cell_pseudobulk_metaintegrator_obj, "Single_Cell_Paired", "sc_pseudobulk"))
curated_sc_pseudobulk_gene_aucs <- sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc < 0.3 | sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]
curated_single_cell_pseudobulk_genes <- curated_sc_pseudobulk_gene_aucs$gene_name
sc_magical_gene_aucs <- curated_sc_pseudobulk_gene_aucs[curated_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
curated_single_cell_magical_genes <- sc_magical_gene_aucs$gene_name

multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(multiome_pseudobulk_genes, multiome_pseudobulk_metaintegrator_obj, "Multiome_Paired", "multiome_pseudobulk"))
curated_multiome_pseudobulk_gene_aucs <- multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc < 0.3 | multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc > 0.7,]
curated_multiome_pseudobulk_genes <- curated_multiome_pseudobulk_gene_aucs$gene_name
multiome_magical_gene_aucs <- curated_multiome_pseudobulk_gene_aucs[curated_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
curated_multiome_magical_genes <- multiome_magical_gene_aucs$gene_name

# Next, let's test our gene lists on the actual bulk RNA-seq data!
# Test on days 2, 5, 8, and 28 for both pseudobulk gene lists and MAGICAL gene lists
# Should I include the extra samples for days that have them?
# Create MetaIntegrator objects for all days (high and low)
# Remove 0 pcr sample
all_D2_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", placebo_counts, placebo_metadata, "2_D2", "2_D_minus_1")
all_D5_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", placebo_counts, placebo_metadata, "2_D5", "2_D_minus_1")
all_D8_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", placebo_counts, placebo_metadata, "2_D8", "2_D_minus_1")
all_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", placebo_counts, placebo_metadata, "2_D28", "2_D_minus_1")

high_D2_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D2", "2_D_minus_1")
high_D5_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D5", "2_D_minus_1")
high_D8_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D8", "2_D_minus_1")
high_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_placebo_counts, high_placebo_metadata, "2_D28", "2_D_minus_1")

low_D2_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_placebo_counts, low_placebo_metadata, "2_D2", "2_D_minus_1")
low_D5_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_placebo_counts, low_placebo_metadata, "2_D5", "2_D_minus_1")
low_D8_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_placebo_counts, low_placebo_metadata, "2_D8", "2_D_minus_1")
low_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_placebo_counts, low_placebo_metadata, "2_D28", "2_D_minus_1")
# Calculate gene AUCs for pseudobulk filtered genes (high and low) - single cell
all_bulk_D2_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, all_D2_bulk_metaintegrator_obj, "Single_Cell_Paired", "all_bulk_D2"))
all_bulk_D5_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, all_D5_bulk_metaintegrator_obj, "Single_Cell_Paired", "all_bulk_D5"))
all_bulk_D8_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, all_D8_bulk_metaintegrator_obj, "Single_Cell_Paired", "all_bulk_D8"))
all_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, all_D28_bulk_metaintegrator_obj, "Single_Cell_Paired", "all_bulk_D28"))

high_bulk_D2_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, high_D2_bulk_metaintegrator_obj, "Single_Cell_Paired", "high_bulk_D2"))
high_bulk_D5_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, high_D5_bulk_metaintegrator_obj, "Single_Cell_Paired", "high_bulk_D5"))
high_bulk_D8_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, high_D8_bulk_metaintegrator_obj, "Single_Cell_Paired", "high_bulk_D8"))
high_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, high_D28_bulk_metaintegrator_obj, "Single_Cell_Paired", "high_bulk_D28"))

low_bulk_D2_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, low_D2_bulk_metaintegrator_obj, "Single_Cell_Paired", "low_bulk_D2"))
low_bulk_D5_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, low_D5_bulk_metaintegrator_obj, "Single_Cell_Paired", "low_bulk_D5"))
low_bulk_D8_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, low_D8_bulk_metaintegrator_obj, "Single_Cell_Paired", "low_bulk_D8"))
low_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, low_D28_bulk_metaintegrator_obj, "Single_Cell_Paired", "low_bulk_D28"))
# Calculate gene AUCs for MAGICAL filtered genes (high and low) - single cell
all_bulk_D2_sc_magical_gene_aucs <- all_bulk_D2_sc_pseudobulk_gene_aucs[all_bulk_D2_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
all_bulk_D5_sc_magical_gene_aucs <- all_bulk_D5_sc_pseudobulk_gene_aucs[all_bulk_D5_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
all_bulk_D8_sc_magical_gene_aucs <- all_bulk_D8_sc_pseudobulk_gene_aucs[all_bulk_D8_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
all_bulk_D28_sc_magical_gene_aucs <- all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]

high_bulk_D2_sc_magical_gene_aucs <- high_bulk_D2_sc_pseudobulk_gene_aucs[high_bulk_D2_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
high_bulk_D5_sc_magical_gene_aucs <- high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
high_bulk_D8_sc_magical_gene_aucs <- high_bulk_D8_sc_pseudobulk_gene_aucs[high_bulk_D8_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
high_bulk_D28_sc_magical_gene_aucs <- high_bulk_D28_sc_pseudobulk_gene_aucs[high_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]

low_bulk_D2_sc_magical_gene_aucs <- low_bulk_D2_sc_pseudobulk_gene_aucs[low_bulk_D2_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
low_bulk_D5_sc_magical_gene_aucs <- low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
low_bulk_D8_sc_magical_gene_aucs <- low_bulk_D8_sc_pseudobulk_gene_aucs[low_bulk_D8_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
low_bulk_D28_sc_magical_gene_aucs <- low_bulk_D28_sc_pseudobulk_gene_aucs[low_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
# Calculate gene AUCs for pseudobulk filtered genes (high and low) - multiome
all_bulk_D2_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, all_D2_bulk_metaintegrator_obj, "Multiome_Paired", "all_bulk_D2"))
all_bulk_D5_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, all_D5_bulk_metaintegrator_obj, "Multiome_Paired", "all_bulk_D5"))
all_bulk_D8_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, all_D8_bulk_metaintegrator_obj, "Multiome_Paired", "all_bulk_D8"))
all_bulk_D28_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, all_D28_bulk_metaintegrator_obj, "Multiome_Paired", "all_bulk_D28"))

high_bulk_D2_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, high_D2_bulk_metaintegrator_obj, "Multiome_Paired", "high_bulk_D2"))
high_bulk_D5_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, high_D5_bulk_metaintegrator_obj, "Multiome_Paired", "high_bulk_D5"))
high_bulk_D8_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, high_D8_bulk_metaintegrator_obj, "Multiome_Paired", "high_bulk_D8"))
high_bulk_D28_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, high_D28_bulk_metaintegrator_obj, "Multiome_Paired", "high_bulk_D28"))

low_bulk_D2_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, low_D2_bulk_metaintegrator_obj, "Multiome_Paired", "low_bulk_D2"))
low_bulk_D5_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, low_D5_bulk_metaintegrator_obj, "Multiome_Paired", "low_bulk_D5"))
low_bulk_D8_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, low_D8_bulk_metaintegrator_obj, "Multiome_Paired", "low_bulk_D8"))
low_bulk_D28_multiome_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, low_D28_bulk_metaintegrator_obj, "Multiome_Paired", "low_bulk_D28"))
# Calculate gene AUCs for MAGICAL filtered genes (high and low) - multiome
all_bulk_D2_multiome_magical_gene_aucs <- all_bulk_D2_multiome_pseudobulk_gene_aucs[all_bulk_D2_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
all_bulk_D5_multiome_magical_gene_aucs <- all_bulk_D5_multiome_pseudobulk_gene_aucs[all_bulk_D5_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
all_bulk_D8_multiome_magical_gene_aucs <- all_bulk_D8_multiome_pseudobulk_gene_aucs[all_bulk_D8_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
all_bulk_D28_multiome_magical_gene_aucs <- all_bulk_D28_multiome_pseudobulk_gene_aucs[all_bulk_D28_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]

high_bulk_D2_multiome_magical_gene_aucs <- high_bulk_D2_multiome_pseudobulk_gene_aucs[high_bulk_D2_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
high_bulk_D5_multiome_magical_gene_aucs <- high_bulk_D5_multiome_pseudobulk_gene_aucs[high_bulk_D5_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
high_bulk_D8_multiome_magical_gene_aucs <- high_bulk_D8_multiome_pseudobulk_gene_aucs[high_bulk_D8_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
high_bulk_D28_multiome_magical_gene_aucs <- high_bulk_D28_multiome_pseudobulk_gene_aucs[high_bulk_D28_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]

low_bulk_D2_multiome_magical_gene_aucs <- low_bulk_D2_multiome_pseudobulk_gene_aucs[low_bulk_D2_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
low_bulk_D5_multiome_magical_gene_aucs <- low_bulk_D5_multiome_pseudobulk_gene_aucs[low_bulk_D5_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
low_bulk_D8_multiome_magical_gene_aucs <- low_bulk_D8_multiome_pseudobulk_gene_aucs[low_bulk_D8_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
low_bulk_D28_multiome_magical_gene_aucs <- low_bulk_D28_multiome_pseudobulk_gene_aucs[low_bulk_D28_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]

# Because of the way we calculate AUC, in our case an AUC of under 0.3 is equally valuable as AUC of over 0.7
# Above, we test genes one at a time as positive genes in our gene set signature
# If we get an AUC of under 0.3, that means that the same gene would score an AUC of over 0.7 as a negative gene
# in our gene set signature
# Why are there more genes that pass AUC in D28 vs days 2/5/8 for low? Must be because I'm training on low D28 data, right?
# Are the genes that pass AUC 0.7 (or AUC 0.3) differentially expressed in bulk data?
auc_df <- data.frame(Filtering_Assay = character(), Filtering_Method = character(), Discovery_Assay = character(), 
                     Discovery_Dataset = character(), Pos_Genes = integer(), Neg_Genes = integer(), Total_Passing_Genes = integer(), 
                     Total_Genes = integer(), Percentage_of_Passing_Genes = double(), stringsAsFactors = FALSE)
auc_names <- c("Filtering_Assay", "Filtering_Method", "Discovery_Assay", "Discovery_Dataset", "Pos_Genes", "Neg_Genes", "Total_Passing_Genes", "Total_Genes", "Percentage_of_Passing_Genes")
# Single cell, cell type pseudobulk filtering
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Single Cell", "Total Pseudobulk", curated_sc_pseudobulk_gene_aucs, "sc_pseudobulk_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "All Bulk D2", all_bulk_D2_sc_pseudobulk_gene_aucs, "all_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "All Bulk D5", all_bulk_D5_sc_pseudobulk_gene_aucs, "all_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "All Bulk D8", all_bulk_D8_sc_pseudobulk_gene_aucs, "all_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "All Bulk D28", all_bulk_D28_sc_pseudobulk_gene_aucs, "all_bulk_D28_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "High Bulk D2", high_bulk_D2_sc_pseudobulk_gene_aucs, "high_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "High Bulk D5", high_bulk_D5_sc_pseudobulk_gene_aucs, "high_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "High Bulk D8", high_bulk_D8_sc_pseudobulk_gene_aucs, "high_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "High Bulk D28", high_bulk_D28_sc_pseudobulk_gene_aucs, "high_bulk_D28_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "Low Bulk D2", low_bulk_D2_sc_pseudobulk_gene_aucs, "low_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "Low Bulk D5", low_bulk_D5_sc_pseudobulk_gene_aucs, "low_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "Low Bulk D8", low_bulk_D8_sc_pseudobulk_gene_aucs, "low_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk", "Bulk RNA-Seq", "Low Bulk D28", low_bulk_D28_sc_pseudobulk_gene_aucs, "low_bulk_D28_gene_auc")
# Single cell, MAGICAL filtering
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Single Cell", "Total Pseudobulk", curated_sc_magical_gene_aucs, "sc_pseudobulk_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "All Bulk D2", all_bulk_D2_sc_magical_gene_aucs, "all_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "All Bulk D5", all_bulk_D5_sc_magical_gene_aucs, "all_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "All Bulk D8", all_bulk_D8_sc_magical_gene_aucs, "all_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "All Bulk D28", all_bulk_D28_sc_magical_gene_aucs, "all_bulk_D28_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "High Bulk D2", high_bulk_D2_sc_magical_gene_aucs, "high_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "High Bulk D5", high_bulk_D5_sc_magical_gene_aucs, "high_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "High Bulk D8", high_bulk_D8_sc_magical_gene_aucs, "high_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "High Bulk D28", high_bulk_D28_sc_magical_gene_aucs, "high_bulk_D28_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "Low Bulk D2", low_bulk_D2_sc_magical_gene_aucs, "low_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "Low Bulk D5", low_bulk_D5_sc_magical_gene_aucs, "low_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "Low Bulk D8", low_bulk_D8_sc_magical_gene_aucs, "low_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Single Cell", "MAGICAL", "Bulk RNA-Seq", "Low Bulk D28", low_bulk_D28_sc_magical_gene_aucs, "low_bulk_D28_gene_auc")
# Multiome, cell type pseudobulk filtering
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Multiome", "Total Pseudobulk", curated_multiome_pseudobulk_gene_aucs, "multiome_pseudobulk_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "All Bulk D2", all_bulk_D2_multiome_pseudobulk_gene_aucs, "all_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "All Bulk D5", all_bulk_D5_multiome_pseudobulk_gene_aucs, "all_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "All Bulk D8", all_bulk_D8_multiome_pseudobulk_gene_aucs, "all_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "All Bulk D28", all_bulk_D28_multiome_pseudobulk_gene_aucs, "all_bulk_D28_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "High Bulk D2", high_bulk_D2_multiome_pseudobulk_gene_aucs, "high_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "High Bulk D5", high_bulk_D5_multiome_pseudobulk_gene_aucs, "high_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "High Bulk D8", high_bulk_D8_multiome_pseudobulk_gene_aucs, "high_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "High Bulk D28", high_bulk_D28_multiome_pseudobulk_gene_aucs, "high_bulk_D28_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "Low Bulk D2", low_bulk_D2_multiome_pseudobulk_gene_aucs, "low_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "Low Bulk D5", low_bulk_D5_multiome_pseudobulk_gene_aucs, "low_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "Low Bulk D8", low_bulk_D8_multiome_pseudobulk_gene_aucs, "low_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "Cell Type Pseudobulk", "Bulk RNA-Seq", "Low Bulk D28", low_bulk_D28_multiome_pseudobulk_gene_aucs, "low_bulk_D28_gene_auc")
# Multiome, MAGICAL filtering
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Multiome", "Total Pseudobulk", curated_multiome_magical_gene_aucs, "multiome_pseudobulk_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "All Bulk D2", all_bulk_D2_multiome_magical_gene_aucs, "all_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "All Bulk D5", all_bulk_D5_multiome_magical_gene_aucs, "all_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "All Bulk D8", all_bulk_D8_multiome_magical_gene_aucs, "all_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "All Bulk D28", all_bulk_D28_multiome_magical_gene_aucs, "all_bulk_D28_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "High Bulk D2", high_bulk_D2_multiome_magical_gene_aucs, "high_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "High Bulk D5", high_bulk_D5_multiome_magical_gene_aucs, "high_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "High Bulk D8", high_bulk_D8_multiome_magical_gene_aucs, "high_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "High Bulk D28", high_bulk_D28_multiome_magical_gene_aucs, "high_bulk_D28_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "Low Bulk D2", low_bulk_D2_multiome_magical_gene_aucs, "low_bulk_D2_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "Low Bulk D5", low_bulk_D5_multiome_magical_gene_aucs, "low_bulk_D5_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "Low Bulk D8", low_bulk_D8_multiome_magical_gene_aucs, "low_bulk_D8_gene_auc")
auc_df <- add_auc_row(auc_df, auc_names, "Multiome", "MAGICAL", "Bulk RNA-Seq", "Low Bulk D28", low_bulk_D28_multiome_magical_gene_aucs, "low_bulk_D28_gene_auc")

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

# Grab list of genes that have AUC > 0.7 in single cell pseudobulk
high_sc_pos_genes <- sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]$gene_name

# See which of these genes are significant in LRT data
high_sc_pos_genes_LRT_pass <- c()
for(gene in high_sc_pos_genes) {
  if(gene %in% rownames(high_placebo_period_2_LRT_analysis_results)) {
    high_sc_pos_genes_LRT_pass <- c(high_sc_pos_genes_LRT_pass, gene)
  }
}

high_placebo_period_2_LRT_analysis_betas <- coef(high_placebo_period_2_LRT_analysis)
high_placebo_period_2_LRT_analysis_betas <- high_placebo_period_2_LRT_analysis_betas[, -c(1:13)]
high_placebo_period_2_LRT_analysis_betas <- high_placebo_period_2_LRT_analysis_betas[rownames(high_placebo_period_2_LRT_analysis_betas) %in% high_sc_pos_genes_LRT_pass,]
high_placebo_period_2_LRT_analysis_thr <- 2.5
colnames(high_placebo_period_2_LRT_analysis_betas) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(high_placebo_period_2_LRT_analysis_betas, breaks=seq(from=-high_placebo_period_2_LRT_analysis_thr, to=high_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14)

# Grab list of genes that have AUC < 0.3 in single cell pseudobulk
high_sc_neg_genes <- sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc < 0.3,]$gene_name
# See which of these genes are significant in LRT data
high_sc_neg_genes_LRT_pass <- c()
for(gene in high_sc_neg_genes) {
  if(gene %in% rownames(high_placebo_period_2_LRT_analysis_results)) {
    high_sc_neg_genes_LRT_pass <- c(high_sc_neg_genes_LRT_pass, gene)
  }
}

high_placebo_period_2_LRT_analysis_betas_neg <- coef(high_placebo_period_2_LRT_analysis)
high_placebo_period_2_LRT_analysis_betas_neg <- high_placebo_period_2_LRT_analysis_betas_neg[, -c(1:13)]
high_placebo_period_2_LRT_analysis_betas_neg <- high_placebo_period_2_LRT_analysis_betas_neg[rownames(high_placebo_period_2_LRT_analysis_betas_neg) %in% high_sc_neg_genes_LRT_pass,]
colnames(high_placebo_period_2_LRT_analysis_betas_neg) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(high_placebo_period_2_LRT_analysis_betas_neg, breaks=seq(from=-1, to=1, length=101),
         cluster_col=FALSE, fontsize_col=14)




# Look at how data compendium papers use convalescent data

# Find intersecting genes between bulk and pseudobulk
# What about genes that have AUC < 0.3 or AUC > 0.7 for ALL time points? How many of those are there?
# Are those genes interesting or not interesting?
# Should we include subsets of times?
# See prediction in bulk datasets using pos and neg signatures, I guess?



# Genes that passed pseudobulk and all high bulk RNA-seq (AUC > 0.7) was IFI27 - interferon based gene. WHY IS THIS NOT WORKING ANYMORE?
high_sc_pos_genes <-sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]$gene_name
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D2_sc_pseudobulk_gene_aucs[high_bulk_D2_sc_pseudobulk_gene_aucs$high_bulk_D2_gene_auc > 0.7,]$gene_name)
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]$gene_name)
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D8_sc_pseudobulk_gene_aucs[high_bulk_D8_sc_pseudobulk_gene_aucs$high_bulk_D8_gene_auc > 0.7,]$gene_name)
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D28_sc_pseudobulk_gene_aucs[high_bulk_D28_sc_pseudobulk_gene_aucs$high_bulk_D28_gene_auc > 0.7,]$gene_name)
# Genes that passed pseudobulk and all high bulk RNA-seq (AUC < 0.3) were TUBA1A and BAG1. WHY IS BAG1 NOW POSITIVE GENE?
high_sc_neg_genes <- sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc < 0.3,]$gene_name
high_sc_neg_genes <- intersect(high_sc_neg_genes, high_bulk_D2_sc_pseudobulk_gene_aucs[high_bulk_D2_sc_pseudobulk_gene_aucs$high_bulk_D2_gene_auc < 0.3,]$gene_name)
high_sc_neg_genes <- intersect(high_sc_neg_genes, high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc < 0.3,]$gene_name)
high_sc_neg_genes <- intersect(high_sc_neg_genes, high_bulk_D8_sc_pseudobulk_gene_aucs[high_bulk_D8_sc_pseudobulk_gene_aucs$high_bulk_D8_gene_auc < 0.3,]$gene_name)
high_sc_neg_genes <- intersect(high_sc_neg_genes, high_bulk_D28_sc_pseudobulk_gene_aucs[high_bulk_D28_sc_pseudobulk_gene_aucs$high_bulk_D28_gene_auc < 0.3,]$gene_name)
# No genes passed pseudobulk and all low bulk RNA-seq (AUC > 0.7)
low_sc_pos_genes <- sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]$gene_name
low_sc_pos_genes <- intersect(low_sc_pos_genes, low_bulk_D2_sc_pseudobulk_gene_aucs[low_bulk_D2_sc_pseudobulk_gene_aucs$low_bulk_D2_gene_auc > 0.7,]$gene_name)
low_sc_pos_genes <- intersect(low_sc_pos_genes, low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc > 0.7,]$gene_name)
low_sc_pos_genes <- intersect(low_sc_pos_genes, low_bulk_D8_sc_pseudobulk_gene_aucs[low_bulk_D8_sc_pseudobulk_gene_aucs$low_bulk_D8_gene_auc > 0.7,]$gene_name)
low_sc_pos_genes <- intersect(low_sc_pos_genes, low_bulk_D28_sc_pseudobulk_gene_aucs[low_bulk_D28_sc_pseudobulk_gene_aucs$low_bulk_D28_gene_auc > 0.7,]$gene_name)
# No genes passed pseudobulk and all low bulk RNA-seq (AUC < 0.3)
low_sc_neg_genes <- sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc < 0.3,]$gene_name
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