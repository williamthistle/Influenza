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
# First, let's remove the 0 pcr sample from low because it's questionable
removed_low_value_aliquots <- rownames(placebo_metadata[placebo_metadata$subject_id == "f18c54d93cef4a4e",])
placebo_metadata <- placebo_metadata[!(placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
placebo_counts <- placebo_counts[,!(colnames(placebo_counts) %in% removed_low_value_aliquots)]
low_placebo_metadata <- low_placebo_metadata[!(low_placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
low_placebo_counts <- low_placebo_counts[,!(colnames(low_placebo_counts) %in% removed_low_value_aliquots)]
# Find bulk RNA-seq associated with the specific samples that we used for single-cell / multiome
# Should be interesting to compare performance in pseudobulk vs actual bulk
single_cell_aliquots <- c("91910a04711aa3dd", "3731a6247ae23831", "2232300b0e3a3d06", "76ea83ff9293871a", "5fdfdbaeb3c8eee8", "981520e7e138d460", "bb3d7b309cb5fc58", "8338411dc3e181e9", "da4fe21a89c8f7f4", "41d248a6ec3b87e2", "e3e01c75894ef461", "4534496c580cb408") # 12 samples - 6 paired
single_cell_subjects <- as.character(unique(all_metadata[all_metadata$aliquot_id %in% single_cell_aliquots,]$subject_id))
multiome_aliquots <- c("717579a2ae2fb6c2", "dde63f8ca98af665", "a464019298ae6682", "d554be0e36e4d789", "e43db0f72b9c2e31", "b82bb7c75d47dac1", "6f609a68dca1261f", "9c6ec1b704700c7d", "7b54cfac7e67b0fa", "575d74707585856a", "c1eb160d7bd1f29f", "8832fff8247b18b9", "abf6d19ee03be1e8", "216bb226181591dd") # 14 samples - 7 paired
multiome_subjects <- as.character(unique(all_metadata[all_metadata$aliquot_id %in% multiome_aliquots,]$subject_id))

single_cell_placebo_metadata <- placebo_metadata[(placebo_metadata$subject_id %in% single_cell_subjects),]
single_cell_placebo_counts <- placebo_counts[,(colnames(placebo_counts) %in% rownames(single_cell_placebo_metadata))]
multiome_placebo_metadata <- placebo_metadata[(placebo_metadata$subject_id %in% multiome_subjects),]
multiome_placebo_counts <- placebo_counts[,(colnames(placebo_counts) %in% rownames(multiome_placebo_metadata))]

high_single_cell_placebo_metadata <- high_placebo_metadata[(high_placebo_metadata$subject_id %in% single_cell_subjects),]
high_single_cell_placebo_counts <- high_placebo_counts[,(colnames(high_placebo_counts) %in% rownames(high_single_cell_placebo_metadata))]
high_multiome_placebo_metadata <- high_placebo_metadata[(high_placebo_metadata$subject_id %in% multiome_subjects),]
high_multiome_placebo_counts <- high_placebo_counts[,(colnames(high_placebo_counts) %in% rownames(high_multiome_placebo_metadata))]

low_single_cell_placebo_metadata <- low_placebo_metadata[(low_placebo_metadata$subject_id %in% single_cell_subjects),]
low_single_cell_placebo_counts <- low_placebo_counts[,(colnames(low_placebo_counts) %in% rownames(low_single_cell_placebo_metadata))]
low_multiome_placebo_metadata <- low_placebo_metadata[(low_placebo_metadata$subject_id %in% multiome_subjects),]
low_multiome_placebo_counts <- low_placebo_counts[,(colnames(low_placebo_counts) %in% rownames(low_multiome_placebo_metadata))]

# Grab DEGs for 2_D28 vs 2_D_minus_1 for bulk (Wald test)
placebo_period_2_D28_vs_D_minus_1_bulk_results <- run_deseq_bulk_analysis_time_series("placebo", placebo_counts,placebo_metadata,
                                                                                      "2_D28", "2_D_minus_1", "~/")

# Create MetaIntegrator objects for D28 specifically for samples that we processed using scRNA-seq or multiome
# We can compare to pseudobulk AUCs for same samples
sc_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", single_cell_placebo_counts, single_cell_placebo_metadata, "2_D28", "2_D_minus_1")
multiome_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", multiome_placebo_counts, multiome_placebo_metadata, "2_D28", "2_D_minus_1")

high_sc_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_single_cell_placebo_counts, high_single_cell_placebo_metadata, "2_D28", "2_D_minus_1")
high_multiome_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_multiome_placebo_counts, high_multiome_placebo_metadata, "2_D28", "2_D_minus_1")

low_sc_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_single_cell_placebo_counts, low_single_cell_placebo_metadata, "2_D28", "2_D_minus_1")
low_multiome_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_multiome_placebo_counts, low_multiome_placebo_metadata, "2_D28", "2_D_minus_1")

# Create MetaIntegrator objects for all days (high and low)
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

# CURRENT PLAN: Use all_D28 auc of >0.7 as filtering on original pseudobulk list (don't include pseudobulk AUC > 0.7)
# alternatively, we could do both. But why would pseudobulk auc > 0.7 (with less data) be useful if we have all_D28 auc of >0.7 (includes that data + more)?
# We also look at MAGICAL gene subset.
# Then, we can look at each of these genes, see whether it's differentially expressed in our bulk data, and look at the trend - check out HumanBase, etc.
# Also, we can look at expression in high and low
# We can also do exactly the same for multiome
# Our single cell data gives us single cell granularity on the genes being expressed. This also applies to the ATAC-seq data and our peaks! mintchip can help us verify our choices,
# but the ATAC-seq gives us the cell type!
# I just need a few good stories. Look at specific genes or families of genes and see why they're interesting!
# What is the behavior of these genes throughout infection and why are they so persistent at D28?
# List cell types for genes when reporting

# Calculate gene AUCs for pseudobulk filtered genes (high and low) - single cell (D28, samples we processed)
sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, sc_D28_bulk_metaintegrator_obj, "Single_Cell_Paired", "sc_bulk_D28"))
high_sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, high_sc_D28_bulk_metaintegrator_obj, "Single_Cell_Paired", "high_sc_bulk_D28"))
low_sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, low_sc_D28_bulk_metaintegrator_obj, "Single_Cell_Paired", "low_sc_bulk_D28"))

# Calculate gene AUCs for MAGICAL genes (high and low) - single cell (D28, samples we processed)
sc_bulk_D28_magical_gene_aucs <- sc_bulk_D28_sc_pseudobulk_gene_aucs[sc_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
high_sc_bulk_D28_magical_gene_aucs <- high_sc_bulk_D28_sc_pseudobulk_gene_aucs[high_sc_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]
low_sc_bulk_D28_magical_gene_aucs <- low_sc_bulk_D28_sc_pseudobulk_gene_aucs[low_sc_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]

# Calculate gene AUCs for pseudobulk filtered genes (high and low) - multiome (D28, samples we processed)
multiome_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, multiome_D28_bulk_metaintegrator_obj, "Multiome_Paired", "multiome_bulk_D28"))
high_multiome_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, high_multiome_D28_bulk_metaintegrator_obj, "Multiome_Paired", "high_multiome_bulk_D28"))
low_multiome_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_multiome_pseudobulk_genes, low_multiome_D28_bulk_metaintegrator_obj, "Multiome_Paired", "low_multiome_bulk_D28"))

# Calculate gene AUCs for MAGICAL genes (high and low) - multiome (D28, samples we processed)
multiome_bulk_D28_magical_gene_aucs <- multiome_bulk_D28_sc_pseudobulk_gene_aucs[multiome_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
high_multiome_bulk_D28_magical_gene_aucs <- high_multiome_bulk_D28_sc_pseudobulk_gene_aucs[high_multiome_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]
low_multiome_bulk_D28_magical_gene_aucs <- low_multiome_bulk_D28_sc_pseudobulk_gene_aucs[low_multiome_bulk_D28_sc_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]

# Calculate gene AUCs for pseudobulk filtered genes (high and low) - single cell
all_bulk_D2_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, all_D2_bulk_metaintegrator_obj, "Single_Cell_Paired", "all_bulk_D2"))
all_bulk_D5_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, all_D5_bulk_metaintegrator_obj, "Single_Cell_Paired", "all_bulk_D5"))
all_bulk_D8_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(curated_single_cell_pseudobulk_genes, all_D8_bulk_metaintegrator_obj, "Single_Cell_Paired", "all_bulk_D8"))
# NOTE - using list of all pseudobulk pass genes (not curated)
all_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(single_cell_pseudobulk_genes, all_D28_bulk_metaintegrator_obj, "Single_Cell_Paired", "all_bulk_D28"))

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
high_placebo_period_2_LRT_analysis_thr <- 1.5
colnames(high_placebo_period_2_LRT_analysis_betas) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(high_placebo_period_2_LRT_analysis_betas, breaks=seq(from=-high_placebo_period_2_LRT_analysis_thr, to=high_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14)

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
         cluster_col=FALSE, fontsize_col=14)






# Look at how data compendium papers use convalescent data

# Find intersecting genes between bulk and pseudobulk
# What about genes that have AUC < 0.3 or AUC > 0.7 for ALL time points? How many of those are there?
# Are those genes interesting or not interesting?
# Should we include subsets of times?
# See prediction in bulk datasets using pos and neg signatures, I guess?
# See how many other genes have same trend as BAG1 (pos AUC in pseudo and fold change is generally negative or vice versa)



# Genes that passed D28 bulk and all bulk RNA-seq (AUC > 0.7) was nothing
all_pass_sc_pos_genes <- all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$all_bulk_D28_gene_auc > 0.7,]$gene_name
all_pass_sc_pos_genes <- intersect(all_pass_sc_pos_genes, all_bulk_D2_sc_pseudobulk_gene_aucs[all_bulk_D2_sc_pseudobulk_gene_aucs$all_bulk_D2_gene_auc > 0.7,]$gene_name)
all_pass_sc_pos_genes <- intersect(all_pass_sc_pos_genes, all_bulk_D5_sc_pseudobulk_gene_aucs[all_bulk_D5_sc_pseudobulk_gene_aucs$all_bulk_D5_gene_auc > 0.7,]$gene_name)
all_pass_sc_pos_genes <- intersect(all_pass_sc_pos_genes, all_bulk_D8_sc_pseudobulk_gene_aucs[all_bulk_D8_sc_pseudobulk_gene_aucs$all_bulk_D8_gene_auc > 0.7,]$gene_name)

# Genes that passed D28 bulk and all bulk RNA-seq (AUC < 0.3) were "OST4"    "C9orf78" "TMA7"    "UQCR11"  "BAG1" 
all_pass_sc_neg_genes <- all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$all_bulk_D28_gene_auc < 0.3,]$gene_name
all_pass_sc_neg_genes <- intersect(all_pass_sc_neg_genes, all_bulk_D2_sc_pseudobulk_gene_aucs[all_bulk_D2_sc_pseudobulk_gene_aucs$all_bulk_D2_gene_auc < 0.3,]$gene_name)
all_pass_sc_neg_genes <- intersect(all_pass_sc_neg_genes, all_bulk_D5_sc_pseudobulk_gene_aucs[all_bulk_D5_sc_pseudobulk_gene_aucs$all_bulk_D5_gene_auc < 0.3,]$gene_name)
all_pass_sc_neg_genes <- intersect(all_pass_sc_neg_genes, all_bulk_D8_sc_pseudobulk_gene_aucs[all_bulk_D8_sc_pseudobulk_gene_aucs$all_bulk_D8_gene_auc < 0.3,]$gene_name)

# Genes that passed D28 bulk and all high bulk RNA-seq (AUC > 0.7) was IL10RA - interferon based gene. 
high_sc_pos_genes <- all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$all_bulk_D28_gene_auc > 0.7,]$gene_name
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D2_sc_pseudobulk_gene_aucs[high_bulk_D2_sc_pseudobulk_gene_aucs$high_bulk_D2_gene_auc > 0.7,]$gene_name)
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]$gene_name)
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D8_sc_pseudobulk_gene_aucs[high_bulk_D8_sc_pseudobulk_gene_aucs$high_bulk_D8_gene_auc > 0.7,]$gene_name)
high_sc_pos_genes <- intersect(high_sc_pos_genes, high_bulk_D28_sc_pseudobulk_gene_aucs[high_bulk_D28_sc_pseudobulk_gene_aucs$high_bulk_D28_gene_auc > 0.7,]$gene_name)
# Genes that passed D28 bulk and all high bulk RNA-seq (AUC < 0.3) were TUBA1A and BAG1. WHY IS BAG1 NOW POSITIVE GENE?
high_sc_neg_genes <-  all_bulk_D28_sc_pseudobulk_gene_aucs[all_bulk_D28_sc_pseudobulk_gene_aucs$all_bulk_D28_gene_auc < 0.3,]$gene_name
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



curated_sc_pseudobulk_gene_aucs
sc_bulk_D28_sc_pseudobulk_gene_aucs