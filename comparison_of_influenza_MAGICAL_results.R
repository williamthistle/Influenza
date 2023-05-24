library(data.table)
library(DESeq2)
library(MetaIntegrator)

base_dir <- "~/GitHub/Influenza/"
single_cell_magical_dir <- "C:/Users/willi/OneDrive - Princeton University/Influenza Analysis/Single Cell RNA-Seq/MAGICAL Analyses/Placebo 6 Sample (Run by Aliza)/"
single_cell_pseudobulk_dir <- paste0(single_cell_magical_dir, "scRNA/pseudo_bulk/")
multiome_magical_dir <- "C:/Users/willi/OneDrive - Princeton University/Influenza Analysis/True Multiome/MAGICAL Analyses/14 Placebo Sample (Final)/"
multiome_pseudobulk_dir <- paste0(multiome_magical_dir, "scRNA_pseudobulk/")

# Tables
single_cell_pseudobulk_gene_table <- read.table(paste0(single_cell_magical_dir, "D28_D1_MAGICAL_12_sample_sc_sc_genes.txt"), sep = "\t", header = TRUE)
single_cell_magical_gene_table <- read.table(paste0(single_cell_magical_dir, "D28_D1_MAGICAL_higher_fc_threshold_results.txt"), sep = "\t", header = TRUE)

multiome_pseudobulk_gene_table <- read.table(paste0(multiome_magical_dir, "D28_D1_MAGICAL_14_sample_multiome_sc_genes.txt"), sep = "\t", header = TRUE)
multiome_magical_gene_table <- read.table(paste0(multiome_magical_dir, "D28_D1_MAGICAL_14_sample_multiome.txt"), sep = "\t", header = TRUE)

multiome_pseudobulk_gene_LR_table <- read.table(paste0(multiome_magical_dir, "D28_D1_MAGICAL_14_sample_multiome_sc_genes_LR.txt"), sep = "\t", header = TRUE)
multiome_magical_gene_LR_table <- read.table(paste0(multiome_magical_dir, "D28_D1_MAGICAL_14_sample_multiome_LR.txt"), sep = "\t", header = TRUE)

# Gene lists
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

# Set up 
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
source(paste0(base_dir, "pseudobulk_analysis_helper.R"))
source(paste0(base_dir, "Data Compendium/Compendium_Functions.R"))
setup_bulk_analysis()
sample_metadata <- read.table(paste0(base_dir, "sample_metadata.tsv"), sep = "\t", header = TRUE)
cell_types <- c("CD4_Naive", "CD8_Naive", "CD4_Memory", "CD8_Memory", "cDC", "HSPC", "pDC", "Platelet", "Plasmablast", "Proliferating", "NK", "T_Naive", "CD14_Mono", "CD16_Mono", "MAIT")

# Create log transformed pseudobulk count tables
single_cell_pseudobulk_counts_log_transformed <- grab_transformed_pseudobulk_counts(single_cell_pseudobulk_dir, cell_types)
multiome_14_pseudobulk_counts_log_transformed <- grab_transformed_pseudobulk_counts(multiome_pseudobulk_dir, cell_types)

# Create MetaIntegrator objects using pseudobulk count tables
single_cell_pseudobulk_metaintegrator_obj <- create_metaintegrator_obj("single cell", single_cell_pseudobulk_counts_log_transformed)
multiome_pseudobulk_metaintegrator_obj <- create_metaintegrator_obj("multiome", multiome_14_pseudobulk_counts_log_transformed)

# Calculate individual AUCs for our gene lists on their respective pseudobulk data
sc_pseudobulk_gene_aucs <- test_individual_genes_on_datasets(single_cell_pseudobulk_genes, single_cell_pseudobulk_metaintegrator_obj, "Single_Cell_Paired", "sc_pseudobulk")
sc_magical_gene_aucs <- sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]

multiome_pseudobulk_gene_aucs <- test_individual_genes_on_datasets(multiome_pseudobulk_genes, multiome_pseudobulk_metaintegrator_obj, "Multiome_Paired", "multiome_pseudobulk")
multiome_magical_gene_aucs <- multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]

# Next, let's test our gene lists on the actual bulk RNA-seq data!
# First, we will use day 5 because it should have the most active infection
# We could also test on days 2, 8 and day 28 (most similar to our pseudobulk data)
# Should I include the extra samples for days that have them?
high_metadata_subset <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D5" | high_placebo_metadata$time_point == "2_D_minus_1",]
low_metadata_subset <- low_placebo_metadata[low_placebo_metadata$time_point == "2_D5" | low_placebo_metadata$time_point == "2_D_minus_1",]

high_counts_subset <- high_placebo_counts[rownames(high_metadata_subset)]
low_counts_subset <- low_placebo_counts[rownames(low_metadata_subset)]

high_counts_subset <- varianceStabilizingTransformation(as.matrix(high_counts_subset))
low_counts_subset <- varianceStabilizingTransformation(as.matrix(low_counts_subset))

high_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", high_counts_subset, high_metadata_subset)
low_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", low_counts_subset, low_metadata_subset)

high_bulk_D5_sc_pseudobulk_gene_aucs <- test_individual_genes_on_datasets(single_cell_pseudobulk_genes, high_bulk_metaintegrator_obj, "Single_Cell_Paired", "high_bulk_D5")
high_bulk_D5_sc_magical_gene_aucs <- high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]

low_bulk_D5_sc_pseudobulk_gene_aucs <- test_individual_genes_on_datasets(single_cell_pseudobulk_genes, low_bulk_metaintegrator_obj, "Single_Cell_Paired", "low_bulk_D5")
low_bulk_D5_sc_magical_gene_aucs <- low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$gene_name %in% single_cell_magical_genes,]

high_bulk_D5_multiome_pseudobulk_gene_aucs <- test_individual_genes_on_datasets(multiome_pseudobulk_genes, high_bulk_metaintegrator_obj, "Multiome_Paired", "high_bulk_D5")
high_bulk_D5_multiome_magical_gene_aucs <- high_bulk_D5_multiome_pseudobulk_gene_aucs[high_bulk_D5_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]

low_bulk_D5_multiome_pseudobulk_gene_aucs <- test_individual_genes_on_datasets(multiome_pseudobulk_genes, low_bulk_metaintegrator_obj, "Multiome_Paired", "low_bulk_D5")
low_bulk_D5_multiome_magical_gene_aucs <- low_bulk_D5_multiome_pseudobulk_gene_aucs[low_bulk_D5_multiome_pseudobulk_gene_aucs$gene_name %in% multiome_magical_genes,]

# What percent of genes have AUC > 0.7 or AUC < 0.3?
# Because of the way we calculate AUC, in our case an AUC of under 0.3 is equally valuable as AUC of over 0.7
# Above, we test genes one at a time as positive genes in our gene set signature
# If we get an AUC of under 0.3, that means that the same gene would score an AUC of over 0.7 as a negative gene
# in our gene set signature
auc_df <- data.frame(Filtering_Assay = character(), Filtering_Method = character(), Discovery_Assay = character(), 
                     Discovery_Dataset = character(), Pos_Genes = integer(), Neg_Genes = integer(), Total_Passing_Genes = integer(), 
                     Total_Genes = integer(), Percentage_of_Passing_Genes = double(), stringsAsFactors = FALSE)
auc_names <- c("Filtering_Assay", "Filtering_Method", "Discovery_Assay", "Discovery_Dataset", "Pos_Genes", "Neg_Genes", "Total_Passing_Genes", "Total_Genes", "Percentage_of_Passing_Genes")
# Filtering: Single cell - cell type pseudobulk
# Discovery dataset: Single cell - total pseudobulk
current_row <- data.frame("Single Cell", "Cell Type Pseudobulk", "Single Cell", "Total Pseudobulk", nrow(sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]),
                          nrow(sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc < 0.3,]), 
                          nrow(sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]) + nrow(sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc < 0.3,]),
                          nrow(sc_pseudobulk_gene_aucs), (nrow(sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]) + 
                                                            nrow(sc_pseudobulk_gene_aucs[sc_pseudobulk_gene_aucs$sc_pseudobulk_gene_auc < 0.3,])) / nrow(sc_pseudobulk_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Single Cell - Cell Type Pseudobulk
# Discovery dataset: Bulk RNA-seq - High Bulk D5
current_row <- data.frame("Single Cell", "Cell Type Pseudobulk", "Bulk RNA-seq", "High Bulk D5", nrow(high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]),
                          nrow(high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc < 0.3,]), 
                          nrow(high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]) + nrow(high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc < 0.3,]),
                          nrow(high_bulk_D5_sc_pseudobulk_gene_aucs), (nrow(high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]) + 
                                                                         nrow(high_bulk_D5_sc_pseudobulk_gene_aucs[high_bulk_D5_sc_pseudobulk_gene_aucs$high_bulk_D5_gene_auc < 0.3,])) / nrow(high_bulk_D5_sc_pseudobulk_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Single Cell - Cell Type Pseudobulk
# Discovery dataset: Bulk RNA-seq - Low Bulk D5
current_row <- data.frame("Single Cell", "Cell Type Pseudobulk", "Bulk RNA-seq", "Low Bulk D5", nrow(low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc > 0.7,]),
                          nrow(low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc < 0.3,]), 
                          nrow(low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc > 0.7,]) + nrow(low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc < 0.3,]),
                          nrow(low_bulk_D5_sc_pseudobulk_gene_aucs), (nrow(low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc > 0.7,]) + 
                                                                        nrow(low_bulk_D5_sc_pseudobulk_gene_aucs[low_bulk_D5_sc_pseudobulk_gene_aucs$low_bulk_D5_gene_auc < 0.3,])) / nrow(low_bulk_D5_sc_pseudobulk_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Single cell - MAGICAL
# Discovery dataset: Single cell - total pseudobulk
current_row <- data.frame("Single Cell", "MAGICAL", "Single Cell", "Total Pseudobulk", nrow(sc_magical_gene_aucs[sc_magical_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]),
                          nrow(sc_magical_gene_aucs[sc_magical_gene_aucs$sc_pseudobulk_gene_auc < 0.3,]), 
                          nrow(sc_magical_gene_aucs[sc_magical_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]) + nrow(sc_magical_gene_aucs[sc_magical_gene_aucs$sc_pseudobulk_gene_auc < 0.3,]),
                          nrow(sc_magical_gene_aucs), (nrow(sc_magical_gene_aucs[sc_magical_gene_aucs$sc_pseudobulk_gene_auc > 0.7,]) + 
                                                         nrow(sc_magical_gene_aucs[sc_magical_gene_aucs$sc_pseudobulk_gene_auc < 0.3,])) / nrow(sc_magical_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Single Cell - MAGICAL
# Discovery dataset: Bulk RNA-seq - High Bulk D5
current_row <- data.frame("Single Cell", "MAGICAL", "Bulk RNA-seq", "High Bulk D5", nrow(high_bulk_D5_sc_magical_gene_aucs[high_bulk_D5_sc_magical_gene_aucs$high_bulk_D5_gene_auc > 0.7,]),
                          nrow(high_bulk_D5_sc_magical_gene_aucs[high_bulk_D5_sc_magical_gene_aucs$high_bulk_D5_gene_auc < 0.3,]), 
                          nrow(high_bulk_D5_sc_magical_gene_aucs[high_bulk_D5_sc_magical_gene_aucs$high_bulk_D5_gene_auc > 0.7,]) + nrow(high_bulk_D5_sc_magical_gene_aucs[high_bulk_D5_sc_magical_gene_aucs$high_bulk_D5_gene_auc < 0.3,]),
                          nrow(high_bulk_D5_sc_magical_gene_aucs), (nrow(high_bulk_D5_sc_magical_gene_aucs[high_bulk_D5_sc_magical_gene_aucs$high_bulk_D5_gene_auc > 0.7,]) + 
                                                                      nrow(high_bulk_D5_sc_magical_gene_aucs[high_bulk_D5_sc_magical_gene_aucs$high_bulk_D5_gene_auc < 0.3,])) / nrow(high_bulk_D5_sc_magical_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Single Cell - MAGICAL
# Discovery dataset: Bulk RNA-seq - Low Bulk D5
current_row <- data.frame("Single Cell", "MAGICAL", "Bulk RNA-seq", "Low Bulk D5", nrow(low_bulk_D5_sc_magical_gene_aucs[low_bulk_D5_sc_magical_gene_aucs$low_bulk_D5_gene_auc > 0.7,]),
                          nrow(low_bulk_D5_sc_magical_gene_aucs[low_bulk_D5_sc_magical_gene_aucs$low_bulk_D5_gene_auc < 0.3,]), 
                          nrow(low_bulk_D5_sc_magical_gene_aucs[low_bulk_D5_sc_magical_gene_aucs$low_bulk_D5_gene_auc > 0.7,]) + nrow(low_bulk_D5_sc_magical_gene_aucs[low_bulk_D5_sc_magical_gene_aucs$low_bulk_D5_gene_auc < 0.3,]),
                          nrow(low_bulk_D5_sc_magical_gene_aucs), (nrow(low_bulk_D5_sc_magical_gene_aucs[low_bulk_D5_sc_magical_gene_aucs$low_bulk_D5_gene_auc > 0.7,]) + 
                                                                     nrow(low_bulk_D5_sc_magical_gene_aucs[low_bulk_D5_sc_magical_gene_aucs$low_bulk_D5_gene_auc < 0.3,])) / nrow(low_bulk_D5_sc_magical_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Multiome - cell type pseudobulk
# Discovery dataset: Multiome - total pseudobulk
current_row <- data.frame("Multiome", "Cell Type Pseudobulk", "Multiome", "Total Pseudobulk", nrow(multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc > 0.7,]),
                          nrow(multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc < 0.3,]), 
                          nrow(multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc > 0.7,]) + nrow(multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc < 0.3,]),
                          nrow(multiome_pseudobulk_gene_aucs), (nrow(multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc > 0.7,]) + 
                                                                  nrow(multiome_pseudobulk_gene_aucs[multiome_pseudobulk_gene_aucs$multiome_pseudobulk_gene_auc < 0.3,])) / nrow(multiome_pseudobulk_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Multiome - Cell Type Pseudobulk
# Discovery dataset: Bulk RNA-seq - High Bulk D5
current_row <- data.frame("Multiome", "Cell Type Pseudobulk", "Bulk RNA-seq", "High Bulk D5", nrow(high_bulk_D5_multiome_pseudobulk_gene_aucs[high_bulk_D5_multiome_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]),
                          nrow(high_bulk_D5_multiome_pseudobulk_gene_aucs[high_bulk_D5_multiome_pseudobulk_gene_aucs$high_bulk_D5_gene_auc < 0.3,]), 
                          nrow(high_bulk_D5_multiome_pseudobulk_gene_aucs[high_bulk_D5_multiome_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]) + nrow(high_bulk_D5_multiome_pseudobulk_gene_aucs[high_bulk_D5_multiome_pseudobulk_gene_aucs$high_bulk_D5_gene_auc < 0.3,]),
                          nrow(high_bulk_D5_multiome_pseudobulk_gene_aucs), (nrow(high_bulk_D5_multiome_pseudobulk_gene_aucs[high_bulk_D5_multiome_pseudobulk_gene_aucs$high_bulk_D5_gene_auc > 0.7,]) + 
                                                                               nrow(high_bulk_D5_multiome_pseudobulk_gene_aucs[high_bulk_D5_multiome_pseudobulk_gene_aucs$high_bulk_D5_gene_auc < 0.3,])) / nrow(high_bulk_D5_multiome_pseudobulk_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Multiome - Cell Type Pseudobulk
# Discovery dataset: Bulk RNA-seq - Low Bulk D5
current_row <- data.frame("Multiome", "Cell Type Pseudobulk", "Bulk RNA-seq", "Low Bulk D5", nrow(low_bulk_D5_multiome_pseudobulk_gene_aucs[low_bulk_D5_multiome_pseudobulk_gene_aucs$low_bulk_D5_gene_auc > 0.7,]),
                          nrow(low_bulk_D5_multiome_pseudobulk_gene_aucs[low_bulk_D5_multiome_pseudobulk_gene_aucs$low_bulk_D5_gene_auc < 0.3,]), 
                          nrow(low_bulk_D5_multiome_pseudobulk_gene_aucs[low_bulk_D5_multiome_pseudobulk_gene_aucs$low_bulk_D5_gene_auc > 0.7,]) + nrow(low_bulk_D5_multiome_pseudobulk_gene_aucs[low_bulk_D5_multiome_pseudobulk_gene_aucs$low_bulk_D5_gene_auc < 0.3,]),
                          nrow(low_bulk_D5_multiome_pseudobulk_gene_aucs), (nrow(low_bulk_D5_multiome_pseudobulk_gene_aucs[low_bulk_D5_multiome_pseudobulk_gene_aucs$low_bulk_D5_gene_auc > 0.7,]) + 
                                                                              nrow(low_bulk_D5_multiome_pseudobulk_gene_aucs[low_bulk_D5_multiome_pseudobulk_gene_aucs$low_bulk_D5_gene_auc < 0.3,])) / nrow(low_bulk_D5_multiome_pseudobulk_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Multiome - MAGICAL
# Discovery dataset: Multiome - total pseudobulk
current_row <- data.frame("Multiome", "MAGICAL", "Multiome", "Total Pseudobulk", nrow(multiome_magical_gene_aucs[multiome_magical_gene_aucs$multiome_pseudobulk_gene_auc > 0.7,]),
                          nrow(multiome_magical_gene_aucs[multiome_magical_gene_aucs$multiome_pseudobulk_gene_auc < 0.3,]), 
                          nrow(multiome_magical_gene_aucs[multiome_magical_gene_aucs$multiome_pseudobulk_gene_auc > 0.7,]) + nrow(multiome_magical_gene_aucs[multiome_magical_gene_aucs$multiome_pseudobulk_gene_auc < 0.3,]),
                          nrow(multiome_magical_gene_aucs), (nrow(multiome_magical_gene_aucs[multiome_magical_gene_aucs$multiome_pseudobulk_gene_auc > 0.7,]) + 
                                                               nrow(multiome_magical_gene_aucs[multiome_magical_gene_aucs$multiome_pseudobulk_gene_auc < 0.3,])) / nrow(multiome_magical_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Multiome - MAGICAL
# Discovery dataset: Bulk RNA-seq - High Bulk D5
current_row <- data.frame("Multiome", "MAGICAL", "Bulk RNA-seq", "High Bulk D5", nrow(high_bulk_D5_multiome_magical_gene_aucs[high_bulk_D5_multiome_magical_gene_aucs$high_bulk_D5_gene_auc > 0.7,]),
                          nrow(high_bulk_D5_multiome_magical_gene_aucs[high_bulk_D5_multiome_magical_gene_aucs$high_bulk_D5_gene_auc < 0.3,]), 
                          nrow(high_bulk_D5_multiome_magical_gene_aucs[high_bulk_D5_multiome_magical_gene_aucs$high_bulk_D5_gene_auc > 0.7,]) + nrow(high_bulk_D5_multiome_magical_gene_aucs[high_bulk_D5_multiome_magical_gene_aucs$high_bulk_D5_gene_auc < 0.3,]),
                          nrow(high_bulk_D5_multiome_magical_gene_aucs), (nrow(high_bulk_D5_multiome_magical_gene_aucs[high_bulk_D5_multiome_magical_gene_aucs$high_bulk_D5_gene_auc > 0.7,]) + 
                                                                               nrow(high_bulk_D5_multiome_magical_gene_aucs[high_bulk_D5_multiome_magical_gene_aucs$high_bulk_D5_gene_auc < 0.3,])) / nrow(high_bulk_D5_multiome_magical_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)
# Filtering: Multiome - MAGICAL
# Discovery dataset: Bulk RNA-seq - Low Bulk D5
current_row <- data.frame("Multiome", "MAGICAL", "Bulk RNA-seq", "Low Bulk D5", nrow(low_bulk_D5_multiome_magical_gene_aucs[low_bulk_D5_multiome_magical_gene_aucs$low_bulk_D5_gene_auc > 0.7,]),
                          nrow(low_bulk_D5_multiome_magical_gene_aucs[low_bulk_D5_multiome_magical_gene_aucs$low_bulk_D5_gene_auc < 0.3,]), 
                          nrow(low_bulk_D5_multiome_magical_gene_aucs[low_bulk_D5_multiome_magical_gene_aucs$low_bulk_D5_gene_auc > 0.7,]) + nrow(low_bulk_D5_multiome_magical_gene_aucs[low_bulk_D5_multiome_magical_gene_aucs$low_bulk_D5_gene_auc < 0.3,]),
                          nrow(low_bulk_D5_multiome_magical_gene_aucs), (nrow(low_bulk_D5_multiome_magical_gene_aucs[low_bulk_D5_multiome_magical_gene_aucs$low_bulk_D5_gene_auc > 0.7,]) + 
                                                                              nrow(low_bulk_D5_multiome_magical_gene_aucs[low_bulk_D5_multiome_magical_gene_aucs$low_bulk_D5_gene_auc < 0.3,])) / nrow(low_bulk_D5_multiome_magical_gene_aucs))
names(current_row) <- auc_names
auc_df <- rbind(auc_df, current_row)








# Find intersect between bulk and pseudobulk
# See if those are usually found in MAGICAL
# Check to see if single cell RNA-seq data are balanced between high and low viral load individuals (4 HVL and 2 LVL, but how many cells for each?)







# Interesting idea - compare AUCs in our gene list(s) to other datasets (e.g., Data Compendium discovery or our own bulk RNA-seq)
# # IDEA - separate into HVL and LVL and see whether it has higher AUC for higher viral load ppl
comparison_aucs <- sc_pseudobulk_aucs
comparison_aucs$discovery_aucs <- sc_discovery_flu_aucs$flu_discovery_gene_auc
comparison_aucs <- comparison_aucs[,c(1,2,4,3)]