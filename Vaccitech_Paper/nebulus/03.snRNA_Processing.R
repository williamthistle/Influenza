# Load extra RNA functions
source("~/00.setup.R")
home_dir <- "/Genomics/ogtr04/wat2/"

################## SETUP ##################
data_path <- paste0(home_dir, "multiome/data")
reference_dir <- paste0(home_dir, "references/")
reference_file_name <- "pbmc_multimodal.h5seurat"
output_dir <- paste0(home_dir, "multiome/analysis/")
analysis_name <- "primary_analysis_7_subject_14_sample"
reference_tissue <- "pbmc_full"
reference_cell_type_attribute <- "celltype.l2"
species <- "human"
record_doublets <- TRUE
data_type <- "RNA"
# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "all_multiome_paired_minus_0_sample"
# Grab samples that we want to analyze
data_tokens <- read.table("~/flu_data_tokens.tsv", header = TRUE)
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))
samples <- data_tokens[data_tokens$token == data_token,]$samples
samples <-  unlist(strsplit(samples, ","))
sample_id_list <- all_sample_id_list[all_sample_id_list %in% samples]

# Load information about samples
sample_metadata <- read.table("~/all_metadata_sheet.tsv", sep = "\t", header = TRUE)
sample_metadata <- sample_metadata[sample_metadata$aliquot_id %in% sample_id_list,]
sample_metadata_for_SPEEDI_df <- sample_metadata
rownames(sample_metadata_for_SPEEDI_df) <- sample_metadata_for_SPEEDI_df$aliquot_id
sample_metadata_for_SPEEDI_df <- sample_metadata_for_SPEEDI_df[c("subject_id", "time_point", "sex", "viral_load")]
sample_metadata_for_SPEEDI_df$time_point[sample_metadata_for_SPEEDI_df$time_point == 'D-1'] <- 'D_minus_1'

# Break down metadata by category
high_viral_load_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$viral_load == "high",]))
low_viral_load_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$viral_load == "low",]))
all_viral_load_samples <- c(high_viral_load_samples, low_viral_load_samples)
d28_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$time_point == "D28",]))
d_minus_1_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$time_point == "D_minus_1",]))
all_day_samples <- c(d28_samples, d_minus_1_samples)
male_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$sex == "M",]))
female_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$sex == "F",]))
all_sex_samples <- c(male_samples, female_samples)

# Normalize paths (in case user provides relative paths)
data_path <- normalize_dir_path(data_path)
reference_dir <- normalize_dir_path(reference_dir)
output_dir <- normalize_dir_path(output_dir)
# ArchR likes to write some files to the working directory, so we'll set our working directory to output_dir
# and then reset it to the original working directory once we're done running SPEEDI
old_wd <- getwd()
# Create output_dir if it doesn't already exist
if (!dir.exists(output_dir)) {dir.create(output_dir)}
# Add "/" to end of output_dir if not already present
last_char_of_output_dir_path <- substr(output_dir, nchar(output_dir), nchar(output_dir))
if(last_char_of_output_dir_path != "/") {
  output_dir <- paste0(output_dir, "/")
}
# Set analysis name
if(is.null(analysis_name)) {
  analysis_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
  analysis_name <- gsub(":", "-", analysis_name)
}
# Update our output dir to be the specific analysis directory
output_dir <- paste0(output_dir, analysis_name, "/")
if (!dir.exists(output_dir)) {dir.create(output_dir)}
setwd(output_dir)
# Output dirs for RNA
RNA_output_dir <- paste0(output_dir, "RNA", "/")
if (!dir.exists(RNA_output_dir)) {dir.create(RNA_output_dir)}
# Create log file
log_file_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
log_file_name <- gsub(":", "-", log_file_name)
log_file_name <- paste0(output_dir, log_file_name)
log_file <- logr::log_open(log_file_name, logdir = FALSE)
# Load reference
reference <- LoadReferenceSPEEDI(reference_tissue = reference_tissue, species = species, reference_dir = reference_dir,
                                 reference_file_name = reference_file_name, log_flag = TRUE)
idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
reference <- reference[,-idx]

# Read in RNA data, filter data, perform initial processing, infer batches, integrate by batch, and process UMAP of integration
all_sc_exp_matrices <- Read_RNA(data_path = data_path, sample_id_list = sample_id_list, log_flag = TRUE)
sc_obj <- FilterRawData_RNA(all_sc_exp_matrices = all_sc_exp_matrices, species = species,
                            record_doublets = record_doublets, output_dir = RNA_output_dir,
                            log_file_path = log_file_name, log_flag = TRUE)
rm(all_sc_exp_matrices)
sc_obj <- InitialProcessing_RNA(sc_obj = sc_obj, species = species, metadata_df = sample_metadata_for_SPEEDI_df, log_flag = TRUE)
# save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".step.4.RNA.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_7_subject_14_sample.step.4.RNA.rds"))
sc_obj <- InferBatches(sc_obj = sc_obj, log_flag = TRUE)
sc_obj <- IntegrateByBatch_RNA(sc_obj = sc_obj, log_flag = TRUE)
# save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".step.6.RNA.rds"))
sc_obj <- VisualizeIntegration(sc_obj = sc_obj, log_flag = TRUE)
# save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".step.7.RNA.rds"))
sc_obj <- MapCellTypes_RNA(sc_obj = sc_obj, reference = reference,
                           reference_cell_type_attribute = reference_cell_type_attribute,
                           output_dir = RNA_output_dir, data_type = "snRNA", log_flag = TRUE)
sc_obj <- MajorityVote_RNA_alt(sc_obj)
save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".step.8.RNA.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_7_subject_14_sample.step.8.RNA.rds"))

print_UMAP_RNA(sc_obj, file_name = "Initial_RNA_UMAP_by_Majority_Vote_Cell_Type_res_3.png",
               group_by_category = "predicted_celltype_majority_vote", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Initial_RNA_UMAP_by_Cluster_res_3.png",
               group_by_category = "seurat_clusters", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Initial_RNA_UMAP_by_Raw_Predicted_Cell_Type_res_3.png",
               group_by_category = "predicted.id", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Initial_RNA_UMAP_by_Viral_Load_res_3.png",
               group_by_category = "viral_load", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Initial_RNA_UMAP_by_Sample_res_3.png",
               group_by_category = "sample", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Initial_RNA_UMAP_by_Day_res_3.png",
               group_by_category = "time_point", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Initial_RNA_UMAP_by_Sex_res_3.png",
               group_by_category = "sex", output_dir = RNA_output_dir,
               log_flag = log_flag)


sc_obj$old.predicted.id <- sc_obj$predicted.id
Cell_type_combined <- sc_obj$predicted.id
idx <- grep("CD4 T", Cell_type_combined)
Cell_type_combined[idx] <- "CD4 Memory"
idx <- grep("CD8 T", Cell_type_combined)
Cell_type_combined[idx] <- "CD8 Memory"
idx <- grep("cDC", Cell_type_combined)
Cell_type_combined[idx] <- "cDC"
idx <- grep("Proliferating", Cell_type_combined)
Cell_type_combined[idx] <- "Proliferating"
#idx <- grep("CD4 Naive", Cell_type_combined)
#Cell_type_combined[idx] <- "T Naive"
#idx <- grep("CD8 Naive", Cell_type_combined)
#Cell_type_combined[idx] <- "T Naive"
#idx <- grep("Treg", Cell_type_combined)
#Cell_type_combined[idx] <- "T Naive"
sc_obj$predicted.id <- Cell_type_combined
sc_obj <- MajorityVote_RNA_alt(sc_obj, res = 3)

cluster_info <- capture_cluster_info(sc_obj)

messy_clusters <- c(1,2,4,7,12,16,17,22,27,28,47)
idxPass <- which(Idents(sc_obj) %in% messy_clusters)
cellsPass <- names(sc_obj$orig.ident[-idxPass])
sc_obj_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

sc_obj_subset <- remove_cells_based_on_umap(sc_obj_subset, -2.25, 0, 1, 2.5) # 14 sample multiome
sc_obj_subset <- remove_cells_based_on_umap(sc_obj_subset, 1.5, 2.5, -3, -1.5) # 14 sample multiome
sc_obj_subset <- remove_cells_based_on_umap(sc_obj_subset, 2.2, 3.5, -6, -4) # 14 sample multiome

print_UMAP_RNA(sc_obj_subset, file_name = "Final_Combined_Cell_Type_RNA_UMAP_by_Majority_Vote_Cell_Type.png",
               group_by_category = "predicted_celltype_majority_vote", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_subset, file_name = "Final_Combined_Cell_Type_RNA_UMAP_by_Cluster.png",
               group_by_category = "seurat_clusters", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_subset, file_name = "Final_Combined_Cell_Type_RNA_UMAP_by_Raw_Predicted_Cell_Type.png",
               group_by_category = "predicted.id", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_subset, file_name = "Final_Combined_Cell_Type_RNA_UMAP_by_Viral_Load.png",
               group_by_category = "viral_load", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_subset, file_name = "Final_Combined_Cell_Type_RNA_UMAP_by_Sample.png",
               group_by_category = "sample", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_subset, file_name = "Final_Combined_Cell_Type_RNA_UMAP_by_Day.png",
               group_by_category = "time_point", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_subset, file_name = "Final_Combined_Cell_Type_RNA_UMAP_by_Sex.png",
               group_by_category = "sex", output_dir = RNA_output_dir,
               log_flag = log_flag)

# Combine cell types for MAGICAL and other analyses that require ATAC-seq (granularity isn't as good for ATAC-seq)
sc_obj <- combine_cell_types_magical(sc_obj)
# Run differential expression for each cell type within each group of interest
differential_genes_dir <- paste0(RNA_output_dir, "diff_genes/", date, "/final/")
if (!dir.exists(differential_genes_dir)) {dir.create(differential_genes_dir, recursive = TRUE)}
run_differential_expression_group(sc_obj, differential_genes_dir, "time_point")
run_differential_expression_group(sc_obj, differential_genes_dir, "viral_load") # NOT RUN YET
run_differential_expression_group(sc_obj, differential_genes_dir, "sex")  # NOT RUN YET

create_magical_cell_type_proportion_file(sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples)
create_magical_cell_type_proportion_file(sc_obj, RNA_output_dir,"viral_load", high_viral_load_samples, d28_samples, male_samples)  # NOT RUN YET
create_magical_cell_type_proportion_file(sc_obj, RNA_output_dir,"sex", high_viral_load_samples, d28_samples, male_samples)  # NOT RUN YET

pseudobulk_rna_dir <- paste0(RNA_output_dir, "pseudobulk_rna/", date, "/")
if (!dir.exists(pseudobulk_rna_dir)) {dir.create(pseudobulk_rna_dir, recursive = TRUE)}
create_magical_cell_type_pseudobulk_files(sc_obj, pseudobulk_rna_dir)

# HVL WORK
idxPass <- which(sc_obj$viral_load %in% "high")
cellsPass <- names(sc_obj$orig.ident[idxPass])
hvl_sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

print_UMAP_RNA(hvl_sc_obj, file_name = "HVL_Final_Combined_Cell_Type_RNA_UMAP_by_Majority_Vote_Cell_Type.png",
               group_by_category = "predicted_celltype_majority_vote", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(hvl_sc_obj, file_name = "HVL_Final_Combined_Cell_Type_RNA_UMAP_by_Cluster.png",
               group_by_category = "seurat_clusters", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(hvl_sc_obj, file_name = "HVL_Final_Combined_Cell_Type_RNA_UMAP_by_Raw_Predicted_Cell_Type.png",
               group_by_category = "predicted.id", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(hvl_sc_obj, file_name = "HVL_Final_Combined_Cell_Type_RNA_UMAP_by_Sample.png",
               group_by_category = "sample", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(hvl_sc_obj, file_name = "HVL_Final_Combined_Cell_Type_RNA_UMAP_by_Day.png",
               group_by_category = "time_point", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(hvl_sc_obj, file_name = "HVL_Final_Combined_Cell_Type_RNA_UMAP_by_Sex.png",
               group_by_category = "sex", output_dir = RNA_output_dir,
               log_flag = log_flag)

HVL_differential_genes_dir <- paste0(RNA_output_dir, "diff_genes/", date, "/HVL_SCT_controlling_for_subject/")
if (!dir.exists(HVL_differential_genes_dir)) {dir.create(HVL_differential_genes_dir, recursive = TRUE)}
run_differential_expression_group(hvl_sc_obj, HVL_differential_genes_dir, "time_point")

create_magical_cell_type_proportion_file(hvl_sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples)
hvl_pseudobulk_rna_dir <- paste0(RNA_output_dir, "pseudobulk_rna/", date, "/HVL_RNA/")
if (!dir.exists(hvl_pseudobulk_rna_dir)) {dir.create(hvl_pseudobulk_rna_dir, recursive = TRUE)}
create_magical_cell_type_pseudobulk_files(hvl_sc_obj, hvl_pseudobulk_rna_dir)

# Run pseudobulk DE
pseudobulk_de_df <- run_de(hvl_sc_obj, replicate_col = "sample", cell_type_col = "magical_cell_types", label_col = "time_point", de_method = "DESeq2")
pseudobulk_de_df <- na.omit(pseudobulk_de_df)
pseudobulk_de_df <- pseudobulk_de_df[pseudobulk_de_df$p_val < 0.05,]
pseudobulk_de_df <- pseudobulk_de_df[pseudobulk_de_df$avg_logFC < -0.3 | pseudobulk_de_df$avg_logFC > 0.3,]

pseudobulk_cell_types_for_correction <- c("B", "CD4_Memory", "CD8_Memory", "CD14_Mono", "CD16_Mono", "NK_MAGICAL", "MAIT", "T_Naive")
DEG_dir <- "/home/wat2/single_cell/analysis/primary_analysis_6_subject_12_sample/RNA/diff_genes/2023-07-06/HVL_SCT/"
final_list_of_genes <- data.frame(Cell_Type = character(), Gene_Name = character(), sc_pval_adj = character(), sc_log2FC = character(), pseudo_bulk_pval = character(),
                                  pseudo_bulk_log2FC = character())
for(current_cell_type in pseudobulk_cell_types_for_correction) {
  current_DEG_table <- read.table(paste0(DEG_dir, "D28-vs-D_minus_1-degs-", current_cell_type, "-time_point.csv"), sep = ",", header = TRUE)
  current_DEG_table <- current_DEG_table[current_DEG_table$p_val_adj < 0.05,]
  current_DEG_table <- current_DEG_table[abs(current_DEG_table$avg_log2FC) > 0.1,]
  current_DEG_table <- current_DEG_table[current_DEG_table$pct.1 > 0.1 | current_DEG_table$pct.2 > 0.1,]
  if(current_cell_type != "NK_MAGICAL") {
    pseudobulk_de_df_cell_type_subset <- pseudobulk_de_df[pseudobulk_de_df$cell_type == sub("_", " ", current_cell_type),]
  } else {
    pseudobulk_de_df_cell_type_subset <- pseudobulk_de_df[pseudobulk_de_df$cell_type == current_cell_type,]
  }
  final_cell_type_genes <- intersect(current_DEG_table$X, pseudobulk_de_df_cell_type_subset$gene)
  for(current_gene in final_cell_type_genes) {
    current_sc_pval_adj <- current_DEG_table[current_DEG_table$X == current_gene,]$p_val_adj
    current_sc_log2FC <- current_DEG_table[current_DEG_table$X == current_gene,]$avg_log2FC
    current_pseudo_bulk_pval <- pseudobulk_de_df_cell_type_subset[pseudobulk_de_df_cell_type_subset$gene == current_gene,]$p_val
    current_pseudo_bulk_log2FC <- pseudobulk_de_df_cell_type_subset[pseudobulk_de_df_cell_type_subset$gene == current_gene,]$avg_logFC
    current_row <- data.frame(current_cell_type, current_gene, current_sc_pval_adj, current_sc_log2FC, current_pseudo_bulk_pval, current_pseudo_bulk_log2FC)
    names(current_row) <- c("Cell_Type", "Gene_Name", "sc_pval_adj", "sc_log2FC", "pseudo_bulk_pval", "pseudo_bulk_log2FC")
    final_list_of_genes <- rbind(final_list_of_genes, current_row)
  }
}

write.table(final_list_of_genes, paste0(DEG_dir, "D28_D1_DESeq2_pseudobulk_genes.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
# Add printing of positive and negative fold change for ease

