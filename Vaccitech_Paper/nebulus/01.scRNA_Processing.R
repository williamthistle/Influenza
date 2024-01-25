# Load extra RNA functions
source("~/00.setup.R")
home_dir <- "/Genomics/ogtr04/wat2/"

################## SETUP ##################
data_path <- paste0(home_dir, "single_cell/data")
reference_dir <- paste0(home_dir, "references/")
reference_file_name <- "pbmc_multimodal.h5seurat"
output_dir <- paste0(home_dir, "single_cell/analysis/")
analysis_name <- "primary_analysis_6_subject_12_sample"
reference_tissue <- "pbmc_full"
reference_cell_type_attribute <- "celltype.l2"
species <- "human"
record_doublets <- TRUE
data_type <- "RNA"
# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "single_cell_paired_sample"
# Grab samples that we want to analyze
data_tokens <- read.table("~/flu_data_tokens.tsv", header = TRUE)
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))
samples <- data_tokens[data_tokens$token == data_token,]$samples
samples <- unlist(strsplit(samples, ","))
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
                            record_doublets = FALSE, output_dir = RNA_output_dir,
                            log_file_path = log_file_name, log_flag = TRUE)
rm(all_sc_exp_matrices)
sc_obj <- InitialProcessing_RNA(sc_obj = sc_obj, species = species, metadata_df = sample_metadata_for_SPEEDI_df, log_flag = TRUE)
sc_obj <- InferBatches(sc_obj = sc_obj, log_flag = TRUE)
# save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".new.batch.inference.4.RNA.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_6_subject_12_sample.new.batch.inference.4.RNA.rds"))
sc_obj <- IntegrateByBatch_RNA(sc_obj = sc_obj, log_flag = TRUE)
sc_obj <- VisualizeIntegration(sc_obj = sc_obj, output_dir = RNA_output_dir, log_flag = TRUE)
sc_obj <- MapCellTypes_RNA(sc_obj = sc_obj, reference = reference,
                           reference_cell_type_attribute = reference_cell_type_attribute,
                           output_dir = RNA_output_dir, log_flag = TRUE)
# save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".new.batch.inference.final.RNA.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_6_subject_12_sample.new.batch.inference.final.RNA.rds"))
# vincy_obj <- readRDS("~/single_cell/analysis/vincy_analysis/integrated_obj_labeled.rds") # VINCY'S ANALYSIS

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
idx <- grep("Treg", Cell_type_combined)
Cell_type_combined[idx] <- "CD4 Memory"
sc_obj$predicted.id <- Cell_type_combined
sc_obj <- MajorityVote_RNA(sc_obj)

# Capture info about each cluster
#cluster_info <- capture_cluster_info(sc_obj)

# Combine cell types for MAGICAL and other analyses that require ATAC-seq (granularity isn't as good for ATAC-seq)
#sc_obj <- combine_cell_types_magical(sc_obj)

# May be 32, not sure
messy_clusters <- c(33)
idxPass <- which(Idents(sc_obj) %in% messy_clusters)
cellsPass <- names(sc_obj$orig.ident[-idxPass])
sc_obj_minus_clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

# Print UMAPs for all subjects (HVL and LVL)
print_UMAP_RNA(sc_obj_minus_clusters, file_name = "Final_RNA_UMAP_by_Majority_Vote_Cell_Type.png",
               group_by_category = "predicted_celltype_majority_vote", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_minus_clusters, file_name = "Final_RNA_UMAP_by_Cluster_test_2.png",
               group_by_category = "seurat_clusters", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_minus_clusters, file_name = "Final_RNA_UMAP_by_Raw_Predicted_Cell_Type.png",
               group_by_category = "predicted.id", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_minus_clusters, file_name = "Final_RNA_UMAP_by_Viral_Load.png",
               group_by_category = "viral_load", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_minus_clusters, file_name = "Final_RNA_UMAP_by_Sample.png",
               group_by_category = "sample", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_minus_clusters, file_name = "Final_RNA_UMAP_by_Day.png",
               group_by_category = "time_point", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj_minus_clusters, file_name = "Final_RNA_UMAP_by_Sex.png",
               group_by_category = "sex", output_dir = RNA_output_dir,
               log_flag = log_flag)

# HVL WORK
idxPass <- which(sc_obj$viral_load %in% "high")
cellsPass <- names(sc_obj$orig.ident[idxPass])
hvl_sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

hvl_sc_obj <- MajorityVote_RNA(hvl_sc_obj)

hvl_cluster_info <- capture_cluster_info(hvl_sc_obj)

messy_clusters <- c(32)
idxPass <- which(Idents(hvl_sc_obj) %in% messy_clusters)
cellsPass <- names(hvl_sc_obj$orig.ident[-idxPass])
hvl_sc_obj <- subset(x = hvl_sc_obj, subset = cell_name %in% cellsPass)

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

# Combine cell types for MAGICAL and other analyses that require ATAC-seq (granularity isn't as good for ATAC-seq)
hvl_sc_obj <- combine_cell_types_magical(hvl_sc_obj)

# Record cell type proportions
create_magical_cell_type_proportion_file(hvl_sc_obj, "/Genomics/ogtr04/wat2/", "time_point", high_viral_load_samples, d28_samples, male_samples)

# save(hvl_sc_obj, file = paste0(RNA_output_dir, analysis_name, ".hvl.new.batch.inference.final.RNA.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_6_subject_12_sample.hvl.new.batch.inference.final.RNA.rds"))

# Update NK_MAGICAL to be NK for MAGICAL cell types
cell_type_combined <- hvl_sc_obj$magical_cell_types
idx <- grep("NK", cell_type_combined)
cell_type_combined[idx] <- "NK"
hvl_sc_obj$magical_cell_types <- cell_type_combined

MAGICAL_file_dir <- paste0(RNA_output_dir, "MAGICAL/")
if (!dir.exists(MAGICAL_file_dir)) {dir.create(MAGICAL_file_dir)}
MAGICAL_cell_metadata_dir <- paste0(MAGICAL_file_dir, "scRNA_Cell_Metadata/")
if (!dir.exists(MAGICAL_cell_metadata_dir)) {dir.create(MAGICAL_cell_metadata_dir)}
MAGICAL_read_counts_dir <- paste0(MAGICAL_file_dir, "scRNA_Read_Counts/")
if (!dir.exists(MAGICAL_read_counts_dir)) {dir.create(MAGICAL_read_counts_dir)}
MAGICAL_genes_dir <- paste0(MAGICAL_file_dir, "scRNA_Genes/")
if (!dir.exists(MAGICAL_genes_dir)) {dir.create(MAGICAL_genes_dir)}

write.table(hvl_sc_obj@assays$RNA@counts@Dimnames[[1]], file = paste0(MAGICAL_genes_dir, "HVL_RNA_genes.tsv"),
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")

for(cell_type in unique(hvl_sc_obj$magical_cell_types)) {
  print(cell_type)
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  
  #RNA assay cell read counts
  cell_index=which(hvl_sc_obj$magical_cell_types==cell_type)
  cell_type_scRNA_counts = as(hvl_sc_obj@assays$RNA@counts[, cell_index], "dgTMatrix")
  saveRDS(cell_type_scRNA_counts, file= paste0(MAGICAL_read_counts_dir, cell_type_for_file_name, "_HVL_RNA_read_counts.rds"))
  
  # Cell metadata
  cell_type_scRNA_meta = hvl_sc_obj@meta.data[cell_index,]
  write.table(data.frame(rownames(cell_type_scRNA_meta), cell_type_scRNA_meta$magical_cell_types, cell_type_scRNA_meta$sample, cell_type_scRNA_meta$time_point),
              file = paste0(MAGICAL_cell_metadata_dir, cell_type_for_file_name, "_HVL_RNA_cell_metadata.tsv"),
              quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = "\t")
}
  
HVL_differential_genes_dir <- paste0(RNA_output_dir, "diff_genes/", date, "/HVL_controlling_for_subject_id/")
if (!dir.exists(HVL_differential_genes_dir)) {dir.create(HVL_differential_genes_dir, recursive = TRUE)}
run_differential_expression_controlling_for_subject_id(hvl_sc_obj, HVL_differential_genes_dir, sample_metadata_for_SPEEDI_df, "time_point", unique(hvl_sc_obj$predicted_celltype_majority_vote))
HVL_differential_genes_MAGICAL_dir <- paste0(RNA_output_dir, "diff_genes/", date, "/HVL_controlling_for_subject_id_MAGICAL/")
if (!dir.exists(HVL_differential_genes_MAGICAL_dir)) {dir.create(HVL_differential_genes_MAGICAL_dir, recursive = TRUE)}
run_differential_expression_controlling_for_subject_id(hvl_sc_obj, HVL_differential_genes_MAGICAL_dir, sample_metadata_for_SPEEDI_df, "time_point", unique(hvl_sc_obj$magical_cell_types))

create_magical_cell_type_proportion_file(hvl_sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples)
hvl_pseudobulk_rna_dir <- paste0(RNA_output_dir, "pseudobulk_rna/", date, "/HVL_RNA/")
if (!dir.exists(hvl_pseudobulk_rna_dir)) {dir.create(hvl_pseudobulk_rna_dir, recursive = TRUE)}
create_magical_cell_type_pseudobulk_files(hvl_sc_obj, hvl_pseudobulk_rna_dir)

# LVL WORK
idxPass <- which(sc_obj$viral_load %in% "low")
cellsPass <- names(sc_obj$orig.ident[idxPass])
lvl_sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

lvl_sc_obj <- MajorityVote_RNA(lvl_sc_obj)

lvl_cluster_info <- capture_cluster_info(lvl_sc_obj)

messy_clusters <- c(30)
idxPass <- which(Idents(lvl_sc_obj) %in% messy_clusters)
cellsPass <- names(lvl_sc_obj$orig.ident[-idxPass])
lvl_sc_obj <- subset(x = lvl_sc_obj, subset = cell_name %in% cellsPass)

print_UMAP_RNA(lvl_sc_obj, file_name = "lvl_Final_Combined_Cell_Type_RNA_UMAP_by_Majority_Vote_Cell_Type.png",
               group_by_category = "predicted_celltype_majority_vote", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(lvl_sc_obj, file_name = "lvl_Final_Combined_Cell_Type_RNA_UMAP_by_Cluster.png",
               group_by_category = "seurat_clusters", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(lvl_sc_obj, file_name = "lvl_Final_Combined_Cell_Type_RNA_UMAP_by_Raw_Predicted_Cell_Type.png",
               group_by_category = "predicted.id", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(lvl_sc_obj, file_name = "lvl_Final_Combined_Cell_Type_RNA_UMAP_by_Sample.png",
               group_by_category = "sample", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(lvl_sc_obj, file_name = "lvl_Final_Combined_Cell_Type_RNA_UMAP_by_Day.png",
               group_by_category = "time_point", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(lvl_sc_obj, file_name = "lvl_Final_Combined_Cell_Type_RNA_UMAP_by_Sex.png",
               group_by_category = "sex", output_dir = RNA_output_dir,
               log_flag = log_flag)

# Combine cell types for MAGICAL and other analyses that require ATAC-seq (granularity isn't as good for ATAC-seq)
lvl_sc_obj <- combine_cell_types_magical(lvl_sc_obj)

# Record cell type proportions
create_magical_cell_type_proportion_file(lvl_sc_obj, "/Genomics/ogtr04/wat2/", "time_point", high_viral_load_samples, d28_samples, male_samples, token = "LVL")

# save(lvl_sc_obj, file = paste0(RNA_output_dir, analysis_name, ".lvl.new.batch.inference.final.RNA.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_6_subject_12_sample.lvl.new.batch.inference.final.RNA.rds"))

# Update NK_MAGICAL to be NK for MAGICAL cell types
cell_type_combined <- lvl_sc_obj$magical_cell_types
idx <- grep("NK", cell_type_combined)
cell_type_combined[idx] <- "NK"
lvl_sc_obj$magical_cell_types <- cell_type_combined

MAGICAL_file_dir <- paste0(RNA_output_dir, "MAGICAL/")
if (!dir.exists(MAGICAL_file_dir)) {dir.create(MAGICAL_file_dir)}
MAGICAL_cell_metadata_dir <- paste0(MAGICAL_file_dir, "scRNA_Cell_Metadata/")
if (!dir.exists(MAGICAL_cell_metadata_dir)) {dir.create(MAGICAL_cell_metadata_dir)}
MAGICAL_read_counts_dir <- paste0(MAGICAL_file_dir, "scRNA_Read_Counts/")
if (!dir.exists(MAGICAL_read_counts_dir)) {dir.create(MAGICAL_read_counts_dir)}
MAGICAL_genes_dir <- paste0(MAGICAL_file_dir, "scRNA_Genes/")
if (!dir.exists(MAGICAL_genes_dir)) {dir.create(MAGICAL_genes_dir)}

write.table(lvl_sc_obj@assays$RNA@counts@Dimnames[[1]], file = paste0(MAGICAL_genes_dir, "lvl_RNA_genes.tsv"),
            quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")

for(cell_type in unique(lvl_sc_obj$magical_cell_types)) {
  print(cell_type)
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  
  #RNA assay cell read counts
  cell_index=which(lvl_sc_obj$magical_cell_types==cell_type)
  cell_type_scRNA_counts = as(lvl_sc_obj@assays$RNA@counts[, cell_index], "dgTMatrix")
  saveRDS(cell_type_scRNA_counts, file= paste0(MAGICAL_read_counts_dir, cell_type_for_file_name, "_lvl_RNA_read_counts.rds"))
  
  # Cell metadata
  cell_type_scRNA_meta = lvl_sc_obj@meta.data[cell_index,]
  write.table(data.frame(rownames(cell_type_scRNA_meta), cell_type_scRNA_meta$magical_cell_types, cell_type_scRNA_meta$sample, cell_type_scRNA_meta$time_point),
              file = paste0(MAGICAL_cell_metadata_dir, cell_type_for_file_name, "_lvl_RNA_cell_metadata.tsv"),
              quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = "\t")
}

lvl_differential_genes_dir <- paste0(RNA_output_dir, "diff_genes/", date, "/lvl_controlling_for_subject_id/")
if (!dir.exists(lvl_differential_genes_dir)) {dir.create(lvl_differential_genes_dir, recursive = TRUE)}
run_differential_expression_controlling_for_subject_id(lvl_sc_obj, lvl_differential_genes_dir, sample_metadata_for_SPEEDI_df, "time_point", unique(lvl_sc_obj$predicted_celltype_majority_vote))
lvl_differential_genes_MAGICAL_dir <- paste0(RNA_output_dir, "diff_genes/", date, "/lvl_controlling_for_subject_id_MAGICAL/")
if (!dir.exists(lvl_differential_genes_MAGICAL_dir)) {dir.create(lvl_differential_genes_MAGICAL_dir, recursive = TRUE)}
run_differential_expression_controlling_for_subject_id(lvl_sc_obj, lvl_differential_genes_MAGICAL_dir, sample_metadata_for_SPEEDI_df, "time_point", unique(lvl_sc_obj$magical_cell_types))

create_magical_cell_type_proportion_file(lvl_sc_obj, RNA_output_dir, "time_point", low_viral_load_samples, d28_samples, male_samples)
lvl_pseudobulk_rna_dir <- paste0(RNA_output_dir, "pseudobulk_rna/", date, "/lvl_RNA/")
if (!dir.exists(lvl_pseudobulk_rna_dir)) {dir.create(lvl_pseudobulk_rna_dir, recursive = TRUE)}
create_magical_cell_type_pseudobulk_files(lvl_sc_obj, lvl_pseudobulk_rna_dir)

### ETC ###

# Find markers for each cluster
run_differential_expression_cluster(hvl_sc_obj, RNA_output_dir)

# Run differential expression for each cell type within each group of interest
#differential_genes_dir <- paste0(RNA_output_dir, "diff_genes/", date, "/alg_4/")
#if (!dir.exists(differential_genes_dir)) {dir.create(differential_genes_dir, recursive = TRUE)}
#run_differential_expression_group(sc_obj, differential_genes_dir, "time_point")
#run_differential_expression_group(sc_obj, differential_genes_dir, "viral_load") # NOT RUN YET
#run_differential_expression_group(sc_obj, differential_genes_dir, "sex")  # NOT RUN YET

#create_magical_cell_type_proportion_file(sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples)
#create_magical_cell_type_proportion_file(sc_obj, RNA_output_dir,"viral_load", high_viral_load_samples, d28_samples, male_samples)  # NOT RUN YET
#create_magical_cell_type_proportion_file(sc_obj, RNA_output_dir,"sex", high_viral_load_samples, d28_samples, male_samples)  # NOT RUN YET

#pseudobulk_rna_dir <- paste0(RNA_output_dir, "pseudobulk_rna/", date, "/")
#if (!dir.exists(pseudobulk_rna_dir)) {dir.create(pseudobulk_rna_dir, recursive = TRUE)}
#create_magical_cell_type_pseudobulk_files(sc_obj, pseudobulk_rna_dir)

DefaultAssay(hvl_sc_obj) <- "SCT"

idxPass <- which(hvl_sc_obj$time_point %in% "D28")
cellsPass <- names(hvl_sc_obj$orig.ident[idxPass])
sc_hvl_sc_obj_d28 <- subset(x = hvl_sc_obj, subset = cell_name %in% cellsPass)

idxPass <- which(hvl_sc_obj$time_point %in% "D_minus_1")
cellsPass <- names(hvl_sc_obj$orig.ident[idxPass])
sc_hvl_sc_obj_d_minus_1 <- subset(x = hvl_sc_obj, subset = cell_name %in% cellsPass)

d28_IRAK3_plot <- FeaturePlot(sc_hvl_sc_obj_d28, features = c("IRAK3"))
ggsave(filename = paste0(RNA_output_dir, "D28_IRAK3.png"), plot = d28_IRAK3_plot, device = "png")
d_minus_1_IRAK3_plot <- FeaturePlot(sc_hvl_sc_obj_d_minus_1, features = c("IRAK3"))
ggsave(filename = paste0(RNA_output_dir, "D_minus_1_IRAK3.png"), plot = d_minus_1_IRAK3_plot, device = "png")






idxPass <- which(hvl_sc_obj$predicted_celltype_majority_vote %in% "Unknown")
cellsPass <- names(hvl_sc_obj$orig.ident[idxPass])
unknown_hvl_sc_obj <- subset(x = hvl_sc_obj, subset = cell_name %in% cellsPass)

