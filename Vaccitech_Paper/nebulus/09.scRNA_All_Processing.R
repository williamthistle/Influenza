# Load extra RNA functions
source("~/00.setup.R")
home_dir <- "/Genomics/ogtr04/wat2/"

################## SETUP ##################
data_path <- paste0(home_dir, "single_cell/data")
reference_dir <- paste0(home_dir, "references/")
reference_file_name <- "pbmc_multimodal.h5seurat"
output_dir <- paste0(home_dir, "single_cell/analysis/")
analysis_name <- "all_single_cell"
reference_tissue <- "pbmc_full"
reference_cell_type_attribute <- "celltype.l2"
species <- "human"
record_doublets <- TRUE
data_type <- "RNA"
# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "all_single_cell"
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
sample_metadata_for_SPEEDI_df <- sample_metadata_for_SPEEDI_df[c("subject_id", "time_point", "sex", "viral_load_category", "treatment")]
sample_metadata_for_SPEEDI_df$time_point[sample_metadata_for_SPEEDI_df$time_point == 'D-1'] <- 'D_minus_1'

# Break down metadata by category
high_viral_load_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$viral_load == "high",]))
low_viral_load_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$viral_load == "low",]))
all_viral_load_samples <- c(high_viral_load_samples, moderate_viral_load_samples, low_viral_load_samples)
d28_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$time_point == "D28",]))
d_minus_1_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$time_point == "D_minus_1",]))
all_day_samples <- c(d28_samples, d_minus_1_samples)
male_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$sex == "M",]))
female_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$sex == "F",]))
all_sex_samples <- c(male_samples, female_samples)
placebo_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$treatment == "PLACEBO",]))
vaccinated_samples <- sort(rownames(sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$treatment == "MVA-NP+M1",]))
all_treatment_samples <- c(placebo_samples, vaccinated_samples)

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
all_sc_exp_matrices <- Read_RNA(input_dir = data_path, sample_id_list = sample_id_list, log_flag = TRUE)
sc_obj <- FilterRawData_RNA(all_sc_exp_matrices = all_sc_exp_matrices, species = species,
                            record_doublets = FALSE, output_dir = RNA_output_dir,
                            log_file_path = log_file_name, log_flag = TRUE)
rm(all_sc_exp_matrices)
sc_obj <- InitialProcessing_RNA(sc_obj = sc_obj, species = species, metadata_df = sample_metadata_for_SPEEDI_df, log_flag = TRUE)
sc_obj <- InferBatches(sc_obj = sc_obj, log_flag = TRUE)
save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".new.batch.inference.4.RNA.rds"))
# load(paste0(RNA_output_dir, "primary_analysis_6_subject_12_sample.new.batch.inference.4.RNA.rds"))
sc_obj <- IntegrateByBatch_RNA(sc_obj = sc_obj, log_flag = TRUE)
sc_obj <- VisualizeIntegration(sc_obj = sc_obj, output_dir = RNA_output_dir, log_flag = TRUE)
sc_obj <- MapCellTypes_RNA(sc_obj = sc_obj, reference = reference,
                           reference_cell_type_attribute = reference_cell_type_attribute,
                           output_dir = RNA_output_dir, log_flag = TRUE)
# save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".new.batch.inference.final.RNA.rds"))
# load(paste0(RNA_output_dir, "all_single_cell.new.batch.inference.final.RNA.rds"))

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

# save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".new.batch.inference.final.2.RNA.rds"))
load(paste0(RNA_output_dir, "all_single_cell.new.batch.inference.final.2.RNA.rds"))

# Capture info about each cluster
cluster_info <- capture_cluster_info(sc_obj)

# Combine cell types for MAGICAL and other analyses that require ATAC-seq (granularity isn't as good for ATAC-seq)
sc_obj <- combine_cell_types_magical(sc_obj)

# Test
messy_clusters <- c(17)
idxPass <- which(Idents(sc_obj) %in% messy_clusters)
cellsPass <- names(sc_obj$orig.ident[-idxPass])
sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

# Print UMAPs for all subjects (HVL and LVL)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Majority_Vote_Cell_Type.png",
               group_by_category = "predicted_celltype_majority_vote", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Cluster.png",
               group_by_category = "seurat_clusters", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Raw_Predicted_Cell_Type.png",
               group_by_category = "predicted.id", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Viral_Load.png",
               group_by_category = "viral_load", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Sample.png",
               group_by_category = "sample", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Day.png",
               group_by_category = "time_point", output_dir = RNA_output_dir,
               log_flag = log_flag)
print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Sex.png",
               group_by_category = "sex", output_dir = RNA_output_dir,
               log_flag = log_flag)

# Separate into relevant subsets

# HVL
idxPass <- which(sc_obj$viral_load_category %in% "high")
cellsPass <- names(sc_obj$orig.ident[idxPass])
hvl_sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

# HVL PLACEBO
idxPass <- which(hvl_sc_obj$treatment %in% "PLACEBO")
cellsPass <- names(hvl_sc_obj$orig.ident[idxPass])
hvl_placebo_sc_obj <- subset(x = hvl_sc_obj, subset = cell_name %in% cellsPass)

# HVL VACCINATED
idxPass <- which(hvl_sc_obj$treatment %in% "MVA-NP+M1")
cellsPass <- names(hvl_sc_obj$orig.ident[idxPass])
hvl_vaccinated_sc_obj <- subset(x = hvl_sc_obj, subset = cell_name %in% cellsPass)

# LVL
idxPass <- which(sc_obj$viral_load_category %in% "low")
cellsPass <- names(sc_obj$orig.ident[idxPass])
lvl_sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

# LVL PLACEBO
idxPass <- which(lvl_sc_obj$treatment %in% "PLACEBO")
cellsPass <- names(lvl_sc_obj$orig.ident[idxPass])
lvl_placebo_sc_obj <- subset(x = lvl_sc_obj, subset = cell_name %in% cellsPass)

# LVL VACCINATED
idxPass <- which(lvl_sc_obj$treatment %in% "MVA-NP+M1")
cellsPass <- names(lvl_sc_obj$orig.ident[idxPass])
lvl_vaccinated_sc_obj <- subset(x = lvl_sc_obj, subset = cell_name %in% cellsPass)

# Print UMAPs for each subset
print_final_UMAPs(hvl_sc_obj, RNA_output_dir, "HVL")
print_final_UMAPs(lvl_sc_obj, RNA_output_dir, "LVL")
print_final_UMAPs(hvl_placebo_sc_obj, RNA_output_dir, "HVL_PLACEBO")
print_final_UMAPs(hvl_vaccinated_sc_obj, RNA_output_dir, "HVL_VACCINATED")
print_final_UMAPs(lvl_placebo_sc_obj, RNA_output_dir, "LVL_PLACEBO")
print_final_UMAPs(lvl_vaccinated_sc_obj, RNA_output_dir, "LVL_VACCINATED")

# Record cell type proportions for each subset
create_RNA_cell_type_proportion_file(hvl_sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples, token = "HVL")
create_RNA_cell_type_proportion_file(lvl_sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples, token = "LVL")
create_RNA_cell_type_proportion_file(hvl_placebo_sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples, token = "HVL_PLACEBO")
create_RNA_cell_type_proportion_file(lvl_placebo_sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples, token = "LVL_PLACEBO")
create_RNA_cell_type_proportion_file(hvl_vaccinated_sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples, token = "HVL_VACCINATED")
create_RNA_cell_type_proportion_file(lvl_vaccinated_sc_obj, RNA_output_dir, "time_point", high_viral_load_samples, d28_samples, male_samples, token = "LVL_VACCINATED")

# Create MAGICAL input files for each relevant subset
create_magical_input_files(hvl_placebo_sc_obj, paste0(RNA_output_dir, "MAGICAL_HVL_PLACEBO_", date, "/"))
create_magical_input_files(lvl_placebo_sc_obj, paste0(RNA_output_dir, "MAGICAL_LVL_PLACEBO_", date, "/"))
create_magical_input_files(hvl_vaccinated_sc_obj, paste0(RNA_output_dir, "MAGICAL_HVL_VACCINATED_", date, "/"))

# Run differential expression analysis, per cell type, for each relevant subset
run_differential_expression_controlling_for_subject_id(hvl_placebo_sc_obj, paste0(RNA_output_dir, "DE_HVL_PLACEBO_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point", magical_cell_types = FALSE)
run_differential_expression_controlling_for_subject_id(lvl_placebo_sc_obj, paste0(RNA_output_dir, "DE_LVL_PLACEBO_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point", magical_cell_types = FALSE)
run_differential_expression_controlling_for_subject_id(hvl_vaccinated_sc_obj, paste0(RNA_output_dir, "DE_HVL_VACCINATED_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point", magical_cell_types = FALSE)

run_differential_expression_controlling_for_subject_id(hvl_placebo_sc_obj, paste0(RNA_output_dir, "DE_HVL_PLACEBO_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point", magical_cell_types = TRUE)
run_differential_expression_controlling_for_subject_id(lvl_placebo_sc_obj, paste0(RNA_output_dir, "DE_LVL_PLACEBO_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point", magical_cell_types = TRUE)
run_differential_expression_controlling_for_subject_id(hvl_vaccinated_sc_obj, paste0(RNA_output_dir, "DE_HVL_VACCINATED_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point", magical_cell_types = TRUE)

