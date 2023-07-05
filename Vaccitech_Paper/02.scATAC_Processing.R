library(SPEEDI)
library(Seurat)
library(parallel)
library(doMC)
# Load extra RNA functions
home_dir <- "~/"
source(paste0(home_dir, "extra_functions/rna/preprocessing_and_qc.R"))
source(paste0(home_dir, "extra_functions/rna/get_stats.R"))
source(paste0(home_dir, "extra_functions/rna/manipulate_data.R"))
source(paste0(home_dir, "extra_functions/rna/differential_expression.R"))
source(paste0(home_dir, "extra_functions/rna/visualization.R"))

################## SETUP ##################
date <- Sys.Date()
data_path <- "~/single_cell/data"
reference_dir <- "~/references/"
reference_file_name <- "pbmc_multimodal.h5seurat"
output_dir <- "~/single_cell/analysis/"
analysis_name <- "primary_analysis_6_subject_12_sample"
reference_tissue <- "pbmc_full"
reference_cell_type_attribute <- "celltype.l2"
species <- "human"
record_doublets <- TRUE
data_type <- "RNA"
# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "single_cell_paired_sample"
# Grab samples that we want to analyze
data_tokens <- read.table(paste0(home_dir, "flu_data_tokens.tsv"), header = TRUE)
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))
samples <- data_tokens[data_tokens$token == data_token,]$samples
samples <-  unlist(strsplit(samples, ","))
sample_id_list <- all_sample_id_list[all_sample_id_list %in% samples]

# Load information about samples
sample_metadata <- read.table(paste0(home_dir, "/all_metadata_sheet.tsv"), sep = "\t", header = TRUE)
sample_metadata <- sample_metadata[sample_metadata$aliquot_id %in% sample_id_list,]
sample_metadata_for_SPEEDI_df <- sample_metadata
rownames(sample_metadata_for_SPEEDI_df) <- sample_metadata_for_SPEEDI_df$aliquot_id
sample_metadata_for_SPEEDI_df <- sample_metadata_for_SPEEDI_df[c("subject_id", "time_point", "sex", "viral_load")]
sample_metadata_for_SPEEDI_df$time_point[sample_metadata_for_SPEEDI_df$time_point == 'D-1'] <- 'D_minus_1'

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

# Output dir for ATAC
ATAC_output_dir <- paste0(output_dir, "ATAC", "/")
if (!dir.exists(ATAC_output_dir)) {dir.create(ATAC_output_dir)}
setwd(ATAC_output_dir)
# Read in ATAC data, filter data, perform initial processing, infer batches, and integrate by batch
atac_proj <- Read_ATAC(data_path = data_path, sample_id_list = sample_id_list, species = species, log_flag = TRUE)
atac_proj <- FilterRawData_ATAC(proj = atac_proj, log_flag = TRUE)
atac_proj <- InitialProcessing_ATAC(proj = atac_proj, log_flag = TRUE)
atac_proj <- IntegrateByBatch_ATAC_alt(proj = atac_proj, log_flag = TRUE)
atac_proj <- MapCellTypes_ATAC(proj = atac_proj, reference = reference,
                               reference_cell_type_attribute = reference_cell_type_attribute, log_flag = TRUE)
ArchR::saveArchRProject(ArchRProj = atac_proj, load = FALSE)
# load ArchR project: atac_proj <- loadArchRProject(path = ATAC_output_dir)

print_cell_type_distributions(atac_proj)
create_cell_type_proportion_MAGICAL_atac(atac_proj, ATAC_output_dir, c("day"), day_metadata)
addArchRGenome("hg38")
atac_proj <- pseudo_bulk_replicates_and_call_peaks(atac_proj)
# Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
# calculations
atac_proj <- addPeakMatrix(atac_proj)
differential_peaks_dir <- paste0(ATAC_output_dir, "diff_peaks/", date, "/")
if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir, recursive = TRUE)}
calculate_daps_for_each_cell_type(atac_proj, differential_peaks_dir)
# Create Peaks.txt file for MAGICAL
peak_txt_file <- create_peaks_file(atac_proj, ATAC_output_dir)
# Create peak_motif_matches.txt file for MAGICAL
create_peak_motif_matches_file(atac_proj, ATAC_output_dir, peak_txt_file)
# Create pseudobulk counts for peaks for each cell type
pseudo_bulk_dir <- paste0(ATAC_output_dir, "pseudo_bulk_atac/", date, "/")
if (!dir.exists(pseudo_bulk_dir)) {dir.create(pseudo_bulk_dir, recursive = TRUE)}
create_pseudobulk_atac(atac_proj, pseudo_bulk_dir)
