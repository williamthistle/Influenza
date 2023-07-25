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
# Output dir for ATAC
ATAC_output_dir <- paste0(output_dir, "ATAC", "/")
if (!dir.exists(ATAC_output_dir)) {dir.create(ATAC_output_dir)}
setwd(ATAC_output_dir)

# Load reference
reference <- LoadReferenceSPEEDI(reference_tissue = reference_tissue, species = species, reference_dir = reference_dir,
                                 reference_file_name = reference_file_name, log_flag = TRUE)
idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
reference <- reference[,-idx]
# Read in ATAC data, filter data, perform initial processing, infer batches, and integrate by batch
atac_proj <- Read_ATAC(data_path = data_path, sample_id_list = sample_id_list, species = species, log_flag = TRUE)
atac_proj <- FilterRawData_ATAC(proj = atac_proj, log_flag = TRUE)
atac_proj <- InitialProcessing_ATAC(proj = atac_proj, log_flag = TRUE)
# Note - when I don't use batch correction (and same numFeatures / resolution as Vincy), my plot looks very similar to Vincy's!
atac_proj <- IntegrateByBatch_ATAC_alt(proj = atac_proj, log_flag = TRUE)
atac_proj <- MapCellTypes_ATAC(proj = atac_proj, reference = reference,
                               reference_cell_type_attribute = reference_cell_type_attribute, log_flag = TRUE)
# save ArchR project: ArchR::saveArchRProject(ArchRProj = atac_proj, load = FALSE)
# load ArchR project: atac_proj <- loadArchRProject(path = paste0(ATAC_output_dir, "ArchROutput"))

atac_proj <- combine_cell_types_atac(atac_proj)

pal <- paletteDiscrete(values = atac_proj$Cell_type_voting)
p1 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p2 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", pal = pal, force = TRUE, keepAxis = TRUE)
p3 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Harmony_clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p4 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p5 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p6 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
ArchR::plotPDF(p1,p2,p3,p4,p5,p6, name = "UMAP_after_Final_Cell_Type_Majority_Voting_plots_combined_cell_types", ArchRProj = atac_proj, addDOC = FALSE, width = 5, height = 5)

cluster_info <- get_cluster_info(atac_proj)

# Remove messy clusters
idxPass <- which(atac_proj$Harmony_clusters %in% c("C4", "C11", "C17", "C21", "C23", "C24", "C25", "C27", "C28", "C29", "C30", "C31", "C32"))
cellsPass <- atac_proj$cellNames[-idxPass]
final_proj <- atac_proj[cellsPass, ]

# Override label
final_proj <- override_cluster_label_atac(final_proj, "C8", "Proliferating")

pal <- paletteDiscrete(values = final_proj$Cell_type_voting)
p1 <- ArchR::plotEmbedding(ArchRProj = final_proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", pal = pal, force = TRUE, keepAxis = TRUE)
p2 <- ArchR::plotEmbedding(ArchRProj = final_proj, colorBy = "cellColData", name = "Harmony_clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p3 <- ArchR::plotEmbedding(ArchRProj = final_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p4 <- ArchR::plotEmbedding(ArchRProj = final_proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
ArchR::plotPDF(p1,p2,p3,p4, name = "UMAP_after_Final_Cell_Type_Majority_Voting_plots_combined_cell_types_minus_messy_clusters", ArchRProj = atac_proj, addDOC = FALSE, width = 5, height = 5)

# Add sample metadata
final_proj <- add_sample_metadata_atac(final_proj, high_viral_load_samples, low_viral_load_samples,
                                       d28_samples, d_minus_1_samples, male_samples, female_samples)
viral_load_metadata <- parse_metadata_for_samples(final_proj, "viral_load", high_viral_load_samples, low_viral_load_samples,
                                                  d28_samples, d_minus_1_samples, male_samples, female_samples)
day_metadata <- parse_metadata_for_samples(final_proj, "time_point", high_viral_load_samples, low_viral_load_samples,
                                           d28_samples, d_minus_1_samples, male_samples, female_samples)
sex_metadata <- parse_metadata_for_samples(final_proj, "sex", high_viral_load_samples, low_viral_load_samples,
                                           d28_samples, d_minus_1_samples, male_samples, female_samples)

# Print distributions for each cell type and create cell type proportions file for MAGICAL
print_cell_type_distributions(final_proj)
create_cell_type_proportion_MAGICAL_atac(final_proj, ATAC_output_dir, c("time_point"), day_metadata)

# Call peaks
addArchRGenome("hg38")
final_proj <- pseudo_bulk_replicates_and_call_peaks(final_proj)
# Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
# calculations
final_proj <- addPeakMatrix(final_proj)
# TODO: Make it NK_MAGICAL instead of NK? So it's synced with DEGs
differential_peaks_dir <- paste0(ATAC_output_dir, "diff_peaks/", date, "/")
if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir, recursive = TRUE)}
calculate_daps_for_each_cell_type(final_proj, differential_peaks_dir)
# Create Peaks.txt file for MAGICAL
peak_txt_file <- create_peaks_file(final_proj, ATAC_output_dir)
# Create peak_motif_matches.txt file for MAGICAL
create_peak_motif_matches_file(final_proj, ATAC_output_dir, peak_txt_file)
# Create pseudobulk counts for peaks for each cell type
pseudo_bulk_dir <- paste0(ATAC_output_dir, "pseudo_bulk_atac/", date, "/")
if (!dir.exists(pseudo_bulk_dir)) {dir.create(pseudo_bulk_dir, recursive = TRUE)}
create_pseudobulk_atac(final_proj, pseudo_bulk_dir)

# Do HVL work
idxPass <- which(final_proj$viral_load %in% c("high"))
cellsPass <- final_proj$cellNames[idxPass]
HVL_final_proj <- final_proj[cellsPass, ]

pal <- paletteDiscrete(values = HVL_final_proj$Cell_type_voting)
p1 <- ArchR::plotEmbedding(ArchRProj = HVL_final_proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", pal = pal, force = TRUE, keepAxis = TRUE)
p2 <- ArchR::plotEmbedding(ArchRProj = HVL_final_proj, colorBy = "cellColData", name = "Harmony_clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p3 <- ArchR::plotEmbedding(ArchRProj = HVL_final_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p4 <- ArchR::plotEmbedding(ArchRProj = HVL_final_proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
ArchR::plotPDF(p1,p2,p3,p4, name = "HVL_UMAP_after_Final_Cell_Type_Majority_Voting_plots_combined_cell_types_minus_messy_clusters", ArchRProj = HVL_final_proj, addDOC = FALSE, width = 5, height = 5)

hvl_day_metadata <- parse_metadata_for_samples(HVL_final_proj, "time_point", high_viral_load_samples, low_viral_load_samples,
                                               d28_samples, d_minus_1_samples, male_samples, female_samples)
hvl_sex_metadata <- parse_metadata_for_samples(HVL_final_proj, "sex", high_viral_load_samples, low_viral_load_samples,
                                               d28_samples, d_minus_1_samples, male_samples, female_samples)

create_cell_type_proportion_MAGICAL_atac(HVL_final_proj, ATAC_output_dir, c("time_point"), hvl_day_metadata)

# Call peaks
addArchRGenome("hg38")
HVL_final_proj <- pseudo_bulk_replicates_and_call_peaks(HVL_final_proj)
# Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
# calculations
HVL_final_proj <- addPeakMatrix(HVL_final_proj)
# TODO: Make it NK_MAGICAL instead of NK? So it's synced with DEGs
HVL_differential_peaks_dir <- paste0(ATAC_output_dir, "diff_peaks/", date, "/HVL/")
if (!dir.exists(HVL_differential_peaks_dir)) {dir.create(HVL_differential_peaks_dir, recursive = TRUE)}
calculate_daps_for_each_cell_type(HVL_final_proj, HVL_differential_peaks_dir)
# Create Peaks.txt file for MAGICAL
HVL_peak_txt_file <- create_peaks_file(HVL_final_proj, ATAC_output_dir)
# Create peak_motif_matches.txt file for MAGICAL
create_peak_motif_matches_file(HVL_final_proj, ATAC_output_dir, HVL_peak_txt_file)
# Create pseudobulk counts for peaks for each cell type
HVL_pseudo_bulk_dir <- paste0(ATAC_output_dir, "pseudo_bulk_atac/", date, "/HVL/")
if (!dir.exists(HVL_pseudo_bulk_dir)) {dir.create(HVL_pseudo_bulk_dir, recursive = TRUE)}
create_pseudobulk_atac(HVL_final_proj, HVL_pseudo_bulk_dir)

