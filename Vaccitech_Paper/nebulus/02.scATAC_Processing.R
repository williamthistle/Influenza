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
# Output dir for ATAC
ATAC_output_dir <- paste0(output_dir, "ATAC", "/")
if (!dir.exists(ATAC_output_dir)) {dir.create(ATAC_output_dir)}
setwd(ATAC_output_dir)
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
# Read in ATAC data, filter data, perform initial processing, infer batches, and integrate by batch
atac_proj <- Read_ATAC(data_path = data_path, sample_id_list = sample_id_list, species = species, log_flag = TRUE)
atac_proj <- FilterRawData_ATAC(proj = atac_proj, log_flag = TRUE)
atac_proj <- InitialProcessing_ATAC(proj = atac_proj, log_flag = TRUE)
# Note - when I don't use batch correction (and same numFeatures / resolution as Vincy), my plot looks very similar to Vincy's!
atac_proj <- IntegrateByBatch_ATAC(proj = atac_proj, output_dir = ATAC_output_dir, log_flag = TRUE)
atac_proj <- MapCellTypes_ATAC(proj = atac_proj, reference = reference, output_dir = ATAC_output_dir,
                               reference_cell_type_attribute = reference_cell_type_attribute, log_flag = TRUE)
# save ArchR project: ArchR::saveArchRProject(ArchRProj = atac_proj, load = FALSE)
# load ArchR project: atac_proj <- loadArchRProject(path = paste0(ATAC_output_dir, "ArchROutput"))

atac_proj <- combine_cell_types_atac(atac_proj)
atac_proj <- MajorityVote_ATAC(proj = atac_proj)

num_cells <- length(atac_proj$cellNames)
num_samples <- length(unique(atac_proj$Sample))
sample_text <- paste0("(", num_samples, " Samples, ", 
                      num_cells, " Cells)")

pal <- paletteDiscrete(values = atac_proj$Cell_type_voting)
p1 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", 
                           name = "Cell_type_voting", embedding = "UMAP", 
                           pal = pal, force = TRUE, keepAxis = TRUE) + 
  ggplot2::ggtitle(paste0("ATAC Data Integration\n(By Majority Vote Cell Type)\n", 
                          sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size = 18), 
                                                         legend.text = ggplot2::element_text(size = 10))
ggplot2::ggsave(filename = paste0(ATAC_output_dir, "Final_ATAC_UMAP_by_Majority_Vote_Cell_Type.png"), 
                plot = p1, device = "png", width = 8, height = 8, 
                units = "in")

atac_proj <- add_sample_metadata_atac(atac_proj, high_viral_load_samples, low_viral_load_samples,
                                       d28_samples, d_minus_1_samples, male_samples, female_samples)
viral_load_metadata <- parse_metadata_for_samples(atac_proj, "viral_load", high_viral_load_samples, low_viral_load_samples,
                                                  d28_samples, d_minus_1_samples, male_samples, female_samples)
day_metadata <- parse_metadata_for_samples(atac_proj, "time_point", high_viral_load_samples, low_viral_load_samples,
                                           d28_samples, d_minus_1_samples, male_samples, female_samples)
sex_metadata <- parse_metadata_for_samples(atac_proj, "sex", high_viral_load_samples, low_viral_load_samples,
                                           d28_samples, d_minus_1_samples, male_samples, female_samples)

cluster_info <- get_cluster_info(atac_proj)

# Remove messy clusters
# C4 looks like unintegrated cells but it comes from LVL so it won't be relevant for our HVL analysis
idxPass <- which(atac_proj$seurat_clusters %in% c("4", "7", "13", "16", "20", "24", "28", "29"))
cellsPass <- atac_proj$cellNames[-idxPass]
atac_proj_minus_clusters <- atac_proj[cellsPass, ]

num_cells <- length(atac_proj_minus_clusters$cellNames)
num_samples <- length(unique(atac_proj_minus_clusters$Sample))
sample_text <- paste0("(", num_samples, " Samples, ", 
                      num_cells, " Cells)")

pal <- paletteDiscrete(values = atac_proj_minus_clusters$Cell_type_voting)
p1 <- ArchR::plotEmbedding(ArchRProj = atac_proj_minus_clusters, colorBy = "cellColData", 
                           name = "Cell_type_voting", embedding = "UMAP", 
                           pal = pal, force = TRUE, keepAxis = TRUE) + 
  ggplot2::ggtitle(paste0("ATAC Data Integration\n(By Majority Vote Cell Type)\n", 
                          sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size = 18), 
                                                         legend.text = ggplot2::element_text(size = 10))
ggplot2::ggsave(filename = paste0(ATAC_output_dir, "Final_ATAC_UMAP_by_Majority_Vote_Cell_Type_Minus_Clusters.png"), 
                plot = p1, device = "png", width = 8, height = 8, 
                units = "in")

# Subset to HVL
idxPass <- which(atac_proj_minus_clusters$viral_load %in% c("high"))
cellsPass <- atac_proj_minus_clusters$cellNames[idxPass]
HVL_proj_minus_clusters <- atac_proj_minus_clusters[cellsPass, ]

num_cells <- length(HVL_proj_minus_clusters$cellNames)
num_samples <- length(unique(HVL_proj_minus_clusters$Sample))
sample_text <- paste0("(", num_samples, " Samples, ", 
                      num_cells, " Cells)")

pal <- paletteDiscrete(values = HVL_proj_minus_clusters$Cell_type_voting)
p1 <- ArchR::plotEmbedding(ArchRProj = HVL_proj_minus_clusters, colorBy = "cellColData", 
                           name = "Cell_type_voting", embedding = "UMAP", 
                           pal = pal, force = TRUE, keepAxis = TRUE) + 
  ggplot2::ggtitle(paste0("ATAC Data Integration\n(By Majority Vote Cell Type, Reclustered)\n", 
                          sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size = 18), 
                                                         legend.text = ggplot2::element_text(size = 10))
ggplot2::ggsave(filename = paste0(ATAC_output_dir, "Final_ATAC_UMAP_by_Majority_Vote_Cell_Type_Minus_Clusters_HVL.png"), 
                plot = p1, device = "png", width = 8, height = 8, 
                units = "in")


# Call peaks
addArchRGenome("hg38")
HVL_proj_minus_clusters <- pseudo_bulk_replicates_and_call_peaks(HVL_proj_minus_clusters)
# Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
# calculations
HVL_proj_minus_clusters <- addPeakMatrix(HVL_proj_minus_clusters)

# Save peak metadata
HVL_peaks <- getPeakSet(HVL_proj_minus_clusters)
HVL_peaks_df <- as.data.frame(HVL_peaks@seqnames)
HVL_peaks_df <- cbind(HVL_peaks_df, as.data.frame(HVL_peaks@ranges))
HVL_peaks_df <- cbind(HVL_peaks_df, as.data.frame(HVL_peaks@elementMetadata))
write.table(x = HVL_peaks_df, file = paste0(ATAC_output_dir, "HVL_peaks_info.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Create Peaks.txt file
peak_txt_file <- create_peaks_file(HVL_proj_minus_clusters, ATAC_output_dir)
# Create peak_motif_matches.txt file
HVL_proj_minus_clusters <- create_peak_motif_matches_file(HVL_proj_minus_clusters, ATAC_output_dir, peak_txt_file)
# save ArchR project: ArchR::saveArchRProject(ArchRProj = HVL_proj_minus_clusters, outputDirectory = paste0(ATAC_output_dir, "HVL"), load = FALSE)
# load ArchR project: HVL_proj_minus_clusters <- loadArchRProject(path = paste0(ATAC_output_dir, "HVL"))
# Create pseudobulk counts for peaks for each cell type
# pseudo_bulk_dir <- paste0(ATAC_output_dir, "pseudo_bulk_atac/", date, "/")
# if (!dir.exists(pseudo_bulk_dir)) {dir.create(pseudo_bulk_dir, recursive = TRUE)}
# create_pseudobulk_atac(HVL_proj_minus_clusters, pseudo_bulk_dir)
# Find DASs
sample_metadata <- HVL_proj_minus_clusters$Sample
for(j in 1:nrow(sample_metadata_for_SPEEDI_df)) {
  sample_metadata <- gsub(rownames(sample_metadata_for_SPEEDI_df)[j], sample_metadata_for_SPEEDI_df[j,1], sample_metadata)
}
HVL_proj_minus_clusters <- addCellColData(ArchRProj = HVL_proj_minus_clusters, data = sample_metadata, cells = HVL_proj_minus_clusters$cellNames, name = "subject_id", force = TRUE)


differential_peaks_dir <- paste0(ATAC_output_dir, "diff_peaks/", date, "/motifs_without_tss_enrichment_bias/")
if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir, recursive = TRUE)}
calculate_daps_for_each_cell_type(HVL_proj_minus_clusters, differential_peaks_dir, sample_metadata_for_SPEEDI_df)


### ETC ###
# Print distributions for each cell type and create cell type proportions file for MAGICAL
print_cell_type_distributions(final_proj)
create_cell_type_proportion_MAGICAL_atac(final_proj, ATAC_output_dir, c("time_point"), day_metadata)


pal <- paletteDiscrete(values = atac_proj$Cell_type_voting)
p1 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p2 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", pal = pal, force = TRUE, keepAxis = TRUE)
p3 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "seurat_clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p4 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p5 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
p6 <- ArchR::plotEmbedding(ArchRProj = atac_proj, colorBy = "cellColData", name = "Batch", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
ArchR::plotPDF(p1,p2,p3,p4,p5,p6, name = "UMAP_new_batch_correction_and_majority_vote", ArchRProj = atac_proj, addDOC = FALSE, width = 5, height = 5)

