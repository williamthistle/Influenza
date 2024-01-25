# Load extra RNA functions
source("~/00.setup.R")
home_dir <- "/Genomics/ogtr04/wat2/"

### SHORTCUT ###
HVL_proj_minus_clusters <- loadArchRProject(path = paste0(ATAC_output_dir, "HVL"))
seurat_atac <- readRDS(file = paste0(ATAC_output_dir, "HVL_seurat.RDS"))

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
idxPass <- which(atac_proj$seurat_clusters %in% c("4", "7", "10", "13", "15", "16", "18", "20", "24", "26", "28", "29"))
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
# save ArchR project: ArchR::saveArchRProject(ArchRProj = HVL_proj_minus_clusters, outputDirectory = paste0(ATAC_output_dir, "HVL"), load = FALSE, overwrite = TRUE)
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

MAGICAL_file_dir <- paste0(ATAC_output_dir, "MAGICAL/")
if (!dir.exists(MAGICAL_file_dir)) {dir.create(MAGICAL_file_dir)}
MAGICAL_cell_metadata_dir <- paste0(MAGICAL_file_dir, "scATAC_Cell_Metadata/")
if (!dir.exists(MAGICAL_cell_metadata_dir)) {dir.create(MAGICAL_cell_metadata_dir)}
MAGICAL_peak_counts_dir <- paste0(MAGICAL_file_dir, "scATAC_Read_Counts/")
if (!dir.exists(MAGICAL_peak_counts_dir)) {dir.create(MAGICAL_peak_counts_dir)}
MAGICAL_candidate_peaks_dir <- paste0(MAGICAL_file_dir, "scATAC_Peak_Coordinates/")
if (!dir.exists(MAGICAL_candidate_peaks_dir)) {dir.create(MAGICAL_candidate_peaks_dir)}
MAGICAL_motif_mapping_prior_dir <- paste0(MAGICAL_file_dir, "scATAC_Motif_Mapping_Prior/")
if (!dir.exists(MAGICAL_motif_mapping_prior_dir)) {dir.create(MAGICAL_motif_mapping_prior_dir)}
for(cell_type in unique(HVL_proj_minus_clusters$Cell_type_voting)) {
  print(cell_type)
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  # Grab cells associated with cell type
  idxPass <- which(HVL_proj_minus_clusters$Cell_type_voting %in% cell_type)
  cellsPass <- HVL_proj_minus_clusters$cellNames[idxPass]
  # Subset ArchR project and peak matrix to associated cells
  HVL_proj_cell_type_subset <- HVL_proj_minus_clusters[cellsPass,]
  # 1) Cell metadata
  HVL_proj_metadata_df <- data.frame(cell_index = seq(length(HVL_proj_cell_type_subset$cellNames)), 
                                     cell_barcode = HVL_proj_cell_type_subset$cellNames, 
                                     cell_type = HVL_proj_cell_type_subset$Cell_type_voting, 
                                     sample = HVL_proj_cell_type_subset$Sample, 
                                     condition = HVL_proj_cell_type_subset$time_point)
  write.table(HVL_proj_metadata_df, file = paste0(MAGICAL_cell_metadata_dir, cell_type_for_file_name, "_HVL_ATAC_cell_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  # 2) ATAC assay cell count
  HVL_peak_matrix_cell_type_subset <- getMatrixFromProject(ArchRProj = HVL_proj_cell_type_subset, useMatrix = "PeakMatrix", useSeqnames = NULL,
                                                verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),
                                                logFile = createLogFile("getMatrixFromProject"))
  
  final_peak_matrix <- HVL_peak_matrix_cell_type_subset@assays@data@listData$PeakMatrix
  final_peak_matrix <- final_peak_matrix[, HVL_proj_cell_type_subset$cellNames]
  saveRDS(final_peak_matrix, file= paste0(MAGICAL_peak_counts_dir, cell_type_for_file_name, "_HVL_ATAC_read_counts.rds"))
}

# 2) Peak set
current_peaks <- getPeakSet(HVL_proj_minus_clusters)
write.table(current_peaks[,1:13], file = paste0(MAGICAL_candidate_peaks_dir, "HVL_ATAC_peak_coordinates.tsv"), quote = FALSE, col.names = FALSE,  sep = "\t")

HVL_proj_minus_clusters <- addMotifAnnotations(ArchRProj = HVL_proj_minus_clusters, motifSet = "cisbp", name = "Motif",
                                                 force = TRUE)
# A Motif-Matches-In-Peaks.rds file will be created under the Annotations folder
peak_motif_mapping <- readRDS(file = paste0(ATAC_output_dir, "HVL/Annotations/Motif-Matches-In-Peaks.rds"))
write.table(lapply(summary(peak_motif_mapping@assays@data@listData$matches), as.numeric), 
            file=paste0(MAGICAL_motif_mapping_prior_dir, "HVL_ATAC_motif_mapping_prior.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = "\t")


differential_peaks_dir <- paste0(ATAC_output_dir, "diff_peaks/", date, "/new_peaks_with_fewer_tnaive/")
if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir, recursive = TRUE)}
calculate_daps_for_each_cell_type(HVL_proj_minus_clusters, differential_peaks_dir, sample_metadata_for_SPEEDI_df)

sample_pseudobulk_peak_counts <- create_pseudobulk_atac_sample(HVL_proj_minus_clusters)
pseudobulk_metadata <- sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$subject_id %in% HVL_proj_minus_clusters$subject_id,]
pseudobulk_metadata$aliquots <- rownames(pseudobulk_metadata)
pseudobulk_metadata <- pseudobulk_metadata[match(colnames(sample_pseudobulk_peak_counts), pseudobulk_metadata$aliquots),]
pseudobulk_analysis <- DESeq2::DESeqDataSetFromMatrix(countData = sample_pseudobulk_peak_counts, colData = pseudobulk_metadata, design = stats::formula(paste("~ subject_id + time_point")))
pseudobulk_analysis <- DESeq2::DESeq(pseudobulk_analysis)
pseudobulk_analysis_results_contrast <- utils::tail(DESeq2::resultsNames(pseudobulk_analysis), n=1)
print(pseudobulk_analysis_results_contrast)
pseudobulk_analysis_results <- DESeq2::results(pseudobulk_analysis, name=pseudobulk_analysis_results_contrast)
pseudobulk_analysis_results <- pseudobulk_analysis_results[rowSums(is.na(pseudobulk_analysis_results)) == 0, ] # Remove NAs
# 4505 DASs 
pseudobulk_analysis_results <- pseudobulk_analysis_results[pseudobulk_analysis_results$pvalue < 0.05,]

# Theoretically, I could create my own hg38 refseq file, but maybe not worth it

HVL_sample_metadata_for_SPEEDI_df <- sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$viral_load == "high",]
create_cell_type_proportion_MAGICAL_atac(HVL_proj_minus_clusters, "/Genomics/ogtr04/wat2/", HVL_sample_metadata_for_SPEEDI_df)


# Convert ArchR to Signac
HVL_peak_matrix <- getPeakMatrix(HVL_proj_minus_clusters)
annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38")

seurat_atac <- ArchR2Signac(
  ArchRProject = HVL_proj_minus_clusters,
  refversion = "hg38",
  #samples = samplelist, # list of samples in the ArchRProject (default will use ArchRProject@cellColData$Sample but another list can be provided)
  fragments_dir = data_path,
  pm = HVL_peak_matrix, # peak matrix from getPeakMatrix()
  fragments_fromcellranger = "Yes", # fragments_fromcellranger This is an Yes or No selection ("NO" | "N" | "No" or "YES" | "Y" | "Yes")
  fragments_file_extension = NULL, # Default - NULL: File_Extension for fragments files (typically they should be '.tsv.gz' or '.fragments.tsv.gz')
  annotation = annotations # annotation from getAnnotation()
)


# Returns TRUE as expected - samples are in correct order
all.equal(colnames(seurat_atac), gsub('#', '_', rownames(HVL_proj_minus_clusters@cellColData)))

# Add gene score matrix from ArchR to Seurat object (could recompute it using Signac, but we'll stick with ArchR matrix for now)
gsm <- getGeneScoreMatrix(ArchRProject = HVL_proj_minus_clusters, SeuratObject = seurat_atac)
seurat_atac[['RNA']] <- CreateAssayObject(counts = gsm)

# This doesn't work currently because the LSI matrix was calculated with the full set of cells, whereas this HVL project is just
# the HVL cells. So I can either recompute the LSI matrix for the HVL proj or transfer the full object from ArchR and then subdivide it to HVL cells
# in the Seurat object.
seurat_atac <- addDimRed(
  ArchRProject = HVL_proj_minus_clusters,
  SeuratObject = seurat_atac,
  addUMAPs = "UMAP",
  reducedDims = "Harmony"
)

seurat_atac$predicted_celltype_majority_vote <- HVL_proj_minus_clusters$Cell_type_voting
seurat_atac$subject_id <- HVL_proj_minus_clusters$subject_id
cell_names <- rownames(seurat_atac@meta.data)
seurat_atac <- Seurat::AddMetaData(seurat_atac, metadata = cell_names, col.name = "cell_name")

differential_peaks_dir <- paste0(ATAC_output_dir, "diff_peaks/", date, "/seurat_peaks/")
if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir, recursive = TRUE)}

cell_types <- unique(seurat_atac$predicted_celltype_majority_vote)

run_differential_expression_controlling_for_subject_id_atac(seurat_atac, differential_peaks_dir, sample_metadata_for_SPEEDI_df, "time_point", cell_types)

# Add GC content
seurat_atac <- RegionStats(seurat_atac, genome = BSgenome.Hsapiens.UCSC.hg38)
# Add motifs
seurat_atac <- AddMotifs(
  object = seurat_atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_pwms_v2
)

# saveRDS(seurat_atac, file = paste0(ATAC_output_dir, "HVL_seurat.RDS"))
# seurat_atac <- readRDS(file = paste0(ATAC_output_dir, "HVL_seurat.RDS"))


motif_output_dir <- paste0(ATAC_output_dir, "motifs/", date, "/")
if (!dir.exists(motif_output_dir)) {dir.create(motif_output_dir, recursive = TRUE)}

motif_input_dir <- "/Genomics/function/pentacon/wat2/single_cell/analysis/primary_analysis_6_subject_12_sample/ATAC/diff_peaks/2023-11-30/seurat_peaks/"

generate_motifs_with_signac(seurat_atac, motif_input_dir, motif_output_dir)

                         
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

