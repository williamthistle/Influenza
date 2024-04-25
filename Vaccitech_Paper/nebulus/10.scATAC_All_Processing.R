# Load extra RNA functions
source("~/00.setup.R")
home_dir <- "/Genomics/function/pentacon/wat2/"

### SHORTCUT ###
# HVL_proj_minus_clusters <- loadArchRProject(path = paste0(ATAC_output_dir, "minus_clusters"))
# seurat_atac <- readRDS(file = paste0(ATAC_output_dir, "seurat_minus_clusters.RDS"))

################## SETUP ##################
data_path <- paste0(home_dir, "single_cell/data")
reference_dir <- paste0(home_dir, "references/")
reference_file_name <- "pbmc_multimodal.h5seurat"
output_dir <- paste0(home_dir, "single_cell/analysis/")
analysis_name <- "all_single_cell"
reference_tissue <- "pbmc_full"
reference_cell_type_attribute <- "celltype.l2"
species <- "human"
# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "all_single_cell"
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
sample_metadata_for_SPEEDI_df <- sample_metadata_for_SPEEDI_df[c("subject_id", "time_point", "sex", "viral_load_category", "treatment")]
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
atac_proj <- Read_ATAC(input_dir = data_path, sample_id_list = sample_id_list, species = species, log_flag = TRUE)
atac_proj <- FilterRawData_ATAC(proj = atac_proj, output_dir = ATAC_output_dir, log_flag = TRUE)

addArchRThreads(threads = 8)
atac_proj <- addIterativeLSI(ArchRProj = atac_proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                        clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                        varFeatures = 20000, dims = 1:30)
atac_proj <- addUMAP(ArchRProj = atac_proj, reducedDims = "IterativeLSI", force = TRUE)
atac_proj <- addClusters(input = atac_proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 2, knnAssign = 30, maxClusters = NULL, force = TRUE)

scRNA <- LoadReferenceSPEEDI(reference_tissue = reference_tissue, species = species, reference_dir = reference_dir,
                             reference_file_name = reference_file_name, log_flag = TRUE)
idx <- which(scRNA$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
scRNA <- scRNA[,-idx]

addArchRThreads(threads = 16)
atac_proj <- addGeneIntegrationMatrix_SPEEDI(
  ArchRProj = atac_proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",#Harmony
  seRNA = scRNA,
  addToArrow = FALSE,
  groupRNA = "celltype.l2",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  normalization.method = "SCT",
  force = TRUE
)

# save ArchR project: ArchR::saveArchRProject(ArchRProj = atac_proj, load = FALSE)
# load ArchR project: atac_proj <- loadArchRProject(path = paste0(ATAC_output_dir, "ArchROutput"))
# load ArchR project: atac_proj <- loadArchRProject(path = paste0(ATAC_output_dir, "no_batch_correction"))

atac_proj <- add_sample_metadata_atac(atac_proj, high_viral_load_samples, low_viral_load_samples,
                                      d28_samples, d_minus_1_samples, male_samples, female_samples, placebo_samples, vaccinated_samples)
viral_load_metadata <- parse_metadata_for_samples(atac_proj, "viral_load", high_viral_load_samples, low_viral_load_samples,
                                                  d28_samples, d_minus_1_samples, male_samples, female_samples, placebo_samples, vaccinated_samples)
day_metadata <- parse_metadata_for_samples(atac_proj, "time_point", high_viral_load_samples, low_viral_load_samples,
                                           d28_samples, d_minus_1_samples, male_samples, female_samples, placebo_samples, vaccinated_samples)
sex_metadata <- parse_metadata_for_samples(atac_proj, "sex", high_viral_load_samples, low_viral_load_samples,
                                           d28_samples, d_minus_1_samples, male_samples, female_samples, placebo_samples, vaccinated_samples)
treatment_metadata <- parse_metadata_for_samples(atac_proj, "treatment", high_viral_load_samples, low_viral_load_samples,
                                                 d28_samples, d_minus_1_samples, male_samples, female_samples, placebo_samples, vaccinated_samples)

atac_proj <- combine_cell_types_atac(atac_proj)
#atac_proj <- MajorityVote_ATAC(proj = atac_proj) - use code from some_random_atac_code

cluster_info <- get_cluster_info(atac_proj) # Make sure you replace function with one from some_random_atac_code

# Print with all clusters
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
ggplot2::ggsave(filename = paste0(ATAC_output_dir, "Final_ATAC_UMAP_by_Majority_Vote_Cell_Type_Full.png"), 
                plot = p1, device = "png", width = 8, height = 8, 
                units = "in")


# Remove messy clusters
idxPass <- which(atac_proj$Clusters %in% c("C1", "C2", "C3", "C8", "C11", "C12", "C13", "C14", "C24", "C28", "C36"))
cellsPass <- atac_proj$cellNames[-idxPass]
atac_proj_minus_clusters <- atac_proj[cellsPass, ]

# Make sure you use code from some_random_atac_code
atac_proj_minus_clusters <- override_cluster_label_atac(atac_proj_minus_clusters, c("C16"), "CD8 Memory")
atac_proj_minus_clusters <- override_cluster_label_atac(atac_proj_minus_clusters, c("C20"), "Proliferating")
atac_proj_minus_clusters <- override_cluster_label_atac(atac_proj_minus_clusters, c("C25", "C34"), "CD4 Naive")

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

# Call peaks
addArchRGenome("hg38")
atac_proj_minus_clusters <- pseudo_bulk_replicates_and_call_peaks(atac_proj_minus_clusters)
# Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
# calculations
atac_proj_minus_clusters <- addPeakMatrix(atac_proj_minus_clusters)

# Save peak metadata
peaks <- getPeakSet(atac_proj_minus_clusters)
peaks_df <- as.data.frame(peaks@seqnames)
peaks_df <- cbind(peaks_df, as.data.frame(peaks@ranges))
peaks_df <- cbind(peaks_df, as.data.frame(peaks@elementMetadata))
write.table(x = peaks_df, file = paste0(ATAC_output_dir, "peaks_info.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# Create Peaks.txt file
peak_txt_file <- create_peaks_file(atac_proj_minus_clusters, ATAC_output_dir)
# Create peak_motif_matches.txt file
atac_proj_minus_clusters <- create_peak_motif_matches_file(atac_proj_minus_clusters, ATAC_output_dir, peak_txt_file)
# save ArchR project: ArchR::saveArchRProject(ArchRProj = atac_proj_minus_clusters, outputDirectory = paste0(ATAC_output_dir, "minus_clusters_2"), load = FALSE, overwrite = TRUE)
# load ArchR project: atac_proj_minus_clusters <- loadArchRProject(path = paste0(ATAC_output_dir, "minus_clusters_2"))
# Find DASs
sample_metadata <- atac_proj_minus_clusters$Sample
for(j in 1:nrow(sample_metadata_for_SPEEDI_df)) {
  sample_metadata <- gsub(rownames(sample_metadata_for_SPEEDI_df)[j], sample_metadata_for_SPEEDI_df[j,1], sample_metadata)
}
atac_proj_minus_clusters <- addCellColData(ArchRProj = atac_proj_minus_clusters, data = sample_metadata, cells = atac_proj_minus_clusters$cellNames, name = "subject_id", force = TRUE)

# Create MAGICAL files for HVL placebo
idxPass <- which(atac_proj_minus_clusters$treatment %in% "PLACEBO")
cellsPass <- atac_proj_minus_clusters$cellNames[idxPass]
placebo_atac_proj <- atac_proj_minus_clusters[cellsPass,]

idxPass <- which(placebo_atac_proj$viral_load %in% "high")
cellsPass <- placebo_atac_proj$cellNames[idxPass]
placebo_hvl_atac_proj <- placebo_atac_proj[cellsPass,]

create_magical_input_files_atac(placebo_hvl_atac_proj, paste0(ATAC_output_dir, "MAGICAL_HVL_PLACEBO_", date, "/"))

# Record cell type proportions for each subset - TODO
# HVL_sample_metadata_for_SPEEDI_df <- sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$viral_load == "low",]
# create_cell_type_proportion_MAGICAL_atac(HVL_proj_minus_clusters, "/Genomics/ogtr04/wat2/", HVL_sample_metadata_for_SPEEDI_df)

# Convert ArchR to Signac
peak_matrix <- getPeakMatrix(atac_proj_minus_clusters)
annotations <- ArchRtoSignac::getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38")

seurat_atac <- ArchRtoSignac::ArchR2Signac(
  ArchRProject = atac_proj_minus_clusters,
  refversion = "hg38",
  #samples = samplelist, # list of samples in the ArchRProject (default will use ArchRProject@cellColData$Sample but another list can be provided)
  fragments_dir = data_path,
  pm = peak_matrix, # peak matrix from getPeakMatrix()
  fragments_fromcellranger = "Yes", # fragments_fromcellranger This is an Yes or No selection ("NO" | "N" | "No" or "YES" | "Y" | "Yes")
  fragments_file_extension = NULL, # Default - NULL: File_Extension for fragments files (typically they should be '.tsv.gz' or '.fragments.tsv.gz')
  annotation = annotations # annotation from getAnnotation()
)


# Returns TRUE as expected - samples are in correct order
all.equal(colnames(seurat_atac), gsub('#', '_', rownames(atac_proj_minus_clusters@cellColData)))

# Add gene score matrix from ArchR to Seurat object (could recompute it using Signac, but we'll stick with ArchR matrix for now)
gsm <- getGeneScoreMatrix(ArchRProject = atac_proj_minus_clusters, SeuratObject = seurat_atac)
seurat_atac[['RNA']] <- CreateAssayObject(counts = gsm)

# This doesn't work currently because the LSI matrix was calculated with the full set of cells, whereas this HVL project is just
# the HVL cells. So I can either recompute the LSI matrix for the HVL proj or transfer the full object from ArchR and then subdivide it to HVL cells
# in the Seurat object.
#seurat_atac <- addDimRed(
#  ArchRProject = HVL_proj_minus_clusters,
#  SeuratObject = seurat_atac,
#  addUMAPs = "UMAP",
#  reducedDims = "Harmony"
#)

seurat_atac$predicted_celltype_majority_vote <- atac_proj_minus_clusters$Cell_type_voting
seurat_atac$subject_id <- atac_proj_minus_clusters$subject_id
cell_names <- rownames(seurat_atac@meta.data)
seurat_atac <- Seurat::AddMetaData(seurat_atac, metadata = cell_names, col.name = "cell_name")

# saveRDS(seurat_atac, file = paste0(ATAC_output_dir, "seurat_minus_clusters.RDS"))
# seurat_atac <- readRDS(file =  paste0(ATAC_output_dir, "seurat_minus_clusters.RDS"))

# Separate into relevant subsets

sc_obj <- seurat_atac

# HVL
idxPass <- which(sc_obj$viral_load %in% "high")
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
idxPass <- which(sc_obj$viral_load %in% "low")
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

run_differential_expression_controlling_for_subject_id_atac(hvl_placebo_sc_obj, paste0(ATAC_output_dir, "DE_HVL_PLACEBO_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point")
run_differential_expression_controlling_for_subject_id_atac(lvl_placebo_sc_obj, paste0(ATAC_output_dir, "DE_LVL_PLACEBO_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point")
run_differential_expression_controlling_for_subject_id_atac(hvl_vaccinated_sc_obj, paste0(ATAC_output_dir, "DE_HVL_VACCINATED_", date, "/"), sample_metadata_for_SPEEDI_df, "time_point")

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
