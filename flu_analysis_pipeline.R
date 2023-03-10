library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(clustree)
library(dplyr)
library(ggplot2)
library(hexbin)
library(mclust)
library(openxlsx)
library(outliers)
library(pheatmap)
library(scDblFinder)
library(SeuratDisk)
library(BiocParallel)
library(stringr)
library(writexl)
library(SoupX)

################## SETUP ##################
date <- Sys.Date()
home_dir <- "~/"

# Load SPEEDI (for RNA-seq analyses)
SPEEDI_dir <- paste0(home_dir, "SPEEDI")
source(paste0(SPEEDI_dir, "/prototype_API.R"))

# Load information about samples
sample_metadata <- read.table(paste0(SPEEDI_dir, "/sample_metadata.tsv"), sep = "\t", header = TRUE)

# Declare data and analysis type
data_type <- "multiome" # Can be multiome or single_cell
analysis_type <- "RNA-seq" # Can be RNA-seq or ATAC-seq

# Directory where all data are for data_type specified above
data_path <- paste0(home_dir, data_type, "/data/")
# List of samples that can potentially be processed for the data type
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))

# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "female_hvl_vs_lvl_multiome"
# Create directory for this particular data in the analysis directory if it doesn't exist
base_analysis_dir <- paste0(home_dir, data_type, "/analysis/", data_token, "/")
if (!dir.exists(base_analysis_dir)) {dir.create(base_analysis_dir)}
# Grab samples that we want to analyze
data_tokens <- read.table(paste0(home_dir, "flu_data_tokens.tsv"), header = TRUE)
samples <- data_tokens[data_tokens$token == data_token,]$samples
samples <-  unlist(strsplit(samples, ","))
sample_id_list <- all_sample_id_list[all_sample_id_list %in% samples]
if(length(sample_id_list) != length(samples)) {
  print(paste0("Requested samples: ", samples))
  print(paste0("Available samples (from your list): ", sample_id_list))
  stop("You have at least one invalid sample in your list")
}
sample_count <- length(sample_id_list)

# analysis_token is used to define a specific analysis
analysis_token <- "default"
analysis_dir <- paste0(home_dir, data_type, "/analysis/", data_token, "/", analysis_token, "/")
if (!dir.exists(analysis_dir)) {dir.create(analysis_dir)}
# Directory for analysis plots
plot_dir <- paste0(analysis_dir, "plots/")
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}

# Identify high / low viral load samples, D28 / D-1 samples, and M / F samples
# Subset metadata to only contain our sample_ids
# We will add these metadata to our Seurat object later
sample_metadata <- sample_metadata[sample_metadata$aliquot %in% sample_id_list,]
high_viral_load_samples <- sort(sample_metadata[sample_metadata$viral_load == "high",]$aliquot)
low_viral_load_samples <- sort(sample_metadata[sample_metadata$viral_load == "low",]$aliquot)
all_viral_load_samples <- c(high_viral_load_samples, low_viral_load_samples)
d28_samples <- sort(sample_metadata[sample_metadata$viral_load == "2_D28",]$aliquot)
d_minus_1_samples <- sort(sample_metadata[sample_metadata$viral_load == "2_D_minus_1",]$aliquot)
all_day_samples <- c(d28_samples, d_minus_1_samples)
male_samples <- sort(sample_metadata[sample_metadata$sex == "M",]$aliquot)
female_samples <- sort(sample_metadata[sample_metadata$sex == "F",]$aliquot)
all_sex_samples <- c(male_samples, female_samples)

# Organize sample list by high viral load followed by low viral load
sample_id_list <- sample_id_list[order(match(sample_id_list, all_viral_load_samples))]

# Parameters for processing both RNA-seq and ATAC data
# save_progress: If you want to save your progress
save_progress <- FALSE
# Parameters for processing RNA-seq
# record_doublets: If you want to run scDblFinder and record which cells are doublets
# run_qc: If you want to plot QC metrics (and not run the rest of the pipeline)
# run_soup: If you want to remove ambient RNA using SoupX
record_doublets = FALSE
run_qc = FALSE
run_soup = FALSE
# Parameters for processing ATAC-seq data
# TODO
################## ANALYSIS ##################
if(analysis_type == "RNA_seq") {
  # Step 1 - grab matrices
  if(run_soup) {
    # TODO: Can I make this parallel?
    all_sc_exp_matrices <- process_matrices_through_soup(data_path, sample_id_list)
  } else {
    all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
  }
  # Run QC (optional)
  if(run_qc) {
    sc_obj <- generate_qc_plots(all_sc_exp_matrices, plot_dir, run_soup, date, high_viral_load_samples, low_viral_load_samples,
                      d28_samples, d_minus_1_samples, male_samples, female_samples)
    # If you run QC, you probably don't want to proceed with the rest of the pipeline, so we stop here
  } else {
    # Step 2 - filter raw data (either with adaptive thresholds or strict thresholds)
    sc_obj <- FilterRawData(all_sc_exp_matrices, human = TRUE, record_doublets = FALSE, adaptive_QC_thresholds = FALSE)
    rm(all_sc_exp_matrices)
    sc_obj <- add_sample_metadata(sc_obj, high_viral_load_samples, low_viral_load_samples,
                        d28_samples, d_minus_1_samples, male_samples, female_samples)
    # Step 3 - normalize data
    sc_obj <- InitialProcessing(sc_obj, human = TRUE)
    if(save_progress) {
      save(sc_obj, file = paste0(analysis_dir, "3_sc_obj.rds"))
    }
    # Step 4 - find batches in data
    sc_obj <- InferBatches(sc_obj)
    # Step 5 - integrate data by batch
    sc_obj <- IntegrateByBatch(sc_obj)
    if(save_progress) {
      save(sc_obj, file = paste0(analysis_dir, "5_sc_obj.rds"))
    }
    # Step 6 - process integrated assay and potentially prepare to find markers via PrepSCTFindMarkers
    sc_obj <- VisualizeIntegration(sc_obj, prep_sct_find_markers = FALSE)
    if(save_progress) {
      save(sc_obj, file = paste0(analysis_dir, "6_sc_obj.rds"))
    }
    # Load reference and remove cell types we don't like
    reference <- LoadReference("PBMC", human = TRUE)
    idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
    reference <- reference[,-idx]
    # Map cell types from reference to our data
    sc_obj <- MapCellTypes(sc_obj, reference, data_type = "snRNA")
    rm(reference)
    # Combine cell types and re-do majority vote
    sc_obj <- combine_cell_types(sc_obj, resolution = 1.5)
    # Print UMAP by cell type (majority vote) and by cluster number - it will currently be messy
    print_UMAP(sc_obj, sample_count, "predicted_celltype_majority_vote", plot_dir, paste0("clusters_by_cell_type_majority_vote_", date, ".png"))
    print_UMAP(sc_obj, sample_count, "seurat_clusters", plot_dir, paste0("clusters_by_cluster_num_", date, ".png"))
    print_UMAP(sc_obj, sample_count, "predicted.id", plot_dir, paste0("clusters_by_cell_type_", date, ".png"))
    # We always want to save our final sc_obj the end
    save(sc_obj, file = paste0(analysis_dir, "7_sc_obj.rds"))
  }
} else if(analysis_type == "ATAC_seq") {
  
  
} else {
  stop("Invalid analysis type")
}

