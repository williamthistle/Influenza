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

# Load information about viral load
viral_load_info <- read.table(paste0(home_dir_SPEEDI, "/viral_load_info.tsv"), sep = "\t", header = TRUE)

# Declare data and analysis type
data_type <- "multiome" # Can be multiome or single_cell
analysis_type <- "RNA-seq" # Can be RNA-seq or ATAC-seq

# Directory where all data are
data_path <- paste0(home_dir, data_type, "/data/")
# List of samples that can potentially be processed for the data type
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))

# data_token is used to choose subset of data that we want to analyze (pre-defined)
data_token <- "all_multiome"
# Grab samples that we want to analyze
data_tokens <- read.table(paste0(home_dir, "flu_data_tokens.tsv"), header = TRUE)
samples <- data_tokens[data_tokens$token == data_token,]$samples
sample_id_list <- all_sample_id_list %in% samples

# analysis_token is used to define a specific analysis
analysis_token <- "default"
analysis_dir <- paste0(home_dir, data_type, "/analysis/", analysis_token, "/")
if (!dir.exists(analysis_dir)) {dir.create(analysis_dir)}
plot_dir <- paste0(analysis_dir, "plots/")
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}

# Identify high / low viral load samples
high_viral_load <- c()
low_viral_load <- c()
for(sample_id in sample_id_list) {
  if(sample_id %in% viral_load_info$aliquot) {
    current_info <- viral_load_info[viral_load_info$aliquot == sample_id,]
    if(current_info$viral_load == "high") {
      high_viral_load <- c(high_viral_load, sample_id)
    } else {
      low_viral_load <- c(low_viral_load, sample_id)
    }
  } else {
    print(sample_id)
    stop("You have a sample that is not labeled as high or low viral load (printed above). Not currently supported.")
  }
}

high_viral_load <- sort(high_viral_load)
low_viral_load <- sort(low_viral_load)
all_viral_load <- c(high_viral_load, low_viral_load)
# Organize sample list by high viral load followed by low viral load
sample_id_list <- sample_id_list[order(match(sample_id_list, all_viral_load))]

# Parameters for processing both RNA-seq and ATAC data
# rerun: If files are already present, don't overwrite them (just move on) 
# save_progress: If you want to save your progress
rerun <- FALSE
save_progress <- FALSE
# Parameters for processing RNA-seq
# record_doublets: If you want to run scDblFinder and record which cells are doublets
# run_qc: If you want to plot QC metrics
# run_soup: If you want to remove ambient RNA using SoupX
record_doublets = FALSE
run_qc = FALSE
run_soup = FALSE
# Parameters for processing ATAC-seq data
# TODO
if(analysis_type == "RNA_seq") {
  # Step 1 - grab matrices
  if(run_soup) {
    all_sc_exp_matrices <- list()
    for(sample_id in sample_id_list) {
      # Load data and estimate soup profile
      sc = load10X(paste0(data_path, sample_id, "/outs/"))
      # Estimate rho
      sc = autoEstCont(sc)
      # Clean the data
      sc_exp_matrix = adjustCounts(sc)
      # Add sample name
      if (grepl("_|\\.", i)) {
        prefix <- paste0(strsplit(sample_id, "_")[[1]][1], "#")
      } else {
        prefix <- paste0(sample_id, "#")
      }
      colnames(sc_exp_matrix) <- paste0(prefix, colnames(sc_exp_matrix))
      # Add to list of all matrices
      all_sc_exp_matrices <- c(all_sc_exp_matrices, sc_exp_matrix)
    }
  } else {
    all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
  }
  # Run QC (optional)
  if(run_qc) {
    sc_obj <- CreateSeuratObject(counts = all_sc_exp_matrices,
                                 assay = "RNA",
                                 min.cells = 3,
                                 min.features = 3,
                                 project = "flu")
    sc_obj$sample <- as.vector(sapply(strsplit(colnames(sc_obj), "#"), "[", 1))
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^MT-",
                                   col.name = "percent.mt") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RPS",
                                   col.name = "percent.rps") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RPL",
                                   col.name = "percent.rpl") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^HB[A|B]",
                                   col.name = "percent.hb")
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RP[SL]",
                                   col.name = "percent.rp")
    # Adding cell_names as metadata is useful (e.g., for subsetting)
    cell_names <- rownames(sc_obj@meta.data)
    sc_obj <- AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")
    # Add viral load metadata to Seurat object
    viral_load_vec <- c()
    for(current_sample in sc_obj$sample) {
      if(current_sample %in% high_viral_load) {
        viral_load_vec <- c(viral_load_vec, "HIGH")
      } else {
        viral_load_vec <- c(viral_load_vec, "LOW")
      }
    }
    sc_obj$viral_load <- viral_load_vec
    # Make sure that sample names are in correct order (all highs together then all lows together)
    sc_obj$sample <- factor(sc_obj$sample, levels = all_viral_load)
    # Set up plotting colors
    plotting_colors <- c(rep("FC4E07", length(high_viral_load)), rep("2E9FDF", length(low_viral_load)))
    # Create violin plots for nFeature_RNA, nCount_RNA, and percent.mt
    p <- VlnPlot(sc_obj, features = c("nFeature_RNA"), split.by = "sample", group.by = "viral_load", raster = FALSE) + scale_fill_manual(values = plotting_colors) +
                                                                                                                                                               xlab("Viral Load")
    ggsave(paste0(plot_dir, "nFeature_violin_plots_", date, "_soup_", run_soup, ".png"), plot = p, device = "png", width = 10, height = 10, units = "in")
    p <- VlnPlot(sc_obj, features = c("nCount_RNA"), split.by = "sample", group.by = "viral_load", raster = FALSE) + scale_fill_manual(values = plotting_colors) +
                                                                                                                                                             xlab("Viral Load")
    
    ggsave(paste0(plot_dir, "nCount_violin_plots_", date, "_soup_", run_soup, ".png"), plot = p, device = "png", width = 10, height = 10, units = "in")
    p <- VlnPlot(sc_obj, features = c("percent.mt"), split.by = "sample", group.by = "viral_load", raster = FALSE) + scale_fill_manual(values = plotting_colors) +
                                                                                                                                                             xlab("Viral Load")
    ggsave(paste0(plot_dir, "percentMT_violin_plots_", date, "_soup_", run_soup, ".png"), plot = p, device = "png", width = 10, height = 10, units = "in")
    # Create nCount vs nFeature scatter plot for each sample
    individual_samples <- SplitObject(sc_obj, split.by = "sample")
    # Visualize nCount vs nFeature via scatter plot for each sample
    nCount_vs_nFeature_plots <- list()
    for (i in 1:length(individual_samples)) {
      p <- FeatureScatter(individual_samples[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
      nCount_vs_nFeature_plots[[i]] <- p
    }
    n <- length(nCount_vs_nFeature_plots)
    nCol <- floor(sqrt(n))
    nCount_vs_nFeature_plots <- do.call("grid.arrange", c(nCount_vs_nFeature_plots, ncol=nCol))
    ggsave(paste0(plot_dir, "nCount_vs_nFeature_plots_", date, "_soup_", run_soup, ".png"), plot = nCount_vs_nFeature_plots, device = "png", width = 15, height = 15, units = "in")
    # If you wanted to run QC, you don't want to proceed with the rest of the pipeline
  } else {
    # Step 2 - filter raw data (either with adaptive thresholds or strict thresholds)
    assign("sc_obj", FilterRawData(all_sc_exp_matrices, human = TRUE, record_doublets = FALSE, adaptive_QC_thresholds = FALSE), envir = .GlobalEnv)
    rm(all_sc_exp_matrices)
    # Adding cell_names as metadata is useful (e.g., for subsetting)
    cell_names <- rownames(sc_obj@meta.data)
    sc_obj <- AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")
    # Step 3 - normalize data
    assign("sc_obj", InitialProcessing(sc_obj, human = TRUE), envir = .GlobalEnv)
    if(save_progress) {
      save(sc_obj, file = paste0(output_dir, "3_", naming_token, "_sc_obj_full.rds"))
    }
    # Step 4 - find batches in data
    assign("sc_obj", InferBatches(sc_obj), envir = .GlobalEnv)
    # Step 5 - integrate data by batch
    assign("sc_obj", IntegrateByBatch(sc_obj), envir = .GlobalEnv)
    if(save_progress) {
      saveRDS(sc_obj, file = paste0(output_dir, "5_", naming_token, "_sc_obj.rds"))
    }
    # Step 6 - process integrated assay and potentially prepare to find markers via PrepSCTFindMarkers
    assign("sc_obj", VisualizeIntegration(sc_obj, prep_sct_find_markers = FALSE), envir = .GlobalEnv)
    if(save_progress) {
      save(sc_obj, file = paste0(output_dir, "6_", naming_token, "_sc_obj.rds"))
    }
    # Load reference and remove cell types we don't like
    reference <- LoadReference("PBMC", human = TRUE)
    idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
    reference <- reference[,-idx]
    # Map cell types from reference to our data
    assign("sc_obj", MapCellTypes(sc_obj, reference, data_type = "snRNA"), envir = .GlobalEnv)
    rm(reference)
    # We always want to save at the end
    save(sc_obj, file = paste0(output_dir, "7_", naming_token, "_sc_obj.rds"))
  }
} else if(analysis_type == "ATAC_seq") {
  
  
} else {
  stop("Invalid analysis type")
}

