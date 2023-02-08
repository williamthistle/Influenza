library(clustree)
library(gridExtra)
library(scDblFinder)

home_dir <- "~/SPEEDI"
source(paste0(home_dir, "/prototype_API.R"))

# naming_token is used to name output files
naming_token <- "high_vs_low_viral_load_D28"
#naming_token <- "all_high_vs_low_viral_load_D28_reduced_references"

# data_path is where input data are stored
data_path <- paste0("/data/home/wat2/", naming_token, "/RNA_seq_data/")
# Get list of samples that will be processed
sample_id_list <- list.dirs(data_path, recursive = FALSE)
sample_id_list <- strsplit(sample_id_list, "/")
sample_id_list <- unlist(lapply(sample_id_list, tail, n = 1L))
sample_count <- length(sample_id_list)
output_dir <- paste0("/data/home/wat2/", naming_token, "/RNA_seq_data_output/")

# Read in files and get necessary stats for QC visualization
all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
sc_obj <- CreateSeuratObject(counts = all_sc_exp_matrices,
                             assay = "RNA",
                             min.cells = 3,
                             min.features = 3,
                             project = "unbias")
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
# Creating cell_names column is useful for subsetting sc_obj
cell_names <- rownames(sc_obj@meta.data)
sc_obj <- AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")
# Assign viral load (high or low) to each cell - these IDs need to be updated depending on your analysis
high_viral_load <- c("b82bb7c75d47dac1", "3c4540710e55f7b1", "6f609a68dca1261f", "7b54cfac7e67b0fa")
low_viral_load <- c("abf6d19ee03be1e8", "216bb226181591dd", "d360f89cf9585dfe")
all_viral_load <- c(high_viral_load, low_viral_load)
viral_load_label <- c(rep("HIGH", length(high_viral_load)), rep("LOW", length(low_viral_load)))
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
sc_obj$sample <- factor(sc_obj$sample, levels = c("3c4540710e55f7b1", "6f609a68dca1261f", "7b54cfac7e67b0fa", "b82bb7c75d47dac1", "216bb226181591dd", "abf6d19ee03be1e8", "d360f89cf9585dfe"))
# Create violin plots for nFeature_RNA, nCount_RNA, and percent.mt
p <- VlnPlot(sc_obj, features = c("nFeature_RNA"), split.by = "sample", group.by = "viral_load", raster = FALSE) + scale_fill_manual(values = c("#FC4E07", "#FC4E07", "#FC4E07", "#FC4E07", 
                                                                                                                                                          "#2E9FDF", "#2E9FDF", "#2E9FDF")) +
                                                                                                                                                           xlab("Viral Load")
ggsave(paste0(output_dir, naming_token, "_nFeature_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
p <- VlnPlot(sc_obj, features = c("nCount_RNA"), split.by = "sample", group.by = "viral_load", raster = FALSE) + scale_fill_manual(values = c("#FC4E07", "#FC4E07", "#FC4E07", "#FC4E07", 
                                                                                                                                                         "#2E9FDF", "#2E9FDF", "#2E9FDF")) +
                                                                                                                                                           xlab("Viral Load")
  
ggsave(paste0(output_dir, naming_token, "_nCount_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
p <- VlnPlot(sc_obj, features = c("percent.mt"), split.by = "sample", group.by = "viral_load", raster = FALSE) + scale_fill_manual(values = c("#FC4E07", "#FC4E07", "#FC4E07", "#FC4E07", 
                                                                                                                                                         "#2E9FDF", "#2E9FDF", "#2E9FDF")) +
                                                                                                                                                           xlab("Viral Load")
ggsave(paste0(output_dir, naming_token, "percent.mt_violin_plots.png"), plot = p, device = "png", width = 10, height = 10, units = "in")
# Create nCount vs nFeature scatter plot for each sample
individual_samples <- SplitObject(sc_obj, split.by = "sample")
# Visualize nCount vs nFeature via scatter plot for each sample - THIS WORKS
nCount_vs_nFeature_plots <- list()
for (i in 1:length(individual_samples)) {
  p <- FeatureScatter(individual_samples[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
  nCount_vs_nFeature_plots[[i]] <- p
}
n <- length(nCount_vs_nFeature_plots)
nCol <- floor(sqrt(n))
nCount_vs_nFeature_plots <- do.call("grid.arrange", c(nCount_vs_nFeature_plots, ncol=nCol))
ggsave(paste0(output_dir, naming_token, "_nCount_vs_nFeature_plots.png"), plot = nCount_vs_nFeature_plots, device = "png", width = 15, height = 15, units = "in")

# RUN AFTER STEP 3 - see what unintegrated plots look like
cell_count <- length(sc_obj$cell_name)
current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
DimPlot(sc_obj, reduction = "umap", group.by = "sample", split.by = "sample", ncol = 3, repel = TRUE, raster = FALSE) + 
  labs(title = current_title) + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(output_dir, naming_token, "_preintegrated_clusters_by_sample.png"), device = "png", dpi = 300, width = 8, height = 8, units = "in")

# Find doublets
sce <- scDblFinder(as.SingleCellExperiment(sc_obj), samples = "sample")
sc_obj <- as.Seurat(sce)
# See distribution of doublets in each sample
doublet_sc_obj <- subset(x = sc_obj, subset = scDblFinder.class %in% "doublet")
print(table(doublet_sc_obj$sample))
# Remove doublets
sc_obj <- subset(x = sc_obj, subset = scDblFinder.class %in% "singlet")


