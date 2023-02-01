library(clustree)
library(gri)

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
sc_obj$sample <- factor(sc_obj$sample, levels = c("3c4540710e55f7b1", "6f609a68dca1261f", "7b54cfac7e67b0fa", "b82bb7c75d47dac1", "216bb226181591dd", "abf6d19ee03be1e8", "d360f89cf9585dfe"))
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

individual_samples <- SplitObject(sc_obj, split.by = "sample")
# Visualize nFeature_RNA via violin plot
nFeature_plots <- list()
for (i in 1:length(individual_samples)) {
  p <- VlnPlot(individual_samples[[i]], features = c("nFeature_RNA"))
  nFeature_plots[[i]] <- p
}
n <- length(nFeature_plots)
nCol <- floor(sqrt(n))
nFeature_plots <- do.call("grid.arrange", c(nFeature_plots, ncol=nCol))
ggsave(paste0(output_dir, naming_token, "_nFeature_violin_plots.png"), plot = nFeature_plots, device = "png", width = 10, height = 10, units = "in")




for(current_sample in unique(sc_obj$sample)) {
  idxPass <- which(sc_obj$sample %in% sample)
  cellsPass <- names(sc_obj$orig.ident[idxPass])
  sc_obj.current.sample <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  FeatureScatter(sc_obj.current.sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
  ggsave(paste0(output_dir, naming_token, "_nCount_vs_nFeature_", sample, ".png"), device = "png", dpi = 300)
}
# High viral
idxPass <- which(sc_obj$viral_load %in% "HIGH")
cellsPass <- names(sc_obj$orig.ident[idxPass])
sc_obj.current.viral.load.high <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
FeatureScatter(sc_obj.current.viral.load.high, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
ggsave(paste0(output_dir, naming_token, "_nCount_vs_nFeature_HIGH.png"), device = "png", dpi = 300)

idxPass <- which(sc_obj$viral_load %in% "LOW")
cellsPass <- names(sc_obj$orig.ident[idxPass])
sc_obj.current.viral.load.low <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
FeatureScatter(sc_obj.current.viral.load.low, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
ggsave(paste0(output_dir, naming_token, "_nCount_vs_nFeature_LOW.png"), device = "png", dpi = 300)

# Run after Step 3
FeaturePlot(sc_obj, "nFeature_RNA")
ggsave(paste0(output_dir, naming_token, "_nFeature_QC_test.png"), device = "png", dpi = 300)
FeaturePlot(sc_obj, "nCount_RNA")
ggsave(paste0(output_dir, naming_token, "_nCount_QC_test.png"), device = "png", dpi = 300)
FeaturePlot(sc_obj, "percent.mt")
ggsave(paste0(output_dir, naming_token, "_percent.mt_test.png"), device = "png", dpi = 300)