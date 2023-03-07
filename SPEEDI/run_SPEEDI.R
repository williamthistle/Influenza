library(clustree)
library(pheatmap)
library(scDblFinder)
library(BiocParallel)

home_dir <- "~/SPEEDI"
source(paste0(home_dir, "/prototype_API.R"))

# naming_token is used to select analysis and name output files
naming_token <- "high_vs_low_viral_load_D28_V3"
date <- Sys.Date()

# data_path is where input data are stored
data_path <- paste0("/data/home/wat2/", naming_token, "/RNA_seq_data/")
# Get list of samples that will be processed
sample_id_list <- list.dirs(data_path, recursive = FALSE)
sample_id_list <- strsplit(sample_id_list, "/")
sample_id_list <- unlist(lapply(sample_id_list, tail, n = 1L))
sample_count <- length(sample_id_list)
output_dir <- paste0("/data/home/wat2/", naming_token, "/RNA_seq_data_output/")
viral_load_info <- read.table(paste0(home_dir, "/viral_load_info.tsv"), sep = "\t", header = TRUE)

# Run SPEEDI on our list of samples
save_progress <- TRUE
sc_obj <- run_SPEEDI(data_path, output_dir, sample_id_list, naming_token, save_progress, use_simplified_reference = FALSE)
if(!save_progress) {
  save.image(paste0(output_dir, "7_", naming_token, ".RData"))
}

#load(paste0(output_dir, "7_", naming_token, ".RData"))
#load(paste0(output_dir, "7_", naming_token, "_sc_obj_scaled.rds"))

# Add cell names to Seurat object - only needed if SPEEDI functions are used directly (as opposed to run_SPEEDI wrapper)
#cell_names <- rownames(sc_obj@meta.data)
#sc_obj <- AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")

# Remove cell types which we don't believe in
idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
reference <- reference[,-idx]

# Assign metadata to each cell
high_viral_load <- c()
low_viral_load <- c()
for(sample_id in unique(sc_obj$sample)) {
  if(sample_id %in% viral_load_info$aliquot) {
    current_info <- viral_load_info[viral_load_info$aliquot == sample_id,]
    if(current_info$viral_load == "high") {
      high_viral_load <- c(high_viral_load, sample_id)
    } else {
      low_viral_load <- c(low_viral_load, sample_id)
    }
  }
}

high_viral_load <- sort(high_viral_load)
low_viral_load <- sort(low_viral_load)

all_viral_load <- c(high_viral_load, low_viral_load)
viral_load_label <- c(rep("HIGH", length(high_viral_load)), rep("LOW", length(low_viral_load)))
viral_load_vec <- c()
for(current_sample in sc_obj$sample) {
  if(current_sample %in% high_viral_load) {
    viral_load_vec <- c(viral_load_vec, "HIGH")
  } else if(current_sample %in% low_viral_load) {
    viral_load_vec <- c(viral_load_vec, "LOW")
  }
}
sc_obj$viral_load <- viral_load_vec

# Combine cell types and re-do majority vote
sc_obj$old.predicted.id <- sc_obj$predicted.id
Cell_type_combined = sc_obj$predicted.id
idx <- grep("CD4 T", Cell_type_combined)
Cell_type_combined[idx] <- "CD4 Memory"
idx <- grep("CD8 T", Cell_type_combined)
Cell_type_combined[idx] <- "CD8 Memory"
idx <- grep("cDC", Cell_type_combined)
Cell_type_combined[idx] <- "cDC"
idx <- grep("Proliferating", Cell_type_combined)
Cell_type_combined[idx] <- "Proliferating"
sc_obj$predicted.id <- Cell_type_combined
sc_obj <- MajorityVote(sc_obj, 2.4)



# Print UMAP by cell type (majority vote) and by cluster number - it will currently be messy
print_UMAP(sc_obj, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_majority_vote_doublets_intact_", date, ".png"))
print_UMAP(sc_obj, sample_count, "seurat_clusters", output_dir, naming_token, paste0("_clusters_by_cluster_num_doublets_intact_", date, ".png"))
print_UMAP(sc_obj, sample_count, "predicted.id", output_dir, naming_token, paste0("_clusters_by_cell_type_doublets_intact_", date, ".png"))

# Split by cluster
cell_count <- length(sc_obj$cell_name)
current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
DimPlot(sc_obj, reduction = "umap", group.by = "predicted.id", split.by = "seurat_clusters", ncol = 5, repel = TRUE, raster = FALSE) +
  labs(title = current_title) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(output_dir, naming_token, "_clusters_by_cluster_", date, ".png"), device = "png", dpi = 300, width = 20, height = 20, units = "in")

# Print viral load plot
print_UMAP(sc_obj, sample_count, "viral_load", output_dir, naming_token, paste0("_clusters_by_viral_load_", date, ".png"))

# See overlap between ATAC and RNA-seq
filtered_ATAC_cells <- read.table(paste0(output_dir, "filtered_ATAC_cells.txt"), comment.char = "")
filtered_ATAC_cells <- filtered_ATAC_cells$V1


# Let's see what the cell type proportions are for our raw data (using predicted.id)
# We aren't that concerned about what the majority vote says right now because we still need to clean up our clusters
raw_cell_tables <- calculate_props_and_counts(sc_obj, viral_load_label, all_viral_load)
write.csv(raw_cell_tables[1], file = paste0(output_dir, naming_token, "_raw_RNA_cell_type_proportion_", date, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(raw_cell_tables[2], file = paste0(output_dir, naming_token, "_raw_RNA_cell_counts_", date, ".csv"), quote = FALSE, row.names = FALSE)

# Let's also use Vincy's code to create a nice heatmap of cell type proportions
sc_obj$sample <- factor(sc_obj$sample, levels = all_viral_load)
voting_cell_type_proportion  <- as.matrix(table(sc_obj$sample, sc_obj$predicted_celltype_majority_vote))
voting_cell_type_proportion <- apply(voting_cell_type_proportion, 1, function(x){x/sum(x)})
# Combine B Cells
rownames(raw_cell_type_proportion)[2] <- "B Cells"
raw_cell_type_proportion[2,] <- raw_cell_type_proportion[2,] + raw_cell_type_proportion[3,] + raw_cell_type_proportion[4,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(3, 4), ]
# Combine Monocytes
rownames(raw_cell_type_proportion)[3] <- "Monocytes"
raw_cell_type_proportion[3,] <- raw_cell_type_proportion[3,] + raw_cell_type_proportion[4,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(4), ]
# Combine T Cells
rownames(raw_cell_type_proportion)[4] <- "T Cells"
raw_cell_type_proportion[4,] <- raw_cell_type_proportion[4,] + raw_cell_type_proportion[6,] + raw_cell_type_proportion[7,] + raw_cell_type_proportion[8,] + raw_cell_type_proportion[10,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(6, 7, 8, 10), ]
# Combine Rare T Cells
rownames(raw_cell_type_proportion)[5] <- "Rare T Cells"
raw_cell_type_proportion[5,] <- raw_cell_type_proportion[5,] + raw_cell_type_proportion[6,] + raw_cell_type_proportion[8,] + raw_cell_type_proportion[10,] + raw_cell_type_proportion[19,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(6, 8, 10, 19), ]
# Combine NK / ILC
rownames(raw_cell_type_proportion)[9] <- "NK / ILC"
raw_cell_type_proportion[9,] <- raw_cell_type_proportion[9,] + raw_cell_type_proportion[11,] + raw_cell_type_proportion[12,]
raw_cell_type_proportion <- raw_cell_type_proportion[-c(11, 12), ]

my_sample_col <- data.frame(sample = rep(c("HIGH", "LOW"), c(length(high_viral_load),length(low_viral_load)))) #The study group by their order
row.names(my_sample_col) <- colnames(voting_cell_type_proportion)

output.plot <- pheatmap(voting_cell_type_proportion, annotation_col = my_sample_col, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.3f")
ggsave(paste0(output_dir, naming_token, "_voting_cell_type_proportions_heatmap_", date, ".png"), plot = output.plot, device = "png", width = 8, height = 8, units = "in")

# Let's use clustree to try to figure out the best clustering resolution
for (res in seq(0, 6, 0.3)) {
  sc_obj <- FindClusters(sc_obj, resolution = res)
}
clustree(sc_obj, prefix = "integrated_snn_res.")
ggsave(paste0(output_dir, naming_token, "_cluster.trees_", date, ".png"), device = "png", width = 20, height = 20, units = "in")

# Now, we should re-run our majority vote with the best resolution
best_res <- 2.4
sc_obj <- MajorityVote(sc_obj, best_res)

# To decide which clusters we need to remove, we will use two tactics:
# 1. Capture information about cell type distributions in each cluster (and cluster mean S score, G2M score, and CC diff for cell cycle info)
#    High S score and high G2M score seem to indicate Proliferating cluster
# 1
raw_cluster_info <- capture_cluster_info(sc_obj)

# Remove messy clusters
#messy_clusters <- c(3, 21, 27) # For multiome F - previous non-scaled
messy_clusters <- c(4, 20, 29, 36) # For multiome F - scaled
#messy_clusters <- c(0, 2, 21, 22, 34) # for multiome F - doublets intact
#messy_clusters <- c(9, 29, 40, 45, 47, 49, 50, 58, 60) # For multiome + scRNA-seq
idxPass <- which(Idents(sc_obj) %in% messy_clusters)
cellsPass <- names(sc_obj$orig.ident[-idxPass])
sc_obj.minus.messy.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_without_messy_clusters_doublets_intact_", date, ".png"))

# Get list of cells and predicted cell types to use for snATAC-seq
# First, we write the cells and predicted cell types from the uncurated (full set of cells)
load("~/high_vs_low_viral_load_D28_V3/RNA_seq_data_output/7_high_vs_low_viral_load_D28_V3_sc_obj_full.rds")
cells_for_ATAC <- data.frame("cells" = sc_obj$cell_name, voted_type = sc_obj$predicted.id)
write.csv(cells_for_ATAC, file = paste0(output_dir, naming_token, "_final_cell_names_uncurated_", date, ".csv"), quote = FALSE, row.names = FALSE)
# Next, we write the cells and predicted cell types from the curated (filtered snRNA-seq set of cells)
load("~/high_vs_low_viral_load_D28_V3/RNA_seq_data_output/7_high_vs_low_viral_load_D28_V3_sc_obj_scaled_minus_messy_clusters.rds")
cells_for_ATAC <- data.frame("cells" = sc_obj$cell_name, voted_type = sc_obj$predicted_celltype_majority_vote)
write.csv(cells_for_ATAC, file = paste0(output_dir, naming_token, "_final_cell_names_curated_", date, ".csv"), quote = FALSE, row.names = FALSE)





print("Performing differential expression between groups (HIGH and LOW) for each cell type")
for (cell_type in unique(sc_obj$predicted_celltype_majority_vote)) {
  print(cell_type)
  idxPass <- which(sc_obj$predicted_celltype_majority_vote %in% cell_type)
  print(length(idxPass))
  cellsPass <- names(sc_obj$orig.ident[idxPass])
  cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  print(table(cells_subset$viral_load))
  DefaultAssay(cells_subset) <- "SCT"
  Idents(cells_subset) <- "viral_load"
  diff_markers <- FindMarkers(cells_subset, ident.1 = "HIGH", ident.2 = "LOW", assay = "SCT", recorrect_umi = FALSE, logfc.threshold = 0, min.pct = 0)
  #diff_markers$p_val_adj = p.adjust(diff_markers$p_val, method='fdr')
  #diff_markers <- diff_markers[diff_markers$avg_log2FC > 0.1 | diff_markers$avg_log2FC < -0.1,]
  #diff_markers <- diff_markers[diff_markers$p_val_adj < 0.05,]
  print(nrow(diff_markers))
  cell_type <- sub(" ", "_", cell_type)
  write.csv(diff_markers, paste0(output_dir, "HIGH-vs-LOW-degs-", cell_type, ".csv"), quote = FALSE)
}




### RANDOM CODE WHICH MAY NOT BE NECESSARY TO USE ###

# Demo of how to get more UMAP axis ticks
test <- DimPlot(sc_obj.minus.messy.clusters, reduction = "umap", group.by = "predicted_celltype_majority_vote", label = TRUE,
                label.size = 3, repel = TRUE, raster = FALSE)
test <- test + scale_x_continuous(breaks = scales::pretty_breaks(n = 20)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))
ggsave(paste0(output_dir, "test.png"), device = "png", dpi = 300)

# To confirm the cell type associated with a given cluster, we can also find the cluster markers (example below)
cluster15.markers <- FindMarkers(sc_obj, ident.1 = 15, min.pct = 0.25)

# Code to compare weird cluster to normal cluster
idxPass <- which(Idents(sc_obj) %in% c(15))
cellsPass <- names(sc_obj$orig.ident[idxPass])
sc_obj.cluster.15 <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
idxPass <- which(Idents(sc_obj) %in% c(1))
cellsPass <- names(sc_obj$orig.ident[idxPass])
sc_obj.cluster.1 <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

# Manually override cluster labels
NK_cluster <- c(17)
for(cluster_id in NK_cluster) {
  print(cluster_id)
  majority_vote <- sc_obj.minus.messy.clusters$predicted_celltype_majority_vote
  idxPass <- which(Idents(sc_obj.minus.messy.clusters) %in% c(cluster_id))
  majority_vote[idxPass] <- "NK"
  sc_obj.minus.messy.clusters$predicted_celltype_majority_vote <- majority_vote
}

CD8_memory_cluster <- c(9)
for(cluster_id in CD8_memory_cluster) {
  print(cluster_id)
  majority_vote <- sc_obj.minus.messy.clusters$predicted_celltype_majority_vote
  idxPass <- which(Idents(sc_obj.minus.messy.clusters) %in% c(cluster_id))
  majority_vote[idxPass] <- "CD8 Memory"
  sc_obj.minus.messy.clusters$predicted_celltype_majority_vote <- majority_vote
}


# Select cells based on UMAP coordinates
umap.coords <- as.data.frame(sc_obj[["umap"]]@cell.embeddings)
umap.coords <- umap.coords[umap.coords$"UMAP_1" > 0 & umap.coords$"UMAP_1" < 2,]
umap.coords <- umap.coords[umap.coords$"UMAP_2" > -3 & umap.coords$"UMAP_2" < 0,]
cellsPass <- rownames(umap.coords)
sc_obj.bridge.cluster <- subset(x = sc_obj, subset = cell_name %in% cellsPass)

# Run SingleR cell annotation
gene_expression <- sc_obj@assays$RNA@data[-c(grep("^MT-", row.names(sc_obj@assays$RNA@data)), grep("^RP[SL]", row.names(sc_obj@assays$RNA@data))),]
singler.fine <- SingleR(test = gene_expression,
                        ref = reference@assays$SCT@data, 
                        labels = reference$celltype.l2,
                        aggr.ref=TRUE,
                        BPPARAM=MulticoreParam(7))
labels <- singler.fine$pruned.labels
names(labels) <- row.names(singler.fine)
sc_obj$SingleR_labels <- labels
save(sc_obj, file = paste0(output_dir, "7_", naming_token, "_sc_obj_with_singleR.rds"))

# Look at platelets
idxPass <- which(sc_obj$predicted.id.with.platelets %in% "Platelet")
cellsPass <- names(sc_obj$orig.ident[idxPass])
cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
table(cells_subset$SingleR_labels)
# Can also remove platelets from reference, and we see vast majority are labeled as CD14 monocytes

sc_obj$predicted.id <- sc_obj$SingleR_labels
sc_obj$predicted.id <- ifelse(is.na(sc_obj$predicted.id), "NA", sc_obj$predicted.id)


for(cell_type in unique(sc_obj$predicted.id)) {
  print(cell_type)
  idxPass <- which(sc_obj$predicted.id %in% cell_type)
  cellsPass <- names(sc_obj$orig.ident[idxPass])
  cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  print(median(cells_subset$predicted.id.score))
}


ggsave(paste0(output_dir, naming_token, "_query_and_reference_overlayed_", date, ".png"), device = "png", dpi = 300, width = 20, height = 20, units = "in")

# Find markers for each cluster
cluster_ids <- unique(sc_obj$seurat_clusters)
for(cluster_id in cluster_ids) {
  print(cluster_id)
  cluster.markers <- FindMarkers(sc_obj, ident.1 = cluster_id, assay = "SCT")
  write.table(cluster.markers, paste0(output_dir, "7_", naming_token, "_", cluster_id, ".txt"), quote = FALSE, sep = "\t")
}

# Alternative
cluster.markers <- FindAllMarkers(sc_obj, assay = "SCT")
write.table(cluster.markers, paste0(output_dir, "7_", naming_token, "_cluster_markers_scaled.txt"), quote = FALSE, sep = "\t")





# The ideal resolution may have changed. Let's check now!
for (res in seq(0, 3, 0.3)) {
  sc_obj.minus.messy.clusters <- FindClusters(sc_obj.minus.messy.clusters, resolution = res)
}
clustree(sc_obj.minus.messy.clusters, prefix = "integrated_snn_res.")
ggsave(paste0(output_dir, naming_token, "_cluster.trees.minus.messy.clusters.PNG"), device = "png", width = 20, height = 20, units = "in")

# Redo majority vote with updated resolution - here, 1.2 looks the best
sc_obj.minus.messy.clusters <- MajorityVote(sc_obj.minus.messy.clusters)
# Print cell type plot without messy clusters
cell_count <- length(sc_obj.minus.messy.clusters$cell_name)
current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, "_clusters_by_cell_type_without_messy_clusters.png")

# Update cluster predictions and produce new plots
processed_cluster_info <- capture_cluster_info(sc_obj.minus.messy.clusters)
#2
for(cluster_id in unique(sc_obj.minus.messy.clusters$seurat_clusters)) {
  print(cluster_id)
  # Remove current cluster
  idxPass <- which(Idents(sc_obj.minus.messy.clusters) %in% c(cluster_id))
  cellsPass <- names(sc_obj.minus.messy.clusters$orig.ident[-idxPass])
  sc_obj.minus.current.cluster <- subset(x = sc_obj.minus.messy.clusters, subset = cell_name %in% cellsPass)
  cell_count <- length(sc_obj.minus.current.cluster$cell_name)
  current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  # Print cell type plot without cluster
  print_UMAP(sc_obj.minus.current.cluster, sample_count, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, paste0("_clusters_by_cell_type_minus_messy_clusters_without_cluster_", cluster_id, ".png"))
}

# Remove some more messy clusters, I guess
messy_clusters_round_2 <- c(1, 11)
idxPass <- which(Idents(sc_obj.minus.messy.clusters) %in% messy_clusters_round_2)
cellsPass <- names(sc_obj.minus.messy.clusters$orig.ident[-idxPass])
sc_obj.minus.messy.clusters <- subset(x = sc_obj.minus.messy.clusters, subset = cell_name %in% cellsPass)

# Plot again
cell_count <- length(sc_obj.minus.messy.clusters$cell_name)
current_title <- paste0("scRNA-seq and/or snRNA-seq Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", current_title, output_dir, naming_token, "_clusters_by_cell_type_without_messy_clusters_round_2.png")

# Print viral load plot
print_UMAP(sc_obj.minus.messy.clusters, "viral_load", current_title, output_dir, naming_token, "_clusters_by_viral_load_without_messy_clusters.png")

# Calculate cell type proportions and cell counts for each cell type (split by sample)
processed_cell_tables <- calculate_props_and_counts(sc_obj.minus.messy.clusters, viral_load_label, all_viral_load)
write.csv(processed_cell_tables[1], file = paste0(output_dir, naming_token, "_RNA_cell_type_proportion.csv"), quote = FALSE, row.names = FALSE)
write.csv(processed_cell_tables[2], file = paste0(output_dir, naming_token, "_RNA_cell_counts.csv"), quote = FALSE, row.names = FALSE)

# 2. Remove one cluster at a time and see how plot looks after removing each cluster
for(cluster_id in unique(sc_obj$seurat_clusters)) {
  print(cluster_id)
  # Remove current cluster
  idxPass <- which(Idents(sc_obj) %in% c(cluster_id))
  cellsPass <- names(sc_obj$orig.ident[-idxPass])
  sc_obj.minus.current.cluster <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  # Print cell type plot without cluster
  print_UMAP(sc_obj.minus.current.cluster, sample_count, "predicted_celltype_majority_vote", output_dir, naming_token, paste0("_clusters_by_cell_type_without_cluster_", cluster_id, "_", date, ".png"))
}


