library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)
library(SeuratDisk)
library(ArchR) # Just used for confusion matrix
output_dir <- "~/scRNA_seq_data_output/"
load(paste0(output_dir, "integrated_obj_after_final_plot.RData"))


library(ArchR)
library(stringr)
library(pheatmap)
library(mclust)
library(writexl)
library(ggplot2)
library(hexbin)
library(SeuratDisk)
library(dplyr)
library(openxlsx)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Matrix)

set.seed(1)
output_dir <- "~/scATAC_seq_data_output/"
load(paste0(output_dir, "atac_after_peak_matrix.RData"))


for (res in seq(0, 1, 0.1)) {
  flu.combined.sct <- FindClusters(flu.combined.sct, resolution = res)
}

clustree(flu.combined.sct, prefix = "integrated_snn_res.")
ggsave(paste0(output_dir, "flu.combined.sct.cluster.trees.PDF"), device = "pdf")

proj <- saveArchRProject()
proj <- loadArchRProject(path = paste0(output_dir, "ArchR"))

load(paste0(output_dir, "atac_after_cell_type_voting_pre.RData"))

cluster_predictions <- vector()
cluster_distributions <- list()
idx <- 1
for (cluster in unique(proj.filtered$Clusters)) {
  idxPass <- which(proj.filtered$Clusters %in% cluster)
  cellsPass <- proj.filtered$cellNames[idxPass]
  filtered_cluster <-proj.filtered[cellsPass,]
  cluster_predictions <- append(cluster_predictions, table(filtered_cluster$Cell_type_voting))
  cluster_distributions[[idx]] <- table(filtered_cluster$Cell_type_combined)
  idx <- idx + 1
}

names(cluster_predictions) <- paste(unique(proj.filtered$Clusters), "-", names(cluster_predictions))
# Way to see memory of all R objects
sort( sapply(ls(),function(x){object.size(get(x))})) 

# Perform majority cell-type voting
# 1) If one type is vast majority, then convert all cells to that type
# 2) If there is an even mix of two cell types, then potentially keep both and feed all other cells into one of the two types
#    The idea behind this was to preserve borders between clusters, but in reality, it turns out that just converting all cells to the majority works better.
performCellMajorityVoteConversion <- function(seurat_obj, cluster_ids, keep_second, fold_into_second) {
  for (cluster_id in cluster_ids) {
    print(cluster_id)
    idxPass <- which(Idents(seurat_obj) %in% cluster_id)
    cellsPass <- names(seurat_obj$orig.ident[idxPass])
    cluster_distribution <- table(seurat_obj$Cell_type_combined[idxPass])
    cluster_distribution <- sort(cluster_distribution, TRUE)
    if(length(cluster_distribution) == 1) {
      print("No need to update this cluster because it only has one cell type")
    } else {
      # If keep_second, then we keep the cell type with second largest count (cluster border)
      if (!missing(keep_second)) {
        # If fold_into is missing, then fold everything else into majority class
        if(missing(fold_into_second)) {
          seurat_obj$Cell_type_combined[idxPass][seurat_obj$Cell_type_combined[idxPass] != names(cluster_distribution)[1] &
                                                   seurat_obj$Cell_type_combined[idxPass] != names(cluster_distribution)[2]] <- names(cluster_distribution)[1] 
        } else {
          seurat_obj$Cell_type_combined[idxPass][seurat_obj$Cell_type_combined[idxPass] != names(cluster_distribution)[1] &
                                                   seurat_obj$Cell_type_combined[idxPass] != names(cluster_distribution)[2]] <- names(cluster_distribution)[2]
        }
      } else {
        # If no keep_second, then just make everything the majority class
        seurat_obj$Cell_type_combined[idxPass][seurat_obj$Cell_type_combined[idxPass] != names(cluster_distribution)[1]] <- names(cluster_distribution)[1]
      }
    }
  }
  seurat_obj
}

plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE) + 
  labs(title = "scATAC-seq Data Integration \n (12 Samples, 64K Cells)") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(output_dir, "test.png"), device = "png", dpi = 300)

