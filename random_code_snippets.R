library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)
library(SeuratDisk)
library(ArchR) # Just used for confusion matrix
project_dir <- "~/POST_VACCINE_VS_PLACEBO/"
output_dir <- paste0(project_dir, "scRNA_seq_data_output/")
image_dir <- paste0(output_dir, "images/")
load(paste0(image_dir, "integrated_obj_after_predictions.RData"))


library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)
library(SeuratDisk)
library(ArchR) # Just used for confusion matrix
output_dir <- paste0("scRNA_seq_data_output/")
load(paste0(output_dir, "integrated_obj_after_predictions.RData"))




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
load(paste0(output_dir, "atac_after_lsi_umap_clusters_2.RData"))


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

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

dap_df <- read.table("MAGICAL_daps.txt", sep = "\t", header = TRUE)
cell_type_order <- c("T Naive", "CD4 Memory", "NK", "CD14 Mono", "CD16 Mono", "MAIT", "B", "CD8 Memory")
cell_type_order <- rev(cell_type_order)
ggplot(dap_df, aes(x=Cell.Type, y=DAPs, fill=Cell.Type, )) +
  geom_bar(stat="identity")+theme_minimal()+coord_flip() + scale_x_discrete(limits = cell_type_order) +
  xlab("Cell Type") + ylab("Number of DAPs") + ggtitle("DAPs per Cell Type") +
  labs(fill = "Cell Type") +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=14)) + 
  theme(plot.title = element_text(face="bold", size=20)) +
  theme(axis.title=element_text(size=18)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.title=element_text(size=16)) +
  scale_fill_discrete(name = "Cell Type") +
  scale_fill_manual(values=cbPalette)
ggsave("DAPs.png", device = "png", dpi = 300)

deg_df <- read.table("MAGICAL_degs.txt", sep = "\t", header = TRUE)
cell_type_order <- c("B", "T Naive", "CD4 Memory", "MAIT", "CD8 Memory", "NK", "CD16 Mono", "CD14 Mono")
ggplot(deg_df, aes(x=Cell.Type, y=DEGs, fill=Cell.Type, )) +
  geom_bar(stat="identity")+theme_minimal()+coord_flip() + scale_x_discrete(limits = cell_type_order) +
  xlab("Cell Type") + ylab("Number of DEGs") + ggtitle("DEGs per Cell Type") +
  labs(fill = "Cell Type") +
  theme(plot.title = element_text(hjust = 0.5), legend.text=element_text(size=14)) + 
  theme(plot.title = element_text(face="bold", size=20)) +
  theme(axis.title=element_text(size=18)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.title=element_text(size=16)) +
  scale_fill_discrete(name = "Cell Type") +
  scale_fill_manual(values=cbPalette)
ggsave("DEGs.png", device = "png", dpi = 300)

mag_associations_df <- read.table("MAGICAL_associations.txt", sep = "\t", header = TRUE)
cell_type_order <- c("CD8 Memory", "B", "MAIT", "CD4 Memory", "CD16 Mono", "CD14 Mono", "NK", "T Naive")
ggplot(mag_associations_df, aes(x=Cell.Type, y=MAGICAL.Associations, fill=Cell.Type)) +
  geom_bar(stat="identity")+theme_minimal()+coord_flip() + scale_x_discrete(limits = cell_type_order) +
  xlab("Cell Type") + ylab("Number of MAGICAL Associations") + ggtitle("MAGICAL Associations per Cell Type") +
  labs(fill = "Cell Type") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(face="bold", size=20)) +
  theme(axis.title=element_text(size=18)) +
  theme(axis.text=element_text(size=14)) +
  theme(legend.position="none") +
  scale_fill_discrete(name = "Cell Type") +
  scale_fill_manual(values=cbPalette)
ggsave("MAGICAL_associations.png", device = "png", dpi = 300)

cd4_memory_df <- read.table("CD4_Memory_Genes.txt")
unique_genes <- unique(cd4_memory_df$V1)
write.table(unique_genes, "CD4_Memory_Genes_Unique.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

cd8_memory_df <- read.table("CD8_Memory_Genes.txt")
unique_genes <- unique(cd8_memory_df$V1)
write.table(unique_genes, "CD8_Memory_Genes_Unique.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

t_naive_df <- read.table("T_Naive_Genes.txt")
unique_genes <- unique(t_naive_df$V1)
write.table(unique_genes, "T_Naive_Genes_Unique.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

mait_df <- read.table("MAIT_Genes.txt")
unique_genes <- unique(mait_df$V1)
write.table(unique_genes, "MAIT_Genes_Unique.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

cd14_mono_df <- read.table("CD14_Mono_Genes.txt")
unique_genes <- unique(cd14_mono_df$V1)
write.table(unique_genes, "CD14_Mono_Genes_Unique.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

cd16_mono_df <- read.table("CD16_Mono_Genes.txt")
unique_genes <- unique(cd16_mono_df$V1)
write.table(unique_genes, "CD16_Mono_Genes_Unique.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

nk_df <- read.table("NK_Genes.txt")
unique_genes <- unique(nk_df$V1)
write.table(unique_genes, "NK_Genes_Unique.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

b_df <- read.table("B_Genes.txt")
unique_genes <- unique(b_df$V1)
write.table(unique_genes, "B_Genes_Unique.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)