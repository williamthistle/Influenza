# Code to look at removing subjects and seeing if it improves DEG signal
subject_3_groupings <- t(combn(unique(hvl_sc_obj$subject_id),3))
diff_results_3 <- list()
for(row_index in 1:nrow(subject_3_groupings)) {
  current_grouping <- subject_3_groupings[row_index,]
  idxPass <- which(hvl_sc_obj$subject_id %in% current_grouping)
  cellsPass <- names(hvl_sc_obj$orig.ident[idxPass])
  current_hvl_sc_obj <- subset(x = hvl_sc_obj, subset = cell_name %in% cellsPass)
  DefaultAssay(current_hvl_sc_obj) <- "SCT"
  HVL_differential_genes_dir <- paste0(RNA_output_dir, "diff_genes/", date, "/HVL_SCT_minus_", row_index, "/")
  if (!dir.exists(HVL_differential_genes_dir)) {dir.create(HVL_differential_genes_dir, recursive = TRUE)}
  run_differential_expression_group(current_hvl_sc_obj, HVL_differential_genes_dir, "time_point")
  DefaultAssay(current_hvl_sc_obj) <- "RNA"
  pseudobulk_de_df <- run_de(current_hvl_sc_obj, replicate_col = "sample", cell_type_col = "magical_cell_types", label_col = "time_point", de_method = "DESeq2")
  pseudobulk_de_df <- na.omit(pseudobulk_de_df)
  pseudobulk_de_df <- pseudobulk_de_df[pseudobulk_de_df$p_val < 0.05,]
  pseudobulk_de_df <- pseudobulk_de_df[pseudobulk_de_df$avg_logFC < -0.3 | pseudobulk_de_df$avg_logFC > 0.3,]
  pseudobulk_cell_types_for_correction <- c("B", "CD4_Memory", "CD8_Memory", "CD14_Mono", "CD16_Mono", "NK_MAGICAL", "T_Naive")
  final_list_of_genes <- data.frame(Cell_Type = character(), Gene_Name = character(), sc_pval_adj = character(), sc_log2FC = character(), pseudo_bulk_pval = character(),
                                    pseudo_bulk_log2FC = character())
  for(current_cell_type in pseudobulk_cell_types_for_correction) {
    current_DEG_table <- read.table(paste0(HVL_differential_genes_dir, "D28-vs-D_minus_1-degs-", current_cell_type, "-time_point.csv"), sep = ",", header = TRUE)
    current_DEG_table <- current_DEG_table[current_DEG_table$p_val_adj < 0.05,]
    current_DEG_table <- current_DEG_table[abs(current_DEG_table$avg_log2FC) > 0.1,]
    current_DEG_table <- current_DEG_table[current_DEG_table$pct.1 > 0.1 | current_DEG_table$pct.2 > 0.1,]
    if(current_cell_type != "NK_MAGICAL") {
      pseudobulk_de_df_cell_type_subset <- pseudobulk_de_df[pseudobulk_de_df$cell_type == sub("_", " ", current_cell_type),]
    } else {
      pseudobulk_de_df_cell_type_subset <- pseudobulk_de_df[pseudobulk_de_df$cell_type == current_cell_type,]
    }
    final_cell_type_genes <- intersect(current_DEG_table$X, pseudobulk_de_df_cell_type_subset$gene)
    for(current_gene in final_cell_type_genes) {
      current_sc_pval_adj <- current_DEG_table[current_DEG_table$X == current_gene,]$p_val_adj
      current_sc_log2FC <- current_DEG_table[current_DEG_table$X == current_gene,]$avg_log2FC
      current_pseudo_bulk_pval <- pseudobulk_de_df_cell_type_subset[pseudobulk_de_df_cell_type_subset$gene == current_gene,]$p_val
      current_pseudo_bulk_log2FC <- pseudobulk_de_df_cell_type_subset[pseudobulk_de_df_cell_type_subset$gene == current_gene,]$avg_logFC
      current_row <- data.frame(current_cell_type, current_gene, current_sc_pval_adj, current_sc_log2FC, current_pseudo_bulk_pval, current_pseudo_bulk_log2FC)
      names(current_row) <- c("Cell_Type", "Gene_Name", "sc_pval_adj", "sc_log2FC", "pseudo_bulk_pval", "pseudo_bulk_log2FC")
      final_list_of_genes <- rbind(final_list_of_genes, current_row)
    }
  }
  diff_results_3[[row_index]] <- final_list_of_genes
}



idxPass <- which(hvl_sc_obj_d28$predicted_celltype_majority_vote %in% "CD16 Mono")
cellsPass <- names(hvl_sc_obj_d28$orig.ident[idxPass])
cd14_mono_hvl_sc_obj_d28 <- subset(x = hvl_sc_obj_d28, subset = cell_name %in% cellsPass)
table(cd14_mono_hvl_sc_obj_d28$sample)







# Code to look at how GeneIntegration labels look for RNA cells
atac_proj <- loadArchRProject(path = paste0(ATAC_output_dir, "ArchROutput"))

atac_proj <- add_sample_metadata_atac(atac_proj, high_viral_load_samples, low_viral_load_samples,
                                      d28_samples, d_minus_1_samples, male_samples, female_samples)
viral_load_metadata <- parse_metadata_for_samples(atac_proj, "viral_load", high_viral_load_samples, low_viral_load_samples,
                                                  d28_samples, d_minus_1_samples, male_samples, female_samples)
day_metadata <- parse_metadata_for_samples(atac_proj, "time_point", high_viral_load_samples, low_viral_load_samples,
                                           d28_samples, d_minus_1_samples, male_samples, female_samples)
sex_metadata <- parse_metadata_for_samples(atac_proj, "sex", high_viral_load_samples, low_viral_load_samples,
                                           d28_samples, d_minus_1_samples, male_samples, female_samples)
atac_proj <- combine_cell_types_atac(atac_proj)
atac_proj <- MajorityVote_ATAC(proj = atac_proj)

idxPass <- which(atac_proj$viral_load %in% c("high"))
cellsPass <- atac_proj$cellNames[idxPass]
HVL_atac_proj <- atac_proj[cellsPass, ]

atac_proj_2 <- loadArchRProject(path = paste0(ATAC_output_dir, "ArchROutput"))
atac_proj_2 <- add_sample_metadata_atac(atac_proj_2, high_viral_load_samples, low_viral_load_samples,
                                        d28_samples, d_minus_1_samples, male_samples, female_samples)
atac_proj_2 <- add_rna_labels_for_atac_data(atac_proj_2, ATAC_output_dir, source_rna_file = "rna_seq_labeled_cells-14_final.csv", subset_to_rna = TRUE)
atac_proj_2 <- combine_cell_types_atac(atac_proj_2)
atac_proj_2$Cell_type_voting <- atac_proj_2$predictedGroup
atac_proj_2 <- remove_cell_types(atac_proj_2, c("HSPC", "Plasmablast", "Proliferating", "MAIT"))
idxPass <- which(atac_proj_2$viral_load %in% c("high"))
cellsPass <- atac_proj_2$cellNames[idxPass]
HVL_atac_proj_2 <- atac_proj_2[cellsPass, ]



idxPass <- which(HVL_atac_proj$cellNames %in% HVL_atac_proj_2$cellNames)
cellsPass <- HVL_atac_proj$cellNames[idxPass]
HVL_atac_proj <- HVL_atac_proj[cellsPass, ]

proj <- HVL_atac_proj
log_flag <- FALSE
reducedDims_param <- "Harmony"

tile_reduc <- ArchR::getReducedDims(ArchRProj = proj, 
                                    reducedDims = reducedDims_param, returnMatrix = TRUE)
tmp <- matrix(stats::rnorm(nrow(tile_reduc) * 3, 10), 
              ncol = nrow(tile_reduc), nrow = 3)
colnames(tmp) <- rownames(tile_reduc)
rownames(tmp) <- paste0("t", seq_len(nrow(tmp)))
obj <- Seurat::CreateSeuratObject(tmp, project = "scATAC", 
                                  min.cells = 0, min.features = 0)
obj[[reducedDims_param]] <- Seurat::CreateDimReducObject(embeddings = tile_reduc, 
                                                         key = paste0(reducedDims_param, "_"), assay = "RNA")
obj <- Seurat::FindNeighbors(obj, reduction = reducedDims_param, 
                             dims = 1:29)
obj <- find_clusters_SPEEDI(obj, resolution = 2, method = "Leiden", 
                            log_flag = log_flag)
obj <- Seurat::RunUMAP(obj, reduction = reducedDims_param, 
                       dims = 1:29, seed.use = get_speedi_seed())
proj <- ArchR::addCellColData(ArchRProj = proj, cells = names(obj$seurat_clusters), 
                              data = as.character(obj$seurat_clusters), name = "seurat_clusters", 
                              force = TRUE)
rm(obj)

HVL_atac_proj <- proj
HVL_atac_proj <- MajorityVote_ATAC(proj = HVL_atac_proj)

num_cells <- length(HVL_atac_proj$cellNames)
num_samples <- length(unique(HVL_atac_proj$Sample))
sample_text <- paste0("(", num_samples, " Samples, ", 
                      num_cells, " Cells)")

pal <- paletteDiscrete(values = HVL_atac_proj$Cell_type_voting)
p1 <- ArchR::plotEmbedding(ArchRProj = HVL_atac_proj, colorBy = "cellColData", 
                           name = "Cell_type_voting", embedding = "UMAP", 
                           pal = pal, force = TRUE, keepAxis = TRUE) + 
  ggplot2::ggtitle(paste0("ATAC Data Integration\n(By Majority Vote Cell Type, Reclustered)\n", 
                          sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size = 18), 
                                                         legend.text = ggplot2::element_text(size = 10))
ggplot2::ggsave(filename = paste0(ATAC_output_dir, "HVL_Final_ATAC_UMAP_by_Majority_Vote_Cell_Type_RNA_Cells_Reclustered.png"), 
                plot = p1, device = "png", width = 8, height = 8, 
                units = "in")

cluster_info <- get_cluster_info(HVL_atac_proj)



