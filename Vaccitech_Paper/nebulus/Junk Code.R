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


