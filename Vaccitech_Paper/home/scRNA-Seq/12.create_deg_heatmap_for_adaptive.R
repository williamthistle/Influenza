# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

adaptive_deg_heatmap_table <- data.frame("CD4 Memory" = numeric(), "CD8 Memory" = numeric(), "MAIT" = numeric(), "T Naive" = numeric(), "B" = numeric(),check.names = FALSE)

# Positive
for(cell_type in colnames(adaptive_deg_heatmap_table)) {
  cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table <- adaptive_sc_pseudobulk_deg_combined_cell_types_table[adaptive_sc_pseudobulk_deg_combined_cell_types_table$Cell_Type == cell_type,]
  other_cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table <- adaptive_sc_pseudobulk_deg_combined_cell_types_table[adaptive_sc_pseudobulk_deg_combined_cell_types_table$Cell_Type != cell_type,]
  ordered_positive_fc_genes <- cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table[order(cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table$sc_log2FC, decreasing = TRUE),]
  total_added_pos_genes <- 0
  index <- 1
  while(total_added_pos_genes < 5) {
    current_pos_gene <- ordered_positive_fc_genes[index,]
    if(!(current_pos_gene$Gene_Name %in% rownames(adaptive_deg_heatmap_table))) {
      total_added_pos_genes <- total_added_pos_genes + 1
      new_row <- c(NA, NA, NA, NA, NA)
      cell_type_index <- grep(current_pos_gene$Cell_Type, colnames(adaptive_deg_heatmap_table))
      new_row[cell_type_index] <- current_pos_gene$sc_log2FC
      other_cell_types_rows <- other_cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table[other_cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table$Gene_Name == current_pos_gene$Gene_Name,]
      if(nrow(other_cell_types_rows) > 0) {
        for(other_cell_type_row_index in 1:nrow(other_cell_types_rows)) {
          other_cell_type <- other_cell_types_rows[other_cell_type_row_index,]$Cell_Type
          other_fc_value <- other_cell_types_rows[other_cell_type_row_index,]$sc_log2FC
          cell_type_index <- grep(other_cell_type, colnames(adaptive_deg_heatmap_table))
          new_row[cell_type_index] <- other_fc_value
        }
      }
      new_row <- t(as.data.frame(new_row))
      colnames(new_row) <- colnames(adaptive_deg_heatmap_table)
      rownames(new_row) <- current_pos_gene$Gene_Name
      adaptive_deg_heatmap_table <- rbind(adaptive_deg_heatmap_table, new_row)
    }
    index <- index + 1
  }
}

# Negative
for(cell_type in colnames(adaptive_deg_heatmap_table)) {
  cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table <- adaptive_sc_pseudobulk_deg_combined_cell_types_table[adaptive_sc_pseudobulk_deg_combined_cell_types_table$Cell_Type == cell_type,]
  other_cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table <- adaptive_sc_pseudobulk_deg_combined_cell_types_table[adaptive_sc_pseudobulk_deg_combined_cell_types_table$Cell_Type != cell_type,]
  ordered_negative_fc_genes <- cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table[order(cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table$sc_log2FC),]
  total_added_neg_genes <- 0
  index <- 1
  while(total_added_neg_genes < 5) {
    current_neg_gene <- ordered_negative_fc_genes[index,]
    if(!(current_neg_gene$Gene_Name %in% rownames(adaptive_deg_heatmap_table))) {
      total_added_neg_genes <- total_added_neg_genes + 1
      new_row <- c(NA, NA, NA, NA, NA)
      cell_type_index <- grep(current_neg_gene$Cell_Type, colnames(adaptive_deg_heatmap_table))
      new_row[cell_type_index] <- current_neg_gene$sc_log2FC
      other_cell_types_rows <- other_cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table[other_cell_type_adaptive_sc_pseudobulk_deg_combined_cell_types_table$Gene_Name == current_neg_gene$Gene_Name,]
      if(nrow(other_cell_types_rows) > 0) {
        for(other_cell_type_row_index in 1:nrow(other_cell_types_rows)) {
          other_cell_type <- other_cell_types_rows[other_cell_type_row_index,]$Cell_Type
          other_fc_value <- other_cell_types_rows[other_cell_type_row_index,]$sc_log2FC
          cell_type_index <- grep(other_cell_type, colnames(adaptive_deg_heatmap_table))
          new_row[cell_type_index] <- other_fc_value
        }
      }
      new_row <- t(as.data.frame(new_row))
      colnames(new_row) <- colnames(adaptive_deg_heatmap_table)
      rownames(new_row) <- current_neg_gene$Gene_Name
      adaptive_deg_heatmap_table <- rbind(adaptive_deg_heatmap_table, new_row)
    }
    index <- index + 1
  }
}

pheatmap::pheatmap(as.matrix(adaptive_deg_heatmap_table), cluster_row=FALSE, cluster_col=FALSE, fontsize_col=14, 
                   main = "Fold Change for Top DEGs in Adaptive Immune Cell Types", filename = "C:/Users/willi/Desktop/adaptive_scRNA_deg_heatmap.png")

mat_df <- as.data.frame(as.table(as.matrix(adaptive_deg_heatmap_table)))

ggplot(mat_df, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "grey80") +
  scale_y_discrete(limits = rev(levels(mat_df$Var1))) +
  theme_minimal() +
  labs(title = "Fold Change for Top DEGs in Adaptive Immune Cell Types",
       x = "Cell Type",
       y = "Gene", fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=14)) + coord_fixed(ratio = 0.1)
