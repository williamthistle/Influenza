# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

innate_deg_heatmap_table <- data.frame("CD14 Mono" = numeric(), "CD16 Mono" = numeric(), "NK" = numeric(), "cDC" = numeric(), "pDC" = numeric(), check.names = FALSE)

# Positive
for(cell_type in colnames(innate_deg_heatmap_table)) {
  cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table <- innate_sc_pseudobulk_deg_combined_cell_types_table[innate_sc_pseudobulk_deg_combined_cell_types_table$Cell_Type == cell_type,]
  other_cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table <- innate_sc_pseudobulk_deg_combined_cell_types_table[innate_sc_pseudobulk_deg_combined_cell_types_table$Cell_Type != cell_type,]
  ordered_positive_fc_genes <- cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table[order(cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table$sc_log2FC, decreasing = TRUE),]
  total_added_pos_genes <- 0
  index <- 1
  while(total_added_pos_genes < 5) {
    current_pos_gene <- ordered_positive_fc_genes[index,]
    if(!(current_pos_gene$Gene_Name %in% rownames(innate_deg_heatmap_table))) {
      total_added_pos_genes <- total_added_pos_genes + 1
      new_row <- c(NA, NA, NA, NA, NA)
      cell_type_index <- grep(current_pos_gene$Cell_Type, colnames(innate_deg_heatmap_table))
      new_row[cell_type_index] <- current_pos_gene$sc_log2FC
      other_cell_types_rows <- other_cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table[other_cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table$Gene_Name == current_pos_gene$Gene_Name,]
      if(nrow(other_cell_types_rows) > 0) {
        for(other_cell_type_row_index in 1:nrow(other_cell_types_rows)) {
          other_cell_type <- other_cell_types_rows[other_cell_type_row_index,]$Cell_Type
          other_fc_value <- other_cell_types_rows[other_cell_type_row_index,]$sc_log2FC
          cell_type_index <- grep(other_cell_type, colnames(innate_deg_heatmap_table))
          new_row[cell_type_index] <- other_fc_value
        }
      }
      new_row <- t(as.data.frame(new_row))
      colnames(new_row) <- colnames(innate_deg_heatmap_table)
      rownames(new_row) <- current_pos_gene$Gene_Name
      innate_deg_heatmap_table <- rbind(innate_deg_heatmap_table, new_row)
    }
    index <- index + 1
  }
}

# Negative
for(cell_type in colnames(innate_deg_heatmap_table)) {
  cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table <- innate_sc_pseudobulk_deg_combined_cell_types_table[innate_sc_pseudobulk_deg_combined_cell_types_table$Cell_Type == cell_type,]
  other_cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table <- innate_sc_pseudobulk_deg_combined_cell_types_table[innate_sc_pseudobulk_deg_combined_cell_types_table$Cell_Type != cell_type,]
  ordered_negative_fc_genes <- cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table[order(cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table$sc_log2FC),]
  total_added_neg_genes <- 0
  index <- 1
  while(total_added_neg_genes < 5) {
    current_neg_gene <- ordered_negative_fc_genes[index,]
    if(!(current_neg_gene$Gene_Name %in% rownames(innate_deg_heatmap_table))) {
      total_added_neg_genes <- total_added_neg_genes + 1
      new_row <- c(NA, NA, NA, NA, NA)
      cell_type_index <- grep(current_neg_gene$Cell_Type, colnames(innate_deg_heatmap_table))
      new_row[cell_type_index] <- current_neg_gene$sc_log2FC
      other_cell_types_rows <- other_cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table[other_cell_type_innate_sc_pseudobulk_deg_combined_cell_types_table$Gene_Name == current_neg_gene$Gene_Name,]
      if(nrow(other_cell_types_rows) > 0) {
        for(other_cell_type_row_index in 1:nrow(other_cell_types_rows)) {
          other_cell_type <- other_cell_types_rows[other_cell_type_row_index,]$Cell_Type
          other_fc_value <- other_cell_types_rows[other_cell_type_row_index,]$sc_log2FC
          cell_type_index <- grep(other_cell_type, colnames(innate_deg_heatmap_table))
          new_row[cell_type_index] <- other_fc_value
        }
      }
      new_row <- t(as.data.frame(new_row))
      colnames(new_row) <- colnames(innate_deg_heatmap_table)
      rownames(new_row) <- current_neg_gene$Gene_Name
      innate_deg_heatmap_table <- rbind(innate_deg_heatmap_table, new_row)
    }
    index <- index + 1
  }
}

#pheatmap::pheatmap(as.matrix(innate_deg_heatmap_table), cluster_row=FALSE, cluster_col=FALSE, fontsize_col=14, 
#                   main = "Fold Change for Top DEGs in Innate Immune Cell Types", filename = "C:/Users/willi/Desktop/innate_scRNA_deg_heatmap.png")

mat_df <- as.data.frame(as.table(as.matrix(innate_deg_heatmap_table)))

ggplot() +
  geom_tile(data = mat_df, aes(x = Var2, y = Var1, fill = Freq, color = "")) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "grey80") +
  scale_y_discrete(limits = rev(levels(mat_df$Var1))) +
  scale_color_manual(values=NA) +              
  guides(color=guide_legend("Not Significant", override.aes=list(fill = "grey80", color="black"), order = 999)) +
  theme_minimal() +
  labs(title = "Fold Change for Top DEGs in Innate Immune Cell Types",
       x = "Cell Type",
       y = "Gene", fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=14)) + coord_fixed(ratio = 0.1) 

