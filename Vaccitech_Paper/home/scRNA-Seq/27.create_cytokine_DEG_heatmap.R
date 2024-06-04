# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

innate_scRNA_hvl_placebo_degs <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type %in% innate_cell_types,]

heatmap_genes <- c("CCL3", "CX3CR1", "CXCL16", "NFIL3", "IL32", "IRAK3", "IL1RAP", "IRF2", "IRF7", "IFNGR1", "IFNG", "JUN", "JUNB", 
                   "FOSL2")

cell_type_vector <- c()
gene_name_vector <- c()
fold_change_vector <- c()
significance_vector <- c()

for(innate_cell_type in innate_cell_types) {
  innate_cell_type_for_file_name <- sub(" ", "_", innate_cell_type)
  current_scRNA_hvl_placebo_degs <- innate_scRNA_hvl_placebo_degs[innate_scRNA_hvl_placebo_degs$Cell_Type == innate_cell_type,]
  current_unfiltered_hvl_placebo_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", innate_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                           sep = "\t", header = TRUE)
  for(heatmap_gene in heatmap_genes) {
    current_heatmap_gene_unfiltered_data <- current_unfiltered_hvl_placebo_degs[rownames(current_unfiltered_hvl_placebo_degs) == heatmap_gene,]
    cell_type_vector <- c(cell_type_vector, innate_cell_type)
    gene_name_vector <- c(gene_name_vector, heatmap_gene)
    fold_change_vector <- c(fold_change_vector, current_heatmap_gene_unfiltered_data$avg_log2FC)
    if(heatmap_gene %in% current_scRNA_hvl_placebo_degs$Gene_Name) {
      current_p_value <- current_scRNA_hvl_placebo_degs[current_scRNA_hvl_placebo_degs$Gene_Name == heatmap_gene,]$sc_pval_adj
      if(current_p_value < 0.001) {
        significance_vector <- c(significance_vector, "***")
      } else if(current_p_value < 0.01) {
        significance_vector <- c(significance_vector, "**")
      } else if(current_p_value < 0.05) {
        significance_vector <- c(significance_vector, "*")
      }
    } else {
      significance_vector <- c(significance_vector, NA)
    }
  }
}

cytokine_heatmap_df <- data.frame(Cell_Type = cell_type_vector, Gene_Name = gene_name_vector, fold_change = fold_change_vector, significant = significance_vector)

cytokine_heatmap_df$Gene_Name <- factor(cytokine_heatmap_df$Gene_Name, levels = unique(cytokine_heatmap_df$Gene_Name))
cytokine_heatmap_df$Cell_Type <- factor(cytokine_heatmap_df$Cell_Type, levels = c("CD14 Mono", "CD16 Mono", "NK", "NK_CD56bright", "cDC", "pDC"))


ggplot() + 
  geom_raster(data = cytokine_heatmap_df, aes(x = Cell_Type, y = Gene_Name, fill = fold_change)) +
  geom_text(data = cytokine_heatmap_df, aes(x = Cell_Type, y = Gene_Name, label = significant), nudge_y = 0.15, nudge_x = 0.40, size = 6) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic() + labs(title = "Fold Change for Cytokine-Related DEG in Innate Immune Cell Types",
                      x = "Cell Type",
                      y = "Gene", fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5))