# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

innate_scRNA_hvl_placebo_degs <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type %in% innate_cell_types,]

heatmap_genes <- c("KAT6A", "HDAC5", "HDAC7", "HDAC9", "ASH1L", "PRDM2", "SETD2", "KDM2B", "KDM3A")

heatmap_gene_types <- list()
heatmap_gene_types[["KAT6A"]] <- "Histone Acetyltransferase"
heatmap_gene_types[["HDAC5"]] <- "Histone Deacetylase"
heatmap_gene_types[["HDAC7"]] <- "Histone Deacetylase"
heatmap_gene_types[["HDAC9"]] <- "Histone Deacetylase"
heatmap_gene_types[["ASH1L"]] <- "Histone Methyltransferase"
heatmap_gene_types[["PRDM2"]] <- "Histone Methyltransferase"
heatmap_gene_types[["SETD2"]] <- "Histone Methyltransferase"
heatmap_gene_types[["KDM2B"]] <- "Histone Demethylase"
heatmap_gene_types[["KDM3A"]] <- "Histone Demethylase"

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

epigenetic_heatmap_df <- data.frame(Cell_Type = cell_type_vector, Gene_Name = gene_name_vector, fold_change = fold_change_vector, significant = significance_vector)

epigenetic_heatmap_df$Gene_Name <- factor(epigenetic_heatmap_df$Gene_Name, levels = unique(epigenetic_heatmap_df$Gene_Name))
epigenetic_heatmap_df$Cell_Type <- factor(epigenetic_heatmap_df$Cell_Type, levels = c("CD14 Mono", "CD16 Mono", "NK", "NK_CD56bright", "cDC", "pDC"))

epigenetic_type <- c()
for(current_row_index in 1:nrow(epigenetic_heatmap_df)) {
  current_row <- epigenetic_heatmap_df[current_row_index,]
  epigenetic_type <- c(epigenetic_type, heatmap_gene_types[[current_row$Gene_Name]])
}

epigenetic_heatmap_df$Epigenetic_Type <- epigenetic_type

category_colors <- epigenetic_heatmap_df %>% 
  distinct(Epigenetic_Type, Gene_Name) %>% 
  mutate(color = case_when(
    Epigenetic_Type == "Histone Acetyltransferase" ~ "red",
    Epigenetic_Type == "Histone Deacetylase" ~ "green",
    Epigenetic_Type == "Histone Methyltransferase" ~ "blue",
    Epigenetic_Type == "Histone Demethylase" ~ "orange",
    TRUE ~ "black" # default color
  )) %>% 
  pull(color, name = Gene_Name)

ggplot() + 
  geom_raster(data = epigenetic_heatmap_df, aes(x = Cell_Type, y = Gene_Name, fill = fold_change)) +
  geom_text(data = epigenetic_heatmap_df, aes(x = Cell_Type, y = Gene_Name, label = significant), nudge_y = 0.15, nudge_x = 0.20, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = "Fold Change for Epigenetic Remodeling Genes in Innate Immune Cell Types",
                      x = "Cell Type",
                      y = "Gene", fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12, color = category_colors[levels(factor(epigenetic_heatmap_df$Gene_Name))])) + coord_fixed(ratio = 0.5)