# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

innate_scRNA_hvl_vaccinated_degs <- scRNA_hvl_vaccinated_degs[scRNA_hvl_vaccinated_degs$Cell_Type %in% innate_cell_types,]

# Find epigenetic remodeling DEGs
search_terms <- c("histone", "methyltransferase", "acetyltransferase", "demethylase")
gene_terms <- c("HIST1H1C", "HIST1H1D", "HIST1H1E", "H2AFZ", "KAT6A", "BRD1", "MSL2", "HDAC5", "HDAC7", "HDAC9", "METTL23", "ASH1L", "PRDM2", "SETD2","KDM2B", "KDM3A",
                              "KAT", "HDAC", "KDM")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c('hgnc_symbol', 'description', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'hgnc_symbol',
  values = unique(innate_scRNA_hvl_vaccinated_degs$Gene_Name),
  mart = ensembl
)

annotations_final <- annotations %>%
  filter(
    str_detect(hgnc_symbol, paste(gene_terms, collapse = "|")) |
      str_detect(description, paste(search_terms, collapse = "|"))
  )

# Final heatmap genes
heatmap_genes <- c("HDAC2", "HDAC9", "SETD3", "KAT6B")

heatmap_gene_types <- list()
heatmap_gene_types[["HDAC2"]] <- "Lysine Deacetylase"
heatmap_gene_types[["HDAC9"]] <- "Lysine Deacetylase"
heatmap_gene_types[["KAT6B"]] <- "Lysine Acetyltransferase"
heatmap_gene_types[["SETD3"]] <- "Lysine Methyltransferase"

cell_type_vector <- c()
gene_name_vector <- c()
fold_change_vector <- c()
significance_vector <- c()

for(innate_cell_type in innate_cell_types) {
  innate_cell_type_for_file_name <- sub(" ", "_", innate_cell_type)
  current_scRNA_hvl_vaccinated_degs <- innate_scRNA_hvl_vaccinated_degs[innate_scRNA_hvl_vaccinated_degs$Cell_Type == innate_cell_type,]
  current_unfiltered_hvl_vaccinated_degs <- read.table(paste0(scRNA_hvl_vaccinated_deg_dir, "D28-vs-D_minus_1-degs-", innate_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                    sep = "\t", header = TRUE)
  for(heatmap_gene in heatmap_genes) {
    current_heatmap_gene_unfiltered_data <- current_unfiltered_hvl_vaccinated_degs[rownames(current_unfiltered_hvl_vaccinated_degs) == heatmap_gene,]
    cell_type_vector <- c(cell_type_vector, innate_cell_type)
    gene_name_vector <- c(gene_name_vector, heatmap_gene)
    fold_change_vector <- c(fold_change_vector, current_heatmap_gene_unfiltered_data$avg_log2FC)
    if(heatmap_gene %in% current_scRNA_hvl_vaccinated_degs$Gene_Name) {
      current_p_value <- current_scRNA_hvl_vaccinated_degs[current_scRNA_hvl_vaccinated_degs$Gene_Name == heatmap_gene,]$sc_pval_adj
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
epigenetic_heatmap_df$Cell_Type <- factor(epigenetic_heatmap_df$Cell_Type, levels = c("CD14 Mono", "CD16 Mono", "cDC", "pDC", "NK", "NK_CD56bright"))

epigenetic_type <- c()
for(current_row_index in 1:nrow(epigenetic_heatmap_df)) {
  current_row <- epigenetic_heatmap_df[current_row_index,]
  epigenetic_type <- c(epigenetic_type, heatmap_gene_types[[current_row$Gene_Name]])
}

epigenetic_heatmap_df$Epigenetic_Type <- epigenetic_type

category_colors <- epigenetic_heatmap_df %>% 
  distinct(Epigenetic_Type, Gene_Name) %>% 
  mutate(color = case_when(
    Epigenetic_Type == "Histone" ~ "#555555",
    Epigenetic_Type == "Lysine Acetyltransferase" ~ "#6a3d9a",
    Epigenetic_Type == "Lysine Deacetylase" ~ "#0b6623",
    Epigenetic_Type == "Lysine Methyltransferase" ~ "#b22222",
    Epigenetic_Type == "Lysine Demethylase" ~ "#ff7f00",
    
    TRUE ~ "black" # default color
  )) %>% 
  pull(color, name = Gene_Name)

epigenetic_remodeling_heatmap_plot <- ggplot() + 
  geom_raster(data = epigenetic_heatmap_df, aes(x = Cell_Type, y = Gene_Name, fill = fold_change)) +
  geom_text(data = epigenetic_heatmap_df, aes(x = Cell_Type, y = Gene_Name, label = significant), nudge_y = 0.15, nudge_x = 0.20, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                                       x = NULL,
                                       y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16, color = category_colors[levels(factor(epigenetic_heatmap_df$Gene_Name))]))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "epigenetic_remodeling_deg_hvl_vaccinated_heatmap.png"), plot = epigenetic_remodeling_heatmap_plot, device='png', dpi=300, width = 5, height = 5, units = "in")

# Flip it!
epigenetic_remodeling_heatmap_plot <- ggplot() + 
  geom_raster(data = epigenetic_heatmap_df, aes(x = Gene_Name, y = Cell_Type, fill = fold_change)) +
  geom_text(data = epigenetic_heatmap_df, aes(x = Gene_Name, y = Cell_Type, label = significant), nudge_y = 0.15, nudge_x = 0.10, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                                       x = NULL,
                                       y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16, color = category_colors[levels(factor(epigenetic_heatmap_df$Gene_Name))]),
        axis.text.y = element_text(size = 16)) # + coord_fixed(ratio = 0.5)

ggsave(filename = paste0("C:/Users/willi/Desktop/", "epigenetic_remodeling_deg_heatmap_hvl_vaccinated_flipped.png"), plot = epigenetic_remodeling_heatmap_plot, device='png', dpi=300, width = 8, height = 3, units = "in")
