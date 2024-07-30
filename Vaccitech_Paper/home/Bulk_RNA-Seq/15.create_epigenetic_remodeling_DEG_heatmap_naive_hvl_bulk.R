# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

bulk_hvl_placebo_degs <- list()
bulk_hvl_placebo_degs[["Day 2"]] <- hvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results
bulk_hvl_placebo_degs[["Day 5"]] <- hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results
bulk_hvl_placebo_degs[["Day 8"]] <- hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results
bulk_hvl_placebo_degs[["Day 28"]] <- hvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results

heatmap_genes <- c("H1-2", "H1-3", "H1-4", "H2AZ1", "KAT6A", "BRD1", "MSL2", "HDAC5", "HDAC7", "HDAC9", "METTL23", "ASH1L", "PRDM2", "SETD2","KDM2B", "KDM3A")

# H1-2 -> HIST1H1C
# H1-3 -> HIST1H1D
# H1-4 -> HIST1H1E
# H2AZ1 -> H2AFZ
heatmap_gene_types <- list()
heatmap_gene_types[["H1-2 (HIST1H1C)"]] <- "Histone"
heatmap_gene_types[["H1-3 (HIST1H1D)"]] <- "Histone"
heatmap_gene_types[["H1-4 (HIST1H1E)"]] <- "Histone"
heatmap_gene_types[["H2AZ1 (H2AFZ)"]] <- "Histone"
heatmap_gene_types[["KAT6A"]] <- "Lysine Acetyltransferase"
heatmap_gene_types[["BRD1"]] <- "Lysine Acetyltransferase"
heatmap_gene_types[["MSL2"]] <- "Lysine Acetyltransferase"
heatmap_gene_types[["HDAC5"]] <- "Lysine Deacetylase"
heatmap_gene_types[["HDAC7"]] <- "Lysine Deacetylase"
heatmap_gene_types[["HDAC9"]] <- "Lysine Deacetylase"
heatmap_gene_types[["ASH1L"]] <- "Lysine Methyltransferase"
heatmap_gene_types[["PRDM2"]] <- "Lysine Methyltransferase"
heatmap_gene_types[["SETD2"]] <- "Lysine Methyltransferase"
heatmap_gene_types[["METTL23"]] <- "Lysine Methyltransferase"
heatmap_gene_types[["KDM2B"]] <- "Lysine Demethylase"
heatmap_gene_types[["KDM3A"]] <- "Lysine Demethylase"

bulk_day_vector <- c()
gene_name_vector <- c()
fold_change_vector <- c()
significance_vector <- c()

for(current_day in names(bulk_hvl_placebo_degs)) {
  current_bulk_degs <- bulk_hvl_placebo_degs[[current_day]][[2]]
  current_unfiltered_bulk_degs <- bulk_hvl_placebo_degs[[current_day]][[8]]
  for(heatmap_gene in heatmap_genes) {
    current_heatmap_gene_unfiltered_data <- current_unfiltered_bulk_degs[rownames(current_unfiltered_bulk_degs) == heatmap_gene,]
    bulk_day_vector <- c(bulk_day_vector, current_day)
    gene_name_vector <- c(gene_name_vector, heatmap_gene)
    fold_change_vector <- c(fold_change_vector, current_heatmap_gene_unfiltered_data$log2FoldChange)
    if(heatmap_gene %in% rownames(current_bulk_degs)) {
      current_p_value <- current_bulk_degs[rownames(current_bulk_degs) == heatmap_gene,]$padj
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

epigenetic_remodeling_heatmap_df <- data.frame(Day = bulk_day_vector, Gene_Name = gene_name_vector, fold_change = fold_change_vector, significant = significance_vector)

epigenetic_remodeling_heatmap_df[epigenetic_remodeling_heatmap_df$Gene_Name == "H1-2",]$Gene_Name <- "H1-2 (HIST1H1C)"
epigenetic_remodeling_heatmap_df[epigenetic_remodeling_heatmap_df$Gene_Name == "H1-3",]$Gene_Name <- "H1-3 (HIST1H1D)"
epigenetic_remodeling_heatmap_df[epigenetic_remodeling_heatmap_df$Gene_Name == "H1-4",]$Gene_Name <- "H1-4 (HIST1H1E)"
epigenetic_remodeling_heatmap_df[epigenetic_remodeling_heatmap_df$Gene_Name == "H2AZ1",]$Gene_Name <- "H2AZ1 (H2AFZ)"

epigenetic_remodeling_heatmap_df$Gene_Name <- factor(epigenetic_remodeling_heatmap_df$Gene_Name, levels = unique(epigenetic_remodeling_heatmap_df$Gene_Name))
epigenetic_remodeling_heatmap_df$Day <- factor(epigenetic_remodeling_heatmap_df$Day, levels = c("Day 2", "Day 5", "Day 8", "Day 28"))

epigenetic_remodeling_type <- c()
for(current_row_index in 1:nrow(epigenetic_remodeling_heatmap_df)) {
  current_row <- epigenetic_remodeling_heatmap_df[current_row_index,]
  epigenetic_remodeling_type <- c(epigenetic_remodeling_type, heatmap_gene_types[[as.character(current_row$Gene_Name)]])
}

epigenetic_remodeling_heatmap_df$epigenetic_remodeling_Type <- epigenetic_remodeling_type

category_colors <- epigenetic_remodeling_heatmap_df %>% 
  distinct(epigenetic_remodeling_Type, Gene_Name) %>% 
  mutate(color = case_when(
    epigenetic_remodeling_Type == "Histone" ~ "#555555",
    epigenetic_remodeling_Type == "Lysine Acetyltransferase" ~ "#6a3d9a",
    epigenetic_remodeling_Type == "Lysine Deacetylase" ~ "#0b6623",
    epigenetic_remodeling_Type == "Lysine Methyltransferase" ~ "#b22222",
    epigenetic_remodeling_Type == "Lysine Demethylase" ~ "#ff7f00",
    
    TRUE ~ "black" # default color
  )) %>% 
  pull(color, name = Gene_Name)


epigenetic_remodeling_heatmap_plot <- ggplot() + 
  geom_raster(data = epigenetic_remodeling_heatmap_df, aes(x = Day, y = Gene_Name, fill = fold_change)) +
  geom_text(data = epigenetic_remodeling_heatmap_df, aes(x = Day, y = Gene_Name, label = significant), nudge_y = 0.15, nudge_x = 0.20, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                      x = NULL,
                      y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16, color = category_colors[levels(factor(epigenetic_remodeling_heatmap_df$Gene_Name))])) # + coord_fixed(ratio = 0.5)

ggsave(filename = paste0("C:/Users/willi/Desktop/", "epigenetic_remodeling_deg_heatmap_naive_hvl_bulk.png"), plot = epigenetic_remodeling_heatmap_plot, device='png', dpi=300, width = 5, height = 8, units = "in")

# Flip it!
epigenetic_remodeling_heatmap_plot <- ggplot() + 
  geom_raster(data = epigenetic_remodeling_heatmap_df, aes(x = Gene_Name, y = Day, fill = fold_change)) +
  geom_text(data = epigenetic_remodeling_heatmap_df, aes(x = Gene_Name, y = Day, label = significant), nudge_y = 0.15, nudge_x = 0.10, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                                       x = NULL,
                                       y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16, color = category_colors[levels(factor(epigenetic_remodeling_heatmap_df$Gene_Name))]),
        axis.text.y = element_text(size = 16)) # + coord_fixed(ratio = 0.5)

ggsave(filename = paste0("C:/Users/willi/Desktop/", "epigenetic_remodeling_deg_heatmap_naive_hvl_bulk_flipped.png"), plot = epigenetic_remodeling_heatmap_plot, device='png', dpi=300, width = 8, height = 3, units = "in")

