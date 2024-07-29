# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

bulk_lvl_placebo_degs <- list()
bulk_lvl_placebo_degs[["Day 2"]] <- lvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results
bulk_lvl_placebo_degs[["Day 5"]] <- lvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results
bulk_lvl_placebo_degs[["Day 8"]] <- lvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results
bulk_lvl_placebo_degs[["Day 28"]] <- lvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results

heatmap_genes <- c("CCL3", "CCL4", "CCL4L2", "CXCL16", "CXCL2",
                   "IL1B",
                   "IFI44L", "IFIT2", "IFIT3", "IFIT5", "IFITM1", "IFITM10", "IFRD1", 
                   "FOS", "FOSB", "JUN", "JUND",
                   "DUSP1",
                   "NFKBIA")

heatmap_gene_types <- list()
heatmap_gene_types[["CCL3"]] <- "Chemokine"
heatmap_gene_types[["CCL4"]] <- "Chemokine"
heatmap_gene_types[["CCL4L2"]] <- "Chemokine"
heatmap_gene_types[["CXCL16"]] <- "Chemokine"
heatmap_gene_types[["CXCL2"]] <- "Chemokine"

heatmap_gene_types[["IL1B"]] <- "Interleukin"

heatmap_gene_types[["IFI44L"]] <- "Interferon"
heatmap_gene_types[["IFIT2"]] <- "Interferon"
heatmap_gene_types[["IFIT3"]] <- "Interferon"
heatmap_gene_types[["IFIT5"]] <- "Interferon"
heatmap_gene_types[["IFITM1"]] <- "Interferon"
heatmap_gene_types[["IFITM10"]] <- "Interferon"
heatmap_gene_types[["IFRD1"]] <- "Interferon"

heatmap_gene_types[["FOS"]] <- "AP-1"
heatmap_gene_types[["FOSB"]] <- "AP-1"
heatmap_gene_types[["JUN"]] <- "AP-1"
heatmap_gene_types[["JUND"]] <- "AP-1"

heatmap_gene_types[["DUSP1"]] <- "MAP Kinase"

heatmap_gene_types[["NFKBIA"]] <- "NF-κB"

bulk_day_vector <- c()
gene_name_vector <- c()
fold_change_vector <- c()
significance_vector <- c()

for(current_day in names(bulk_lvl_placebo_degs)) {
  current_bulk_degs <- bulk_lvl_placebo_degs[[current_day]][[2]]
  current_unfiltered_bulk_degs <- bulk_lvl_placebo_degs[[current_day]][[8]]
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

cytokine_heatmap_df <- data.frame(Day = bulk_day_vector, Gene_Name = gene_name_vector, fold_change = fold_change_vector, significant = significance_vector)

cytokine_heatmap_df$Gene_Name <- factor(cytokine_heatmap_df$Gene_Name, levels = unique(cytokine_heatmap_df$Gene_Name))
cytokine_heatmap_df$Day <- factor(cytokine_heatmap_df$Day, levels = c("Day 2", "Day 5", "Day 8", "Day 28"))

cytokine_type <- c()
for(current_row_index in 1:nrow(cytokine_heatmap_df)) {
  current_row <- cytokine_heatmap_df[current_row_index,]
  cytokine_type <- c(cytokine_type, heatmap_gene_types[[as.character(current_row$Gene_Name)]])
}

cytokine_heatmap_df$Cytokine_Type <- cytokine_type

category_colors <- cytokine_heatmap_df %>% 
  distinct(Cytokine_Type, Gene_Name) %>% 
  mutate(color = case_when(
    Cytokine_Type == "Interferon" ~ "#e31a1c",
    Cytokine_Type == "Interleukin" ~ "#1f78b4",
    Cytokine_Type == "Chemokine" ~ "#33a02c",
    Cytokine_Type == "AP-1" ~ "orange",
    Cytokine_Type == "MAP Kinase" ~ "purple",
    Cytokine_Type == "JAK-STAT" ~ "#8e3563",
    Cytokine_Type == "NF-κB" ~ "#CC79A7",
    TRUE ~ "black" # default color
  )) %>% 
  pull(color, name = Gene_Name)


cytokine_heatmap_plot <- ggplot() + 
  geom_raster(data = cytokine_heatmap_df, aes(x = Day, y = Gene_Name, fill = fold_change)) +
  geom_text(data = cytokine_heatmap_df, aes(x = Day, y = Gene_Name, label = significant), nudge_y = 0.15, nudge_x = 0.20, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                                       x = NULL,
                                       y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16, color = category_colors[levels(factor(cytokine_heatmap_df$Gene_Name))])) # + coord_fixed(ratio = 0.5)

ggsave(filename = paste0("C:/Users/willi/Desktop/", "cytokine_deg_heatmap_naive_lvl_bulk.png"), plot = cytokine_heatmap_plot, device='png', dpi=300, width = 5, height = 8, units = "in")

# Flip it!
cytokine_heatmap_plot <- ggplot() + 
  geom_raster(data = cytokine_heatmap_df, aes(x = Gene_Name, y = Day, fill = fold_change)) +
  geom_text(data = cytokine_heatmap_df, aes(x = Gene_Name, y = Day, label = significant), nudge_y = 0.15, nudge_x = 0.10, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                                       x = NULL,
                                       y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16, color = category_colors[levels(factor(cytokine_heatmap_df$Gene_Name))]),
        axis.text.y = element_text(size = 16)) # + coord_fixed(ratio = 0.5)

ggsave(filename = paste0("C:/Users/willi/Desktop/", "cytokine_deg_heatmap_naive_lvl_bulk_flipped.png"), plot = cytokine_heatmap_plot, device='png', dpi=300, width = 15, height = 3, units = "in")

