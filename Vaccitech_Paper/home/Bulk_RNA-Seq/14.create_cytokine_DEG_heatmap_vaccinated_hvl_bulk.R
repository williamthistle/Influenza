# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

bulk_hvl_vaccinated_degs <- list()
bulk_hvl_vaccinated_degs[["Day 2"]] <- hvl_full_time_series_vaccinated_period_2_D2_vs_D_minus_1_results
bulk_hvl_vaccinated_degs[["Day 5"]] <- hvl_full_time_series_vaccinated_period_2_D5_vs_D_minus_1_results
bulk_hvl_vaccinated_degs[["Day 8"]] <- hvl_full_time_series_vaccinated_period_2_D8_vs_D_minus_1_results
bulk_hvl_vaccinated_degs[["Day 28"]] <- hvl_full_time_series_vaccinated_period_2_D28_vs_D_minus_1_results

heatmap_genes <- c("CX3CR1", "CCR2", "CCL3", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR4",
                   "CSF1R", "NFIL3", "IL10", "IL10RA", "IL1R1", "IL1R2", "IL1RAP",
                   "IFNG", "IRF1", "IRF2", "IRF8", "OAS1", "IFI16", "IFIH1", "IFITM1", "IFNAR1", "IFNGR1", "IRF2BP2", 
                   "FOS", "JUN", "FOSL2", "JDP2",
                   "DUSP1", "MAPK14", "MAP2K3", "MAP3K8", "MAPK6", "MAPK8", "MAPKAPK2", "DUSP4",
                   "NFKBIA", "NFKBIZ",
                   "JAK2", "STAT1")

heatmap_gene_types <- list()
heatmap_gene_types[["CCL3"]] <- "Chemokine"
heatmap_gene_types[["CCR2"]] <- "Chemokine"
heatmap_gene_types[["CX3CR1"]] <- "Chemokine"
heatmap_gene_types[["CXCL16"]] <- "Chemokine"
heatmap_gene_types[["CXCL2"]] <- "Chemokine"
heatmap_gene_types[["CXCL3"]] <- "Chemokine"
heatmap_gene_types[["CXCL8"]] <- "Chemokine"
heatmap_gene_types[["CXCR4"]] <- "Chemokine"

heatmap_gene_types[["CSF1R"]] <- "Interleukin"
heatmap_gene_types[["NFIL3"]] <- "Interleukin"
heatmap_gene_types[["IL10"]] <- "Interleukin"
heatmap_gene_types[["IL10RA"]] <- "Interleukin"
heatmap_gene_types[["IL1R1"]] <- "Interleukin"
heatmap_gene_types[["IL1R2"]] <- "Interleukin"
heatmap_gene_types[["IL1RAP"]] <- "Interleukin"

heatmap_gene_types[["OAS1"]] <- "Interferon"
heatmap_gene_types[["IFI16"]] <- "Interferon"
heatmap_gene_types[["IFIH1"]] <- "Interferon"
heatmap_gene_types[["IFITM1"]] <- "Interferon"
heatmap_gene_types[["IFNAR1"]] <- "Interferon"
heatmap_gene_types[["IFNG"]] <- "Interferon"
heatmap_gene_types[["IFNGR1"]] <- "Interferon"
heatmap_gene_types[["IRF1"]] <- "Interferon"
heatmap_gene_types[["IRF2"]] <- "Interferon"
heatmap_gene_types[["IRF2BP2"]] <- "Interferon"
heatmap_gene_types[["IRF8"]] <- "Interferon"

heatmap_gene_types[["FOS"]] <- "AP-1"
heatmap_gene_types[["FOSL2"]] <- "AP-1"
heatmap_gene_types[["JUN"]] <- "AP-1"
heatmap_gene_types[["JDP2"]] <- "AP-1"

heatmap_gene_types[["DUSP1"]] <- "MAP Kinase"
heatmap_gene_types[["DUSP4"]] <- "MAP Kinase"
heatmap_gene_types[["MAP2K3"]] <- "MAP Kinase"
heatmap_gene_types[["MAP3K8"]] <- "MAP Kinase"
heatmap_gene_types[["MAPK14"]] <- "MAP Kinase"
heatmap_gene_types[["MAPK6"]] <- "MAP Kinase"
heatmap_gene_types[["MAPK8"]] <- "MAP Kinase"
heatmap_gene_types[["MAPKAPK2"]] <- "MAP Kinase"

heatmap_gene_types[["NFKBIA"]] <- "NF-κB"
heatmap_gene_types[["NFKBIZ"]] <- "NF-κB"

heatmap_gene_types[["JAK2"]] <- "JAK-STAT"
heatmap_gene_types[["STAT1"]] <- "JAK-STAT"

bulk_day_vector <- c()
gene_name_vector <- c()
fold_change_vector <- c()
significance_vector <- c()

for(current_day in names(bulk_hvl_vaccinated_degs)) {
  current_bulk_degs <- bulk_hvl_vaccinated_degs[[current_day]][[2]]
  current_unfiltered_bulk_degs <- bulk_hvl_vaccinated_degs[[current_day]][[8]]
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

ggsave(filename = paste0("C:/Users/wat2/Desktop/", "cytokine_deg_heatmap_vaccinated_hvl_bulk.png"), plot = cytokine_heatmap_plot, device='png', dpi=300, width = 5, height = 8, units = "in")

# Flip it!
cytokine_heatmap_plot <- ggplot() + 
  geom_raster(data = cytokine_heatmap_df, aes(x = Gene_Name, y = Day, fill = fold_change)) +
  geom_text(data = cytokine_heatmap_df, aes(x = Gene_Name, y = Day, label = significant), nudge_y = 0.15, nudge_x = 0.10, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                                       x = NULL,
                                       y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16, color = category_colors[levels(factor(cytokine_heatmap_df$Gene_Name))]),
        axis.text.y = element_text(size = 16)) # + coord_fixed(ratio = 0.5)

ggsave(filename = paste0("C:/Users/wat2/Desktop/", "cytokine_deg_heatmap_vaccinated_hvl_bulk_flipped.png"), plot = cytokine_heatmap_plot, device='png', dpi=300, width = 15, height = 3, units = "in")

