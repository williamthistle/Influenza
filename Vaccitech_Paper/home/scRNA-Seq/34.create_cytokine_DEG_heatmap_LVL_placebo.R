# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

innate_scRNA_lvl_placebo_degs <- scRNA_lvl_placebo_degs[scRNA_lvl_placebo_degs$Cell_Type %in% innate_cell_types,]

# Find cytokine DEGs
search_terms <- c("interferon", "interleukin", "chemokine", "AP-1")
gene_terms <- c("CCL3", "CX3CR1", "CCL3L1", "CXCL16", 
                "IL32", "CASP1", "CSF1R", "NFIL3", "IRAK3", "IL1RAP", "PTGES",   
                "IRF2", "IRF7", "IFNG", "OAS1", "MNDA", "PSMB9", "IFNGR1", "DNAJC3", "GBP5", "USP38",  
                "JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "JDP2", 
                "MAP3K11", "CSK", "DUSP1", "DUSP2", "DUSP6", "ABHD17A", "RIPK1", "MAPK7", "MAPK8", "MAP2K1", "MAP3K8", "MAP3K20", "MAP4K3", "MAPKAPK2",
                "PTK2B", "RELL1", "MINK1", "BRAF",
                "NMRAL1", "NFKBIA", "NFKBIZ", "NFKB1",
                "JAK1", "STAT3", "STAT4",
                "IRF", "STAT", "JAK", "MAP", "DUSP", "IRAK")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c('hgnc_symbol', 'description', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'hgnc_symbol',
  values = unique(innate_scRNA_lvl_placebo_degs$Gene_Name),
  mart = ensembl
)

annotations_final <- annotations %>%
  filter(
    str_detect(hgnc_symbol, paste(gene_terms, collapse = "|")) |
      str_detect(description, paste(search_terms, collapse = "|"))
  )

# Final heatmap genes
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

cell_type_vector <- c()
gene_name_vector <- c()
fold_change_vector <- c()
significance_vector <- c()

for(innate_cell_type in innate_cell_types) {
  innate_cell_type_for_file_name <- sub(" ", "_", innate_cell_type)
  current_scRNA_lvl_placebo_degs <- innate_scRNA_lvl_placebo_degs[innate_scRNA_lvl_placebo_degs$Cell_Type == innate_cell_type,]
  current_unfiltered_lvl_placebo_degs <- read.table(paste0(scRNA_lvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", innate_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                       sep = "\t", header = TRUE)
  for(heatmap_gene in heatmap_genes) {
    current_heatmap_gene_unfiltered_data <- current_unfiltered_lvl_placebo_degs[rownames(current_unfiltered_lvl_placebo_degs) == heatmap_gene,]
    cell_type_vector <- c(cell_type_vector, innate_cell_type)
    gene_name_vector <- c(gene_name_vector, heatmap_gene)
    fold_change_vector <- c(fold_change_vector, current_heatmap_gene_unfiltered_data$avg_log2FC)
    if(heatmap_gene %in% current_scRNA_lvl_placebo_degs$Gene_Name) {
      current_p_value <- current_scRNA_lvl_placebo_degs[current_scRNA_lvl_placebo_degs$Gene_Name == heatmap_gene,]$sc_pval_adj
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
cytokine_heatmap_df$Cell_Type <- factor(cytokine_heatmap_df$Cell_Type, levels = c("CD14 Mono", "CD16 Mono", "cDC", "pDC", "NK", "NK_CD56bright"))

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
  geom_raster(data = cytokine_heatmap_df, aes(x = Cell_Type, y = Gene_Name, fill = fold_change)) +
  geom_text(data = cytokine_heatmap_df, aes(x = Cell_Type, y = Gene_Name, label = significant), nudge_y = 0.15, nudge_x = 0.20, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                                       x = NULL,
                                       y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16, color = category_colors[levels(factor(cytokine_heatmap_df$Gene_Name))])) # + coord_fixed(ratio = 0.5)

ggsave(filename = paste0("C:/Users/willi/Desktop/", "lvl_placebo_cytokine_deg_heatmap.png"), plot = cytokine_heatmap_plot, device='png', dpi=300, width = 5, height = 8, units = "in")

# Flip it!
cytokine_heatmap_plot <- ggplot() + 
  geom_raster(data = cytokine_heatmap_df, aes(x = Gene_Name, y = Cell_Type, fill = fold_change)) +
  geom_text(data = cytokine_heatmap_df, aes(x = Gene_Name, y = Cell_Type, label = significant), nudge_y = 0.15, nudge_x = 0.10, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = NULL,
                                       x = NULL,
                                       y = NULL, fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16, color = category_colors[levels(factor(cytokine_heatmap_df$Gene_Name))]),
        axis.text.y = element_text(size = 16)) # + coord_fixed(ratio = 0.5)

ggsave(filename = paste0("C:/Users/willi/Desktop/", "lvl_placebo_cytokine_deg_heatmap_flipped.png"), plot = cytokine_heatmap_plot, device='png', dpi=300, width = 15, height = 3, units = "in")

