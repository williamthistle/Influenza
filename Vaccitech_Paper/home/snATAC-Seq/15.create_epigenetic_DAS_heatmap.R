# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

epigenetic_remodeling_das_coords <- read.table(paste0(scATAC_hvl_placebo_das_dir, "epigenetic_remodeling_DAS_coordinates.tsv"), sep = "\t", header = TRUE)
innate_cell_types_atac <- c("CD14 Mono", "CD16 Mono", "cDC", "pDC", "NK")

cell_type_vector <- c()
gene_name_vector <- c()
chr_vector <- c()
start_vector <- c()
end_vector <- c()
fold_change_vector <- c()
significance_vector <- c()
epigenetic_remodeling_type_vector <- c()

for(current_row in 1:nrow(epigenetic_remodeling_das_coords)) {
  # Grab info associated with current site
  current_coords <- epigenetic_remodeling_das_coords[current_row,]
  current_chr <- current_coords$chr
  current_start <- current_coords$start
  current_end <- current_coords$end
  current_gene <- current_coords$associated_gene
  current_epigenetic_remodeling_type <- current_coords$epigenetic_remodeling_type
  for(innate_cell_type in innate_cell_types_atac) {
    # Read in relevant files
    innate_cell_type_for_file_name <- sub(" ", "_", innate_cell_type)
    current_scATAC_hvl_placebo_upregulated_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "annotated/D28-vs-D_minus_1-degs-", innate_cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_0.01_fc_0.1_upregulated_promoter_subset_annotated.tsv"),
                                                 sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    current_scATAC_hvl_placebo_downregulated_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "annotated/D28-vs-D_minus_1-degs-", innate_cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_0.01_fc_0.1_downregulated_promoter_subset_annotated.tsv"),
                                                             sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    current_unfiltered_hvl_placebo_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "annotated/D28-vs-D_minus_1-degs-", innate_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc_unfiltered_annotated.tsv"),
                                                     sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    # Add relevant info for plotting
    cell_type_vector <- c(cell_type_vector, innate_cell_type)
    gene_name_vector <- c(gene_name_vector, current_gene)
    chr_vector <- c(chr_vector, current_chr)
    start_vector <- c(start_vector, current_start)
    end_vector <- c(end_vector, current_end)
    epigenetic_remodeling_type_vector <- c(epigenetic_remodeling_type_vector, current_epigenetic_remodeling_type)
    unfiltered_info <- current_unfiltered_hvl_placebo_das[current_unfiltered_hvl_placebo_das$seqnames == current_chr & current_unfiltered_hvl_placebo_das$start == current_start & current_unfiltered_hvl_placebo_das$end == current_end,]
    unfiltered_info <- unfiltered_info[unfiltered_info$pct.1 > 0 & unfiltered_info$pct.2 > 0,]
    if(nrow(unfiltered_info) > 0) {
      fold_change_vector <- c(fold_change_vector, unfiltered_info$avg_log2FC)
    } else {
      fold_change_vector <- c(fold_change_vector, 0)
    }
    if(current_chr %in% current_scATAC_hvl_placebo_upregulated_das$seqnames && current_start %in% current_scATAC_hvl_placebo_upregulated_das$start && current_end %in% current_scATAC_hvl_placebo_upregulated_das$end) {
      current_info <- current_scATAC_hvl_placebo_upregulated_das[current_scATAC_hvl_placebo_upregulated_das$seqnames == current_chr & current_scATAC_hvl_placebo_upregulated_das$start == current_start & current_scATAC_hvl_placebo_upregulated_das$end == current_end,]
      current_p_value <- current_info$sc_pval
      if(current_p_value < 0.001) {
        significance_vector <- c(significance_vector, "***")
      } else if(current_p_value < 0.01) {
        significance_vector <- c(significance_vector, "**")
      } else if(current_p_value < 0.05) {
        significance_vector <- c(significance_vector, "*")
      }
    } else if(current_chr %in% current_scATAC_hvl_placebo_downregulated_das$seqnames && current_start %in% current_scATAC_hvl_placebo_downregulated_das$start && current_end %in% current_scATAC_hvl_placebo_downregulated_das$end) {
      current_info <- current_scATAC_hvl_placebo_downregulated_das[current_scATAC_hvl_placebo_downregulated_das$seqnames == current_chr & current_scATAC_hvl_placebo_downregulated_das$start == current_start & current_scATAC_hvl_placebo_downregulated_das$end == current_end,]
      current_p_value <- current_info$sc_pval
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

epigenetic_remodeling_das_heatmap_df <- data.frame(Cell_Type = cell_type_vector, Gene_Name = gene_name_vector, chr = chr_vector,
                                  start = start_vector, end = end_vector, fold_change = fold_change_vector, 
                                  significant = significance_vector, epigenetic_remodeling_type = epigenetic_remodeling_type_vector)

category_colors <- epigenetic_remodeling_das_heatmap_df %>% 
  distinct(epigenetic_remodeling_type, Gene_Name) %>% 
  mutate(color = case_when(
    epigenetic_remodeling_type == "Histone" ~ "#e31a1c",
    epigenetic_remodeling_type == "Lysine Demethylase" ~ "#1f78b4",
    epigenetic_remodeling_type == "Lysine Methyltransferase" ~ "#33a02c",
    epigenetic_remodeling_type == "AT-Rich Interaction Domain" ~ "#a6cee3",
    epigenetic_remodeling_type == "Lysine Acetyltransferase" ~ "#6a3d9a",
    epigenetic_remodeling_type == "DNA Methyltransferase" ~ "#cab2d6",
    epigenetic_remodeling_type == "Sirtuin" ~ "turquoise",
    epigenetic_remodeling_type == "Tet Methylcytosine Dioxygenase" ~ "maroon",
    TRUE ~ "black" # default color
  )) %>% 
  pull(color, name = Gene_Name)

# Approach 1: Plot barplot for significant DAS for each cell type 
significant_epigenetic_remodeling_das_heatmap_df <- epigenetic_remodeling_das_heatmap_df[!is.na(epigenetic_remodeling_das_heatmap_df$significant),]
# Remove pDC
significant_epigenetic_remodeling_das_heatmap_df <- significant_epigenetic_remodeling_das_heatmap_df[significant_epigenetic_remodeling_das_heatmap_df$Cell_Type != "pDC",]
# Remove genes
removed_genes <- c("SIRT2", "SETD7", "SETD3", "ARID3B", "ARID5B", "ARID1B")
significant_epigenetic_remodeling_das_heatmap_df <- significant_epigenetic_remodeling_das_heatmap_df[!(significant_epigenetic_remodeling_das_heatmap_df$Gene_Name %in% removed_genes),]
# Create list to store barplots
epigenetic_remodeling_das_plots <- list()

for(cell_type in unique(significant_epigenetic_remodeling_das_heatmap_df$Cell_Type)) {
  current_epigenetic_remodeling_das_heatmap_df <- significant_epigenetic_remodeling_das_heatmap_df[significant_epigenetic_remodeling_das_heatmap_df$Cell_Type == cell_type,]
  current_epigenetic_remodeling_das_heatmap_df <- current_epigenetic_remodeling_das_heatmap_df[rev(order(current_epigenetic_remodeling_das_heatmap_df$fold_change)),]
  current_epigenetic_remodeling_das_heatmap_df$Gene_Name <- factor(current_epigenetic_remodeling_das_heatmap_df$Gene_Name, levels = current_epigenetic_remodeling_das_heatmap_df$Gene_Name)
  current_epigenetic_remodeling_das_heatmap_df$epigenetic_remodeling_type <- factor(current_epigenetic_remodeling_das_heatmap_df$epigenetic_remodeling_type, levels = c("Interferon", "Interleukin", "Chemokine"))
  current_plot <- ggplot(current_epigenetic_remodeling_das_heatmap_df, aes(x = Gene_Name, y = fold_change, fill = Gene_Name)) +
    geom_col() + geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
    geom_text(aes(label = significant,
                  vjust = ifelse(fold_change >= 0, 0.4, 1)),
              size = 10) + theme_bw() + ylim(-5, 5) +
    scale_fill_manual(values = category_colors) + 
    labs(title = cell_type,
         x = NULL,
         y = "Fold Change", fill = "Epigenetic Remodeling Type") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1), text=element_text(size=32))
  
  epigenetic_remodeling_das_plots[[cell_type]] <- current_plot
}

ggsave("C:/Users/willi/Desktop/epigenetic_remodeling_barplot_without_legend.png", plot = patchwork::wrap_plots(epigenetic_remodeling_das_plots, ncol = 2, nrow = 2), height = 10, width = 14)

# Approach 2: Plot heatmap for all DAS across all cell types
cytokine_das_heatmap_df$Gene_Name <- factor(cytokine_das_heatmap_df$Gene_Name, levels = c("CCL17", "CCR7", "CXCL12", "CXCR3", "XCR1", "IRF1", "IRF4", "IRF5", "IRF8", "IRF2BPL", "IFNGR1", "IFI6", "IFI27L1", "IFI35", "IFI44L", "IFITM3", "IL1B", "IL20", "IL2RB", "IL5RA", "IL6R", "IL15RA", "IL17RA", "IL21R", "IL1RN", "IRAK1BP1"))
cytokine_das_heatmap_df$Cell_Type <- factor(cytokine_das_heatmap_df$Cell_Type, levels = c("CD14 Mono", "CD16 Mono", "NK", "cDC", "pDC"))

cytokine_das_heatmap <- ggplot() + 
  geom_raster(data = cytokine_das_heatmap_df, aes(x = Cell_Type, y = Gene_Name, fill = fold_change)) +
  geom_text(data = cytokine_das_heatmap_df, aes(x = Cell_Type, y = Gene_Name, label = significant), nudge_y = 0.15, nudge_x = 0.30, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = "Fold Change for Cytokine-Related DAS in Innate Immune Cell Types",
                                       x = "Cell Type",
                                       y = "Gene", fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12, color = category_colors[levels(factor(cytokine_das_heatmap_df$Gene_Name))])) + coord_fixed(ratio = 0.4)

ggsave("C:/Users/willi/Desktop/cytokine_das_heatmap.png", plot = cytokine_das_heatmap, height = 10, width = 10)
