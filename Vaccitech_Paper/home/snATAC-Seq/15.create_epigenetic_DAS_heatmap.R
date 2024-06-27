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





# Alternative - all of the sites
epigenetic_remodeling_das_heatmap <- ggplot() + 
  geom_raster(data = epigenetic_remodeling_das_heatmap_df, aes(x = Cell_Type, y = Gene_Name, fill = fold_change)) +
  geom_text(data = epigenetic_remodeling_das_heatmap_df, aes(x = Cell_Type, y = Gene_Name, label = significant), nudge_y = 0.15, nudge_x = 0.30, size = 4) + scale_fill_gradient2(low="navy", mid="white", high="red") +
  theme_classic(base_size = 14) + labs(title = "Fold Change for Promoter Sites Found in Epigenetic Remodeling Genes in Innate Immune Cell Types",
                                       x = "Cell Type",
                                       y = "Gene", fill = "Fold Change") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12)) + coord_fixed(ratio = 0.40)

ggsave("C:/Users/willi/Desktop/epigenetic_remodeling_das_heatmap.png", plot = epigenetic_remodeling_das_heatmap, height = 10, width = 10)



# epigenetic_remodeling_das_heatmap_df <- epigenetic_remodeling_das_heatmap_df[!is.na(epigenetic_remodeling_das_heatmap_df$significant),]

# epigenetic_remodeling_das_heatmap_df$epigenetic_remodeling_type <- c("Histone", "Histone", "Histone Modification",
#                                            "Histone Modification", "DNA Methylation", "Histone Modification",
#                                            "Histone Modification", "Histone Modification", "Chromatin Remodeling",
#                                            "Chromatin Remodeling", "Histone", "Histone Modification",
#                                            "Histone Modification", "Histone Modification", "Chromatin Remodeling", "Histone Modification",
#                                            "Histone", "Histone", "Histone Modification",
#                                            "Histone Modification", "Histone Modification", "DNA Methylation",
#                                            "Chromatin Remodeling", "Histone", "Histone Modification",
#                                            "Histone", "Histone", "Histone Modification",
#                                            "Histone Modification", "Histone Modification", "Histone Modification",
#                                            "DNA Methylation", "Chromatin Remodeling", "Histone Modification",
#                                            "Histone", "Histone Modification", "Histone Modification",
#                                            "Histone Modification",  "Histone Modification", "Histone",
#                                            "Histone", "Histone Modification")
# 
# # Temporary removal of one row
# # epigenetic_remodeling_das_heatmap_df <- epigenetic_remodeling_das_heatmap_df[-c(31),]
# 
# epigenetic_remodeling_das_plots <- list()
# 
# for(cell_type in unique(epigenetic_remodeling_das_heatmap_df$Cell_Type)) {
#   current_epigenetic_remodeling_das_heatmap_df <- epigenetic_remodeling_das_heatmap_df[epigenetic_remodeling_das_heatmap_df$Cell_Type == cell_type,]
#   current_epigenetic_remodeling_das_heatmap_df <- current_epigenetic_remodeling_das_heatmap_df[rev(order(current_epigenetic_remodeling_das_heatmap_df$fold_change)),]
#   current_epigenetic_remodeling_das_heatmap_df$Gene_Name <- factor(current_epigenetic_remodeling_das_heatmap_df$Gene_Name, levels = current_epigenetic_remodeling_das_heatmap_df$Gene_Name)
#   current_epigenetic_remodeling_das_heatmap_df$epigenetic_remodeling_type <- factor(current_epigenetic_remodeling_das_heatmap_df$epigenetic_remodeling_type, levels = c("Histone", "Histone Modification", "DNA Methylation", "Chromatin Remodeling"))
#   if(length(unique(current_epigenetic_remodeling_das_heatmap_df$epigenetic_remodeling_type)) == 2) {
#     barplot_colors <- c("#F8766D", "#00BA38")
#   } else if(length(unique(current_epigenetic_remodeling_das_heatmap_df$epigenetic_remodeling_type)) == 3) {
#     barplot_colors <- c("#F8766D", "#00BA38", "#619CFF")
#   } else {
#     barplot_colors <- c("#F8766D", "#00BA38", "#619CFF", "#7cae00")
#   }
#   
#   current_plot <- ggplot(current_epigenetic_remodeling_das_heatmap_df, aes(x = Gene_Name, y = fold_change, fill = epigenetic_remodeling_type)) +
#     geom_col() + geom_text(aes(label = significant,
#                           vjust = ifelse(fold_change >= 0, 0, 1)),
#                            size = 6) + theme_bw() + ylim(-5, 5) + theme(text = element_text(size = 16)) +
#     scale_fill_manual(values = barplot_colors) + 
#     labs(title = cell_type,
#          x = "Gene Name",
#          y = "Fold Change", fill = "Epigenetic Remodeling Gene Type") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position="none")
#   
#   epigenetic_remodeling_das_plots[[cell_type]] <- current_plot
# }
# 
# ggsave("C:/Users/willi/Desktop/test_with_legend.png", plot = patchwork::wrap_plots(epigenetic_remodeling_das_plots, ncol = 2, nrow = 3), height = 10, width = 20)
#   
