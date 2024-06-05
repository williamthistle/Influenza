### THIS FUNCTION:
### RUNS MAGICAL AND CREATES 1) OVERALL MAGICAL DF, AND 2) TF-CENTRIC MAGICAL DF

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))
source(paste0(base_dir, "extra_functions/my_MAGICAL_helper_functions.R"))

overall_magical_df <- read.table(paste0(MAGICAL_hvl_placebo_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"), sep = "\t", header = TRUE)
magical_gene_overlap_df <- create_magical_gene_overlap_df(overall_magical_df, hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results, hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results)
magical_tf_vs_gene_df <- create_tf_targets_df(overall_magical_df)

# MAP3K11
current_cell_type <- "CD14 Mono"
current_cell_type_for_file_name <- gsub(" ", "_", current_cell_type)
MAP3K11_magical_info <- overall_magical_df[overall_magical_df$Cell_Type == current_cell_type_for_file_name,]
MAP3K11_magical_info <- MAP3K11_magical_info[MAP3K11_magical_info$Gene_symbol == "MAP3K11",]
site_plots <- list()
for(current_row_index in 1:nrow(MAP3K11_magical_info)) {
  current_row <- MAP3K11_magical_info[current_row_index,]
  # Parse out list of TFs
  current_row <- current_row %>%
    mutate(TFs.binding.prob. = strsplit(as.character(TFs.binding.prob.), ",")) %>%
    unnest(TFs.binding.prob.) %>%
    mutate(TFs.binding.prob. = trimws(TFs.binding.prob.)) %>%
    filter(!is.na(TFs.binding.prob.))
  current_row <- subset(current_row, !grepl("^\\(", TFs.binding.prob.))
  list_of_current_tfs <- current_row$TFs.binding.prob.
  gene_name_vector <- c()
  fold_change_vector <- c()
  significance_vector <- c()
  current_scRNA_hvl_placebo_degs <- innate_scRNA_hvl_placebo_degs[innate_scRNA_hvl_placebo_degs$Cell_Type == current_cell_type,]
  current_unfiltered_hvl_placebo_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", current_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                      sep = "\t", header = TRUE)
  for(current_tf in list_of_current_tfs) {
    current_gene_unfiltered_data <- current_unfiltered_hvl_placebo_degs[rownames(current_unfiltered_hvl_placebo_degs) == current_tf,]
    gene_name_vector <- c(gene_name_vector, current_tf)
    fold_change_vector <- c(fold_change_vector, current_gene_unfiltered_data$avg_log2FC)
    if(current_tf %in% current_scRNA_hvl_placebo_degs$Gene_Name) {
      current_p_value <- current_scRNA_hvl_placebo_degs[current_scRNA_hvl_placebo_degs$Gene_Name == current_tf,]$sc_pval_adj
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
  current_site_df <- data.frame(TF = gene_name_vector, Fold_Change = fold_change_vector, significant = significance_vector)
  
  current_site_df <- current_site_df[rev(order(current_site_df$Fold_Change)),]
  current_site_df$TF <- factor(current_site_df$TF, levels = current_site_df$TF)
  current_site_df$Fold_Change_Direction <- ifelse(current_site_df$Fold_Change > 0, "Positive", "Negative")
  current_site_df$Fold_Change_Direction <- factor(current_site_df$Fold_Change_Direction, levels = c("Positive", "Negative"))
  
  current_plot <- ggplot(current_site_df, aes(x = TF, y = Fold_Change, fill = Fold_Change_Direction)) +
    geom_col() + geom_text(aes(label = significant,
                               vjust = ifelse(Fold_Change >= 0, 0, 1)),
                           size = 6) + theme_bw() + ylim(-1, 1) + theme(text = element_text(size = 16)) +
    labs(title = paste0("Site ", current_row_index),
         x = "Gene Name",
         y = "Fold Change", fill = "Fold Change Direction") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="none")
  
  site_plots[[current_row_index]] <- current_plot
}

ggsave("C:/Users/willi/Desktop/MAP3K11_TF_plots.png", plot = patchwork::wrap_plots(site_plots, ncol = 2, nrow = 2), height = 10, width = 14)



  
 
