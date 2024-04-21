# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

find_sc_and_bulk_correlation <- function(first_deg_subset, second_deg_subset, cell_type) {
  correlation_return_vals <- list()
  print(paste0("Number of sc DEGs in first set of data is: ", nrow(first_deg_subset)))
  compare_first_df <- first_deg_subset
  compare_second_df <- second_deg_subset[rownames(second_deg_subset) %in% compare_first_df$Gene_Name,]
  compare_first_df <- compare_first_df[compare_first_df$Gene_Name %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% compare_first_df$Gene_Name,]
  compare_first_df <- compare_first_df[order(compare_first_df$Gene_Name),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$pseudo_bulk_log2FC,
                                             second_fc = compare_second_df$log2FoldChange)
  # Calculate correlation - Pearson
  print("Calculating Pearson correlation")
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  correlation_return_vals[[1]] <- correlation_val$estimate
  print("Calculating Spearman correlation")
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc, method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  correlation_return_vals[[2]] <- correlation_val$estimate
  
  # Plot correlation
  correlation_plot <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("Single Cell FC") + ylab("Bulk FC") + labs(title = cell_type) + xlim(-1.5, 1.5) + ylim(-1.5, 1.5)
  correlation_return_vals[[3]] <- correlation_plot
  
  return(correlation_return_vals)
}


find_sc_and_bulk_correlations <- function(sc_degs, bulk_degs, cell_type) {
  # Unfiltered bulk
  bulk_degs_unfiltered <- bulk_degs[[8]]
  # Find correlation
  sc_degs_vs_bulk_degs_unfiltered_corr <- find_sc_and_bulk_correlation(sc_degs, bulk_degs_unfiltered, cell_type)
  
  return(sc_degs_vs_bulk_degs_unfiltered_corr)
}

correlation_cell_types <- c("CD14 Mono", "CD16 Mono", "cDC", "NK", "CD4 Memory", "CD8 Memory", "MAIT", "B naive", "B memory")

# D28
sc_correlations_D28 <- list()
for(correlation_cell_type in correlation_cell_types) {
  current_corr <- find_sc_and_bulk_correlations(scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == correlation_cell_type,], 
                                                hvl_placebo_period_2_D28_vs_D_minus_1_results, correlation_cell_type)
  sc_correlations_D28[[correlation_cell_type]] <- current_corr
}

sc_correlation_plots <- lapply(sc_correlations_D28, function(x) x[[3]])
ggsave("C:/Users/willi/Desktop/sc_D28.png", plot = patchwork::wrap_plots(sc_correlation_plots, ncol = 3, nrow = 3), height = 10, width = 10)

# D8
sc_correlations_D8 <- list()
for(correlation_cell_type in correlation_cell_types) {
  current_corr <- find_sc_and_bulk_correlations(scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == correlation_cell_type,], 
                                                hvl_placebo_period_2_D8_vs_D_minus_1_results, correlation_cell_type)
  sc_correlations_D8[[correlation_cell_type]] <- current_corr
}

sc_correlation_plots <- lapply(sc_correlations_D8, function(x) x[[3]])
ggsave("C:/Users/willi/Desktop/sc_D8.png", plot = patchwork::wrap_plots(sc_correlation_plots, ncol = 3, nrow = 3), height = 10, width = 10)

# D5
sc_correlations_D5 <- list()
for(correlation_cell_type in correlation_cell_types) {
  current_corr <- find_sc_and_bulk_correlations(scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == correlation_cell_type,], 
                                                hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results, correlation_cell_type)
  sc_correlations_D5[[correlation_cell_type]] <- current_corr
}

sc_correlation_plots <- lapply(sc_correlations_D5, function(x) x[[3]])
ggsave("C:/Users/willi/Desktop/sc_D5.png", plot = patchwork::wrap_plots(sc_correlation_plots, ncol = 3, nrow = 3), height = 10, width = 10)

# D2
sc_correlations_D2 <- list()
for(correlation_cell_type in correlation_cell_types) {
  current_corr <- find_sc_and_bulk_correlations(scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == correlation_cell_type,], 
                                                hvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results, correlation_cell_type)
  sc_correlations_D2[[correlation_cell_type]] <- current_corr
}

sc_correlation_plots <- lapply(sc_correlations_D2, function(x) x[[3]])
ggsave("C:/Users/willi/Desktop/sc_D2.png", plot = patchwork::wrap_plots(sc_correlation_plots, ncol = 3, nrow = 3), height = 10, width = 10)

# D28 MVL
sc_correlations_D28_MVL <- list()
for(correlation_cell_type in correlation_cell_types) {
  current_corr <- find_sc_and_bulk_correlations(scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == correlation_cell_type,], 
                                                mvl_placebo_period_2_D28_vs_D_minus_1_results, correlation_cell_type)
  sc_correlations_D28_MVL[[correlation_cell_type]] <- current_corr
}

sc_correlation_plots <- lapply(sc_correlations_D28_MVL, function(x) x[[3]])
ggsave("C:/Users/willi/Desktop/sc_D28_MVL.png", plot = patchwork::wrap_plots(sc_correlation_plots, ncol = 3, nrow = 3), height = 10, width = 10)

# D28 LVL
sc_correlations_D28_LVL <- list()
for(correlation_cell_type in correlation_cell_types) {
  current_corr <- find_sc_and_bulk_correlations(scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == correlation_cell_type,], 
                                                lvl_placebo_period_2_D28_vs_D_minus_1_results, correlation_cell_type)
  sc_correlations_D28_LVL[[correlation_cell_type]] <- current_corr
}

sc_correlation_plots <- lapply(sc_correlations_D28_LVL, function(x) x[[3]])
ggsave("C:/Users/willi/Desktop/sc_D28_LVL.png", plot = patchwork::wrap_plots(sc_correlation_plots, ncol = 3, nrow = 3), height = 10, width = 10)

# D28 HVL without SC subjects
sc_correlations_D28_without_sc_subjects <- list()
for(correlation_cell_type in correlation_cell_types) {
  current_corr <- find_sc_and_bulk_correlations(scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == correlation_cell_type,], 
                                                hvl_placebo_period_2_D28_vs_D_minus_1_results_without_sc_subjects, correlation_cell_type)
  sc_correlations_D28[[correlation_cell_type]] <- current_corr
}

sc_correlation_plots <- lapply(sc_correlations_D28, function(x) x[[3]])
ggsave("C:/Users/willi/Desktop/sc_D28_without_sc_subjects.png", plot = patchwork::wrap_plots(sc_correlation_plots, ncol = 3, nrow = 3), height = 10, width = 10)