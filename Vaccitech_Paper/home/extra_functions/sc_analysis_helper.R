combine_sc_deg_results <- function(input_dir) {
  file_paths <- list.files(input_dir, pattern = "subject_id_final\\.tsv$", full.names = TRUE)
  files <- list()
  for(file_path in file_paths) {
    current_table <- read.table(file_path, sep = "\t", header = TRUE)
    files[[file_path]] <- current_table
  }
  final_table <- do.call(rbind, files)
  rownames(final_table) <- 1:nrow(final_table)
  write.table(final_table, file = paste0(input_dir, "D28-vs-D_minus_1-degs-time_point.final.list.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(final_table[final_table$sc_log2FC > 0,], file = paste0(input_dir, "D28-vs-D_minus_1-degs-time_point.final.pos.list.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(final_table[final_table$sc_log2FC < 0,], file = paste0(input_dir, "D28-vs-D_minus_1-degs-time_point.final.neg.list.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  return(final_table)
}

# Evaluate cell type proportion changes in scRNA-seq data
evaluate_scRNA_cell_type_proportion_changes <- function(cell_type_prop_table) {
  cell_type_prop_table$Condition <- factor(cell_type_prop_table$Condition, levels = c("D_minus_1", "D28"))
  # Find unadjusted p-values for each cell type of interest
  # We use coin::wilcox_test because it handles ties (which happen when there are 0s in the cell type proportions)
  # Unfortunately, we can't use loops because of the way the function call works (I think)
  cell_type_proportion_p_values <- c()
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD16_Mono ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD8_Naive ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD14_Mono ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD4_Naive ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD4_Memory ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(B ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD8_Memory ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(NK ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(MAIT ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(cDC ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(pDC ~ Condition, data = cell_type_prop_table)))
  # Adjust for multiple hypothesis testing
  cell_type_proportion_p_values_adjusted <- p.adjust(cell_type_proportion_p_values, method = "BH")
  names(cell_type_proportion_p_values_adjusted) <- c("CD16 Mono", "CD8 Naive", "CD14 Mono", "CD4 Naive", "CD4 Memory", "B", "CD8 Memory",
                                                     "NK", "MAIT", "cDC", "pDC")
  return(cell_type_proportion_p_values_adjusted)
}

# Evaluate cell type proportion changes in scATAC-seq data
evaluate_scATAC_cell_type_proportion_changes <- function(cell_type_prop_table) {
  cell_type_prop_table$Condition <- factor(cell_type_prop_table$Condition, levels = c("D_minus_1", "D28"))
  # Find unadjusted p-values for each cell type of interest
  # We use coin::wilcox_test because it handles ties (which happen when there are 0s in the cell type proportions)
  # Unfortunately, we can't use loops because of the way the function call works (I think)
  cell_type_proportion_p_values <- c()
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD16_Mono ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD8_Naive ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD14_Mono ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD4_Naive ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD4_Memory ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(B ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(CD8_Memory ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(NK ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(MAIT ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(cDC ~ Condition, data = cell_type_prop_table)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(pDC ~ Condition, data = cell_type_prop_table)))
  # Adjust for multiple hypothesis testing
  cell_type_proportion_p_values_adjusted <- p.adjust(cell_type_proportion_p_values, method = "BH")
  names(cell_type_proportion_p_values_adjusted) <- c("CD16 Mono", "CD8 Naive", "CD14 Mono", "CD4 Naive", "CD4 Memory", "B", "CD8 Memory",
                                                     "NK", "MAIT", "cDC", "pDC")
  return(cell_type_proportion_p_values_adjusted)
}
