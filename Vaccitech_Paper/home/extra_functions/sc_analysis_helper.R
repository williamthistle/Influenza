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