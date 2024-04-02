combine_sc_deg_results <- function(input_dir) {
  file_paths <- list.files(input_dir, pattern = "subject_id_final\\.tsv$", full.names = TRUE)
  files <- list()
  for(file_path in file_paths) {
    files <- append(files, read.table(file_path, sep = "\t", header = TRUE))
  }
  final_table <- do.call(rbind, files)
  return(final_table)
}