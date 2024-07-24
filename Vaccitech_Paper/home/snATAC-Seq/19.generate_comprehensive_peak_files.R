# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

ATAC_result_dirs <- c(scATAC_hvl_placebo_das_dir, scATAC_hvl_vaccinated_das_dir, scATAC_lvl_placebo_das_dir)
ATAC_innate_cell_types <- c("CD14_Mono", "CD16_Mono", "cDC", "pDC", "NK")

# Generate pseudobulk corrected table
for(ATAC_result_dir in ATAC_result_dirs) {
  current_ATAC_result_files <- list()
  for(ATAC_innate_cell_type in ATAC_innate_cell_types) {
    current_ATAC_result_files[[ATAC_innate_cell_type]] <- read.table(paste0(ATAC_result_dir, "D28-vs-D_minus_1-degs-", ATAC_innate_cell_type, "-time_point-controlling_for_subject_id_final_pct_0.01.tsv"),
                                                                     sep = "\t", header = TRUE)
  }
  combined_ATAC_results <- do.call(rbind, current_ATAC_result_files)
  write.table(combined_ATAC_results, file = paste0(ATAC_result_dir, "combined_final_results.tsv"), sep = "\t", row.names = FALSE,
              quote = FALSE)
}

# Generate SC table
for(ATAC_result_dir in ATAC_result_dirs) {
  current_ATAC_result_files <- list()
  for(ATAC_innate_cell_type in ATAC_innate_cell_types) {
    current_ATAC_result_files[[ATAC_innate_cell_type]] <- read.table(paste0(ATAC_result_dir, "D28-vs-D_minus_1-degs-", ATAC_innate_cell_type, "-time_point-controlling_for_subject_id_sc_pct_0.01.tsv"),
                                                                     sep = "\t", header = TRUE)
  }
  combined_ATAC_results <- do.call(rbind, current_ATAC_result_files)
  write.table(combined_ATAC_results, file = paste0(ATAC_result_dir, "combined_sc_results.tsv"), sep = "\t",
              quote = FALSE)
}
