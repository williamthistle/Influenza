base_dir <- "~/"
scRNA_data_dir <- paste0(base_dir, "scRNA/")
scATAC_data_dir <- paste0(base_dir, "scATAC/")
multiome_data_dir <- paste0(base_dir, "multiome/")
scRNA_qc_file <- paste0(base_dir, "ECHO_FLU_Vaccitech_PBMC_scrnaseq_coded_qc_report_WT.csv")
multiome_qc_file <- paste0(base_dir, "Stanford_FLU_combined_qc_metric_coded_09015022_qc_data.csv")
overall_metadata_file <- paste0(base_dir, "20220609_metadataECHO_Vaccitech_Coded.csv")
# Grab scRNA aliquot names
scRNA_aliquot_list <- list.dirs(scRNA_data_dir, recursive = FALSE)
scRNA_aliquot_list <- strsplit(scRNA_aliquot_list, "/")
scRNA_aliquot_list <- unlist(lapply(scRNA_aliquot_list, tail, n = 1L))
# Grab scATAC aliquot names
scATAC_aliquot_list <- list.dirs(scATAC_data_dir, recursive = FALSE)
scATAC_aliquot_list <- strsplit(scATAC_aliquot_list, "/")
scATAC_aliquot_list <- unlist(lapply(scATAC_aliquot_list, tail, n = 1L))
# Grab multiome aliquot names
multiome_aliquot_list <- list.dirs(multiome_data_dir, recursive = FALSE)
multiome_aliquot_list <- strsplit(multiome_aliquot_list, "/")
multiome_aliquot_list <- unlist(lapply(multiome_aliquot_list, tail, n = 1L))
# Read in tables
scRNA_qc <- read.csv(scRNA_qc_file)
multiome_qc <- read.csv(multiome_qc_file)
overall_metadata <- read.csv(overall_metadata_file)
