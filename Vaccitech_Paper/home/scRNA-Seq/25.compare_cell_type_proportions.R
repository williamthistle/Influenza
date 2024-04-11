# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

hvl_placebo_cell_type_props_unfiltered <- read.table(paste0(scRNA_hvl_dir, "PLACEBO/RNA_cell_type_proportion_time_point_HVL_PLACEBO_unfiltered.csv"), sep = ",", header = TRUE)
hvl_placebo_cell_type_props <- read.table(paste0(scRNA_hvl_dir, "PLACEBO/RNA_cell_type_proportion_time_point_HVL_PLACEBO.csv"), sep = ",", header = TRUE)

hvl_placebo_cell_type_props_unfiltered_results <- evaluate_scRNA_cell_type_proportion_changes(hvl_placebo_cell_type_props_unfiltered)
hvl_placebo_cell_type_props_results <- evaluate_scRNA_cell_type_proportion_changes(hvl_placebo_cell_type_props)

