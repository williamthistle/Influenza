# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Read in data
hvl_placebo_cell_type_props <- read.table(paste0(scATAC_hvl_dir, "Placebo/HVL_placebo_ATAC_cell_type_proportion_time_point.csv"), sep = ",", header = TRUE)
hvl_vaccinated_cell_type_props <- read.table(paste0(scATAC_hvl_dir, "Vaccinated/HVL_vaccinated_ATAC_cell_type_proportion_time_point.csv"), sep = ",", header = TRUE)
lvl_placebo_cell_type_props <- read.table(paste0(scATAC_lvl_dir, "Placebo/LVL_placebo_ATAC_cell_type_proportion_time_point.csv"), sep = ",", header = TRUE)

# Find adjusted p-values
hvl_placebo_cell_type_props_results <- evaluate_scATAC_cell_type_proportion_changes(hvl_placebo_cell_type_props)
hvl_vaccinated_cell_type_props_results <- evaluate_scATAC_cell_type_proportion_changes(hvl_vaccinated_cell_type_props)
lvl_placebo_cell_type_props_results <- evaluate_scATAC_cell_type_proportion_changes(lvl_placebo_cell_type_props)