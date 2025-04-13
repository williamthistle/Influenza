# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Read in data
# hvl_placebo_cell_type_props_unfiltered <- read.table(paste0(scRNA_hvl_dir, "Placebo/RNA_cell_type_proportion_time_point_HVL_PLACEBO_unfiltered.csv"), sep = ",", header = TRUE)
hvl_placebo_cell_type_props <- read.table(paste0(scRNA_hvl_dir, "Placebo/RNA_cell_type_proportion_time_point_HVL_PLACEBO.csv"), sep = ",", header = TRUE)

# hvl_vaccinated_cell_type_props_unfiltered <- read.table(paste0(scRNA_hvl_dir, "Vaccinated/RNA_cell_type_proportion_time_point_HVL_VACCINATED_unfiltered.csv"), sep = ",", header = TRUE)
hvl_vaccinated_cell_type_props <- read.table(paste0(scRNA_hvl_dir, "Vaccinated/RNA_cell_type_proportion_time_point_HVL_VACCINATED.csv"), sep = ",", header = TRUE)

# lvl_placebo_cell_type_props_unfiltered <- read.table(paste0(scRNA_lvl_dir, "Placebo/RNA_cell_type_proportion_time_point_LVL_PLACEBO_unfiltered.csv"), sep = ",", header = TRUE)
lvl_placebo_cell_type_props <- read.table(paste0(scRNA_lvl_dir, "Placebo/RNA_cell_type_proportion_time_point_LVL_PLACEBO.csv"), sep = ",", header = TRUE)

# Find adjusted p-values
# hvl_placebo_cell_type_props_unfiltered_results <- evaluate_scRNA_cell_type_proportion_changes(hvl_placebo_cell_type_props_unfiltered)
hvl_placebo_cell_type_props_results <- evaluate_scRNA_cell_type_proportion_changes(hvl_placebo_cell_type_props)
hvl_placebo_cell_type_props_results[[1]]
hvl_placebo_cell_type_props_results[[2]]

# hvl_vaccinated_cell_type_props_unfiltered_results <- evaluate_scRNA_cell_type_proportion_changes(hvl_vaccinated_cell_type_props_unfiltered)
hvl_vaccinated_cell_type_props_results <- evaluate_scRNA_cell_type_proportion_changes(hvl_vaccinated_cell_type_props)
hvl_vaccinated_cell_type_props_results[[1]]
hvl_vaccinated_cell_type_props_results[[2]]

# lvl_placebo_cell_type_props_unfiltered_results <- evaluate_scRNA_cell_type_proportion_changes(lvl_placebo_cell_type_props_unfiltered)
lvl_placebo_cell_type_props_results <- evaluate_scRNA_cell_type_proportion_changes(lvl_placebo_cell_type_props)
lvl_placebo_cell_type_props_results[[1]]
lvl_placebo_cell_type_props_results[[2]]