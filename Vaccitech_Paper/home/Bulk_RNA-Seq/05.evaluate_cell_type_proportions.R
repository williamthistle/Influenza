# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# HVL (Placebo)

# Period 1 (D2 vs D minus 1) - No significant changes
high_placebo_period_1_D2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(high_placebo_metadata, bulk_cell_types, "1_D2", "1_D_minus_1", absolute_score = TRUE)

# Period 1 (D8 vs D minus 1) - No significant changes
high_placebo_period_1_D8_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(high_placebo_metadata, bulk_cell_types, "1_D8", "1_D_minus_1", absolute_score = TRUE)

# Period 1 (D28 vs D minus 1) - No significant changes
high_placebo_period_1_D28_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(high_placebo_metadata, bulk_cell_types, "1_D28", "1_D_minus_1", absolute_score = TRUE)

# Period 2 (D minus 1 vs D minus 2) - No significant changes
high_placebo_period_2_D_minus_2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(high_placebo_metadata, bulk_cell_types, "2_D_minus_1", "2_D_minus_2", absolute_score = TRUE)

# Period 2 (D2 vs D minus 1) - No significant changes
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(high_placebo_metadata, bulk_cell_types, "2_D2", "2_D_minus_1")

# Period 2 (D5 vs D minus 1) - T.cells.regulatory..Tregs., NK.cells.resting, and Monocytes are all significant
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(high_placebo_metadata, bulk_cell_types, "2_D5", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D8 vs D minus 1) - T.cells.CD4.memory.resting, NK.cells.resting, and Neutrophils are all significant
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(high_placebo_metadata, bulk_cell_types, "2_D8", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D28 vs D minus 1) - No significant changes
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(high_placebo_metadata, bulk_cell_types, "2_D28", "2_D_minus_1", absolute_score = TRUE)

# LVL (Placebo)

# Period 1 (D2 vs D minus 1) - No significant changes
low_placebo_period_1_D2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(low_placebo_metadata, bulk_cell_types, "1_D2", "1_D_minus_1", absolute_score = TRUE)

# Period 1 (D8 vs D minus 1) - No significant changes
low_placebo_period_1_D8_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(low_placebo_metadata, bulk_cell_types, "1_D8", "1_D_minus_1", absolute_score = TRUE)

# Period 1 (D28 vs D minus 1) - No significant changes
low_placebo_period_1_D28_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(low_placebo_metadata, bulk_cell_types, "1_D28", "1_D_minus_1", absolute_score = TRUE)

# Period 2 (D minus 1 vs D minus 2) - No significant changes
low_placebo_period_2_D_minus_2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(low_placebo_metadata, bulk_cell_types, "2_D_minus_1", "2_D_minus_2", absolute_score = TRUE)

# Period 2 (D2 vs D minus 1) - No significant changes
low_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(low_placebo_metadata, bulk_cell_types, "2_D2", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D5 vs D minus 1) - No significant changes
low_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(low_placebo_metadata, bulk_cell_types, "2_D5", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D8 vs D minus 1) - No significant changes
low_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(low_placebo_metadata, bulk_cell_types, "2_D8", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D28 vs D minus 1) - No significant changes
low_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(low_placebo_metadata, bulk_cell_types, "2_D28", "2_D_minus_1", absolute_score = TRUE)

# HVL (Vaccinated)

# Period 1 (D2 vs D minus 1) - CD8 T cells, Tregs, NK cells, and neutrophils are all significant
hvl_vaccinated_period_1_D2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(hvl_vaccinated_metadata, bulk_cell_types, "1_D2", "1_D_minus_1", absolute_score = TRUE)

# Period 1 (D8 vs D minus 1) - No significant changes
hvl_vaccinated_period_1_D8_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(hvl_vaccinated_metadata, bulk_cell_types, "1_D8", "1_D_minus_1", absolute_score = TRUE)

# Period 1 (D28 vs D minus 1) - No significant changes
hvl_vaccinated_period_1_D28_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(hvl_vaccinated_metadata, bulk_cell_types, "1_D28", "1_D_minus_1", absolute_score = TRUE)

# Period 2 (D minus 1 vs D minus 2) - No significant changes
hvl_vaccinated_period_2_D_minus_2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(hvl_vaccinated_metadata, bulk_cell_types, "2_D_minus_1", "2_D_minus_2", absolute_score = TRUE)

# Period 2 (D2 vs D minus 1) - No significant changes
hvl_vaccinated_period_2_D2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(hvl_vaccinated_metadata, bulk_cell_types, "2_D2", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D5 vs D minus 1) - Monocytes are significant
hvl_vaccinated_period_2_D5_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(hvl_vaccinated_metadata, bulk_cell_types, "2_D5", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D8 vs D minus 1) - NK cells, Monocytes, and Neutrophils are all significant
hvl_vaccinated_period_2_D8_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(hvl_vaccinated_metadata, bulk_cell_types, "2_D8", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D28 vs D minus 1) - CD4 memory and Tregs are significant
hvl_vaccinated_period_2_D28_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(hvl_vaccinated_metadata, bulk_cell_types, "2_D28", "2_D_minus_1", absolute_score = TRUE)

# LVL (Vaccinated)

# Period 1 (D2 vs D minus 1) - CD8 T cells, CD4 Naive T cells, and Neutrophils are all significant
lvl_vaccinated_period_1_D2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(lvl_vaccinated_metadata, bulk_cell_types, "1_D2", "1_D_minus_1", absolute_score = TRUE)

# Period 1 (D8 vs D minus 1) - CD4 Memory T cells, Tregs, NK cells, and Neutrophils are all significant
lvl_vaccinated_period_1_D8_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(lvl_vaccinated_metadata, bulk_cell_types, "1_D8", "1_D_minus_1", absolute_score = TRUE)

# Period 1 (D28 vs D minus 1) - No significant changes
lvl_vaccinated_period_1_D28_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(lvl_vaccinated_metadata, bulk_cell_types, "1_D28", "1_D_minus_1", absolute_score = TRUE)

# Period 2 (D minus 1 vs D minus 2) - No significant changes
lvl_vaccinated_period_2_D_minus_2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(lvl_vaccinated_metadata, bulk_cell_types, "2_D_minus_1", "2_D_minus_2", absolute_score = TRUE)

# Period 2 (D2 vs D minus 1) - No significant changes
lvl_vaccinated_period_2_D2_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(lvl_vaccinated_metadata, bulk_cell_types, "2_D2", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D5 vs D minus 1) - No significant changes
lvl_vaccinated_period_2_D5_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(lvl_vaccinated_metadata, bulk_cell_types, "2_D5", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D8 vs D minus 1) - No significant changes
lvl_vaccinated_period_2_D8_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(lvl_vaccinated_metadata, bulk_cell_types, "2_D8", "2_D_minus_1", absolute_score = TRUE)

# Period 2 (D28 vs D minus 1) - B Cells naive, CD4 Memory T cells, and Tregs are all significant
lvl_vaccinated_period_2_D28_vs_D_minus_1_cell_type_proportion_changes <- evaluate_bulk_cell_type_proportion_changes(lvl_vaccinated_metadata, bulk_cell_types, "2_D28", "2_D_minus_1", absolute_score = TRUE)
