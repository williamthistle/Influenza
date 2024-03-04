# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))



bulk_cell_types <- c("B.cells.naive", "T.cells.CD8", "T.cells.CD4.naive", "T.cells.CD4.memory.resting", "T.cells.regulatory..Tregs.", "NK.cells.resting",
                     "Monocytes", "Mast.cells.resting", "Neutrophils")





### EVALUATE CELL TYPE PROPORTIONS FOR 2 D2 vs 2 D MINUS 1
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D2" | high_placebo_metadata$time_point == "2_D_minus_1",]
# Remove subjects that only have one time point (not both)
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata <- high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata[high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata$subject_id  
                                                                                                                                           %in% names(table(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata$subject_id)[table(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata$subject_id) == 2]),]
# Find unadjusted p-values for each cell type of interest
# We use coin::wilcox_test because it handles ties (which happen when there are 0s in the cell type proportions)
# Unfortunately, we can't use loops because of the way the function call works (I think)
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c()

high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(B.cells.naive ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD8 ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD4.naive ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD4.memory.resting ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.regulatory..Tregs. ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(NK.cells.resting ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Monocytes ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Mast.cells.resting ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Neutrophils ~ time_point, data = high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_metadata)))

# Nothing significant
high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values_adjusted <- p.adjust(high_placebo_period_2_D2_vs_D_minus_1_cell_type_proportion_p_values, method = "BH")

### EVALUATE CELL TYPE PROPORTIONS FOR 2 D5 vs 2 D MINUS 1
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D5" | high_placebo_metadata$time_point == "2_D_minus_1",]
# Remove subjects that only have one time point (not both)
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata <- high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata[high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata$subject_id  
                                                                                                                                           %in% names(table(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata$subject_id)[table(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata$subject_id) == 2]),]
# Find unadjusted p-values for each cell type of interest
# We use coin::wilcox_test because it handles ties (which happen when there are 0s in the cell type proportions)
# Unfortunately, we can't use loops because of the way the function call works (I think)
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c()

high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(B.cells.naive ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD8 ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD4.naive ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD4.memory.resting ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.regulatory..Tregs. ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(NK.cells.resting ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Monocytes ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Mast.cells.resting ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Neutrophils ~ time_point, data = high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_metadata)))

# T.cells.regulatory..Tregs., NK.cells.resting, and Monocytes are all significant
high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values_adjusted <- p.adjust(high_placebo_period_2_D5_vs_D_minus_1_cell_type_proportion_p_values, method = "BH")

### EVALUATE CELL TYPE PROPORTIONS FOR 2 D8 vs 2 D MINUS 1
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D8" | high_placebo_metadata$time_point == "2_D_minus_1",]
# Remove subjects that only have one time point (not both)
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata <- high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata[high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata$subject_id  
                                                                                                                                           %in% names(table(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata$subject_id)[table(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata$subject_id) == 2]),]
# Find unadjusted p-values for each cell type of interest
# We use coin::wilcox_test because it handles ties (which happen when there are 0s in the cell type proportions)
# Unfortunately, we can't use loops because of the way the function call works (I think)
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c()

high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(B.cells.naive ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD8 ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD4.naive ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD4.memory.resting ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.regulatory..Tregs. ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(NK.cells.resting ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Monocytes ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Mast.cells.resting ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Neutrophils ~ time_point, data = high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_metadata)))

# T.cells.CD4.memory.resting, NK.cells.resting, and Neutrophils are all significant
high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values_adjusted <- p.adjust(high_placebo_period_2_D8_vs_D_minus_1_cell_type_proportion_p_values, method = "BH")

### EVALUATE CELL TYPE PROPORTIONS FOR 2 D28 vs 2 D MINUS 1
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D28" | high_placebo_metadata$time_point == "2_D_minus_1",]
# Remove subjects that only have one time point (not both)
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata <- high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata[high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata$subject_id  
                                                                                                                                           %in% names(table(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata$subject_id)[table(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata$subject_id) == 2]),]
# Find unadjusted p-values for each cell type of interest
# We use coin::wilcox_test because it handles ties (which happen when there are 0s in the cell type proportions)
# Unfortunately, we can't use loops because of the way the function call works (I think)
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c()

high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(B.cells.naive ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD8 ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD4.naive ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.CD4.memory.resting ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(T.cells.regulatory..Tregs. ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(NK.cells.resting ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Monocytes ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Mast.cells.resting ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values <- c(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values,
                                                                         coin::pvalue(coin::wilcox_test(Neutrophils ~ time_point, data = high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_metadata)))

# Nothing significant
high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values_adjusted <- p.adjust(high_placebo_period_2_D28_vs_D_minus_1_cell_type_proportion_p_values, method = "BH")



