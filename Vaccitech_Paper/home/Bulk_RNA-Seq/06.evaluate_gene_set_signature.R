# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Create metaintegrator obj for data
high_D28_bulk_metaintegrator_obj <- create_metaintegrator_obj("bulk", larger_list_high_placebo_counts, larger_list_high_placebo_metadata, "2_D28", "2_D_minus_1")


high_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", larger_list_high_placebo_counts, larger_list_high_placebo_metadata,
                                                                                      "2_D28", "2_D_minus_1", paste0(bulk_results_dir, "hvl_bulk_placebo_period_2_D28_vs_D_minus_1/"), "high")