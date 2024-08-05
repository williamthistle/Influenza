### THIS FUNCTION:
### RUNS MAGICAL AND CREATES 1) OVERALL MAGICAL DF, AND 2) TF-CENTRIC MAGICAL DF

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))
source(paste0(base_dir, "extra_functions/my_MAGICAL_helper_functions.R"))

hvl_naive_magical_df <- read.table(paste0(MAGICAL_hvl_placebo_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"), sep = "\t", header = TRUE)
lvl_naive_magical_df <- read.table(paste0(MAGICAL_lvl_placebo_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"), sep = "\t", header = TRUE)
hvl_vaccinated_magical_df <- read.table(paste0(MAGICAL_hvl_vaccinated_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"), sep = "\t", header = TRUE)

hvl_naive_magical_df_tf <- create_tf_vs_cell_type_df(hvl_naive_magical_df)
lvl_naive_magical_df_tf <- create_tf_vs_cell_type_df(lvl_naive_magical_df)
hvl_vaccinated_magical_df_tf <- create_tf_vs_cell_type_df(hvl_vaccinated_magical_df)

hvl_naive_magical_tf_heatmap_plot <- create_tf_heatmap_plot(hvl_naive_magical_df_tf)
lvl_naive_magical_tf_heatmap_plot <- create_tf_heatmap_plot(lvl_naive_magical_df_tf)
hvl_vaccinated_magical_tf_heatmap_plot <- create_tf_heatmap_plot(hvl_vaccinated_magical_df_tf)

ggsave("C:/Users/wat2/Desktop/lvl_naive_MAGICAL_TF_plot.png", plot = lvl_naive_magical_tf_heatmap_plot, width = 10, height = 10, units = "in", dpi = 300)
ggsave("C:/Users/wat2/Desktop/hvl_vaccinated_MAGICAL_TF_plot.png", plot = hvl_vaccinated_magical_tf_heatmap_plot, width = 10, height = 10, units = "in", dpi = 300)
