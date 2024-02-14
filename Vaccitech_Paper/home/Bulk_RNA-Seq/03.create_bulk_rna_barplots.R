# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

day_vector <- c("Day 2", "Day 2", "Day 5", "Day 5", "Day 8", "Day 8", "Day 28", "Day 28")
fc_dir_vector <- rep(c("Upregulated", "Downregulated"), 4)
count_vector <- c(nrow(raw_high_placebo_period_2_D2_vs_D_minus_1_results[raw_high_placebo_period_2_D2_vs_D_minus_1_results$log2FoldChange > 0,]),
                  nrow(raw_high_placebo_period_2_D2_vs_D_minus_1_results[raw_high_placebo_period_2_D2_vs_D_minus_1_results$log2FoldChange < 0,]),
                  nrow(raw_high_placebo_period_2_D5_vs_D_minus_1_results[raw_high_placebo_period_2_D5_vs_D_minus_1_results$log2FoldChange > 0,]),
                  nrow(raw_high_placebo_period_2_D5_vs_D_minus_1_results[raw_high_placebo_period_2_D5_vs_D_minus_1_results$log2FoldChange < 0,]),
                  nrow(raw_high_placebo_period_2_D8_vs_D_minus_1_results[raw_high_placebo_period_2_D8_vs_D_minus_1_results$log2FoldChange > 0,]),
                  nrow(raw_high_placebo_period_2_D8_vs_D_minus_1_results[raw_high_placebo_period_2_D8_vs_D_minus_1_results$log2FoldChange < 0,]),
                  nrow(raw_high_placebo_period_2_D28_vs_D_minus_1_results[raw_high_placebo_period_2_D28_vs_D_minus_1_results$log2FoldChange > 0,]),
                  nrow(raw_high_placebo_period_2_D28_vs_D_minus_1_results[raw_high_placebo_period_2_D28_vs_D_minus_1_results$log2FoldChange < 0,]))

bulk_rna_deg_plot_df <- data.frame(day = day_vector, fc_direction = fc_dir_vector, count = count_vector)

bulk_rna_deg_plot_df$day <- factor(bulk_rna_deg_plot_df$day, levels = c("Day 2", "Day 5", "Day 8", "Day 28"))
bulk_rna_deg_plot_df$fc_direction <- factor(bulk_rna_deg_plot_df$fc_direction, levels = c("Upregulated", "Downregulated"))

ggplot(bulk_rna_deg_plot_df, aes(x = factor(day), y = count, fill = fc_direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Upregulated" = "maroon", "Downregulated" = "steelblue")) +
  labs(title = "RNA-Seq DEGs (Post-Exposure vs. Pre-Exposure to H3N2)",
       x = "Post-Exposure Day",
       y = "Count") +
  theme_minimal(base_size = 18) + guides(fill=guide_legend(title="Fold Change Direction")) + 
  theme(plot.title = element_text(hjust = 0.5))
