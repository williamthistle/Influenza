# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

bulk_degs_lineplot_table <- data.frame(Day = character(), Fold_Change_Threshold = character(), Fold_Change_Direction = character(), Count = numeric())

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 2", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Upregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results[[2]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 2", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Downregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results[[2]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 2", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Upregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results[[5]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 2", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Downregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results[[5]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 5", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Upregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[2]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 5", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Downregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[2]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 5", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Upregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[5]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 5", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Downregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[5]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 8", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Upregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results[[2]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 8", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Downregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results[[2]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 8", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Upregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results[[5]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 8", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Downregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results[[5]]$log2FoldChange < 0)))
bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 28", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Upregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results[[2]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 28", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Downregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results[[2]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 28", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Upregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results[[5]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 28", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Downregulated",
                                             Count = sum(hvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results[[5]]$log2FoldChange < 0)))

bulk_degs_lineplot_table$Day <- factor(bulk_degs_lineplot_table$Day, levels = c("Day 2", "Day 5", "Day 8", "Day 28"))
bulk_degs_lineplot_table$Fold_Change_Threshold <- factor(bulk_degs_lineplot_table$Fold_Change_Threshold, levels = c("0.1", "0.585"))

bulk_degs_lineplot_table$Count <- bulk_degs_lineplot_table$Count + 1

ggplot(data = bulk_degs_lineplot_table, aes(x = Day, y = Count, color = paste(Fold_Change_Threshold, Fold_Change_Direction, sep = "_"))) +
  geom_line() +  # Plot lines
  geom_point() +  # Add points along the lines
  labs(
    title = "Count over Days",
    x = "Day",
    y = "Count (log2 scale)",
    color = "Threshold and Direction"
  ) +
  scale_y_continuous(
    trans = "log2"  # Use log2 scale for y-axis
  ) +
  theme_minimal()

ggsave(filename = paste0("C:/Users/willi/Desktop/", "bulk_deg_jitterplot.tiff"), plot = bulk_degs_jitter_plot, device='tiff', dpi=300)