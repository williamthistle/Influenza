# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# HVL
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
bulk_degs_lineplot_table$Fold_Change_Direction <- factor(bulk_degs_lineplot_table$Fold_Change_Direction, levels = c("Upregulated", "Downregulated"))

# bulk_degs_lineplot_table$Count[bulk_degs_lineplot_table$Count == 0] <- bulk_degs_lineplot_table$Count[bulk_degs_lineplot_table$Count == 0] + 1

ggplot(data = bulk_degs_lineplot_table, aes(x = Day, y = Count, color = Fold_Change_Direction, shape = Fold_Change_Threshold)) +
  geom_line(aes(group = interaction(Fold_Change_Threshold, Fold_Change_Direction))) +  
  geom_point(size = 3) +  # Add points along the lines
  labs(
    title = "RNA-Seq DEGs by Day Post-Exposure",
    x = "Day",
    y = "DEG Count",
    color = "Fold Change Direction",
    shape = "Fold Change Threshold"
  ) +   theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(
    breaks = c(0, 500, 1000, 1500),
    limits = c(0, 1500)
  ) + guides(color = guide_legend(order = 1), 
          shape = guide_legend(order = 2))

# LVL
bulk_degs_lineplot_table <- data.frame(Day = character(), Fold_Change_Threshold = character(), Fold_Change_Direction = character(), Count = numeric())

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 2", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Upregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results[[2]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 2", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Downregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results[[2]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 2", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Upregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results[[5]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 2", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Downregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results[[5]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 5", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Upregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[2]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 5", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Downregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[2]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 5", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Upregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[5]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 5", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Downregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[5]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 8", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Upregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results[[2]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 8", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Downregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results[[2]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 8", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Upregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results[[5]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 8", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Downregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results[[5]]$log2FoldChange < 0)))
bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 28", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Upregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results[[2]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 28", Fold_Change_Threshold = "0.1", Fold_Change_Direction = "Downregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results[[2]]$log2FoldChange < 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 28", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Upregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results[[5]]$log2FoldChange > 0)))

bulk_degs_lineplot_table <- rbind(bulk_degs_lineplot_table,
                                  data.frame(Day = "Day 28", Fold_Change_Threshold = "0.585", Fold_Change_Direction = "Downregulated",
                                             Count = sum(lvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results[[5]]$log2FoldChange < 0)))

bulk_degs_lineplot_table$Day <- factor(bulk_degs_lineplot_table$Day, levels = c("Day 2", "Day 5", "Day 8", "Day 28"))
bulk_degs_lineplot_table$Fold_Change_Threshold <- factor(bulk_degs_lineplot_table$Fold_Change_Threshold, levels = c("0.1", "0.585"))
bulk_degs_lineplot_table$Fold_Change_Direction <- factor(bulk_degs_lineplot_table$Fold_Change_Direction, levels = c("Upregulated", "Downregulated"))

# bulk_degs_lineplot_table$Count[bulk_degs_lineplot_table$Count == 0] <- bulk_degs_lineplot_table$Count[bulk_degs_lineplot_table$Count == 0] + 1

ggplot(data = bulk_degs_lineplot_table, aes(x = Day, y = Count, color = Fold_Change_Direction, shape = Fold_Change_Threshold)) +
  geom_line(aes(group = interaction(Fold_Change_Threshold, Fold_Change_Direction))) +  
  geom_point(size = 3) +  # Add points along the lines
  labs(
    title = "RNA-Seq DEGs by Day Post-Exposure",
    x = "Day",
    y = "DEG Count",
    color = "Fold Change Direction",
    shape = "Fold Change Threshold"
  ) +   theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5)) + guides(color = guide_legend(order = 1), 
             shape = guide_legend(order = 2))
