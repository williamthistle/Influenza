# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

find_bulk_correlation <- function(first_deg_subset, second_deg_subset) {
  correlation_return_vals <- list()
  primary_bulk_degs <- rownames(first_deg_subset)
  print(paste0("Number of bulk DEGs in first set of data is: ", length(primary_bulk_degs)))
  compare_first_df <- first_deg_subset
  compare_second_df <- second_deg_subset[rownames(second_deg_subset) %in% primary_bulk_degs,]
  compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$log2FoldChange,
                                             second_fc = compare_second_df$log2FoldChange)
  # Calculate correlation - Pearson
  print("Calculating Pearson correlation")
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  correlation_return_vals[[1]] <- correlation_val$estimate
  print("Calculating Spearman correlation")
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc, method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  correlation_return_vals[[2]] <- correlation_val$estimate
  
  # Plot correlation
  correlation_plot <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr() + xlab("Vaccinated FC") + ylab("Naive FC")
  correlation_return_vals[[3]] <- correlation_plot
  
  return(correlation_return_vals)
}


find_bulk_correlations <- function(first_degs, second_degs) {
  # Unfiltered
  first_degs_unfiltered <- first_degs[[8]]
  second_degs_unfiltered <- second_degs[[8]]
  
  # Filtered (0 logFC threshold)
  first_degs_0_filtered <- first_degs[[1]]
  second_degs_0_filtered <- second_degs[[1]]
  
  # Filtered (0.1 logFC threshold)
  first_degs_0.1_filtered <- first_degs[[2]]
  second_degs_0.1_filtered <- second_degs[[2]]
  
  # Filtered (0.585 logFC threshold)
  first_degs_0.585_filtered <- first_degs[[5]]
  second_degs_0.585_filtered <- second_degs[[5]]
  
  if(nrow(first_degs_0_filtered) > 10) {
    first_degs_0_filtered_vs_second_degs_unfiltered_corr <- find_bulk_correlation(first_degs_0_filtered, second_degs_unfiltered)
  } else {
    first_degs_0_filtered_vs_second_degs_unfiltered_corr <- NULL
  }
  if(nrow(first_degs_0.1_filtered) > 10) {
    first_degs_0.1_filtered_vs_second_degs_unfiltered_corr <- find_bulk_correlation(first_degs_0.1_filtered, second_degs_unfiltered)
  } else {
    first_degs_0.1_filtered_vs_second_degs_unfiltered_corr <- NULL
  }
  if(nrow(first_degs_0.585_filtered) > 10) {
    first_degs_0.585_filtered_vs_second_degs_unfiltered_corr <- find_bulk_correlation(first_degs_0.585_filtered, second_degs_unfiltered)
  } else {
    first_degs_0.585_filtered_vs_second_degs_unfiltered_corr <- NULL
  }
  return(list(first_degs_0_filtered_vs_second_degs_unfiltered_corr, first_degs_0.1_filtered_vs_second_degs_unfiltered_corr, first_degs_0.585_filtered_vs_second_degs_unfiltered_corr))
}

# Period 1

# Day 2 HVL vs LVL
# Doesn't work - too few DEGs for either side
period_1_day_2_hvl_vs_lvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_1_D2_vs_D_minus_1_results, lvl_full_time_series_placebo_period_1_D2_vs_D_minus_1_results)

# Day 8 HVL vs LVL
# Doesn't work - too few DEGs for either side
period_1_day_8_hvl_vs_lvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_1_D8_vs_D_minus_1_results, lvl_full_time_series_placebo_period_1_D8_vs_D_minus_1_results)

# Day 28 HVL vs LVL
period_1_day_28_hvl_vs_lvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_1_D28_vs_D_minus_1_results, lvl_full_time_series_placebo_period_1_D28_vs_D_minus_1_results)

# Period 2

# Day 2 HVL vs LVL
# Doesn't work - too few DEGs for either side
period_2_day_2_hvl_vs_lvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results, lvl_full_time_series_placebo_period_2_D2_vs_D_minus_1_results)

# Day 5 HVL vs LVL
period_2_day_5_hvl_vs_lvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results, lvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results)

# Day 8 HVL vs LVL
period_2_day_8_hvl_vs_lvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results, lvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results)

# Day 8 HVL vs LVL (full)
period_2_day_8_hvl_vs_lvl_full_bulk_correlations <- find_bulk_correlations(hvl_placebo_period_2_D8_vs_D_minus_1_results, lvl_placebo_period_2_D8_vs_D_minus_1_results)

# Day 28 HVL vs LVL
period_2_day_28_hvl_vs_lvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results, lvl_full_time_series_placebo_period_2_D28_vs_D_minus_1_results)

# Day 28 HVL vs LVL (full)
period_2_day_28_hvl_vs_lvl_full_bulk_correlations <- find_bulk_correlations(hvl_placebo_period_2_D28_vs_D_minus_1_results, lvl_placebo_period_2_D28_vs_D_minus_1_results)

# Day 28 LVL vs HVL (full)
period_2_day_28_lvl_vs_hvl_full_bulk_correlations <- find_bulk_correlations(lvl_placebo_period_2_D28_vs_D_minus_1_results, hvl_placebo_period_2_D28_vs_D_minus_1_results)

# Day 5 HVL vs Day 28 HVL
period_2_day_5_hvl_vs_day_28_hvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results, hvl_placebo_period_2_D28_vs_D_minus_1_results)

# Day 8 HVL vs Day 28 HVL
period_2_day_8_hvl_vs_day_28_hvl_bulk_correlations <- find_bulk_correlations(hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results, hvl_placebo_period_2_D28_vs_D_minus_1_results)

### Placebo vs Vaccinated

# Period 2

# Day 8 Placebo HVL vs Day 8 Vaccinated HVL
period_2_day_8_placebo_hvl_vs_vaccinated_hvl_bulk_correlations <- find_bulk_correlations(hvl_placebo_period_2_D8_vs_D_minus_1_results, hvl_vaccinated_period_2_D8_vs_D_minus_1_results)

# Day 28 Vaccinated HVL vs Day 28 Placebo HVL
period_2_day_28_vaccinated_hvl_vs_placebo_hvl_bulk_correlations <- find_bulk_correlations(hvl_vaccinated_period_2_D28_vs_D_minus_1_results, hvl_placebo_period_2_D28_vs_D_minus_1_results)


