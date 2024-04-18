# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

bulk_correlations <- list()
bulk_correlation_plots <- list()

# Comparing 2_D8
# Unfiltered
hvl_placebo_period_2_D8_vs_D_minus_1_results_unfiltered <- hvl_placebo_period_2_D8_vs_D_minus_1_results[[8]]
hvl_vaccinated_period_2_D8_vs_D_minus_1_results_unfiltered <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results[[8]]

# Filtered (0 logFC threshold)
hvl_placebo_period_2_D8_vs_D_minus_1_results_0_filtered <- hvl_placebo_period_2_D8_vs_D_minus_1_results[[1]]
hvl_vaccinated_period_2_D8_vs_D_minus_1_results_0_filtered <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results[[1]]

# Filtered (0.1 logFC threshold)
hvl_placebo_period_2_D8_vs_D_minus_1_results_0.1_filtered <- hvl_placebo_period_2_D8_vs_D_minus_1_results[[2]]
hvl_vaccinated_period_2_D8_vs_D_minus_1_results_0.1_filtered <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results[[2]]

# Filtered (0.585 logFC threshold)
hvl_placebo_period_2_D8_vs_D_minus_1_results_0.585_filtered <- hvl_placebo_period_2_D8_vs_D_minus_1_results[[5]]
hvl_vaccinated_period_2_D8_vs_D_minus_1_results_0.585_filtered <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results[[5]]


# Check correlation between filtered (0 logFC) in vaccinated and unfiltered placebo
primary_bulk_degs <- rownames(hvl_vaccinated_period_2_D8_vs_D_minus_1_results_0_filtered)
print(paste0("Number of bulk DEGs in vaccinated data is: ", length(primary_bulk_degs)))
compare_first_df <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results_0_filtered
compare_second_df <- hvl_placebo_period_2_D8_vs_D_minus_1_results_unfiltered[rownames(hvl_placebo_period_2_D8_vs_D_minus_1_results_unfiltered) %in% primary_bulk_degs,]
compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$log2FoldChange,
                                             second_fc = compare_second_df$log2FoldChange)
# Calculate correlation
correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
print(correlation_val$estimate)
print(correlation_val$p.value)
  
# Plot correlation
ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
  geom_point(size = 2) +
  sm_statCorr(corr_method = "spearman") + xlab("Vaccinated FC") + ylab("Naive FC") + labs(title = "D8 (0 logFC threshold)")
  
# Check correlation between filtered (0 logFC) in placebo and unfiltered vaccinated
primary_bulk_degs <- rownames(hvl_placebo_period_2_D8_vs_D_minus_1_results_0.1_filtered)
print(paste0("Number of bulk DEGs in primary data is: ", length(primary_bulk_degs)))
compare_first_df <- hvl_placebo_period_2_D8_vs_D_minus_1_results_0.1_filtered
compare_second_df <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results_unfiltered[rownames(hvl_vaccinated_period_2_D8_vs_D_minus_1_results_unfiltered) %in% primary_bulk_degs,]
compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]

comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$log2FoldChange,
                                           second_fc = compare_second_df$log2FoldChange)
# Calculate correlation
correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
print(correlation_val$estimate)
print(correlation_val$p.value)

# Plot correlation
ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Naive FC") + ylab("Vaccinated FC") + labs(title = "D8 (0.1 logFC threshold)")

# Check correlation between filtered (0.585 logFC) in placebo and unfiltered vaccinated
primary_bulk_degs <- rownames(hvl_placebo_period_2_D8_vs_D_minus_1_results_0.585_filtered)
print(paste0("Number of bulk DEGs in primary data is: ", length(primary_bulk_degs)))
compare_first_df <- hvl_placebo_period_2_D8_vs_D_minus_1_results_0.585_filtered
compare_second_df <- hvl_vaccinated_period_2_D8_vs_D_minus_1_results_unfiltered[rownames(hvl_vaccinated_period_2_D8_vs_D_minus_1_results_unfiltered) %in% primary_bulk_degs,]
compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]

comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$log2FoldChange,
                                           second_fc = compare_second_df$log2FoldChange)
# Calculate correlation
correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
print(correlation_val$estimate)
print(correlation_val$p.value)

# Plot correlation
ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
  geom_point(size = 2) +
  sm_statCorr() + xlab("Naive FC") + ylab("Vaccinated FC") + labs(title = "D8 (0.585 logFC threshold)")

