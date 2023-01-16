library(DESeq2)
library(data.table)
library(pheatmap)

set.seed(1234)

##### SETUP #####
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
setup_bulk_analysis()

run_deseq_bulk_analysis=function(counts, metadata, test_time, baseline_time, output_dir) {
  # Select the two relevant time points from our metadata
  metadata_subset <- metadata[metadata$time_point == test_time | metadata$time_point == baseline_time,]
  # Remove subjects that only have one time point (not both)
  metadata_subset <- metadata_subset[metadata_subset$subject_id  %in% names(table(metadata_subset$subject_id)[table(metadata_subset$subject_id) == 2]),]
  # Select subset of counts associated with subjects
  counts_subset <- counts[rownames(metadata_subset)]
  # Run DESeq2
  current_analysis <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + time_point)
  current_analysis <- DESeq(current_analysis)
  # Grab results with alpha = 0.05 and lfcThreshold = 0.1
  current_analysis_results <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 0.1)
  current_analysis_results <- current_analysis_results[order(current_analysis_results$padj),]
  current_analysis_results <- subset(current_analysis_results, padj < 0.05)
  write.table(rownames(current_analysis_results), paste0(output_dir, test_time, "_vs_", baseline_time, "_0.1.txt", quote = FALSE, row.names = FALSE, col.names = FALSE))
  # Grab results with alpha = 0.05 and lfcThreshold = 0.585 (1.5 fold increase)
  current_analysis_results_1.5 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 0.585)
  current_analysis_results_1.5 <- current_analysis_results_1.5[order(current_analysis_results_1.5$padj),]
  current_analysis_results_1.5 <- subset(current_analysis_results_1.5, padj < 0.05)
  write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_time, "_vs_", baseline_time, "_0.585.txt", quote = FALSE, row.names = FALSE, col.names = FALSE))
  # Grab results with alpha = 0.05 and lfcThreshold = 1
  current_analysis_results_2 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 1)
  current_analysis_results_2 <- current_analysis_results_2[order(current_analysis_results_2$padj),]
  current_analysis_results_2 <- subset(current_analysis_results_2, padj < 0.05)
  write.table(rownames(current_analysis_results_2), paste0(output_dir, test_time, "_vs_", baseline_time, "_1.txt", quote = FALSE, row.names = FALSE, col.names = FALSE))
  return(list(current_analysis_results, current_analysis_results_1.5, current_analysis_results_2))
}

#### PERIOD 1 ####
# Look at Period 1 DEGs - shouldn't find many (or any) since vaccination didn't occur
# 1 D2 vs 1 D minus 1 - 0/0/0 DEGs found
placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis(placebo_counts, placebo_metadata,
                                                                    "1_D2", "1_D_minus_1")
# 1 D8 vs 1 D minus 1 - 1/0/0 DEGs found
placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis(placebo_counts, placebo_metadata,
                                                                    "1_D8", "1_D_minus_1")
# 1 D28 vs 1 D minus 1 - 2/1/1 DEGs found
placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis(placebo_counts, placebo_metadata,
                                                                    "1_D28", "1_D_minus_1")
#### PERIOD 2 ####
# Look at Period 2 DEGs
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 460/2/2 DEGs found
placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis(placebo_counts, placebo_metadata,
                                                                     "2_D_minus_1", "2_D_minus_2")
# 2 D2 vs 2 D minus 1 - 101/1/1 DEGs
placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis(placebo_counts, placebo_metadata,
                                                                           "2_D2", "2_D_minus_1")
# 2 D5 vs 2 D minus 1 - 2409/98/24 DEGs
placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis(placebo_counts, placebo_metadata,
                                                                    "2_D5", "2_D_minus_1")
# 2 D8 vs 2 D minus 1 - 2712/176/27 DEGs
placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis(placebo_counts, placebo_metadata,
                                                                    "2_D8", "2_D_minus_1")
# 2 D28 vs 2 D minus 1 - 2179/4/1 DEGs
placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis(placebo_counts, placebo_metadata,
                                                                    "2_D28", "2_D_minus_1")
# LRT test for Period 2
# 5757 DEGs found (note that LRT doesn't have a fold change threshold). 
placebo_period_2_LRT_metadata <- placebo_metadata[placebo_metadata$time_point == "2_D28" | placebo_metadata$time_point == "2_D8" | 
                                                    placebo_metadata$time_point == "2_D5" | placebo_metadata$time_point == "2_D2" |
                                                    placebo_metadata$time_point == "2_D_minus_1",]
placebo_period_2_LRT_metadata <- placebo_period_2_LRT_metadata[placebo_period_2_LRT_metadata$subject_id 
                                                                                       %in% names(table(placebo_period_2_LRT_metadata$subject_id)
                                                                                                  [table(placebo_period_2_LRT_metadata$subject_id) == 5]),]
placebo_period_2_LRT_counts <- placebo_counts[rownames(placebo_period_2_LRT_metadata)]
placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = placebo_period_2_LRT_counts,
                                                                     colData = placebo_period_2_LRT_metadata,
                                                                     design = ~ subject_id + time_point)
placebo_period_2_LRT_analysis <- DESeq(placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id)
placebo_period_2_LRT_analysis_results <- results(placebo_period_2_LRT_analysis, alpha = 0.05)
placebo_period_2_LRT_analysis_results <- placebo_period_2_LRT_analysis_results[order(placebo_period_2_LRT_analysis_results$padj),]
placebo_period_2_LRT_analysis_results <- subset(placebo_period_2_LRT_analysis_results, padj < 0.05)
# Plot heatmap with top 20 genes
period_2_without_2_D_minus_2_time_point_analysis_betas <- coef(placebo_period_2_LRT_analysis)
period_2_without_2_D_minus_2_time_point_analysis_results_for_plotting <- results(placebo_period_2_LRT_analysis, alpha = 0.05)
period_2_without_2_D_minus_2_time_point_analysis_20_topGenes <- head(order(period_2_without_2_D_minus_2_time_point_analysis_results_for_plotting$padj),20)
period_2_without_2_D_minus_2_time_point_analysis_20_mat <- period_2_without_2_D_minus_2_time_point_analysis_betas[period_2_without_2_D_minus_2_time_point_analysis_20_topGenes, -c(1:23)]
period_2_without_2_D_minus_2_time_point_analysis_thr <- 1.5
period_2_without_2_D_minus_2_time_point_analysis_20_mat[period_2_without_2_D_minus_2_time_point_analysis_20_mat < -period_2_without_2_D_minus_2_time_point_analysis_thr] <- -period_2_without_2_D_minus_2_time_point_analysis_thr
period_2_without_2_D_minus_2_time_point_analysis_20_mat[period_2_without_2_D_minus_2_time_point_analysis_20_mat > period_2_without_2_D_minus_2_time_point_analysis_thr] <- period_2_without_2_D_minus_2_time_point_analysis_thr
colnames(period_2_without_2_D_minus_2_time_point_analysis_20_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(period_2_without_2_D_minus_2_time_point_analysis_20_mat, breaks=seq(from=-period_2_without_2_D_minus_2_time_point_analysis_thr, to=period_2_without_2_D_minus_2_time_point_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(data_dir, "period_2_without_2_D_minus_2_time_point_analysis_top_20_genes.png"))
# Plot heatmap with top 200 genes (and clean up heatmap plot by removing row names and clustering)
period_2_without_2_D_minus_2_time_point_analysis_200_topGenes <- head(order(period_2_without_2_D_minus_2_time_point_analysis_results_for_plotting$padj),200)
period_2_without_2_D_minus_2_time_point_analysis_200_mat <- period_2_without_2_D_minus_2_time_point_analysis_betas[period_2_without_2_D_minus_2_time_point_analysis_200_topGenes, -c(1:23)]
colnames(period_2_without_2_D_minus_2_time_point_analysis_200_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
period_2_without_2_D_minus_2_time_point_analysis_200_mat[period_2_without_2_D_minus_2_time_point_analysis_200_mat < -period_2_without_2_D_minus_2_time_point_analysis_thr] <- -period_2_without_2_D_minus_2_time_point_analysis_thr
period_2_without_2_D_minus_2_time_point_analysis_200_mat[period_2_without_2_D_minus_2_time_point_analysis_200_mat > period_2_without_2_D_minus_2_time_point_analysis_thr] <- period_2_without_2_D_minus_2_time_point_analysis_thr
pheatmap(period_2_without_2_D_minus_2_time_point_analysis_200_mat, breaks=seq(from=-period_2_without_2_D_minus_2_time_point_analysis_thr, to=period_2_without_2_D_minus_2_time_point_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, show_rownames = FALSE, cluster_row = FALSE, filename = paste0(data_dir, "period_2_without_2_D_minus_2_time_point_analysis_top_200_genes.png"))
save.image(paste0(data_dir, "placebo_bulk_RNA_obj.RData"))
