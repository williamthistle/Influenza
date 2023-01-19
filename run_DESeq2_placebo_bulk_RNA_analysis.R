library(DESeq2)
library(data.table)
library(pheatmap)

set.seed(1234)

##### SETUP #####
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
setup_bulk_analysis()

#### PERIOD 1 ####
# Look at Period 1 DEGs - shouldn't find many (or any) since vaccination didn't occur
# We don't need to separate high and low viral load people here since they're not infected yet?
# 1 D2 vs 1 D minus 1 - 0/0/0 DEGs found
placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
                                                                    "1_D2", "1_D_minus_1", data_dir)
# 1 D8 vs 1 D minus 1 - 1/0/0 DEGs found
placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
                                                                    "1_D8", "1_D_minus_1", data_dir)
# 1 D28 vs 1 D minus 1 - 2/1/1 DEGs found
placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
                                                                    "1_D28", "1_D_minus_1", data_dir)
#### PERIOD 2 ####
# Look at Period 2 DEGs
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 460/2/2 DEGs found - seems weird that 460 are found, even with a low logFC threshold
placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
                                                                     "2_D_minus_1", "2_D_minus_2", data_dir)
# 2 D2 vs 2 D minus 1 - 101/1/1 DEGs
placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
                                                                           "2_D2", "2_D_minus_1", data_dir)
# 2 D5 vs 2 D minus 1 - 2409/98/24 DEGs
placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
                                                                    "2_D5", "2_D_minus_1", data_dir)
# 2 D8 vs 2 D minus 1 - 2712/176/27 DEGs
placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
                                                                    "2_D8", "2_D_minus_1", data_dir)
# 2 D28 vs 2 D minus 1 - 2179/4/1 DEGs
placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
                                                                    "2_D28", "2_D_minus_1", data_dir)
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

#### HIGH VIRAL LOAD ####
# Grab subject IDs with full time series of data
# Note that if a subject is full period 2 then they have full period 1 as well, so that's why we can use
# placebo_period_2_LRT_metadata for this task
full_time_series_placebo <- unique(placebo_period_2_LRT_metadata$subject_id)
# Reorder subject IDs according to viral load (high to low)
full_time_series_placebo <- full_time_series_placebo[order(match(full_time_series_placebo,viral_load_primary$SUBJID))]
# Top 13 will be high viral load and bottom 10 will be low viral load 
high_viral_load_subjects <- full_time_series_placebo[1:13]
low_viral_load_subjects <- tail(full_time_series_placebo, n = 10)
# Grab high viral load placebo counts and metadata
high_placebo_aliquots <- rownames(placebo_metadata[placebo_metadata$subject_id %in% high_viral_load_subjects,])
high_placebo_counts <- placebo_counts[,high_placebo_aliquots]
high_placebo_metadata <- placebo_metadata[high_placebo_aliquots,]
# Grab low viral load placebo counts and metadata
low_placebo_aliquots <- rownames(placebo_metadata[placebo_metadata$subject_id %in% low_viral_load_subjects,])
low_placebo_counts <- placebo_counts[,low_placebo_aliquots]
low_placebo_metadata <- placebo_metadata[low_placebo_aliquots,]
#### PERIOD 2 HIGH VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 609/1/0/0 DEGs found - seems weird that 609 are found, even with a low logFC threshold
# Could there be something that happened shortly before infection that made subjects more susceptible to getting sick?
# Could it be a single patient (or two) messing everything up?
high_placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis("placebo", high_placebo_counts, high_placebo_metadata,
                                                                           "2_D_minus_1", "2_D_minus_2", data_dir, "high")
# 2 D2 vs 2 D minus 1 - 299/0/0/0 DEGs
high_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", high_placebo_counts, high_placebo_metadata,
                                                                    "2_D2", "2_D_minus_1", data_dir, "high")
# 2 D5 vs 2 D minus 1 - 3432/600/228/51 DEGs
high_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", high_placebo_counts, high_placebo_metadata,
                                                                    "2_D5", "2_D_minus_1", data_dir, "high")
# 2 D8 vs 2 D minus 1 - 1940/269/59/8 DEGs
high_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", high_placebo_counts, high_placebo_metadata,
                                                                    "2_D8", "2_D_minus_1", data_dir, "high")
# 2 D28 vs 2 D minus 1 - 34/2/0/0 DEGs
high_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", high_placebo_counts, high_placebo_metadata,
                                                                     "2_D28", "2_D_minus_1", data_dir, "high")


#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 0/0/0/0 DEGs found
low_placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis("placebo", low_placebo_counts, low_placebo_metadata,
                                                                                "2_D_minus_1", "2_D_minus_2", data_dir, "low")
# 2 D2 vs 2 D minus 1 - 1/1/1/1 DEGs
low_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "2_D2", "2_D_minus_1", data_dir, "low")
# 2 D5 vs 2 D minus 1 - 0/0/0/0 DEGs
low_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "2_D5", "2_D_minus_1", data_dir, "low")
# 2 D8 vs 2 D minus 1 - 0/0/0/0 DEGs
low_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "2_D8", "2_D_minus_1", data_dir, "low")
# 2 D28 vs 2 D minus 1 - 51/0/0/0 DEGs
low_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", low_placebo_counts, low_placebo_metadata,
                                                                          "2_D28", "2_D_minus_1", data_dir, "low")

save.image(paste0(data_dir, "placebo_bulk_RNA_obj.RData"))
