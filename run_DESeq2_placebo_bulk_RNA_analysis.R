library(DESeq2)
library(data.table)
library(pheatmap)

set.seed(1234)

##### SETUP #####
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
setup_bulk_analysis()

# Grab the main 23 subjects (13 high viral load, 10 low viral load) that have all 10 time points
placebo_full_time_series_metadata <- placebo_metadata[placebo_metadata$subject_id 
                                                               %in% names(table(placebo_metadata$subject_id)
                                                                          [table(placebo_metadata$subject_id) == 10]),]
# Grab subject IDs for main 23 subjects
placebo_full_time_series_subjects <- unique(placebo_full_time_series_metadata$subject_id)
# Reorder subject IDs according to viral load (high to low)
full_time_series_subjects <- placebo_full_time_series_subjects[order(match(placebo_full_time_series_subjects,viral_load_primary$SUBJID))]
# Top 13 will be high viral load and bottom 10 will be low viral load 
high_viral_load_subjects <- placebo_full_time_series_subjects[1:13]
low_viral_load_subjects <- tail(placebo_full_time_series_subjects, n = 10)
# Grab high viral load placebo counts and metadata
high_placebo_aliquots <- rownames(placebo_metadata[placebo_metadata$subject_id %in% high_viral_load_subjects,])
high_placebo_counts <- placebo_counts[,high_placebo_aliquots]
high_placebo_metadata <- placebo_metadata[high_placebo_aliquots,]
# Grab low viral load placebo counts and metadata
low_placebo_aliquots <- rownames(placebo_metadata[placebo_metadata$subject_id %in% low_viral_load_subjects,])
low_placebo_counts <- placebo_counts[,low_placebo_aliquots]
low_placebo_metadata <- placebo_metadata[low_placebo_aliquots,]
# Now we are ready to run differential expression

#### PERIOD 1 HIGH VIRAL LOAD ####
# We expect ~0 DEGs because placebo was used (no actual vaccination)
# 0/0/0/0 DEGs
high_placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D2", "1_D_minus_1", data_dir, "high")
# 0/0/0/0 DEGs
high_placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D8", "1_D_minus_1", data_dir, "high")
# 1/1/1/1 DEGs
high_placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D28", "1_D_minus_1", data_dir, "high")
#### PERIOD 1 LOW VIRAL LOAD ####
# We expect ~0 DEGs because placebo was used (no actual vaccination)
# 2/0/0/0 DEGs
low_placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "1_D2", "1_D_minus_1", data_dir, "low")
# 0/0/0/0 DEGs
low_placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "1_D8", "1_D_minus_1", data_dir, "low")
# 0/0/0/0 DEGs
low_placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", low_placebo_counts, low_placebo_metadata,
                                                                          "1_D28", "1_D_minus_1", data_dir, "low")
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
# Find length of intersection between D5 (1 FC) and D8 (1 FC)
length(intersect(rownames(high_placebo_period_2_D5_vs_D_minus_1_results[[3]]), rownames(high_placebo_period_2_D8_vs_D_minus_1_results[[3]])))
# Find length of intersection between D5 (2 FC) and D8 (2 FC)
length(intersect(rownames(high_placebo_period_2_D5_vs_D_minus_1_results[[4]]), rownames(high_placebo_period_2_D8_vs_D_minus_1_results[[4]])))
# LRT test
# 6431 DEGs found (note that LRT doesn't have a fold change threshold)
high_placebo_period_2_LRT_metadata <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D28" | high_placebo_metadata$time_point == "2_D8" | 
                                                              high_placebo_metadata$time_point == "2_D5" | high_placebo_metadata$time_point == "2_D2" |
                                                              high_placebo_metadata$time_point == "2_D_minus_1",]
high_placebo_period_2_LRT_metadata <- high_placebo_period_2_LRT_metadata[high_placebo_period_2_LRT_metadata$subject_id 
                                                               %in% names(table(high_placebo_period_2_LRT_metadata$subject_id)
                                                                          [table(high_placebo_period_2_LRT_metadata$subject_id) == 5]),]
high_placebo_period_2_LRT_counts <- high_placebo_counts[rownames(high_placebo_period_2_LRT_metadata)]
high_placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = high_placebo_period_2_LRT_counts,
                                                        colData = high_placebo_period_2_LRT_metadata,
                                                        design = ~ subject_id + time_point)
high_placebo_period_2_LRT_analysis <- DESeq(high_placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id)
high_placebo_period_2_LRT_analysis_results <- results(high_placebo_period_2_LRT_analysis, alpha = 0.05)
high_placebo_period_2_LRT_analysis_results <- high_placebo_period_2_LRT_analysis_results[order(high_placebo_period_2_LRT_analysis_results$padj),]
high_placebo_period_2_LRT_analysis_results <- subset(high_placebo_period_2_LRT_analysis_results, padj < 0.05)
# Plot heatmap with top 20 genes
high_placebo_period_2_LRT_analysis_betas <- coef(high_placebo_period_2_LRT_analysis)
high_placebo_period_2_LRT_analysis_results_for_plotting <- results(high_placebo_period_2_LRT_analysis, alpha = 0.05)
high_placebo_period_2_LRT_analysis_20_topGenes <- head(order(high_placebo_period_2_LRT_analysis_results_for_plotting$padj),20)
high_placebo_period_2_LRT_analysis_20_mat <- high_placebo_period_2_LRT_analysis_betas[high_placebo_period_2_LRT_analysis_20_topGenes, -c(1:13)]
high_placebo_period_2_LRT_analysis_thr <- 6
high_placebo_period_2_LRT_analysis_20_mat[high_placebo_period_2_LRT_analysis_20_mat < -high_placebo_period_2_LRT_analysis_thr] <- -high_placebo_period_2_LRT_analysis_thr
high_placebo_period_2_LRT_analysis_20_mat[high_placebo_period_2_LRT_analysis_20_mat > high_placebo_period_2_LRT_analysis_thr] <- high_placebo_period_2_LRT_analysis_thr
colnames(high_placebo_period_2_LRT_analysis_20_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(high_placebo_period_2_LRT_analysis_20_mat, breaks=seq(from=-high_placebo_period_2_LRT_analysis_thr, to=high_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(data_dir, "high_placebo_period_2_LRT_analysis_top_20_genes.png"))
# Plot heatmap with top 200 genes (and clean up heatmap plot by removing row names and clustering)
high_placebo_period_2_LRT_analysis_200_topGenes <- head(order(high_placebo_period_2_LRT_analysis_results_for_plotting$padj),200)
high_placebo_period_2_LRT_analysis_200_mat <- high_placebo_period_2_LRT_analysis_betas[high_placebo_period_2_LRT_analysis_200_topGenes, -c(1:13)]
colnames(high_placebo_period_2_LRT_analysis_200_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
high_placebo_period_2_LRT_analysis_200_mat[high_placebo_period_2_LRT_analysis_200_mat < -high_placebo_period_2_LRT_analysis_thr] <- -high_placebo_period_2_LRT_analysis_thr
high_placebo_period_2_LRT_analysis_200_mat[high_placebo_period_2_LRT_analysis_200_mat > high_placebo_period_2_LRT_analysis_thr] <- high_placebo_period_2_LRT_analysis_thr
pheatmap(high_placebo_period_2_LRT_analysis_200_mat, breaks=seq(from=-high_placebo_period_2_LRT_analysis_thr, to=high_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, show_rownames = FALSE, cluster_row = FALSE, filename = paste0(data_dir, "high_placebo_period_2_LRT_analysis_top_200_genes.png"))
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
# LRT test
# 0 DEGs found (note that LRT doesn't have a fold change threshold)
low_placebo_period_2_LRT_metadata <- low_placebo_metadata[low_placebo_metadata$time_point == "2_D28" | low_placebo_metadata$time_point == "2_D8" | 
                                                              low_placebo_metadata$time_point == "2_D5" | low_placebo_metadata$time_point == "2_D2" |
                                                              low_placebo_metadata$time_point == "2_D_minus_1",]
low_placebo_period_2_LRT_metadata <- low_placebo_period_2_LRT_metadata[low_placebo_period_2_LRT_metadata$subject_id 
                                                                         %in% names(table(low_placebo_period_2_LRT_metadata$subject_id)
                                                                                    [table(low_placebo_period_2_LRT_metadata$subject_id) == 5]),]
low_placebo_period_2_LRT_counts <- low_placebo_counts[rownames(low_placebo_period_2_LRT_metadata)]
low_placebo_period_2_LRT_analysis <- DESeqDataSetFromMatrix(countData = low_placebo_period_2_LRT_counts,
                                                             colData = low_placebo_period_2_LRT_metadata,
                                                             design = ~ subject_id + time_point)
low_placebo_period_2_LRT_analysis <- DESeq(low_placebo_period_2_LRT_analysis, test = "LRT", reduced = ~ subject_id)
low_placebo_period_2_LRT_analysis_results <- results(low_placebo_period_2_LRT_analysis, alpha = 0.05)
low_placebo_period_2_LRT_analysis_results <- low_placebo_period_2_LRT_analysis_results[order(low_placebo_period_2_LRT_analysis_results$padj),]
low_placebo_period_2_LRT_analysis_results <- subset(low_placebo_period_2_LRT_analysis_results, padj < 0.05)

save.image(paste0(data_dir, "placebo_bulk_RNA_obj.RData"))

# Should I do high vs low viral load differential expression here? 



































#### PERIOD 1 ####
# Look at Period 1 DEGs - shouldn't find many (or any) since vaccination didn't occur
# We don't need to separate high and low viral load people here since they're not infected yet?
# 1 D2 vs 1 D minus 1 - 0/0/0 DEGs found
#placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
#                                                                    "1_D2", "1_D_minus_1", data_dir)
# 1 D8 vs 1 D minus 1 - 1/0/0 DEGs found
#placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
#                                                                    "1_D8", "1_D_minus_1", data_dir)
# 1 D28 vs 1 D minus 1 - 2/1/1 DEGs found
#placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
#                                                                    "1_D28", "1_D_minus_1", data_dir)
#### PERIOD 2 ####
# Look at Period 2 DEGs
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 460/2/2 DEGs found - seems weird that 460 are found, even with a low logFC threshold
#placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
#                                                                     "2_D_minus_1", "2_D_minus_2", data_dir)
# 2 D2 vs 2 D minus 1 - 101/1/1 DEGs
#placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
#                                                                           "2_D2", "2_D_minus_1", data_dir)
# 2 D5 vs 2 D minus 1 - 2409/98/24 DEGs
#placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
#                                                                    "2_D5", "2_D_minus_1", data_dir)
# 2 D8 vs 2 D minus 1 - 2712/176/27 DEGs
#placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
#                                                                    "2_D8", "2_D_minus_1", data_dir)
# 2 D28 vs 2 D minus 1 - 2179/4/1 DEGs
#placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("placebo", placebo_counts, placebo_metadata,
#                                                                    "2_D28", "2_D_minus_1", data_dir)
