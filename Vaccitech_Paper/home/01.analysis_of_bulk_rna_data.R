# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

#### PERIOD 1 HIGH VIRAL LOAD ####
# We expect ~0 DEGs because placebo was used (no actual vaccination)
# 0/0/0/0/0 DEGs
high_placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D2", "1_D_minus_1", paste0(bulk_data_dir, "high_placebo_period_1_D2_vs_D_minus_1"), "high")
raw_high_placebo_period_1_D2_vs_D_minus_1_results <- high_placebo_period_1_D2_vs_D_minus_1_results[[5]]
# 0/0/0/0/0 DEGs
high_placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D8", "1_D_minus_1", paste0(bulk_data_dir, "high_placebo_period_1_D8_vs_D_minus_1"), "high")
raw_high_placebo_period_1_D8_vs_D_minus_1_results <- high_placebo_period_1_D8_vs_D_minus_1_results[[5]]
# 1/1/1/1/16 DEGs
high_placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "1_D28", "1_D_minus_1", paste0(bulk_data_dir, "high_placebo_period_1_D28_vs_D_minus_1"), "high")
raw_high_placebo_period_1_D28_vs_D_minus_1_results <- high_placebo_period_1_D28_vs_D_minus_1_results[[5]]
#### PERIOD 1 LOW VIRAL LOAD ####
# We expect ~0 DEGs because placebo was used (no actual vaccination)
# 1/1/0/0/1 DEGs
low_placebo_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "1_D2", "1_D_minus_1", paste0(bulk_data_dir, "low_placebo_period_1_D2_vs_D_minus_1"), "low")
raw_low_placebo_period_1_D2_vs_D_minus_1_results <- low_placebo_period_1_D2_vs_D_minus_1_results[[5]]
# 0/0/0/0/0 DEGs
low_placebo_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "1_D8", "1_D_minus_1", paste0(bulk_data_dir, "low_placebo_period_1_D8_vs_D_minus_1"), "low")
raw_low_placebo_period_1_D8_vs_D_minus_1_results <- low_placebo_period_1_D8_vs_D_minus_1_results[[5]]
# 19/0/0/0/63 DEGs
low_placebo_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                          "1_D28", "1_D_minus_1", paste0(bulk_data_dir, "low_placebo_period_1_D28_vs_D_minus_1"), "low")
raw_low_placebo_period_1_D28_vs_D_minus_1_results <- low_placebo_period_1_D28_vs_D_minus_1_results[[5]]

#### PERIOD 2 HIGH VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 609/1/0/0/2154 DEGs found - seems weird that so many are found
# Could there be something that happened shortly before infection that made subjects more susceptible to getting sick?
# Could it be a single patient (or two) messing everything up?
# Is it just inflammation from getting blood drawn?
# If so, why would this signal exist for HVL and not LVL?
high_placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                "2_D_minus_1", "2_D_minus_2", paste0(bulk_data_dir, "high_placebo_period_2_D_minus_1_vs_D_minus_2"), "high")
raw_high_placebo_period_2_D_minus_1_vs_D_minus_2_results <- high_placebo_period_2_D_minus_1_vs_D_minus_2_results[[5]]
# 0.585 is good, and -0.585 is good
collected_fmd_high_placebo_period_2_D_minus_1_vs_D_minus_2 <- run_fmd_on_flu_data(raw_high_placebo_period_2_D_minus_1_vs_D_minus_2_results)

# 2 D2 vs 2 D minus 1 - 299/0/0/0/1384 DEGs
high_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                         "2_D2", "2_D_minus_1", paste0(bulk_data_dir, "high_placebo_period_2_D2_vs_D_minus_1"), "high")
raw_high_placebo_period_2_D2_vs_D_minus_1_results <- high_placebo_period_2_D2_vs_D_minus_1_results[[5]]
# 0.585 is good, and -0.585 is good
collected_fmd_high_placebo_period_2_D2_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_placebo_period_2_D2_vs_D_minus_1_results)

# 2 D5 vs 2 D minus 1 - 3432/600/228/51/5276 DEGs
high_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                         "2_D5", "2_D_minus_1", paste0(bulk_data_dir, "high_placebo_period_2_D5_vs_D_minus_1"), "high")
raw_high_placebo_period_2_D5_vs_D_minus_1_results <- high_placebo_period_2_D5_vs_D_minus_1_results[[5]]
# 0.585 or 1 is good, and -0.585 or -1 is good 
collected_fmd_high_placebo_period_2_D5_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_placebo_period_2_D5_vs_D_minus_1_results, c(0.585, 1, 2))

# 2 D8 vs 2 D minus 1 - 1940/269/59/8/3314 DEGs
high_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                         "2_D8", "2_D_minus_1", paste0(bulk_data_dir, "high_placebo_period_2_D8_vs_D_minus_1"), "high")
raw_high_placebo_period_2_D8_vs_D_minus_1_results <- high_placebo_period_2_D8_vs_D_minus_1_results[[5]]
# ???
collected_fmd_high_placebo_period_2_D8_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_placebo_period_2_D8_vs_D_minus_1_results, c(0.585, 1, 2))
# 2 D28 vs 2 D minus 1 - 34/2/0/0/1394 DEGs
high_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                          "2_D28", "2_D_minus_1", paste0(bulk_data_dir, "high_placebo_period_2_D28_vs_D_minus_1"), "high")
raw_high_placebo_period_2_D28_vs_D_minus_1_results <- high_placebo_period_2_D28_vs_D_minus_1_results[[5]]
# 0.3 is good, and -0.3 is good
collected_fmd_high_placebo_period_2_D28_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_placebo_period_2_D28_vs_D_minus_1_results)

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
# Write top 2000 genes to file
high_placebo_period_2_LRT_analysis_2000_topGenes <- rownames(high_placebo_period_2_LRT_analysis_results)[1:2000]
if (!dir.exists(paste0(bulk_data_dir, "high_placebo_period_2_LRT"))) {dir.create(paste0(bulk_data_dir, "high_placebo_period_2_LRT")) }
write.table(high_placebo_period_2_LRT_analysis_2000_topGenes, file = paste0(bulk_data_dir, "high_placebo_period_2_LRT/LRT_top_200_genes.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
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
         cluster_col=FALSE, fontsize_col=14, filename = paste0(bulk_data_dir, "high_placebo_period_2_LRT/high_placebo_period_2_LRT_analysis_top_20_genes.png"))
# Plot heatmap with top 200 genes (and clean up heatmap plot by removing row names and clustering)
high_placebo_period_2_LRT_analysis_200_topGenes <- head(order(high_placebo_period_2_LRT_analysis_results_for_plotting$padj),200)
high_placebo_period_2_LRT_analysis_200_mat <- high_placebo_period_2_LRT_analysis_betas[high_placebo_period_2_LRT_analysis_200_topGenes, -c(1:13)]
colnames(high_placebo_period_2_LRT_analysis_200_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
high_placebo_period_2_LRT_analysis_200_mat[high_placebo_period_2_LRT_analysis_200_mat < -high_placebo_period_2_LRT_analysis_thr] <- -high_placebo_period_2_LRT_analysis_thr
high_placebo_period_2_LRT_analysis_200_mat[high_placebo_period_2_LRT_analysis_200_mat > high_placebo_period_2_LRT_analysis_thr] <- high_placebo_period_2_LRT_analysis_thr
pheatmap(high_placebo_period_2_LRT_analysis_200_mat, breaks=seq(from=-high_placebo_period_2_LRT_analysis_thr, to=high_placebo_period_2_LRT_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, show_rownames = FALSE, cluster_row = FALSE, filename = paste0(bulk_data_dir, "high_placebo_period_2_LRT/high_placebo_period_2_LRT_analysis_top_200_genes.png"))
#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 0/0/0/0/0 DEGs found
low_placebo_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                               "2_D_minus_1", "2_D_minus_2", paste0(bulk_data_dir, "low_placebo_period_2_D_minus_1_vs_D_minus_2"), "low")
raw_low_placebo_period_2_D_minus_1_vs_D_minus_2_results <- low_placebo_period_2_D_minus_1_vs_D_minus_2_results[[5]]
# 2 D2 vs 2 D minus 1 - 1/1/1/1/2 DEGs
low_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                        "2_D2", "2_D_minus_1", paste0(bulk_data_dir, "low_placebo_period_2_D2_vs_D_minus_1"), "low")
raw_low_placebo_period_2_D2_vs_D_minus_1_results <- low_placebo_period_2_D2_vs_D_minus_1_results[[5]]
# 2 D5 vs 2 D minus 1 - 0/0/0/0/1 DEGs
low_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                        "2_D5", "2_D_minus_1", paste0(bulk_data_dir, "low_placebo_period_2_D5_vs_D_minus_1"), "low")
raw_low_placebo_period_2_D5_vs_D_minus_1_results <- low_placebo_period_2_D5_vs_D_minus_1_results[[5]]
# 2 D8 vs 2 D minus 1 - 0/0/0/0/0 DEGs
low_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                        "2_D8", "2_D_minus_1", paste0(bulk_data_dir, "low_placebo_period_2_D8_vs_D_minus_1"), "low")
raw_low_placebo_period_2_D8_vs_D_minus_1_results <- low_placebo_period_2_D8_vs_D_minus_1_results[[5]]
# 2 D28 vs 2 D minus 1 - 51/0/0/0/647 DEGs
low_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", low_placebo_counts, low_placebo_metadata,
                                                                         "2_D28", "2_D_minus_1", paste0(bulk_data_dir, "low_placebo_period_2_D28_vs_D_minus_1"), "low")
raw_low_placebo_period_2_D28_vs_D_minus_1_results <- low_placebo_period_2_D28_vs_D_minus_1_results[[5]]
# 0.3 is good, and -0.3 is good
collected_fmd_low_placebo_period_2_D28_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_placebo_period_2_D28_vs_D_minus_1_results)

save.image(paste0(onedrive_dir, "Influenza Analysis/bulk_RNA_analysis.RData"))



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

# HIGH VS LOW VIRAL LOAD
# 2 D minus 2
# 144/143/142/138/151
# Note: If I don't control for age and sex, I find nothing.
placebo_period_2_D_minus_2_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D_minus_2", "HIGH", "LOW", paste0(bulk_data_dir, "placebo_period_2_D_minus_2_high_vs_low/"), "2_D_minus_2")
# 2 D minus 1
# 167/157/155/151/178
# Note: If I don't control for age and sex, I find nothing.
placebo_period_2_D_minus_1_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D_minus_1", "HIGH", "LOW", paste0(bulk_data_dir, "placebo_period_2_D_minus_1_high_vs_low/"), "2_D_minus_1")
# 2 D2
# 139/138/137/134/139
placebo_period_2_D2_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                                     "2_D2", "HIGH", "LOW", paste0(bulk_data_dir, "placebo_period_2_D2_high_vs_low/"), "2_D2")
# 2 D5
# Note: If I don't control for age and sex, I still find a lot of DEGs.
# https://hb.flatironinstitute.org/module/overview/?body_tag=9a43682013acc9404b08417d64f9e60f2436bfba all upregulated
# 883/269/155/88/1653
placebo_period_2_D5_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D5", "HIGH", "LOW", paste0(bulk_data_dir, "placebo_period_2_D5_high_vs_low/"), "2_D5")
# 2 D8
# https://hb.flatironinstitute.org/module/overview/?body_tag=e7397a8502c94a074b03dedb68c473ef8912a8e9 all upregulated
# 298/207/196/187/365
placebo_period_2_D8_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D8", "HIGH", "LOW", paste0(bulk_data_dir, "placebo_period_2_D8_high_vs_low/"), "2_D8")
# 2 D28
# 193/180/172/170/209
placebo_period_2_D28_high_vs_low_results <- run_deseq_bulk_analysis_viral_load("both", both_placebo_counts, both_placebo_metadata,
                                                                              "2_D28", "HIGH", "LOW", paste0(bulk_data_dir, "placebo_period_2_D8_high_vs_low/"), "2_D28")

# VACCINATED
#### PERIOD 1 HIGH T CELL RESPONSE ####
# 6234/1686/678/163/7799 DEGs
high_vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "1_D2", "1_D_minus_1", paste0(bulk_data_dir, "high_vaccinated_period_1_D2_vs_D_minus_1"), "high")
raw_high_vaccinated_period_1_D2_vs_D_minus_1_results <- high_vaccinated_period_1_D2_vs_D_minus_1_results[[5]]
# 59/15/2/0/96 DEGs
high_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "1_D8", "1_D_minus_1", paste0(bulk_data_dir, "high_vaccinated_period_1_D8_vs_D_minus_1"), "high")
raw_high_vaccinated_period_1_D8_vs_D_minus_1_results <- high_vaccinated_period_1_D8_vs_D_minus_1_results[[5]]
# 0/0/0/0/0 DEGs
high_vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                      "1_D28", "1_D_minus_1", paste0(bulk_data_dir, "high_vaccinated_period_1_D28_vs_D_minus_1"), "high")
raw_high_vaccinated_period_1_D28_vs_D_minus_1_results <- high_vaccinated_period_1_D28_vs_D_minus_1_results[[5]]
#### PERIOD 1 LOW VIRAL LOAD ####
# 2280/312/129/18/3854 DEGs
low_vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                    "1_D2", "1_D_minus_1", paste0(bulk_data_dir, "low_vaccinated_period_1_D2_vs_D_minus_1"), "low")
raw_low_vaccinated_period_1_D2_vs_D_minus_1_results <- low_vaccinated_period_1_D2_vs_D_minus_1_results[[5]]
# 70/0/0/0/521 DEGs
low_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                    "1_D8", "1_D_minus_1", paste0(bulk_data_dir, "low_vaccinated_period_1_D8_vs_D_minus_1"), "low")
raw_low_vaccinated_period_1_D8_vs_D_minus_1_results <- low_vaccinated_period_1_D8_vs_D_minus_1_results[[5]]
# 2/0/0/0/3 DEGs
low_vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                     "1_D28", "1_D_minus_1", paste0(bulk_data_dir, "low_vaccinated_period_1_D28_vs_D_minus_1"), "low")
raw_low_vaccinated_period_1_D28_vs_D_minus_1_results <- low_vaccinated_period_1_D28_vs_D_minus_1_results[[5]]

#### PERIOD 2 HIGH VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 294/0/0/0/1225 DEGs found - seems weird that so many are found
high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                            "2_D_minus_1", "2_D_minus_2", paste0(bulk_data_dir, "high_vaccinated_period_2_D_minus_1_vs_D_minus_2"), "high")
raw_high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_high_vaccinated_period_2_D_minus_1_vs_D_minus_2 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D_minus_1_vs_D_minus_2_results)

# 2 D2 vs 2 D minus 1 - 0/0/0/0/2 DEGs
high_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D2", "2_D_minus_1", paste0(bulk_data_dir, "high_vaccinated_period_2_D2_vs_D_minus_1"), "high")
raw_high_vaccinated_period_2_D2_vs_D_minus_1_results <- high_vaccinated_period_2_D2_vs_D_minus_1_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_high_vaccinated_period_2_D2_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D2_vs_D_minus_1_results)

# 2 D5 vs 2 D minus 1 - 704/20/0/0/1588 DEGs
high_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D5", "2_D_minus_1", paste0(bulk_data_dir, "high_vaccinated_period_2_D5_vs_D_minus_1"), "high")
raw_high_vaccinated_period_2_D5_vs_D_minus_1_results <- high_vaccinated_period_2_D5_vs_D_minus_1_results[[5]]
# 0.585 or 1 is good, and -0.585 or -1 is good 
#collected_fmd_high_vaccinated_period_2_D5_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D5_vs_D_minus_1_results, c(0.585, 1, 2))

# 2 D8 vs 2 D minus 1 - 1359/85/1/0/2441 DEGs
high_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                     "2_D8", "2_D_minus_1", paste0(bulk_data_dir, "high_vaccinated_period_2_D8_vs_D_minus_1"), "high")
raw_high_vaccinated_period_2_D8_vs_D_minus_1_results <- high_vaccinated_period_2_D8_vs_D_minus_1_results[[5]]
# ???
# collected_fmd_high_vaccinated_period_2_D8_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D8_vs_D_minus_1_results, c(0.585, 1, 2))
# 2 D28 vs 2 D minus 1 - 6/4/4/3/21 DEGs
high_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                                      "2_D28", "2_D_minus_1", paste0(bulk_data_dir, "high_vaccinated_period_2_D28_vs_D_minus_1"), "high")
raw_high_vaccinated_period_2_D28_vs_D_minus_1_results <- high_vaccinated_period_2_D28_vs_D_minus_1_results[[5]]
# 0.3 is good, and -0.3 is good
#collected_fmd_high_vaccinated_period_2_D28_vs_D_minus_1 <- run_fmd_on_flu_data(raw_high_vaccinated_period_2_D28_vs_D_minus_1_results)

#### PERIOD 2 LOW VIRAL LOAD ####
# 2 D minus 2 vs 2 D minus 1 - should be virtually zero unless some weird stuff happened between blood draws
# 1/0/0/0/9 DEGs found - seems weird that so many are found
low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                               "2_D_minus_1", "2_D_minus_2", paste0(bulk_data_dir, "low_vaccinated_period_2_D_minus_1_vs_D_minus_2"), "low")
raw_low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results <- low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_low_vaccinated_period_2_D_minus_1_vs_D_minus_2 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D_minus_1_vs_D_minus_2_results)

# 2 D2 vs 2 D minus 1 - 2/0/0/0/2 DEGs
low_vaccinated_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D2", "2_D_minus_1", paste0(bulk_data_dir, "low_vaccinated_period_2_D2_vs_D_minus_1"), "low")
raw_low_vaccinated_period_2_D2_vs_D_minus_1_results <- low_vaccinated_period_2_D2_vs_D_minus_1_results[[5]]
# 0.585 is good, and -0.585 is good
#collected_fmd_low_vaccinated_period_2_D2_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D2_vs_D_minus_1_results)

# 2 D5 vs 2 D minus 1 - 101/20/0/0/366 DEGs
low_vaccinated_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D5", "2_D_minus_1", paste0(bulk_data_dir, "low_vaccinated_period_2_D5_vs_D_minus_1"), "low")
raw_low_vaccinated_period_2_D5_vs_D_minus_1_results <- low_vaccinated_period_2_D5_vs_D_minus_1_results[[5]]
# 0.585 or 1 is good, and -0.585 or -1 is good 
#collected_fmd_low_vaccinated_period_2_D5_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D5_vs_D_minus_1_results, c(0.585, 1, 2))

# 2 D8 vs 2 D minus 1 - 4/0/0/0/15 DEGs
low_vaccinated_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                        "2_D8", "2_D_minus_1", paste0(bulk_data_dir, "low_vaccinated_period_2_D8_vs_D_minus_1"), "low")
#raw_low_vaccinated_period_2_D8_vs_D_minus_1_results <- low_vaccinated_period_2_D8_vs_D_minus_1_results[[5]]
# ???
# collected_fmd_low_vaccinated_period_2_D8_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D8_vs_D_minus_1_results, c(0.585, 1, 2))
# 2 D28 vs 2 D minus 1 - 106/0/0/0/800 DEGs
low_vaccinated_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                                         "2_D28", "2_D_minus_1", paste0(bulk_data_dir, "low_vaccinated_period_2_D28_vs_D_minus_1"), "low")
raw_low_vaccinated_period_2_D28_vs_D_minus_1_results <- low_vaccinated_period_2_D28_vs_D_minus_1_results[[5]]
# 0.3 is good, and -0.3 is good
#collected_fmd_low_vaccinated_period_2_D28_vs_D_minus_1 <- run_fmd_on_flu_data(raw_low_vaccinated_period_2_D28_vs_D_minus_1_results)