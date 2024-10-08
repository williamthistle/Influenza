library(DESeq2)
library(data.table)
library(pheatmap)

set.seed(1234)

##### SETUP #####
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
setup_bulk_analysis()

# Grab the main 14 subjects (7 high T cell response, 7 low T cell response) that have all 10 time points
vaccinated_full_time_series_metadata <- vaccinated_metadata[vaccinated_metadata$subject_id 
                                              %in% names(table(vaccinated_metadata$subject_id)
                                                         [table(vaccinated_metadata$subject_id) == 10]),]

# Grab subject IDs for main 23 subjects
vaccinated_full_time_series_subjects <- unique(vaccinated_full_time_series_metadata$subject_id)
# Reorder subject IDs according to T cell response (high to low) - STILL NEED TO FIGURE OUT WHICH COLUMN TO USE FOR T CELL RESPONSE
vaccinated_full_time_series_subjects <- vaccinated_full_time_series_subjects[order(match(vaccinated_full_time_series_subjects,viral_load_primary$SUBJID))]
# Top 7 will be high T cell response and bottom 7 will be low T cell response 
high_t_cell_subjects <- vaccinated_full_time_series_subjects[1:7]
low_t_cell_subjects <- tail(vaccinated_full_time_series_subjects, n = 7)
# Grab high T cell vaccinated counts and metadata
high_vaccinated_aliquots <- rownames(vaccinated_metadata[vaccinated_metadata$subject_id %in% high_t_cell_subjects,])
high_vaccinated_counts <- vaccinated_counts[,high_vaccinated_aliquots]
high_vaccinated_metadata <- vaccinated_metadata[high_vaccinated_aliquots,]
# Grab low T cell vaccinated counts and metadata
low_vaccinated_aliquots <- rownames(vaccinated_metadata[vaccinated_metadata$subject_id %in% low_t_cell_subjects,])
low_vaccinated_counts <- vaccinated_counts[,low_vaccinated_aliquots]
low_vaccinated_metadata <- vaccinated_metadata[low_vaccinated_aliquots,]
# Now we are ready to run differential expression

#### PERIOD 1 HIGH T CELL RESPONSE ####
# 1 D2 vs 1 D minus 1 -  DEGs found
high_vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                    "1_D2", "1_D_minus_1", data_dir)
# 1 D8 vs 1 D minus 1 -  DEGs found
high_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                    "1_D8", "1_D_minus_1", data_dir)
# 1 D28 vs 1 D minus 1 - DEGs found
high_vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", high_vaccinated_counts, high_vaccinated_metadata,
                                                                     "1_D28", "1_D_minus_1", data_dir)

#### PERIOD 1 LOW T CELL RESPONSE ####
# 1 D2 vs 1 D minus 1 - 7087/1163/409 DEGs found
low_vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                            "1_D2", "1_D_minus_1", data_dir)
# 1 D8 vs 1 D minus 1 - 370/47/7 DEGs found
low_vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                            "1_D8", "1_D_minus_1", data_dir)
# 1 D28 vs 1 D minus 1 - 3/0/0 DEGs found
low_vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", low_vaccinated_counts, low_vaccinated_metadata,
                                                                             "1_D28", "1_D_minus_1", data_dir)
# LRT test for Period 1
# ??? DEGs found (note that LRT doesn't have a fold change threshold). 
vaccinated_period_1_LRT_metadata <- vaccinated_metadata[vaccinated_metadata$time_point == "1_D28" | vaccinated_metadata$time_point == "1_D8" | 
                                                    vaccinated_metadata$time_point == "1_D2" |
                                                    vaccinated_metadata$time_point == "1_D_minus_1",]
vaccinated_period_1_LRT_metadata <- vaccinated_period_1_LRT_metadata[vaccinated_period_1_LRT_metadata$subject_id 
                                                               %in% names(table(vaccinated_period_1_LRT_metadata$subject_id)
                                                                          [table(vaccinated_period_1_LRT_metadata$subject_id) == 4]),]
vaccinated_period_1_LRT_counts <- vaccinated_counts[rownames(vaccinated_period_1_LRT_metadata)]
vaccinated_period_1_LRT_analysis <- DESeqDataSetFromMatrix(countData = vaccinated_period_1_LRT_counts,
                                                        colData = vaccinated_period_1_LRT_metadata,
                                                        design = ~ subject_id + time_point)
vaccinated_period_1_LRT_analysis <- DESeq(vaccinated_period_1_LRT_analysis, test = "LRT", reduced = ~ subject_id)
vaccinated_period_1_LRT_analysis_results <- results(vaccinated_period_1_LRT_analysis, alpha = 0.05)
vaccinated_period_1_LRT_analysis_results <- vaccinated_period_1_LRT_analysis_results[order(vaccinated_period_1_LRT_analysis_results$padj),]
vaccinated_period_1_LRT_analysis_results <- subset(vaccinated_period_1_LRT_analysis_results, padj < 0.05)
# Plot heatmap with top 20 genes
vaccinated_period_1_LRT_time_point_analysis_betas <- coef(vaccinated_period_1_LRT_analysis)
vaccinated_period_1_LRT_time_point_analysis_results_for_plotting <- results(vaccinated_period_1_LRT_analysis, alpha = 0.05)
vaccinated_period_1_LRT_time_point_analysis_20_topGenes <- head(order(vaccinated_period_1_LRT_time_point_analysis_results_for_plotting$padj),20)
vaccinated_period_1_LRT_time_point_analysis_20_mat <- vaccinated_period_1_LRT_time_point_analysis_betas[vaccinated_period_1_LRT_time_point_analysis_20_topGenes, -c(1:14)]
vaccinated_period_1_LRT_time_point_analysis_thr <- 5
vaccinated_period_1_LRT_time_point_analysis_20_mat[vaccinated_period_1_LRT_time_point_analysis_20_mat < -vaccinated_period_1_LRT_time_point_analysis_thr] <- -vaccinated_period_1_LRT_time_point_analysis_thr
vaccinated_period_1_LRT_time_point_analysis_20_mat[vaccinated_period_1_LRT_time_point_analysis_20_mat > vaccinated_period_1_LRT_time_point_analysis_thr] <- vaccinated_period_1_LRT_time_point_analysis_thr
colnames(vaccinated_period_1_LRT_time_point_analysis_20_mat) <- c("Day 2 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(vaccinated_period_1_LRT_time_point_analysis_20_mat, breaks=seq(from=-vaccinated_period_1_LRT_time_point_analysis_thr, to=vaccinated_period_1_LRT_time_point_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(data_dir, "vaccinated_period_1_LRT_time_point_analysis_top_20_genes.png"))
# Plot heatmap with top 200 genes (and clean up heatmap plot by removing row names and clustering)
vaccinated_period_1_LRT_time_point_analysis_200_topGenes <- head(order(vaccinated_period_1_LRT_time_point_analysis_results_for_plotting$padj),200)
vaccinated_period_1_LRT_time_point_analysis_200_mat <- vaccinated_period_1_LRT_time_point_analysis_betas[vaccinated_period_1_LRT_time_point_analysis_200_topGenes, -c(1:14)]
colnames(vaccinated_period_1_LRT_time_point_analysis_200_mat) <- c("Day 2 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
vaccinated_period_1_LRT_time_point_analysis_200_mat[vaccinated_period_1_LRT_time_point_analysis_200_mat < -vaccinated_period_1_LRT_time_point_analysis_thr] <- -vaccinated_period_1_LRT_time_point_analysis_thr
vaccinated_period_1_LRT_time_point_analysis_200_mat[vaccinated_period_1_LRT_time_point_analysis_200_mat > vaccinated_period_1_LRT_time_point_analysis_thr] <- vaccinated_period_1_LRT_time_point_analysis_thr
pheatmap(vaccinated_period_1_LRT_time_point_analysis_200_mat, breaks=seq(from=-vaccinated_period_1_LRT_time_point_analysis_thr, to=vaccinated_period_1_LRT_time_point_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, show_rownames = FALSE, cluster_row = FALSE, filename = paste0(data_dir, "vaccinated_period_1_LRT_time_point_analysis_top_200_genes.png"))
save.image(paste0(data_dir, "vaccinated_bulk_RNA_obj.RData"))
