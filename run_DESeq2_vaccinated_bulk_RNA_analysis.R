library(DESeq2)
library(data.table)
library(pheatmap)

set.seed(1234)

##### SETUP #####
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
setup_bulk_analysis()

#### PERIOD 1 ####
# Look at Period 1 DEGs - should find a lot
# 1 D2 vs 1 D minus 1 - 7087/1163/409 DEGs found
vaccinated_period_1_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", vaccinated_counts, vaccinated_metadata,
                                                                    "1_D2", "1_D_minus_1", data_dir)
# 1 D8 vs 1 D minus 1 - 1/0/0 DEGs found
vaccinated_period_1_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", vaccinated_counts, vaccinated_metadata,
                                                                    "1_D8", "1_D_minus_1", data_dir)
# 1 D28 vs 1 D minus 1 - 2/1/1 DEGs found
vaccinated_period_1_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis("vaccinated", vaccinated_counts, vaccinated_metadata,
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
vaccinated_period_1_LRT_time_point_analysis_betas <- coef(vaccinated_period_2_LRT_analysis)
vaccinated_period_1_LRT_time_point_analysis_results_for_plotting <- results(vaccinated_period_2_LRT_analysis, alpha = 0.05)
vaccinated_period_1_LRT_time_point_analysis_20_topGenes <- head(order(vaccinated_period_1_LRT_time_point_analysis_results_for_plotting$padj),20)
vaccinated_period_1_LRT_time_point_analysis_20_mat <- vaccinated_period_1_LRT_time_point_analysis_betas[vaccinated_period_1_LRT_time_point_analysis_20_topGenes, -c(1:23)]
vaccinated_period_1_LRT_time_point_analysis_thr <- 1.5
vaccinated_period_1_LRT_time_point_analysis_20_mat[vaccinated_period_1_LRT_time_point_analysis_20_mat < -vaccinated_period_1_LRT_time_point_analysis_thr] <- -vaccinated_period_1_LRT_time_point_analysis_thr
vaccinated_period_1_LRT_time_point_analysis_20_mat[vaccinated_period_1_LRT_time_point_analysis_20_mat > vaccinated_period_1_LRT_time_point_analysis_thr] <- vaccinated_period_1_LRT_time_point_analysis_thr
colnames(vaccinated_period_1_LRT_time_point_analysis_20_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(vaccinated_period_1_LRT_time_point_analysis_20_mat, breaks=seq(from=-vaccinated_period_1_LRT_time_point_analysis_thr, to=vaccinated_period_1_LRT_time_point_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(data_dir, "vaccinated_period_1_LRT_time_point_analysis_top_20_genes.png"))
# Plot heatmap with top 200 genes (and clean up heatmap plot by removing row names and clustering)
vaccinated_period_1_LRT_time_point_analysis_200_topGenes <- head(order(vaccinated_period_1_LRT_time_point_analysis_results_for_plotting$padj),200)
vaccinated_period_1_LRT_time_point_analysis_200_mat <- vaccinated_period_1_LRT_time_point_analysis_betas[vaccinated_period_1_LRT_time_point_analysis_200_topGenes, -c(1:23)]
colnames(vaccinated_period_1_LRT_time_point_analysis_200_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
vaccinated_period_1_LRT_time_point_analysis_200_mat[vaccinated_period_1_LRT_time_point_analysis_200_mat < -vaccinated_period_1_LRT_time_point_analysis_thr] <- -vaccinated_period_1_LRT_time_point_analysis_thr
vaccinated_period_1_LRT_time_point_analysis_200_mat[vaccinated_period_1_LRT_time_point_analysis_200_mat > vaccinated_period_1_LRT_time_point_analysis_thr] <- vaccinated_period_1_LRT_time_point_analysis_thr
pheatmap(vaccinated_period_1_LRT_time_point_analysis_200_mat, breaks=seq(from=-vaccinated_period_1_LRT_time_point_analysis_thr, to=vaccinated_period_1_LRT_time_point_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, show_rownames = FALSE, cluster_row = FALSE, filename = paste0(data_dir, "vaccinated_period_1_LRT_time_point_analysis_top_200_genes.png"))
save.image(paste0(data_dir, "vaccinated_bulk_RNA_obj.RData"))
