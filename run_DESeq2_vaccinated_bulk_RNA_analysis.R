library(DESeq2)
library(data.table)
library(pheatmap)

set.seed(1234)

##### SETUP #####
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
setup_bulk_analysis()

#### PERIOD 1 ####
# Look at Period 1 DEGs - should find a bunch since vaccination did occur
# 1 D2 vs 1 D minus 1 - 409 DEGs found
vaccinated_period_1_D2_vs_D_minus_1_metadata <- vaccinated_metadata[vaccinated_metadata$time_point == "1_D2" | vaccinated_metadata$time_point == "1_D_minus_1",]
vaccinated_period_1_D2_vs_D_minus_1_metadata <- vaccinated_period_1_D2_vs_D_minus_1_metadata[vaccinated_period_1_D2_vs_D_minus_1_metadata$subject_id 
                                                                                                     %in% names(table(vaccinated_period_1_D2_vs_D_minus_1_metadata$subject_id)
                                                                                                                [table(vaccinated_period_1_D2_vs_D_minus_1_metadata$subject_id) == 2]),]
vaccinated_period_1_D2_vs_D_minus_1_counts <- vaccinated_counts[rownames(vaccinated_period_1_D2_vs_D_minus_1_metadata)]
vaccinated_period_1_D2_vs_D_minus_1_analysis <- DESeqDataSetFromMatrix(countData = vaccinated_period_1_D2_vs_D_minus_1_counts,
                                                                    colData = vaccinated_period_1_D2_vs_D_minus_1_metadata,
                                                                    design = ~ subject_id + time_point)
vaccinated_period_1_D2_vs_D_minus_1_analysis <- DESeq(vaccinated_period_1_D2_vs_D_minus_1_analysis)
vaccinated_period_1_D2_vs_D_minus_1_analysis_results <- results(vaccinated_period_1_D2_vs_D_minus_1_analysis, contrast = c("time_point", "1_D2", "1_D_minus_1"), alpha = 0.05, lfcThreshold = 1)
vaccinated_period_1_D2_vs_D_minus_1_analysis_results <- vaccinated_period_1_D2_vs_D_minus_1_analysis_results[order(vaccinated_period_1_D2_vs_D_minus_1_analysis_results$padj),]
vaccinated_period_1_D2_vs_D_minus_1_analysis_results <- subset(vaccinated_period_1_D2_vs_D_minus_1_analysis_results, padj < 0.05)

##### ALL 10 TIMEPOINTS (PERIOD 1, PERIOD 2, BOTH PERIODS) #####
# We will begin by using subjects that have all 10 timepoints for our tests
full_time_vaccinated_metadata <- vaccinated_metadata[vaccinated_metadata$subject_id 
                                               %in% names(table(vaccinated_metadata$subject_id)
                                                          [table(vaccinated_metadata$subject_id) == 10]),]
full_time_vaccinated_counts <- vaccinated_counts[rownames(full_time_vaccinated_metadata)]
# Label each patient as low or high viral load?
full_time_vaccinated_metadata$time_point <- as.factor(full_time_vaccinated_metadata$time_point)
levels(full_time_vaccinated_metadata$time_point) <- all_factors
full_time_vaccinated_metadata$sex <- as.factor(full_time_vaccinated_metadata$sex)
full_time_vaccinated_metadata$age <- as.factor(full_time_vaccinated_metadata$age)
full_time_vaccinated_time_point_analysis <- DESeqDataSetFromMatrix(countData = full_time_vaccinated_counts,
                                                                  colData = full_time_vaccinated_metadata,
                                                                  design = ~ time_point + sex + age)
pca_vst <- vst(full_time_vaccinated_time_point_analysis)
# Clearly, we have two groups via PC2 (high and load viral load?)
plotPCA(pca_vst, intgroup = c("time_point"))
pca_vst_nums <- assay(pca_vst)
pcs <- prcomp(t(pca_vst_nums))
pcs <- pcs$x
pcs <- as.data.frame(pcs)



# LRT TESTS
# Period 1 (pre- and post-vaccination)
# Only keep period 1 metadata and vaccinated_counts
period_1_vaccinated_metadata <- full_time_vaccinated_metadata[full_time_vaccinated_metadata$period == 1,]
period_1_vaccinated_counts <- vaccinated_counts[rownames(period_1_vaccinated_metadata)]
# Factorize time point (with associated factor levels) and sex
period_1_vaccinated_metadata$time_point <- as.factor(period_1_vaccinated_metadata$time_point)
levels(period_1_vaccinated_metadata$time_point) <- period_1_factors
period_1_vaccinated_metadata$sex <- as.factor(period_1_vaccinated_metadata$sex)
period_1_vaccinated_metadata$age <- as.factor(period_1_vaccinated_metadata$age)
# Run DESeq2 analysis
period_1_vaccinated_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_1_vaccinated_counts,
                                                       colData = period_1_vaccinated_metadata,
                                                       design = ~ time_point + sex + age)
period_1_vaccinated_time_point_analysis <- DESeq(period_1_vaccinated_time_point_analysis, test="LRT", reduced = ~ sex + age)
# DEBUGGING
period_1_vaccinated_time_point_analysis <- DESeq(period_1_vaccinated_time_point_analysis)
period_1_vaccinated_time_point_analysis_results <- results(period_1_vaccinated_time_point_analysis, contrast = c("time_point", "1_D8", "1_D_minus_1"), 
                                                                        alpha = 0.05, lfcThreshold = 1)
period_1_vaccinated_time_point_analysis_results <- period_1_vaccinated_time_point_analysis_results[order(period_1_vaccinated_time_point_analysis_results$padj),]
period_1_vaccinated_time_point_analysis_results <- subset(period_1_vaccinated_time_point_analysis_results, padj < 0.05)
# DEBUGGING
period_1_vaccinated_time_point_analysis_results <- results(period_1_vaccinated_time_point_analysis, alpha = 0.05)
period_1_vaccinated_time_point_analysis_results <- period_1_vaccinated_time_point_analysis_results[order(period_1_vaccinated_time_point_analysis_results$padj),]
period_1_vaccinated_time_point_analysis_results <- subset(period_1_vaccinated_time_point_analysis_results, padj < 0.05)
#period_1_vaccinated_time_point_analysis_results <- results(period_1_vaccinated_time_point_analysis, contrast = c("time_point", "1_D8", "1_D_minus_1"), 
#                                                           alpha = 0.05, lfcThreshold = 1)

# There are a lot of changes! But do changes persist over time?
# Let's create a heatmap to see!
period_1_vaccinated_betas <- coef(period_1_vaccinated_time_point_analysis)
period_1_vaccinated_time_point_analysis_results_for_plotting <- results(period_1_vaccinated_time_point_analysis, alpha = 0.05)
period_1_vaccinated_topGenes <- head(order(period_1_vaccinated_time_point_analysis_results_for_plotting$padj),20)
period_1_vaccinated_mat <- period_1_vaccinated_betas[period_1_vaccinated_topGenes, -c(1, 5, 6, 7, 8)]
period_1_vaccinated_thr <- 6 
period_1_vaccinated_mat[period_1_vaccinated_mat < -period_1_vaccinated_thr] <- -period_1_vaccinated_thr
period_1_vaccinated_mat[period_1_vaccinated_mat > period_1_vaccinated_thr] <- period_1_vaccinated_thr
colnames(period_1_vaccinated_mat) <- c("Day 2 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(period_1_vaccinated_mat, breaks=seq(from=-period_1_vaccinated_thr, to=period_1_vaccinated_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(data_dir, "period_1_time_point_analysis_vaccinated_top_20_genes.png"))
# It looks like changes are heightened in Day 2 and may not persist
# We can do Wald tests to get a better sense of how each day compares to our baseline (pairwise)
period_1_vaccinated_wald_tests <- c()
period_1_vaccinated_wald_test_names <- resultsNames(period_1_vaccinated_time_point_analysis)
# Remove first and last entries (intercept and sex)
period_1_vaccinated_wald_test_names <- period_1_vaccinated_wald_test_names[c(-1, -5, -6, -7, -8)]
index <- 1
for (current_name in period_1_vaccinated_wald_test_names) {
  current_results <- results(period_1_vaccinated_time_point_analysis, name = current_name, test = "Wald", 
                             alpha = 0.05, lfcThreshold = 1) # NOTE LFC THRESHOLD OF 1
  current_results <- current_results[order(current_results$padj),]
  current_results <- subset(current_results, padj < 0.05)
  period_1_vaccinated_wald_tests[[index]] <- current_results
  index <- index + 1
}
# We see 750 in Day 2, 194 in Day 8, and 209 in Day 28
# Note that without LFC threshold, we see 8458 in Day 2, 216 in Day 8, and 249 in Day 28
###
# Next, since we have more samples for Day -1 and Day 8, let's do a Wald test 
# between the samples in those categories to see if we pick up any differences
# Find subject(s) to exclude (at least one subject doesn't have both 1_D8 and 1_D_minus_1 for some reason)
excluded_subjects <- setdiff(vaccinated_metadata[vaccinated_metadata$time_point == "1_D_minus_1",]$subject_id, vaccinated_metadata[vaccinated_metadata$time_point == "1_D8",]$subject_id)
period_1_vaccinated_more_samples_metadata <- vaccinated_metadata[vaccinated_metadata$period == 1,]
period_1_vaccinated_more_samples_metadata <- period_1_vaccinated_more_samples_metadata[(period_1_vaccinated_more_samples_metadata$time_point == "1_D_minus_1" | 
                                                                                          period_1_vaccinated_more_samples_metadata$time_point == "1_D8"),]
period_1_vaccinated_more_samples_metadata <- subset(period_1_vaccinated_more_samples_metadata,!(subject_id==excluded_subjects))
period_1_vaccinated_more_samples_counts <- vaccinated_counts[rownames(period_1_vaccinated_more_samples_metadata)]
period_1_vaccinated_more_samples_metadata$time_point <- as.factor(period_1_vaccinated_more_samples_metadata$time_point)
levels(period_1_vaccinated_more_samples_metadata$time_point) <- period_1_more_vaccination_data_factors
period_1_vaccinated_more_samples_metadata$sex <- as.factor(period_1_vaccinated_more_samples_metadata$sex)
period_1_vaccinated_more_samples_metadata$age <- as.factor(period_1_vaccinated_more_samples_metadata$age)
period_1_vaccinated_more_samples_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_1_vaccinated_more_samples_counts,
                                                   colData = period_1_vaccinated_more_samples_metadata,
                                                   design = ~ time_point + sex + age)
period_1_vaccinated_more_samples_time_point_analysis <- DESeq(period_1_vaccinated_more_samples_time_point_analysis)
period_1_vaccinated_more_samples_time_point_analysis_results <- results(period_1_vaccinated_more_samples_time_point_analysis, contrast = c("time_point", "1_D8", "1_D_minus_1"), 
                           alpha = 0.05, lfcThreshold = 1)
period_1_vaccinated_more_samples_time_point_analysis_results <- period_1_vaccinated_more_samples_time_point_analysis_results[order(period_1_vaccinated_more_samples_time_point_analysis_results$padj),]
period_1_vaccinated_more_samples_time_point_analysis_results <- subset(period_1_vaccinated_more_samples_time_point_analysis_results, padj < 0.05)
# Why do we only get 1 gene here (1546 without LFC threshold) when our analysis above has 194 with much larger effect size?
# I should run this exact same analysis with LRT followed by Wald test to see if results are the exact same
period_1_vaccinated_more_samples_time_point_analysis_LRT <- DESeqDataSetFromMatrix(countData = period_1_vaccinated_more_samples_counts,
                                                                               colData = period_1_vaccinated_more_samples_metadata,
                                                                               design = ~ time_point + sex + age)
period_1_vaccinated_more_samples_time_point_analysis_LRT <- DESeq(period_1_vaccinated_more_samples_time_point_analysis_LRT, test = "LRT", reduced = ~ sex + age)
period_1_vaccinated_more_samples_time_point_analysis_results_LRT <- results(period_1_vaccinated_more_samples_time_point_analysis_LRT, alpha = 0.05)
period_1_vaccinated_more_samples_time_point_analysis_results_LRT <- period_1_vaccinated_more_samples_time_point_analysis_results_LRT[order(period_1_vaccinated_more_samples_time_point_analysis_results_LRT$padj),]
period_1_vaccinated_more_samples_time_point_analysis_results_LRT <- subset(period_1_vaccinated_more_samples_time_point_analysis_results_LRT, padj < 0.05)
current_results <- results(period_1_vaccinated_more_samples_time_point_analysis_LRT, name = "time_point_1_D8_vs_1_D_minus_1", test = "Wald", 
                           alpha = 0.05, lfcThreshold = 1) # NOTE LFC THRESHOLD OF 1
current_results <- current_results[order(current_results$padj),]
current_results <- subset(current_results, padj < 0.05)
# We still only get 1 gene (slightly different adjusted p-value, but whatever)
# Let's try running LRT + Wald Test on D8 vs D-1 for fewer samples dataset - we get 0 genes
# Test randomly selected subset of 14 subjects to compare to above analysis
random_subjects <- sample(unique(period_1_vaccinated_more_samples_metadata$subject_id), 14)
period_1_vaccinated_more_samples_subset_metadata <- subset(period_1_vaccinated_more_samples_metadata,subject_id %in% random_subjects)
period_1_vaccinated_more_samples_subset_metadata$time_point <- as.factor(period_1_vaccinated_more_samples_subset_metadata$time_point)
period_1_vaccinated_more_samples_subset_counts <- vaccinated_counts[rownames(period_1_vaccinated_more_samples_subset_metadata)]
levels(period_1_vaccinated_more_samples_subset_metadata$time_point) <- period_1_more_vaccination_data_factors
period_1_vaccinated_more_samples_subset_metadata$sex <- as.factor(period_1_vaccinated_more_samples_subset_metadata$sex)
period_1_vaccinated_more_samples_subset_metadata$age <- as.factor(period_1_vaccinated_more_samples_subset_metadata$age)
period_1_vaccinated_more_samples_subset_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_1_vaccinated_more_samples_subset_counts,
                                                                               colData = period_1_vaccinated_more_samples_subset_metadata,
                                                                               design = ~ time_point + sex + age)
period_1_vaccinated_more_samples_subset_time_point_analysis <- DESeq(period_1_vaccinated_more_samples_subset_time_point_analysis)
period_1_vaccinated_more_samples_subset_time_point_analysis_results <- results(period_1_vaccinated_more_samples_subset_time_point_analysis, contrast = c("time_point", "1_D8", "1_D_minus_1"), 
                                                                        alpha = 0.05, lfcThreshold = 1)
period_1_vaccinated_more_samples_subset_time_point_analysis_results <- period_1_vaccinated_more_samples_subset_time_point_analysis_results[order(period_1_vaccinated_more_samples_subset_time_point_analysis_results$padj),]
period_1_vaccinated_more_samples_subset_time_point_analysis_results <- subset(period_1_vaccinated_more_samples_subset_time_point_analysis_results, padj < 0.05)
# 0 genes (24 genes without LFC threshold)

save.image(paste0(data_dir, "vaccinated_bulk_RNA_obj.RData"))
