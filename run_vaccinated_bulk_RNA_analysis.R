library(DESeq2)
library(data.table)
library(pheatmap)

##### SETUP #####
# Read in count and metadata files
base_dir <- "C:/Users/wat2/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
data_dir <- "C:/Users/wat2/Documents/local_data_files/"
#load(paste0(data_dir, "bulk_RNA_obj.RData"))
counts <- fread(paste0(data_dir, "rsem_genes_count.processed.txt"), header = T, sep = "\t")
all_metadata_file <- paste0(base_dir, "all_metadata_sheet.tsv")
all_metadata <- read.table(all_metadata_file, header = TRUE, sep = "\t")
# Take gene_id column from counts and use contents as the rownames of counts
row.names <- as.character(counts$gene_id)
counts <- counts[,2:ncol(counts)]
counts <- as.data.frame(counts)
rownames(counts) <- row.names
# Only keep metadata for bulk RNA-seq aliquots
all_metadata <- all_metadata[all_metadata$bulkRNA_seq == TRUE,]
# Add period into time point (by itself, time point isn't unique - we need period information
# to distinguish between D-1 in period 1 vs D-1 in period 2, for example)
all_metadata$time_point <- paste0(all_metadata$period, "_", all_metadata$time_point)
# Make time point names safe for DESeq2
all_metadata$time_point[all_metadata$time_point == '1_D1 predose'] <- '1_D_minus_1'
all_metadata$time_point[all_metadata$time_point == '2_D-2'] <- '2_D_minus_2'
all_metadata$time_point[all_metadata$time_point == '2_D-1'] <- '2_D_minus_1'
# Divide metadata into vaccinated
vaccinated_metadata <- all_metadata[all_metadata$treatment == "MVA-NP+M1",]
# Find vaccinated and placebo-associated counts
kept_aliquots <- vaccinated_metadata$aliquot_id
vaccinated_counts <- counts[kept_aliquots]
# Sort columns in counts and rows for each so they're in same order (for DESeq2)
colnames(counts) <- sort(colnames(counts))
rownames(all_metadata) <- all_metadata$aliquot_id
rownames(all_metadata) <- sort(rownames(all_metadata))
colnames(vaccinated_counts) <- sort(colnames(vaccinated_counts))
rownames(vaccinated_metadata) <- vaccinated_metadata$aliquot_id
rownames(vaccinated_metadata) <- sort(rownames(vaccinated_metadata))
# Drop aliquot ID column (it's stored in rownames)
all_metadata <- subset(all_metadata, select = -c(aliquot_id))
vaccinated_metadata = subset(vaccinated_metadata, select = -c(aliquot_id))
# Probably OK to round expected counts from RSEM data. DESeq2 expects integers
counts <- round(counts)
vaccinated_counts <- round(vaccinated_counts)
# Currently not filtering sex associated genes
#sex_associated_genes <- find_sex_associated_genes(paste0(data_dir, "sex_associated_genes/"))
# Order factor levels for period 1, period 2, and all time points
period_1_factors <- c("1_D_minus_1", "1_D2", "1_D8", "1_D28")
period_1_more_vaccination_data_factors <- c("1_D_minus_1", "1_D8")
period_2_factors <- c("2_D_minus_2", "2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
period_2_without_2_D_minus_2_factors <- c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
all_factors <- c(period_1_factors, period_2_factors)


########### VACCINATED ########### 
# We will begin by using subjects that have all 10 timepoints for our tests
full_time_vaccinated_metadata <- vaccinated_metadata[vaccinated_metadata$subject_id 
                                               %in% names(table(vaccinated_metadata$subject_id)
                                                          [table(vaccinated_metadata$subject_id) == 10]),]
full_time_vaccinated_counts <- vaccinated_counts[rownames(full_time_vaccinated_metadata)]
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
period_1_vaccinated_time_point_analysis_results <- results(period_1_vaccinated_time_point_analysis, alpha = 0.05)
period_1_vaccinated_time_point_analysis_results <- period_1_vaccinated_time_point_analysis_results[order(period_1_vaccinated_time_point_analysis_results$padj),]
period_1_vaccinated_time_point_analysis_results <- subset(period_1_vaccinated_time_point_analysis_results, padj < 0.05)
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
