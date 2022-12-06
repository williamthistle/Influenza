library(DESeq2)
library(data.table)

##### SETUP #####
# Read in count and metadata files
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
data_dir <- "C:/Users/willi/Documents/local_data_files/"
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
# Divide metadata into vaccinated and placebo
vaccinated_metadata <- all_metadata[all_metadata$treatment == "MVA-NP+M1",]
placebo_metadata <- all_metadata[all_metadata$treatment == "PLACEBO",]
# Find vaccinated and placebo-associated counts
kept_aliquots <- vaccinated_metadata$aliquot_id
vaccinated_counts <- counts[kept_aliquots]
kept_aliquots <- placebo_metadata$aliquot_id
placebo_counts <- counts[kept_aliquots]
# Sort columns in counts and rows for each so they're in same order (for DESeq2)
colnames(vaccinated_counts) <- sort(colnames(vaccinated_counts))
rownames(vaccinated_metadata) <- vaccinated_metadata$aliquot_id
rownames(vaccinated_metadata) <- sort(rownames(vaccinated_metadata))
colnames(placebo_counts) <- sort(colnames(placebo_counts))
rownames(placebo_metadata) <- placebo_metadata$aliquot_id
rownames(placebo_metadata) <- sort(rownames(placebo_metadata))
# Drop aliquot ID column (it's stored in rownames)
vaccinated_metadata = subset(vaccinated_metadata, select = -c(aliquot_id))
placebo_metadata = subset(placebo_metadata, select = -c(aliquot_id))
# Probably OK to round expected counts from RSEM data. DESeq2 expects integers
vaccinated_counts <- round(vaccinated_counts)
placebo_counts <- round(placebo_counts)
# Currently not filtering sex associated genes
#sex_associated_genes <- find_sex_associated_genes(paste0(data_dir, "sex_associated_genes/"))

# Order factor levels for period 1, period 2, and all time points
period_1_factors <- c("1_D_minus_1", "1_D2", "1_D8", "1_D28")
period_2_factors <- c("2_D_minus_2", "2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
period_2_without_2_D_minus_2_factors <- c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
all_factors <- c(period_1_factors, period_2_factors)


########### VACCINATED ########### 
# Basically, the only question we're concerned about currently is - what does differential expression look like
# during period 1 (vaccination)?
# First, we will use subjects that have all 10 timepoints for our tests
full_time_vaccinated_metadata <- vaccinated_metadata[vaccinated_metadata$subject_id 
                                               %in% names(table(vaccinated_metadata$subject_id)
                                                          [table(vaccinated_metadata$subject_id) == 10]),]
full_time_vaccinated_counts <- vaccinated_counts[rownames(full_time_vaccinated_metadata)]
# LRT TESTS
# It may make sense to focus primarily on Period 2 (pre- and post challenge)
# since nothing happens during Period 1 for vaccinated subjects (good check that our differential expression 
# is working properly at least)
# PERIOD 1
# Only keep period 1 metadata and vaccinated_counts
period_1_vaccinated_metadata <- full_time_vaccinated_metadata[full_time_vaccinated_metadata$period == 1,]
period_1_vaccinated_counts <- vaccinated_counts[rownames(period_1_vaccinated_metadata)]
# Factorize time point (with associated factor levels) and sex
period_1_vaccinated_metadata$time_point <- as.factor(period_1_vaccinated_metadata$time_point)
levels(period_1_vaccinated_metadata$time_point) <- period_1_factors
period_1_vaccinated_metadata$sex <- as.factor(period_1_vaccinated_metadata$sex)
# Run DESeq2 analysis
period_1_vaccinated_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_1_vaccinated_counts,
                                                       colData = period_1_vaccinated_metadata,
                                                       design = ~ time_point + sex)
period_1_vaccinated_time_point_analysis <- DESeq(period_1_vaccinated_time_point_analysis, test="LRT", reduced = ~ sex)
period_1_vaccinated_time_point_analysis_results <- results(period_1_vaccinated_time_point_analysis, alpha = 0.05)
period_1_vaccinated_time_point_analysis_results <- period_1_vaccinated_time_point_analysis_results[order(period_1_vaccinated_time_point_analysis_results$padj),]
# Note that log2FoldChange is not part of LRT, so we should just ignore it
period_1_vaccinated_time_point_analysis_results <- subset(period_1_vaccinated_time_point_analysis_results, padj < 0.05)
# There are a lot of changes! Looks like we're wise to avoid using vaccinated subjects

########### PLACEBO ########### 
##### ALL 10 TIMEPOINTS (PERIOD 1, PERIOD 2, BOTH PERIODS) #####
# First, we will use subjects that have all 10 timepoints for our tests
full_time_placebo_metadata <- placebo_metadata[placebo_metadata$subject_id 
                                               %in% names(table(placebo_metadata$subject_id)
                                                          [table(placebo_metadata$subject_id) == 10]),]
full_time_placebo_counts <- placebo_counts[rownames(full_time_placebo_metadata)]
# LRT TESTS
# It may make sense to focus primarily on Period 2 (pre- and post challenge)
# since nothing happens during Period 1 for placebo subjects (good check that our differential expression 
# is working properly at least)
# PERIOD 1
# Only keep period 1 metadata and placebo_counts
period_1_placebo_metadata <- full_time_placebo_metadata[full_time_placebo_metadata$period == 1,]
period_1_placebo_counts <- placebo_counts[rownames(period_1_placebo_metadata)]
# Factorize time point (with associated factor levels) and sex
period_1_placebo_metadata$time_point <- as.factor(period_1_placebo_metadata$time_point)
levels(period_1_placebo_metadata$time_point) <- period_1_factors
period_1_placebo_metadata$sex <- as.factor(period_1_placebo_metadata$sex)
# Run DESeq2 analysis
period_1_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_1_placebo_counts,
                              colData = period_1_placebo_metadata,
                              design = ~ time_point + sex)
period_1_time_point_analysis <- DESeq(period_1_time_point_analysis, test="LRT", reduced = ~ sex)
period_1_time_point_analysis_results <- results(period_1_time_point_analysis, alpha = 0.05)
period_1_time_point_analysis_results <- period_1_time_point_analysis_results[order(period_1_time_point_analysis_results$padj),]
# Note that log2FoldChange is not part of LRT, so we should just ignore it
period_1_time_point_analysis_results <- subset(period_1_time_point_analysis_results, padj < 0.05)
# PERIOD 2
# Only keep period 2 metadata and placebo_counts
period_2_placebo_metadata <- full_time_placebo_metadata[full_time_placebo_metadata$period == 2,]
period_2_placebo_counts <- placebo_counts[rownames(period_2_placebo_metadata)]
# Factorize time point (with associated factor levels) and sex
period_2_placebo_metadata$time_point <- as.factor(period_2_placebo_metadata$time_point)
levels(period_2_placebo_metadata$time_point) <- period_2_factors
period_2_placebo_metadata$sex <- as.factor(period_2_placebo_metadata$sex)
# Run DESeq2 analysis
period_2_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_2_placebo_counts,
                                                       colData = period_2_placebo_metadata,
                                                       design = ~ time_point + sex)
period_2_time_point_analysis <- DESeq(period_2_time_point_analysis, test="LRT", reduced = ~ sex)
period_2_time_point_analysis_results <- results(period_2_time_point_analysis, alpha = 0.05)
period_2_time_point_analysis_results <- period_2_time_point_analysis_results[order(period_2_time_point_analysis_results$padj),]
# Note that log2FoldChange is not part of LRT, so we should just ignore it
period_2_time_point_analysis_results <- subset(period_2_time_point_analysis_results, padj < 0.05)
# PERIOD 2 MINUS 2 D-2
# Because 2 D-2 and 2 D-1 are quite different, for now, let's just use 2 D-1 for Period 2 analysis
# Remove 2 D-2
period_2_without_2_D_minus_2_placebo_metadata <- full_time_placebo_metadata[full_time_placebo_metadata$period == 2,]
period_2_without_2_D_minus_2_placebo_metadata <- period_2_without_2_D_minus_2_placebo_metadata[period_2_without_2_D_minus_2_placebo_metadata$time_point != "2_D_minus_2",]
period_2_without_2_D_minus_2_placebo_counts <- placebo_counts[rownames(period_2_without_2_D_minus_2_placebo_metadata)]
# Factorize time point (with associated factor levels) and sex
period_2_without_2_D_minus_2_placebo_metadata$time_point <- as.factor(period_2_without_2_D_minus_2_placebo_metadata$time_point)
levels(period_2_without_2_D_minus_2_placebo_metadata$time_point) <- period_2_without_2_D_minus_2_factors
period_2_without_2_D_minus_2_placebo_metadata$sex <- as.factor(period_2_without_2_D_minus_2_placebo_metadata$sex)
# Run DESeq2 analysis
period_2_without_2_D_minus_2_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_2_without_2_D_minus_2_placebo_counts,
                                                       colData = period_2_without_2_D_minus_2_placebo_metadata,
                                                       design = ~ time_point + sex)
period_2_without_2_D_minus_2_time_point_analysis <- DESeq(period_2_without_2_D_minus_2_time_point_analysis, test="LRT", reduced = ~ sex)
period_2_without_2_D_minus_2_time_point_analysis_results <- results(period_2_without_2_D_minus_2_time_point_analysis, alpha = 0.05)
period_2_without_2_D_minus_2_time_point_analysis_results <- period_2_without_2_D_minus_2_time_point_analysis_results[order(period_2_without_2_D_minus_2_time_point_analysis_results$padj),]
# Note that log2FoldChange is not part of LRT, so we should just ignore it
period_2_without_2_D_minus_2_time_point_analysis_results <- subset(period_2_without_2_D_minus_2_time_point_analysis_results, padj < 0.05)
# Create heatmap
betas <- coef(period_2_without_2_D_minus_2_time_point_analysis)
period_2_without_2_D_minus_2_time_point_analysis_results_for_plotting <- results(period_2_without_2_D_minus_2_time_point_analysis, alpha = 0.05)
topGenes <- head(order(period_2_without_2_D_minus_2_time_point_analysis_results_for_plotting$padj),20)
mat <- betas[topGenes, -c(1, 6)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
colnames(mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE, fontsize_col=14)
# BOTH PERIODS (ALL TIME POINTS)
# Factorize time point (with associated factor levels) and sex
full_time_placebo_metadata$time_point <- as.factor(full_time_placebo_metadata$time_point)
levels(full_time_placebo_metadata$time_point) <- all_factors
full_time_placebo_metadata$sex <- as.factor(full_time_placebo_metadata$sex)
# Run DESeq2 analysis
time_point_analysis <- DESeqDataSetFromMatrix(countData = full_time_placebo_counts,
                                                       colData = full_time_placebo_metadata,
                                                       design = ~ time_point + sex)
time_point_analysis <- DESeq(time_point_analysis, test="LRT", reduced=~sex)
time_point_analysis_results <- results(time_point_analysis, alpha = 0.05)
time_point_analysis_results <- time_point_analysis_results[order(time_point_analysis_results$padj),]
# Note that log2FoldChange is not part of LRT, so we should just ignore it
time_point_analysis_results <- subset(time_point_analysis_results, padj < 0.05)
# WALD (PAIRWISE) TESTS
wald_time_point_analysis <- DESeqDataSetFromMatrix(countData = full_time_placebo_counts,
                                              colData = full_time_placebo_metadata,
                                              design = ~ time_point + sex)
wald_time_point_analysis <- DESeq(wald_time_point_analysis)
period_1_baseline <- "1_D_minus_1"
period_2_baseline <- "2_D_minus_2"
period_1_contrasts <- c("1_D2", "1_D8", "1_D28")
period_2_contrasts <- c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
period_1_wald_tests <- c()
period_2_wald_tests <- c()
# Period 1 Wald tests
index <- 1
for (current_contrast in period_1_contrasts) {
  current_results <- results(wald_time_point_analysis, contrast = c("time_point", current_contrast, period_1_baseline), 
                                                                    alpha = 0.05)
  current_results <- current_results[order(current_results$padj),]
  current_results <- subset(current_results, padj < 0.05)
  period_1_wald_tests[[index]] <- current_results
  index <- index + 1
}
# Period 2 Wald tests
index <- 1
for (current_contrast in period_2_contrasts) {
  current_results <- results(wald_time_point_analysis, contrast = c("time_point", current_contrast, period_2_baseline), 
                             alpha = 0.05)
  current_results <- current_results[order(current_results$padj),]
  current_results <- subset(current_results, padj < 0.05)
  period_2_wald_tests[[index]] <- current_results
  index <- index + 1
}



##### MORE DATA TIMEPOINTS (PERIOD 2) #####
# We have more data (46 aliquots versus 23) for 2_D-1, 2_D8, and 2_D28, 
# so let's do LRT and Wald tests for that subset specifically
period_2_more_metadata <- placebo_metadata[(placebo_metadata$time_point == "2_D_minus_1" | 
                                          placebo_metadata$time_point == "2_D8" | 
                                          placebo_metadata$time_point == "2_D28"),]
period_2_more_placebo_counts <- placebo_counts[rownames(period_2_more_metadata)]
period_2_more_metadata$time_point <- as.factor(period_2_more_metadata$time_point)
levels(period_2_more_metadata$time_point) <- c("2_D_minus_1", "2_D8", "2_D28")
period_2_more_metadata$sex <- as.factor(period_2_more_metadata$sex)
# Run DESeq2 analysis
period_2_more_data_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_2_more_placebo_counts,
                                                       colData = period_2_more_metadata,
                                                       design = ~ time_point + sex)
period_2_more_data_time_point_analysis <- DESeq(period_2_more_data_time_point_analysis, test="LRT", reduced=~sex)
period_2_more_data_time_point_analysis_results <- results(period_2_more_data_time_point_analysis, alpha = 0.05)
period_2_more_data_time_point_analysis_results <- period_2_more_data_time_point_analysis_results[order(period_2_more_data_time_point_analysis_results$padj),]
# Note that log2FoldChange is not part of LRT, so we should just ignore it
period_2_more_data_time_point_analysis_results <- subset(period_2_more_data_time_point_analysis_results, padj < 0.05)
# Let's compare this result to 23 profiles for same time points and see whether we get more data than overall period 2 analysis
period_2_subset_placebo_metadata <- period_2_placebo_metadata[(period_2_placebo_metadata$time_point == "2_D_minus_1" | 
                                                                 period_2_placebo_metadata$time_point == "2_D8" | 
                                                                 period_2_placebo_metadata$time_point == "2_D28"),]
period_2_subset_placebo_counts <- placebo_counts[rownames(period_2_subset_placebo_metadata)]
# Factorize time point (with associated factor levels) and sex
period_2_subset_placebo_metadata$time_point <- as.character(period_2_subset_placebo_metadata$time_point)
period_2_subset_placebo_metadata$time_point <- factor(period_2_subset_placebo_metadata$time_point)
levels(period_2_subset_placebo_metadata$time_point) <- c("2_D_minus_1", "2_D8", "2_D28")
period_2_subset_placebo_metadata$sex <- as.factor(period_2_subset_placebo_metadata$sex)
# Run DESeq2 analysis
period_2_subset_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_2_subset_placebo_counts,
                                                       colData = period_2_subset_placebo_metadata,
                                                       design = ~ time_point + sex)
period_2_subset_time_point_analysis <- DESeq(period_2_subset_time_point_analysis, test="LRT", reduced=~sex)
period_2_subset_time_point_analysis_results <- results(period_2_subset_time_point_analysis, alpha = 0.05)
period_2_subset_time_point_analysis_results <- period_2_subset_time_point_analysis_results[order(period_2_subset_time_point_analysis_results$padj),]
# Note that log2FoldChange is not part of LRT, so we should just ignore it
period_2_subset_time_point_analysis_results <- subset(period_2_subset_time_point_analysis_results, padj < 0.05)

