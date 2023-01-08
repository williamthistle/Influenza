library(DESeq2)
library(data.table)
library(pheatmap)

##### SETUP #####
# Read in count and metadata files
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
data_dir <- "C:/Users/willi/Documents/local_data_files/"
#load(paste0(data_dir, "placebo_bulk_RNA_obj.RData"))
counts <- fread(paste0(data_dir, "rsem_genes_count.processed.txt"), header = T, sep = "\t")
all_metadata_file <- paste0(base_dir, "all_metadata_sheet.tsv")
all_metadata <- read.table(all_metadata_file, header = TRUE, sep = "\t")
viral_load_file <- paste0(base_dir, "bulk_RNA_viral_load.tsv")
viral_load <- read.table(viral_load_file, header = TRUE, sep = "\t")
viral_load_primary <- viral_load[viral_load$PARAMCD == "QPCRAUC",]
viral_load_primary <- viral_load_primary[viral_load_primary$TRT01A == "PLACEBO",]
viral_load_primary$AVAL <- as.numeric(viral_load_primary$AVAL)
# Organize by viral load (high to low) and grab top 13 subjects - they will be high
viral_load_primary <- viral_load_primary[order(viral_load_primary$AVAL, decreasing = TRUE),]
high_viral_load_subjects <- viral_load_primary$SUBJID[1:13]
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
# Divide metadata into placebo
placebo_metadata <- all_metadata[all_metadata$treatment == "PLACEBO",]
viral_load_vector <- c()
for (subject_id in placebo_metadata$subject_id) {
  if(subject_id %in% relevant_subjects) {
    viral_load_vector <- c(viral_load_vector, "HIGH")
  } else {
    viral_load_vector <- c(viral_load_vector, "LOW")
  }
}
placebo_metadata$viral_load <- viral_load_vector
# Find placebo-associated counts
kept_aliquots <- placebo_metadata$aliquot_id
placebo_counts <- counts[kept_aliquots]
# Sort columns in counts and rows for each so they're in same order (for DESeq2)
colnames(counts) <- sort(colnames(counts))
rownames(all_metadata) <- all_metadata$aliquot_id
rownames(all_metadata) <- sort(rownames(all_metadata))
colnames(placebo_counts) <- sort(colnames(placebo_counts))
rownames(placebo_metadata) <- placebo_metadata$aliquot_id
rownames(placebo_metadata) <- sort(rownames(placebo_metadata))
# Drop aliquot ID column (it's stored in rownames)
all_metadata <- subset(all_metadata, select = -c(aliquot_id))
placebo_metadata = subset(placebo_metadata, select = -c(aliquot_id))
# Probably OK to round expected counts from RSEM data. DESeq2 expects integers
counts <- round(counts)
placebo_counts <- round(placebo_counts)
# Currently not filtering sex associated genes
#sex_associated_genes <- find_sex_associated_genes(paste0(data_dir, "sex_associated_genes/"))
# Order factor levels for period 1, period 2, and all time points
period_1_factors <- c("1_D_minus_1", "1_D2", "1_D8", "1_D28")
period_1_more_vaccination_data_factors <- c("1_D_minus_1", "1_D8")
period_2_factors <- c("2_D_minus_2", "2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
period_2_without_2_D_minus_2_factors <- c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
all_factors <- c(period_1_factors, period_2_factors)


########### PLACEBO ########### 
##### ALL 10 TIMEPOINTS (PERIOD 1, PERIOD 2, BOTH PERIODS) #####
# First, we will use subjects that have all 10 timepoints for our tests
full_time_placebo_metadata <- placebo_metadata[placebo_metadata$subject_id 
                                               %in% names(table(placebo_metadata$subject_id)
                                                          [table(placebo_metadata$subject_id) == 10]),]
full_time_placebo_counts <- placebo_counts[rownames(full_time_placebo_metadata)]
# Label each patient as low or high viral load?
full_time_placebo_metadata$time_point <- as.factor(full_time_placebo_metadata$time_point)
levels(full_time_placebo_metadata$time_point) <- all_factors
full_time_placebo_metadata$sex <- as.factor(full_time_placebo_metadata$sex)
full_time_placebo_metadata$age <- as.factor(full_time_placebo_metadata$age)
full_time_placebo_metadata$viral_load <- as.factor(full_time_placebo_metadata$viral_load)
# DEBUGGING
#full_time_placebo_metadata <- full_time_placebo_metadata[full_time_placebo_metadata$time_point == "2_D8" | full_time_placebo_metadata$time_point == "2_D_minus_1",]
#full_time_placebo_counts <- placebo_counts[rownames(full_time_placebo_metadata)]
# DEBUGGING
full_time_placebo_time_point_analysis <- DESeqDataSetFromMatrix(countData = full_time_placebo_counts,
                                                                   colData = full_time_placebo_metadata,
                                                                   design = ~ time_point + sex + age + viral_load)


pca_vst <- vst(full_time_placebo_time_point_analysis, blind = FALSE)
# Clearly, we have two groups via PC2 (high and load viral load?)
plotPCA(pca_vst, intgroup = c("time_point"))
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
period_1_placebo_metadata$age <- as.factor(period_1_placebo_metadata$age)
# Run DESeq2 analysis
period_1_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_1_placebo_counts,
                              colData = period_1_placebo_metadata,
                              design = ~ time_point + sex + age + viral_load)
period_1_time_point_analysis <- DESeq(period_1_time_point_analysis, test="LRT", reduced = ~ sex + age + viral_load)
period_1_time_point_analysis_results <- results(period_1_time_point_analysis, alpha = 0.05)
period_1_time_point_analysis_results <- period_1_time_point_analysis_results[order(period_1_time_point_analysis_results$padj),]
period_1_time_point_analysis_results <- subset(period_1_time_point_analysis_results, padj < 0.05)
# PERIOD 2
# Only keep period 2 metadata and placebo_counts
period_2_placebo_metadata <- full_time_placebo_metadata[full_time_placebo_metadata$period == 2,]
period_2_placebo_counts <- placebo_counts[rownames(period_2_placebo_metadata)]
# Factorize time point (with associated factor levels) and sex
period_2_placebo_metadata$time_point <- as.factor(period_2_placebo_metadata$time_point)
levels(period_2_placebo_metadata$time_point) <- period_2_factors
period_2_placebo_metadata$sex <- as.factor(period_2_placebo_metadata$sex)
period_2_placebo_metadata$age <- as.factor(period_2_placebo_metadata$age)
# Run DESeq2 analysis
period_2_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_2_placebo_counts,
                                                       colData = period_2_placebo_metadata,
                                                       design = ~ time_point + sex + age + viral_load)
period_2_time_point_analysis <- DESeq(period_2_time_point_analysis, test="LRT", reduced = ~ sex + age + viral_load)
period_2_time_point_analysis_results <- results(period_2_time_point_analysis, alpha = 0.05)
period_2_time_point_analysis_results <- period_2_time_point_analysis_results[order(period_2_time_point_analysis_results$padj),]
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
period_2_without_2_D_minus_2_placebo_metadata$age <- as.factor(period_2_without_2_D_minus_2_placebo_metadata$age)
# Run DESeq2 analysis
period_2_without_2_D_minus_2_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_2_without_2_D_minus_2_placebo_counts,
                                                       colData = period_2_without_2_D_minus_2_placebo_metadata,
                                                       design = ~ time_point + sex + age)
period_2_without_2_D_minus_2_time_point_analysis <- DESeq(period_2_without_2_D_minus_2_time_point_analysis, test="LRT", reduced = ~ sex + age)
period_2_without_2_D_minus_2_time_point_analysis_results <- results(period_2_without_2_D_minus_2_time_point_analysis, alpha = 0.05)
period_2_without_2_D_minus_2_time_point_analysis_results <- period_2_without_2_D_minus_2_time_point_analysis_results[order(period_2_without_2_D_minus_2_time_point_analysis_results$padj),]
period_2_without_2_D_minus_2_time_point_analysis_results <- subset(period_2_without_2_D_minus_2_time_point_analysis_results, padj < 0.05)
# Run Wald tests
period_2_without_2_D_minus_2_placebo_tests <- c()
period_2_without_2_D_minus_2_placebo_wald_test_names <- resultsNames(period_2_without_2_D_minus_2_time_point_analysis)
period_2_without_2_D_minus_2_placebo_wald_test_names <- period_2_without_2_D_minus_2_placebo_wald_test_names[c(-1, -length(period_2_without_2_D_minus_2_placebo_wald_test_names))]
index <- 1
for (current_name in period_2_without_2_D_minus_2_placebo_wald_test_names) {
  current_results <- results(period_2_without_2_D_minus_2_time_point_analysis, name = current_name, test = "Wald", 
                             alpha = 0.05, lfcThreshold = 1)
  current_results <- current_results[order(current_results$padj),]
  current_results <- subset(current_results, padj < 0.05)
  period_2_without_2_D_minus_2_placebo_tests[[index]] <- current_results
  index <- index + 1
}

# Create heatmap
period_2_without_2_D_minus_2_time_point_analysis_betas <- coef(period_2_without_2_D_minus_2_time_point_analysis)
period_2_without_2_D_minus_2_time_point_analysis_results_for_plotting <- results(period_2_without_2_D_minus_2_time_point_analysis, alpha = 0.05)
period_2_without_2_D_minus_2_time_point_analysis_20_topGenes <- head(order(period_2_without_2_D_minus_2_time_point_analysis_results_for_plotting$padj),20)
period_2_without_2_D_minus_2_time_point_analysis_20_mat <- period_2_without_2_D_minus_2_time_point_analysis_betas[period_2_without_2_D_minus_2_time_point_analysis_20_topGenes, -c(1, 6, 7, 8, 9)]
period_2_without_2_D_minus_2_time_point_analysis_thr <- 3 
period_2_without_2_D_minus_2_time_point_analysis_20_mat[period_2_without_2_D_minus_2_time_point_analysis_20_mat < -period_2_without_2_D_minus_2_time_point_analysis_thr] <- -period_2_without_2_D_minus_2_time_point_analysis_thr
period_2_without_2_D_minus_2_time_point_analysis_20_mat[period_2_without_2_D_minus_2_time_point_analysis_20_mat > period_2_without_2_D_minus_2_time_point_analysis_thr] <- period_2_without_2_D_minus_2_time_point_analysis_thr
colnames(period_2_without_2_D_minus_2_time_point_analysis_20_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
pheatmap(period_2_without_2_D_minus_2_time_point_analysis_20_mat, breaks=seq(from=-period_2_without_2_D_minus_2_time_point_analysis_thr, to=period_2_without_2_D_minus_2_time_point_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(data_dir, "period_2_without_2_D_minus_2_time_point_analysis_top_20_genes.png"))
period_2_without_2_D_minus_2_sig_genes <- rownames(period_2_without_2_D_minus_2_time_point_analysis_results)
period_2_without_2_D_minus_2_sig_genes_fcs <- period_2_without_2_D_minus_2_time_point_analysis_betas[period_2_without_2_D_minus_2_sig_genes, -c(1, 6, 7, 8, 9)]
# D8
period_2_without_2_D_minus_2_sig_genes_fcs_D8 <- period_2_without_2_D_minus_2_sig_genes_fcs[,3]
period_2_without_2_D_minus_2_sig_genes_fcs_D8 <- period_2_without_2_D_minus_2_sig_genes_fcs_D8[period_2_without_2_D_minus_2_sig_genes_fcs_D8 <= -1 | period_2_without_2_D_minus_2_sig_genes_fcs_D8 >= 1]
period_2_without_2_D_minus_2_time_point_analysis_D8_gene_ontology_genes <- names(period_2_without_2_D_minus_2_sig_genes_fcs_D8)
write.table(period_2_without_2_D_minus_2_time_point_analysis_D8_gene_ontology_genes, paste0(data_dir, "period_2_without_2_D_minus_2_time_point_analysis_D8_gene_ontology_genes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# D28
period_2_without_2_D_minus_2_sig_genes_fcs_D28 <- period_2_without_2_D_minus_2_sig_genes_fcs[,4]
period_2_without_2_D_minus_2_sig_genes_fcs_D28 <- period_2_without_2_D_minus_2_sig_genes_fcs_D28[period_2_without_2_D_minus_2_sig_genes_fcs_D28 <= -1 | period_2_without_2_D_minus_2_sig_genes_fcs_D28 >= 1]
period_2_without_2_D_minus_2_time_point_analysis_D28_gene_ontology_genes <- names(period_2_without_2_D_minus_2_sig_genes_fcs_D28)
write.table(period_2_without_2_D_minus_2_time_point_analysis_D28_gene_ontology_genes, paste0(data_dir, "period_2_without_2_D_minus_2_time_point_analysis_D28_gene_ontology_genes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# Plot heatmap with top 200 genes (and clean up heatmap plot by removing row names and clustering)
period_2_without_2_D_minus_2_time_point_analysis_200_topGenes <- head(order(period_2_without_2_D_minus_2_time_point_analysis_results_for_plotting$padj),200)
period_2_without_2_D_minus_2_time_point_analysis_200_mat <- period_2_without_2_D_minus_2_time_point_analysis_betas[period_2_without_2_D_minus_2_time_point_analysis_200_topGenes, -c(1, 6, 7, 8, 9)]
colnames(period_2_without_2_D_minus_2_time_point_analysis_200_mat) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
period_2_without_2_D_minus_2_time_point_analysis_200_mat[period_2_without_2_D_minus_2_time_point_analysis_200_mat < -period_2_without_2_D_minus_2_time_point_analysis_thr] <- -period_2_without_2_D_minus_2_time_point_analysis_thr
period_2_without_2_D_minus_2_time_point_analysis_200_mat[period_2_without_2_D_minus_2_time_point_analysis_200_mat > period_2_without_2_D_minus_2_time_point_analysis_thr] <- period_2_without_2_D_minus_2_time_point_analysis_thr
pheatmap(period_2_without_2_D_minus_2_time_point_analysis_200_mat, breaks=seq(from=-period_2_without_2_D_minus_2_time_point_analysis_thr, to=period_2_without_2_D_minus_2_time_point_analysis_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, show_rownames = FALSE, cluster_row = FALSE, filename = paste0(data_dir, "period_2_without_2_D_minus_2_time_point_analysis_top_200_genes.png"))


# BOTH PERIODS (ALL TIME POINTS)
# Factorize time point (with associated factor levels) and sex
full_time_placebo_metadata$time_point <- as.factor(full_time_placebo_metadata$time_point)
levels(full_time_placebo_metadata$time_point) <- all_factors
full_time_placebo_metadata$sex <- as.factor(full_time_placebo_metadata$sex)
# Run DESeq2 analysis
full_time_point_analysis <- DESeqDataSetFromMatrix(countData = full_time_placebo_counts,
                                                       colData = full_time_placebo_metadata,
                                                       design = ~ time_point + sex + age)
full_time_point_analysis <- DESeq(full_time_point_analysis, test="LRT", reduced= ~ sex + age)
full_time_point_analysis_results <- results(full_time_point_analysis, alpha = 0.05)
full_time_point_analysis_results <- full_time_point_analysis_results[order(full_time_point_analysis_results$padj),]
full_time_point_analysis_results <- subset(full_time_point_analysis_results, padj < 0.05)
# WALD (PAIRWISE) TESTS
wald_full_time_point_analysis <- DESeqDataSetFromMatrix(countData = full_time_placebo_counts,
                                              colData = full_time_placebo_metadata,
                                              design = ~ time_point + sex + age)
wald_full_time_point_analysis <- DESeq(wald_full_time_point_analysis)
full_period_1_baseline <- "1_D_minus_1"
full_period_2_baseline <- "2_D_minus_2"
full_period_1_contrasts <- c("1_D2", "1_D8", "1_D28")
full_period_2_contrasts <- c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
full_period_1_wald_tests <- c()
full_period_2_wald_tests <- c()
# Period 1 Wald tests
index <- 1
for (current_contrast in full_period_1_contrasts) {
  current_results <- results(wald_full_time_point_analysis, contrast = c("time_point", current_contrast, full_period_1_baseline), 
                                                                    alpha = 0.05)
  current_results <- current_results[order(current_results$padj),]
  current_results <- subset(current_results, padj < 0.05)
  full_period_1_wald_tests[[index]] <- current_results
  index <- index + 1
}
# Period 2 Wald tests
index <- 1
for (current_contrast in full_period_2_contrasts) {
  current_results <- results(wald_full_time_point_analysis, contrast = c("time_point", current_contrast, full_period_2_baseline), 
                             alpha = 0.05)
  current_results <- current_results[order(current_results$padj),]
  current_results <- subset(current_results, padj < 0.05)
  full_period_2_wald_tests[[index]] <- current_results
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
period_2_more_metadata$age <- as.factor(period_2_more_metadata$age)
# Run DESeq2 analysis
period_2_more_data_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_2_more_placebo_counts,
                                                       colData = period_2_more_metadata,
                                                       design = ~ time_point + sex + age)
period_2_more_data_time_point_analysis <- DESeq(period_2_more_data_time_point_analysis, test="LRT", reduced= ~ sex + age)
period_2_more_data_time_point_analysis_results <- results(period_2_more_data_time_point_analysis, alpha = 0.05)
period_2_more_data_time_point_analysis_results <- period_2_more_data_time_point_analysis_results[order(period_2_more_data_time_point_analysis_results$padj),]
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
period_2_subset_placebo_metadata$age <- as.factor(period_2_subset_placebo_metadata$age)
# Run DESeq2 analysis
period_2_subset_time_point_analysis <- DESeqDataSetFromMatrix(countData = period_2_subset_placebo_counts,
                                                       colData = period_2_subset_placebo_metadata,
                                                       design = ~ time_point + sex + age)
period_2_subset_time_point_analysis <- DESeq(period_2_subset_time_point_analysis, test="LRT", reduced= ~ sex + age)
period_2_subset_time_point_analysis_results <- results(period_2_subset_time_point_analysis, alpha = 0.05)
period_2_subset_time_point_analysis_results <- period_2_subset_time_point_analysis_results[order(period_2_subset_time_point_analysis_results$padj),]
period_2_subset_time_point_analysis_results <- subset(period_2_subset_time_point_analysis_results, padj < 0.05)

save.image(paste0(data_dir, "placebo_bulk_RNA_obj.RData"))
