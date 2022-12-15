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
# Sort columns in counts and rows for each so they're in same order (for DESeq2)
colnames(counts) <- sort(colnames(counts))
rownames(all_metadata) <- all_metadata$aliquot_id
rownames(all_metadata) <- sort(rownames(all_metadata))
# Drop aliquot ID column (it's stored in rownames)
all_metadata <- subset(all_metadata, select = -c(aliquot_id))
# Probably OK to round expected counts from RSEM data. DESeq2 expects integers
counts <- round(counts)
# Currently not filtering sex associated genes
#sex_associated_genes <- find_sex_associated_genes(paste0(data_dir, "sex_associated_genes/"))
# Order factor levels for period 1, period 2, and all time points
period_1_factors <- c("1_D_minus_1", "1_D2", "1_D8", "1_D28")
period_1_more_vaccination_data_factors <- c("1_D_minus_1", "1_D8")
period_2_factors <- c("2_D_minus_2", "2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
period_2_without_2_D_minus_2_factors <- c("2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
all_factors <- c(period_1_factors, period_2_factors)

# VACCINATED + PLACEBO TIME SERIES ANALYSIS - is there interaction between vaccination status and time point?
full_time_all_metadata <- all_metadata[all_metadata$subject_id 
                                           %in% names(table(all_metadata$subject_id)
                                                      [table(all_metadata$subject_id) == 10]),]
full_time_all_counts <- counts[rownames(full_time_all_metadata)]
full_time_all_metadata$time_point <- as.factor(full_time_all_metadata$time_point)
levels(full_time_all_metadata$time_point) <- all_factors
full_time_all_metadata$sex <- as.factor(full_time_all_metadata$sex)
full_time_all_metadata$treatment <- as.factor(full_time_all_metadata$treatment)
# Run DESeq2 analysis
all_time_point_analysis <- DESeqDataSetFromMatrix(countData = full_time_all_counts,
                                              colData = full_time_all_metadata,
                                              design = ~ time_point + sex + treatment + time_point:treatment)
all_time_point_analysis <- DESeq(all_time_point_analysis, test="LRT", reduced= ~ time_point + sex + treatment)
all_time_point_analysis_results <- results(all_time_point_analysis, alpha = 0.05)
all_time_point_analysis_results <- all_time_point_analysis_results[order(all_time_point_analysis_results$padj),]
all_time_point_analysis_results <- subset(all_time_point_analysis_results, padj < 0.05)
all_betas <- coef(all_time_point_analysis)
all_time_point_analysis_results_for_plotting <- results(all_time_point_analysis, alpha = 0.05)
all_topGenes <- head(order(all_time_point_analysis_results_for_plotting$padj),20)
all_mat <- all_betas[all_topGenes, -c(1, 11, 12)]
all_thr <- 3 
all_mat[all_mat < -all_thr] <- -all_thr
all_mat[all_mat > all_thr] <- all_thr
colnames(all_mat) <- c("Pl Day 2 (P1)", "Pl Day 8 (P1)", "Pl Day 28 (P1)", "Pl Day -2 (P2)", "Pl Day -1 (P2)", 
                       "Pl Day 2 (P2)", "Pl Day 5 (P2)", "Pl Day 8 (P2)", "Pl Day 28 (P2)", 
                       "V Day 2 (P1)", "V Day 8 (P1)", "V Day 28 (P1)", "V Day -2 (P2)", "V Day -1 (P2)", 
                       "V Day 2 (P2)", "V Day 5 (P2)", "V Day 8 (P2)", "V Day 28 (P2)")
pheatmap(all_mat, breaks=seq(from=-all_thr, to=all_thr, length=101),
         cluster_col=FALSE, fontsize_col=14, filename = paste0(data_dir, "all_time_point_analysis_top_20_genes.png"))
# Shuffle treatment to see whether we get the same number of DEGs with randomly assigned vaccinated and placebo subjects
shuffled_metadata <- full_time_all_metadata
treatment_shuffled_vector <- shuffled_metadata$treatment
treatment_shuffled_vector <- table(treatment_shuffled_vector)
treatment_shuffled_vector <- c(rep("MVA-NP+M1", treatment_shuffled_vector["MVA-NP+M1"] / 10), rep("PLACEBO", treatment_shuffled_vector["PLACEBO"] / 10))
treatment_shuffled_vector <- sample(treatment_shuffled_vector)
treatment_shuffled_vector <- rep(treatment_shuffled_vector, each = 10)
shuffled_metadata$treatment <- treatment_shuffled_vector
shuffled_time_point_analysis <- DESeqDataSetFromMatrix(countData = full_time_all_counts,
                                                  colData = shuffled_metadata,
                                                  design = ~ time_point + sex + treatment + time_point:treatment)
shuffled_time_point_analysis <- DESeq(shuffled_time_point_analysis, test="LRT", reduced= ~ time_point + sex + treatment)
shuffled_time_point_analysis_results <- results(shuffled_time_point_analysis, alpha = 0.05)
shuffled_time_point_analysis_results <- shuffled_time_point_analysis_results[order(shuffled_time_point_analysis_results$padj),]
shuffled_time_point_analysis_results <- subset(shuffled_time_point_analysis_results, padj < 0.05)

save.image(paste0(data_dir, "both_bulk_RNA_obj.RData"))
