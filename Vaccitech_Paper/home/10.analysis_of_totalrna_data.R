# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

rownames(totalRNA_raw_counts) <- totalRNA_raw_counts$target_id 
totalRNA_raw_counts <- totalRNA_raw_counts[,-c(1)]

sorted_col_names <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13",
                      "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26")
totalRNA_raw_counts <- totalRNA_raw_counts[, sorted_col_names]
colnames(totalRNA_raw_counts) <- totalRNA_metadata_table$Aliquot

associated_subjects <- c()
associated_time_points <- c()
for(current_aliquot in colnames(totalRNA_raw_counts)) {
  associated_subject <- as.character(all_metadata[all_metadata$aliquot_id == current_aliquot,]$subject_id)
  associated_time_point <- as.character(all_metadata[all_metadata$aliquot_id == current_aliquot,]$time_point)
  associated_subjects <- c(associated_subjects, associated_subject)
  associated_time_points <- c(associated_time_points, associated_time_point)
}

totalRNA_metadata_table$Subject <- associated_subjects
totalRNA_metadata_table$time_point <- associated_time_points
rownames(totalRNA_metadata_table) <- totalRNA_metadata_table$Aliquot

# Run DESeq2
current_analysis <- DESeqDataSetFromMatrix(countData = totalRNA_raw_counts, colData = totalRNA_metadata_table, design = ~ Subject + time_point)
current_analysis <- DESeq(current_analysis)
# 24 genes with lfcThreshold = 0
current_analysis_results_0 <- results(current_analysis, contrast = c("time_point", "2_D28", "2_D_minus_1"), alpha = 0.05)
current_analysis_results_0 <- current_analysis_results_0[order(current_analysis_results_0$padj),]
current_analysis_results_0 <- subset(current_analysis_results_0, padj < 0.05)
# 0 genes with lfcThreshold = 0.1
current_analysis_results_0.1 <- results(current_analysis, contrast = c("time_point", "2_D28", "2_D_minus_1"), alpha = 0.05, lfcThreshold = 0.1)
current_analysis_results_0.1 <- current_analysis_results_0.1[order(current_analysis_results_0.1$padj),]
current_analysis_results_0.1 <- subset(current_analysis_results_0.1, padj < 0.05)

# Loaded genemap from Wayne's normalization script
current_genes <- rownames(current_analysis_results_0)
iim=match(current_genes, genemap$ensembl_gene_id)
current_gene_symbols <- genemap[iim,  "hgnc_symbol", drop=F]
current_gene_symbols$original_id <- current_genes

# Don't think there's much here to analyze
