# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

rownames(miRNA_raw_counts) <- miRNA_raw_counts$miRNA.precursor 
miRNA_raw_counts <- miRNA_raw_counts[,-c(1)]

sorted_col_names <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13",
                      "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26")
miRNA_raw_counts <- miRNA_raw_counts[, sorted_col_names]
colnames(miRNA_raw_counts) <- miRNA_metadata_table$Aliquot

associated_subjects <- c()
associated_time_points <- c()
for(current_aliquot in colnames(miRNA_raw_counts)) {
  associated_subject <- as.character(all_metadata[all_metadata$aliquot_id == current_aliquot,]$subject_id)
  associated_time_point <- as.character(all_metadata[all_metadata$aliquot_id == current_aliquot,]$time_point)
  associated_subjects <- c(associated_subjects, associated_subject)
  associated_time_points <- c(associated_time_points, associated_time_point)
}

miRNA_metadata_table$Subject <- associated_subjects
miRNA_metadata_table$time_point <- associated_time_points
rownames(miRNA_metadata_table) <- miRNA_metadata_table$Aliquot

# Run DESeq2
current_analysis <- DESeqDataSetFromMatrix(countData = miRNA_raw_counts, colData = miRNA_metadata_table, design = ~ Subject + time_point)
current_analysis <- DESeq(current_analysis)
current_analysis_results <- results(current_analysis, contrast = c("time_point", "2_D28", "2_D_minus_1"), alpha = 0.05)
current_analysis_results <- current_analysis_results[order(current_analysis_results$padj),]
current_analysis_results <- subset(current_analysis_results, padj < 0.05)

# Nothing here to analyze - we can maybe do a PC plot?