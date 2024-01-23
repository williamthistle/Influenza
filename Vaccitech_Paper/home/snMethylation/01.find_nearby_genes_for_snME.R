# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))
magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

snME_cell_types <- c("B-Mem", "B-Naive", "Monocyte", "NK-cell2", "Tc-Mem", "Tc-Naive", "Th-Mem", "Th-Naive")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

for(snME_cell_type in snME_cell_types) {
  snME_dms_for_current_cell_type <- snME_dms[snME_dms$celltype == snME_cell_type,]
  upregulated_snME_dms_for_current_cell_type <- snME_dms_for_current_cell_type[snME_dms_for_current_cell_type$methylation == "upregulated",]
  downregulated_snME_dms_for_current_cell_type <- snME_dms_for_current_cell_type[snME_dms_for_current_cell_type$methylation == "downregulated",]
  upregulated_differential_analysis_results <- upregulated_snME_dms_for_current_cell_type[,c(2,3,4)]
  upregulated_differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(upregulated_differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
  upregulated_differential_analysis_results_df <- as.data.frame(upregulated_differential_analysis_results)
  upregulated_differential_analysis_results_df <- subset(upregulated_differential_analysis_results_df, distanceToTSS >= -2000 & distanceToTSS <= 500)
  write.table(upregulated_differential_analysis_results_df, file = paste0(upregulated_dir, marker, "_", peak_set, "_FC_", fc, "_upregulated_annotated.tsv"), 
              sep = "\t", quote = FALSE)

  downregulated_differential_analysis_results <- downregulated_snME_dms_for_current_cell_type[,c(2,3,4)]
  downregulated_differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(downregulated_differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
  downregulated_differential_analysis_results_df <- as.data.frame(downregulated_differential_analysis_results)
  downregulated_differential_analysis_results_df <- subset(downregulated_differential_analysis_results_df, distanceToTSS >= -2000 & distanceToTSS <= 500)
  write.table(downregulated_differential_analysis_results_df, file = paste0(downregulated_dir, marker, "_", peak_set, "_FC_", fc, "_downregulated_annotated.tsv"), 
              sep = "\t", quote = FALSE)
}