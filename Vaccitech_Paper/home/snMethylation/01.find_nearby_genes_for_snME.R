# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))
magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

run_fmd_on_snME <- function(gene_list) {
  if(length(gene_list) > 1 & length(gene_list) < 2000) {
    current_fmd_result <- SPEEDI::RunFMD_RNA(gene_list, "blood")
  } else { 
    current_fmd_result <- "EMPTY OR OVER 2000 PEAKS (TOO MANY)"
  }
  return(current_fmd_result)
}

# Find snME within promoter regions of genes for each cell type
snME_cell_types <- c("B-Mem", "B-Naive", "Monocyte", "NK-cell2", "Tc-Mem", "Tc-Naive", "Th-Mem", "Th-Naive")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate genes for methylated promoters
for(snME_cell_type in snME_cell_types) {
  snME_dms_for_current_cell_type <- snME_dms[snME_dms$celltype == snME_cell_type,]
  upregulated_snME_dms_for_current_cell_type <- snME_dms_for_current_cell_type[snME_dms_for_current_cell_type$methylation == "upregulated",]
  downregulated_snME_dms_for_current_cell_type <- snME_dms_for_current_cell_type[snME_dms_for_current_cell_type$methylation == "downregulated",]
  upregulated_differential_analysis_results <- upregulated_snME_dms_for_current_cell_type[,c(2,3,4)]
  upregulated_differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(upregulated_differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
  upregulated_differential_analysis_results_df <- as.data.frame(upregulated_differential_analysis_results)
  write.table(upregulated_differential_analysis_results_df, file = paste0(snME_results_dir, "Closest_Gene_Analysis_No_Promoter_Subset/", snME_cell_type, "_D28_hypermethylated_annotated_genes.tsv"), 
              sep = "\t", quote = FALSE)
  upregulated_differential_analysis_results_df <- subset(upregulated_differential_analysis_results_df, distanceToTSS >= -2000 & distanceToTSS <= 500)
  write.table(upregulated_differential_analysis_results_df, file = paste0(snME_results_dir, "Closest_Gene_Analysis/", snME_cell_type, "_D28_hypermethylated_annotated_genes.tsv"), 
              sep = "\t", quote = FALSE)
  
  downregulated_differential_analysis_results <- downregulated_snME_dms_for_current_cell_type[,c(2,3,4)]
  downregulated_differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(downregulated_differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
  downregulated_differential_analysis_results_df <- as.data.frame(downregulated_differential_analysis_results)
  write.table(downregulated_differential_analysis_results_df, file = paste0(snME_results_dir, "Closest_Gene_Analysis_No_Promoter_Subset/", snME_cell_type, "_D_minus_1_hypermethylated_annotated_genes.tsv"), 
              sep = "\t", quote = FALSE)
  downregulated_differential_analysis_results_df <- subset(downregulated_differential_analysis_results_df, distanceToTSS >= -2000 & distanceToTSS <= 500)
  write.table(downregulated_differential_analysis_results_df, file = paste0(snME_results_dir, "Closest_Gene_Analysis/", snME_cell_type, "_D_minus_1_hypermethylated_annotated_genes.tsv"), 
              sep = "\t", quote = FALSE)
}

# Run HB on positive (hypermethylated in D28) genes
# Note - I manually subset me3 to TSS (-3k to 3k), so I should add code to do that automatically below
pos_fmd_snME_list <- list()
for(snME_cell_type in snME_cell_types) {
  pos_snME_results_file <- paste0(snME_results_dir, "Closest_Gene_Analysis/", snME_cell_type, "_D28_hypermethylated_annotated_genes.tsv")
  pos_snME_cell_type_annotations <- read.table(pos_snME_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  pos_genes <- unique(pos_snME_cell_type_annotations$SYMBOL)
  pos_results <- run_fmd_on_mintchip(pos_genes)
  pos_fmd_snME_list[[snME_cell_type]] <- pos_results
}

neg_fmd_snME_list <- list()
for(snME_cell_type in snME_cell_types) {
  neg_snME_results_file <- paste0(snME_results_dir, "Closest_Gene_Analysis/", snME_cell_type, "_D_minus_1_hypermethylated_annotated_genes.tsv")
  neg_snME_cell_type_annotations <- read.table(neg_snME_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  neg_genes <- unique(neg_snME_cell_type_annotations$SYMBOL)
  neg_results <- run_fmd_on_mintchip(neg_genes)
  neg_fmd_snME_list[[snME_cell_type]] <- neg_results
}

# Create annotation plot
peak_annotation_plots <- list()
innate_snME_cell_types <- c("Monocyte", "NK-cell2")
for(snME_cell_type in innate_snME_cell_types) {
  # Upregulated
  up_differential_analysis_results_file <- paste0(snME_results_dir, "Closest_Gene_Analysis_No_Promoter_Subset/", snME_cell_type, "_D28_hypermethylated_annotated_genes.tsv")
  up_differential_analysis_results <- read.table(up_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  # Downregulated
  down_differential_analysis_results_file <- paste0(snME_results_dir, "Closest_Gene_Analysis_No_Promoter_Subset/", snME_cell_type, "_D_minus_1_hypermethylated_annotated_genes.tsv")
  down_differential_analysis_results <- read.table(down_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  # Combine
  differential_analysis_results <- rbind(up_differential_analysis_results, down_differential_analysis_results)
  differential_analysis_results <- differential_analysis_results[,c(1,2,3,4,5)]
  differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
  if(snME_cell_type == "NK-cell2") {
    snME_cell_type <- "NK"
  }
  peak_annotation_plots[[snME_cell_type]] <- differential_analysis_results
}

plotAnnoBar(peak_annotation_plots, ylab = "Percentage", title = "Distribution of Genomic Features for snME DMR")



