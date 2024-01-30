# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))
magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

run_fmd_on_snATAC <- function(gene_list) {
  if(length(gene_list) > 1 & length(gene_list) < 2000) {
    current_fmd_result <- SPEEDI::RunFMD_RNA(gene_list, "blood")
  } else { 
    current_fmd_result <- "EMPTY OR OVER 2000 PEAKS (TOO MANY)"
  }
  return(current_fmd_result)
}

# Find snME within promoter regions of genes for each cell type
snATAC_cell_types <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "T Naive", "Proliferating")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Create annotated up and downregulated peak files (annotated)
for(snATAC_cell_type in snATAC_cell_types) {
  differential_analysis_results_file <- paste0(mintchip_das_dir, marker, "/", marker, "_", peak_set, "_FC_", fc, ".tsv")
      if(file.size(differential_analysis_results_file) != 1 && file.size(differential_analysis_results_file) != 75) {
        differential_analysis_results <- read.table(differential_analysis_results_file, sep = "\t", header = TRUE)
        upregulated_differential_analysis_results <- differential_analysis_results[differential_analysis_results$Fold > 0,]
        downregulated_differential_analysis_results <- differential_analysis_results[differential_analysis_results$Fold < 0,]
        upregulated_dir <- paste0(marker_dir, "upregulated/")
        downregulated_dir <- paste0(marker_dir, "downregulated/")
        all_dir <- paste0(marker_dir, "all/")
        if (!dir.exists(upregulated_dir)) {dir.create(upregulated_dir)}
        if (!dir.exists(downregulated_dir)) {dir.create(downregulated_dir)}
        if (!dir.exists(all_dir)) {dir.create(all_dir)}
        write.table(differential_analysis_results, file = paste0(all_dir, marker, "_", peak_set, "_FC_", fc, ".tsv"), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
        if(nrow(upregulated_differential_analysis_results) > 0) {
          write.table(upregulated_differential_analysis_results, file = paste0(upregulated_dir, marker, "_", peak_set, "_FC_", fc, "_upregulated.tsv"), 
                      sep = "\t", quote = FALSE, row.names = FALSE)
          upregulated_differential_analysis_results <- upregulated_differential_analysis_results[,c(1,2,3,4,5)]
          upregulated_differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(upregulated_differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
          upregulated_differential_analysis_results <- as.data.frame(upregulated_differential_analysis_results)
          write.table(upregulated_differential_analysis_results, file = paste0(upregulated_dir, marker, "_", peak_set, "_FC_", fc, "_upregulated_annotated.tsv"), 
                      sep = "\t", quote = FALSE, row.names = FALSE)
        }
        if(nrow(downregulated_differential_analysis_results) > 0) {
          write.table(downregulated_differential_analysis_results, file = paste0(downregulated_dir, marker, "_", peak_set, "_FC_", fc, "_downregulated.tsv"), 
                      sep = "\t", quote = FALSE, row.names = FALSE)
          downregulated_differential_analysis_results <- downregulated_differential_analysis_results[,c(1,2,3,4,5)]
          downregulated_differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(downregulated_differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
          downregulated_differential_analysis_results <- as.data.frame(downregulated_differential_analysis_results)
          write.table(downregulated_differential_analysis_results, file = paste0(downregulated_dir, marker, "_", peak_set, "_FC_", fc, "_downregulated_annotated.tsv"), 
                      sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }
    }
  }
}


for(fc in c(0, 0.1, 0.2, 0.3, 0.585, 1, 2)) {



# Run HB on positive (hypermethylated in D28) genes
pos_fmd_snATAC_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  pos_snATAC_results_file <- paste0(snATAC_results_dir, "Closest_Gene_Analysis/", snATAC_cell_type, "_D28_hypermethylated_annotated_genes.tsv")
  pos_snATAC_cell_type_annotations <- read.table(pos_snATAC_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  pos_genes <- unique(pos_snATAC_cell_type_annotations$SYMBOL)
  pos_results <- run_fmd_on_mintchip(pos_genes)
  pos_fmd_snATAC_list[[snATAC_cell_type]] <- pos_results
}

neg_fmd_snATAC_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  neg_snATAC_results_file <- paste0(snATAC_results_dir, "Closest_Gene_Analysis/", snATAC_cell_type, "_D_minus_1_hypermethylated_annotated_genes.tsv")
  neg_snATAC_cell_type_annotations <- read.table(neg_snATAC_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  neg_genes <- unique(neg_snATAC_cell_type_annotations$SYMBOL)
  neg_results <- run_fmd_on_mintchip(neg_genes)
  neg_fmd_snATAC_list[[snATAC_cell_type]] <- neg_results
}
