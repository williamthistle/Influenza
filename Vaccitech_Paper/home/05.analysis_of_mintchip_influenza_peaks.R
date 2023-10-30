# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

run_fmd_on_mintchip <- function(mintchip_marker_das_list_annotated) {
  fmd_results <- list()
  index_1 <- 1
  for(marker in mintchip_marker_das_list_annotated) {
    fmd_marker_results <- list()
    index_2 <- 1
    for(current_list_of_peaks in marker) {
      if(length(current_list_of_peaks) > 1) {
        current_fmd_result <- SPEEDI::RunFMD_RNA(unique(current_list_of_peaks$SYMBOL), "blood")
      } else { 
        current_fmd_result <- "EMPTY"
      }
      fmd_marker_results[[index_2]] <- current_fmd_result
      index_2 <- index_2 + 1
    }
    fmd_results[[index_1]] <- fmd_marker_results
    index_1 <- index_1 + 1
  }
  return(fmd_results)
}


mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")
mintchip_marker_das_list <- list()

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

for(marker in mintchip_markers) {
  differential_analysis_results_0_filtered <- read.table(paste0(mintchip_das_dir, marker, "_D28_D1_diff_filtered.tsv"), sep = "\t", header = TRUE)
  differential_analysis_results_0.1_filtered <- read.table(paste0(mintchip_das_dir, marker, "_D28_D1_diff_filtered_0.1.tsv"), sep = "\t", header = TRUE)
  differential_analysis_results_0.2_filtered <- read.table(paste0(mintchip_das_dir, marker, "_D28_D1_diff_filtered_0.2.tsv"), sep = "\t", header = TRUE)
  differential_analysis_results_0.3_filtered <- read.table(paste0(mintchip_das_dir, marker, "_D28_D1_diff_filtered_0.3.tsv"), sep = "\t", header = TRUE)
  differential_analysis_results_0.585_filtered <- read.table(paste0(mintchip_das_dir, marker, "_D28_D1_diff_filtered_0.585.tsv"), sep = "\t", header = TRUE)
  differential_analysis_results_1_filtered <- read.table(paste0(mintchip_das_dir, marker, "_D28_D1_diff_filtered_1.tsv"), sep = "\t", header = TRUE)
  differential_analysis_results_2_filtered <- read.table(paste0(mintchip_das_dir, marker, "_D28_D1_diff_filtered_2.tsv"), sep = "\t", header = TRUE)
  mintchip_marker_das_list[[marker]] <- list(differential_analysis_results_0_filtered, differential_analysis_results_0.1_filtered, 
                                             differential_analysis_results_0.2_filtered, differential_analysis_results_0.3_filtered,
                                             differential_analysis_results_0.585_filtered, differential_analysis_results_1_filtered, 
                                        differential_analysis_results_2_filtered)
}

pos_mintchip_marker_das_list_annotated <- list()
neg_mintchip_marker_das_list_annotated <- list()

for(marker in mintchip_markers) {
  old_das <- mintchip_marker_das_list[[marker]]
  current_pos_annotated_das <- list()
  current_neg_annotated_das <- list()
  index <- 1
  for(current_das in old_das) {
    if(nrow(current_das) > 0) {
      pos_das <- current_das[current_das$log2FoldChange > 0,]
      if(nrow(pos_das) > 0) {
        pos_das <- pos_das[,c(1,2,3,4,6)]
        pos_das$length <- pos_das$end - pos_das$start
        pos_das <- annotatePeak(makeGRangesFromDataFrame(pos_das), TxDb = txdb, annoDb = "org.Hs.eg.db")
        pos_das <- as.data.frame(pos_das)
        current_pos_annotated_das[[index]] <- pos_das
      } else {
        current_pos_annotated_das[[index]] <- "EMPTY"
      }
      neg_das <- current_das[current_das$log2FoldChange < 0,]
      if(nrow(neg_das) > 0) {
        neg_das <- neg_das[,c(1,2,3,4,6)]
        neg_das$length <- neg_das$end - neg_das$start
        neg_das <- annotatePeak(makeGRangesFromDataFrame(neg_das), TxDb = txdb, annoDb = "org.Hs.eg.db")
        neg_das <- as.data.frame(neg_das)
        current_neg_annotated_das[[index]] <- neg_das
      } else {
        current_neg_annotated_das[[index]] <- "EMPTY"
      }
    } else {
      current_pos_annotated_das[[index]] <- "EMPTY"
      current_neg_annotated_das[[index]] <- "EMPTY"
    }
    index <- index + 1
  }
  pos_mintchip_marker_das_list_annotated[[marker]] <- current_pos_annotated_das
  neg_mintchip_marker_das_list_annotated[[marker]] <- current_neg_annotated_das
}

pos_fmd <- run_fmd_on_mintchip(pos_mintchip_marker_das_list_annotated)
neg_fmd <- run_fmd_on_mintchip(neg_mintchip_marker_das_list_annotated)
  


