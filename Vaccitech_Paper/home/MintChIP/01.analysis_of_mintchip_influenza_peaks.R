### THIS FUNCTION:
### 1) RUNS CHIPSEEKER ON EACH SET OF UPREGULATED AND DOWNREGULATED MINTCHIP PEAKS
### (WITH DIFFERENT FOLDCHANGE THRESHOLDS - 0, 0.1, 0.2, 0.3, 0.585, 1, 2)
### I CAN LOOK AT FUNCTION OF INDIVIDUAL HIGH EFFECT SIZE GENES
### 2) RUNS FMD ON EACH SET OF PEAKS TO GET HUMANBASE PLOTS
### THESE PLOTS COULD BE USED IN MAIN FIGURES OR SUPPLEMENTS
### I CAN ALSO FIND PATHWAY ENRICHMENT IN MODULES
### OR SHOULD I AVOID USING HB AND JUST HAVE PATHWAY ENRICHMENT FOR EACH MARKER?
### LESS GRANULAR AND SPECIFIC BUT WOULD BE WAY MORE CONDENSED
### OTHERWISE, I AM LOOKING AT 4 PLOTS FOR EACH MARKER (HB UP, REACTOME UP,
### HB DOWN, REACTOME DOWN)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

run_fmd_on_mintchip <- function(gene_list) {
  if(length(gene_list) > 1 & length(gene_list) < 2000) {
    current_fmd_result <- SPEEDI::RunFMD_RNA(gene_list, "blood")
  } else { 
    current_fmd_result <- "EMPTY OR OVER 2000 PEAKS (TOO MANY)"
  }
  return(current_fmd_result)
}


mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

promoter_terms <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")

# Create annotated up and downregulated peak files (annotated)
for(marker in mintchip_markers) {
  marker_dir <- paste0(mintchip_das_dir, marker, "/")
  for(peak_set in c("DESeq2", "edgeR", "consensus_peak_set")) {
    for(fc in c(0, 0.1, 0.2, 0.3, 0.585, 1, 2)) {
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

# Run FMD (pos and neg)
# TODO: Try enrichr to see if it works better?
pos_fmd_list <- list()
for(marker in mintchip_markers) {
  marker_dir <- paste0(mintchip_das_dir, marker, "/")
  pos_fmd_list[[marker]] <- list()
  for(fc in c(0, 0.1, 0.2, 0.3, 0.585, 1, 2)) {
    pos_differential_analysis_results_file <- paste0(mintchip_das_dir, marker, "/upregulated/", marker, "_DESeq2_FC_", fc, "_upregulated_annotated.tsv")
    if(file.exists(pos_differential_analysis_results_file) && file.size(pos_differential_analysis_results_file) != 1 && file.size(pos_differential_analysis_results_file) != 75) {
      pos_differential_analysis_results_file <- read.table(pos_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      # pos_differential_analysis_results_file <- subset(pos_differential_analysis_results_file, annotation %in% promoter_terms)
      pos_genes <- unique(pos_differential_analysis_results_file$SYMBOL)
      pos_results <- run_fmd_on_mintchip(pos_genes)
      pos_fmd_list[[marker]][[as.character(fc)]] <- pos_results
    }
  }
}

# saveRDS(pos_fmd_list, file = paste0(mintchip_fmd_dir, "pos_fmd_V2.RDS"))
pos_fmd_list <- readRDS(paste0(mintchip_fmd_dir, "pos_fmd_V2.RDS"))


neg_fmd_list <- list()
for(marker in mintchip_markers) {
  marker_dir <- paste0(mintchip_das_dir, marker, "/")
  neg_fmd_list[[marker]] <- list()
  for(fc in c(0, 0.1, 0.2, 0.3, 0.585, 1, 2)) {
    neg_differential_analysis_results_file <- paste0(mintchip_das_dir, marker, "/downregulated/", marker, "_DESeq2_FC_", fc, "_downregulated_annotated.tsv")
    if(file.exists(neg_differential_analysis_results_file) && file.size(neg_differential_analysis_results_file) != 1 && file.size(neg_differential_analysis_results_file) != 75) {
      neg_differential_analysis_results_file <- read.table(neg_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      neg_genes <- unique(neg_differential_analysis_results_file$SYMBOL)
      neg_results <- run_fmd_on_mintchip(neg_genes)
      neg_fmd_list[[marker]][[as.character(fc)]] <- neg_results
    }
  }
}

# saveRDS(neg_fmd_list, file = paste0(mintchip_fmd_dir, "neg_fmd_V2.RDS"))
neg_fmd_list <- readRDS(paste0(mintchip_fmd_dir, "neg_fmd_V2.RDS"))

# Create peak annotation plot
mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")
mintchip_markers_alt <- c("H3K27me3", "H3K27Ac", "H3K4me3", "H3K36me3", "H3K9me3", "H3K4me1")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peak_annotation_plots <- list()
for(marker in mintchip_markers) {
  differential_analysis_results_file <- paste0(mintchip_das_dir, marker, "/all/", marker, "_DESeq2_FC_0.1.tsv")
  differential_analysis_results <- read.table(differential_analysis_results_file, sep = "\t", header = TRUE)
  differential_analysis_results <- differential_analysis_results[,c(1,2,3,4,5)]
  differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
  print(marker)
  print(differential_analysis_results)
  peak_annotation_plots[[marker]] <- differential_analysis_results
}

mintchip_annotation_barplots <- plotAnnoBar(peak_annotation_plots, ylab = "Percentage", title = NULL) + theme_classic(base_size = 18)
ggsave(filename = paste0("C:/Users/willi/Desktop/", "mintchip_genomic_features_2.png"), plot = mintchip_annotation_barplots, device='png', dpi=300, width = 8, units = "in")

      
