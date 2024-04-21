# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))
# magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

run_fmd_on_snATAC <- function(gene_list) {
  if(length(gene_list) > 1 & length(gene_list) < 2000) {
    current_fmd_result <- SPEEDI::RunFMD_RNA(gene_list, "blood")
  } else { 
    current_fmd_result <- "EMPTY OR OVER 2000 PEAKS (TOO MANY)"
  }
  return(current_fmd_result)
}

snATAC_cell_types <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "CD4 Naive", "CD8 Naive", "cDC", "pDC", "MAIT", "Proliferating")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Create annotated up and downregulated peak files (annotated)
# NOTE: Check histogram of fold change for CD16 Mono.
snATAC_peak_annotated_dir <- paste0(scATAC_hvl_placebo_das_dir, "annotated/")
if (!dir.exists(snATAC_peak_annotated_dir)) {dir.create(snATAC_peak_annotated_dir)}
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  differential_analysis_results_file <- paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_0.01.tsv")
  differential_analysis_results <- read.table(differential_analysis_results_file, sep = "\t", header = TRUE)
  upregulated_differential_analysis_results <- differential_analysis_results[differential_analysis_results$sc_log2FC > 0,]
  downregulated_differential_analysis_results <- differential_analysis_results[differential_analysis_results$sc_log2FC < 0,]
  for(fc in c(0.1, 0.2, 0.3, 0.585, 1, 2)) {
    # Upregulated
    upregulated_differential_analysis_results_fc_subset <- upregulated_differential_analysis_results[upregulated_differential_analysis_results$sc_log2FC >= fc,]
    if(nrow(upregulated_differential_analysis_results_fc_subset) > 0) {
      current_peaks <- upregulated_differential_analysis_results_fc_subset$Peak_Name
      chromosomes <- sapply(strsplit(current_peaks, "-"), `[`, 1)
      start_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 2))
      end_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 3))
      upregulated_differential_analysis_results_fc_subset$seqnames <- chromosomes
      upregulated_differential_analysis_results_fc_subset$start <- start_coords
      upregulated_differential_analysis_results_fc_subset$end <- end_coords
      upregulated_differential_analysis_results_fc_subset <- upregulated_differential_analysis_results_fc_subset[,c(7,8,9)]
      upregulated_differential_analysis_results_fc_subset <- annotatePeak(makeGRangesFromDataFrame(upregulated_differential_analysis_results_fc_subset), TxDb = txdb, annoDb = "org.Hs.eg.db")
      upregulated_differential_analysis_results_fc_subset <- as.data.frame(upregulated_differential_analysis_results_fc_subset)
      write.table(upregulated_differential_analysis_results_fc_subset, file = paste0(snATAC_peak_annotated_dir, snATAC_cell_type_for_file_name, "_FC_", fc, "_upregulated_annotated.tsv"), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
    # Downregulated
    downregulated_differential_analysis_results_fc_subset <- downregulated_differential_analysis_results[downregulated_differential_analysis_results$sc_log2FC <= -fc,]
    if(nrow(downregulated_differential_analysis_results_fc_subset) > 0) {
      current_peaks <- downregulated_differential_analysis_results_fc_subset$Peak_Name
      chromosomes <- sapply(strsplit(current_peaks, "-"), `[`, 1)
      start_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 2))
      end_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 3))
      downregulated_differential_analysis_results_fc_subset$seqnames <- chromosomes
      downregulated_differential_analysis_results_fc_subset$start <- start_coords
      downregulated_differential_analysis_results_fc_subset$end <- end_coords
      downregulated_differential_analysis_results_fc_subset <- downregulated_differential_analysis_results_fc_subset[,c(7,8,9)]
      downregulated_differential_analysis_results_fc_subset <- annotatePeak(makeGRangesFromDataFrame(downregulated_differential_analysis_results_fc_subset), TxDb = txdb, annoDb = "org.Hs.eg.db")
      downregulated_differential_analysis_results_fc_subset <- as.data.frame(downregulated_differential_analysis_results_fc_subset)
      write.table(downregulated_differential_analysis_results_fc_subset, file = paste0(snATAC_peak_annotated_dir, snATAC_cell_type_for_file_name, "_FC_", fc, "_downregulated_annotated.tsv"), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}

# Run FMD (pos and neg)
pos_fmd_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  print(snATAC_cell_type)
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  pos_fmd_list[[snATAC_cell_type_for_file_name]] <- list()
  for(fc in c(0.1, 0.2, 0.3, 0.585, 1, 2)) {
    pos_differential_analysis_results_file <- paste0(snATAC_peak_annotated_dir, snATAC_cell_type_for_file_name, "_FC_", fc, "_upregulated_annotated.tsv")
    if(file.exists(pos_differential_analysis_results_file) && file.size(pos_differential_analysis_results_file) != 1 && file.size(pos_differential_analysis_results_file) != 75) {
      pos_differential_analysis_results_file <- read.table(pos_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      pos_genes <- unique(pos_differential_analysis_results_file$SYMBOL)
      pos_results <- run_fmd_on_snATAC(pos_genes)
      pos_fmd_list[[snATAC_cell_type_for_file_name]][[as.character(fc)]] <- pos_results
    }
  }
}

# saveRDS(pos_fmd_list, file = paste0(snATAC_peak_annotated_dir, "pos_fmd.RDS"))
# pos_fmd_list <- readRDS(paste0(snATAC_peak_annotated_dir, "pos_fmd.RDS"))

neg_fmd_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  neg_fmd_list[[snATAC_cell_type_for_file_name]] <- list()
  for(fc in c(0.1, 0.2, 0.3, 0.585, 1, 2)) {
    neg_differential_analysis_results_file <- paste0(snATAC_peak_annotated_dir, snATAC_cell_type_for_file_name, "_FC_", fc, "_downregulated_annotated.tsv")
    if(file.exists(neg_differential_analysis_results_file) && file.size(neg_differential_analysis_results_file) != 1 && file.size(neg_differential_analysis_results_file) != 75) {
      neg_differential_analysis_results_file <- read.table(neg_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      neg_genes <- unique(neg_differential_analysis_results_file$SYMBOL)
      neg_results <- run_fmd_on_snATAC(neg_genes)
      neg_fmd_list[[snATAC_cell_type_for_file_name]][[as.character(fc)]] <- neg_results
    }
  }
}

# saveRDS(neg_fmd_list, file = paste0(snATAC_peak_annotated_dir, "neg_fmd.RDS"))
# neg_fmd_list <- readRDS(paste0(snATAC_peak_annotated_dir, "neg_fmd.RDS"))

peak_annotation_plots <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  differential_analysis_results_file <- paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01.tsv")
  differential_analysis_results <- read.table(differential_analysis_results_file, sep = "\t", header = TRUE)
  current_peaks <- differential_analysis_results$Peak_Name
  chromosomes <- sapply(strsplit(current_peaks, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 3))
  differential_analysis_results$seqnames <- chromosomes
  differential_analysis_results$start <- start_coords
  differential_analysis_results$end <- end_coords
  differential_analysis_results <- differential_analysis_results[,c(7,8,9)]
  differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
  peak_annotation_plots[[snATAC_cell_type]] <- differential_analysis_results
}

plotAnnoBar(peak_annotation_plots, ylab = "Percentage", title = "Distribution of Genomic Features for snATAC-Seq Data")

# Run FMD (pos and neg, promoter regions)
pos_fmd_promoter_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  pos_fmd_promoter_list[[snATAC_cell_type_for_file_name]] <- list()
  for(fc in c(0.1)) {
    pos_differential_analysis_results_file <- paste0(snATAC_peak_annotated_dir, snATAC_cell_type_for_file_name, "_FC_", fc, "_upregulated_annotated.tsv")
    if(file.exists(pos_differential_analysis_results_file) && file.size(pos_differential_analysis_results_file) != 1 && file.size(pos_differential_analysis_results_file) != 75) {
      pos_differential_analysis_results_file <- read.table(pos_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      pos_differential_analysis_results_file <- subset(pos_differential_analysis_results_file, distanceToTSS >= -2000 & distanceToTSS <= 500)
      pos_genes <- unique(pos_differential_analysis_results_file$SYMBOL)
      pos_results <- run_fmd_on_snATAC(pos_genes)
      pos_fmd_promoter_list[[snATAC_cell_type_for_file_name]][[as.character(fc)]] <- pos_results
    }
  }
}

saveRDS(pos_fmd_promoter_list, file = paste0(snATAC_peak_annotated_dir, "pos_fmd_promoter.RDS"))
# pos_fmd_promoter_list <- readRDS(paste0(snATAC_peak_annotated_dir, "pos_fmd_promoter.RDS"))

neg_fmd_promoter_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  neg_fmd_promoter_list[[snATAC_cell_type_for_file_name]] <- list()
  for(fc in c(0.1)) {
    neg_differential_analysis_results_file <- paste0(snATAC_peak_annotated_dir, snATAC_cell_type_for_file_name, "_FC_", fc, "_downregulated_annotated.tsv")
    if(file.exists(neg_differential_analysis_results_file) && file.size(neg_differential_analysis_results_file) != 1 && file.size(neg_differential_analysis_results_file) != 75) {
      neg_differential_analysis_results_file <- read.table(neg_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      neg_differential_analysis_results_file <- subset(neg_differential_analysis_results_file, distanceToTSS >= -2000 & distanceToTSS <= 500)
      neg_genes <- unique(neg_differential_analysis_results_file$SYMBOL)
      neg_results <- run_fmd_on_snATAC(neg_genes)
      neg_fmd_promoter_list[[snATAC_cell_type_for_file_name]][[as.character(fc)]] <- neg_results
    }
  }
}

saveRDS(neg_fmd_promoter_list, file = paste0(snATAC_peak_annotated_dir, "neg_fmd_promoter.RDS"))
# neg_fmd_promoter_list <- readRDS(paste0(snATAC_peak_annotated_dir, "neg_fmd_promoter.RDS"))
