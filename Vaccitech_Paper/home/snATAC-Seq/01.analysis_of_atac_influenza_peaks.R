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

peak_annotation_plots <- list()

innate_snATAC_cell_types <- c("pDC", "CD16 Mono", "cDC", "NK", "CD14 Mono")
snATAC_peak_annotated_dir <- paste0(scATAC_hvl_placebo_das_dir, "annotated/")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

for(snATAC_cell_type in innate_snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  differential_analysis_results_file <- paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_0.01.tsv")
  differential_analysis_results <- read.table(differential_analysis_results_file, sep = "\t", header = TRUE)
  current_peaks <- differential_analysis_results$Peak_Name
  chromosomes <- sapply(strsplit(current_peaks, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 3))
  differential_analysis_results$seqnames <- chromosomes
  differential_analysis_results$start <- start_coords
  differential_analysis_results$end <- end_coords
  differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
  peak_annotation_plots[[snATAC_cell_type]] <- differential_analysis_results
}

plotAnnoBar(peak_annotation_plots, ylab = "Percentage", title = "Distribution of Genomic Features for scATAC-Seq DAS")

get_enrichr_results <- function(gene_list) {
  # Set up databases
  # Use dbs <- listEnrichrDbs() to get full list of databases
  dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "Reactome_2022")
  
  # Get results
  results <- enrichr(gene_list, dbs)
  results[[1]] <- results[[1]][results[[1]]$Adjusted.P.value < 0.05,]
  results[[2]] <- results[[2]][results[[2]]$Adjusted.P.value < 0.05,]
  results[[3]] <- results[[3]][results[[3]]$Adjusted.P.value < 0.05,]
  results[[4]] <- results[[4]][results[[4]]$Adjusted.P.value < 0.05,]
  return(results)
}

# Run enrichment analysis (pos and neg, promoter regions)
pos_fmd_promoter_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  pos_fmd_promoter_list[[snATAC_cell_type_for_file_name]] <- list()
  for(fc in c(0.1, 0.2, 0.3, 0.585, 1, 2)) {
    pos_differential_analysis_results_file <- paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_0.01_fc_", fc, "_upregulated_annotated.tsv")
    if(file.exists(pos_differential_analysis_results_file) && file.size(pos_differential_analysis_results_file) != 1 && file.size(pos_differential_analysis_results_file) != 75) {
      pos_differential_analysis_results_file <- read.table(pos_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      pos_differential_analysis_results_file <- subset(pos_differential_analysis_results_file, annotation %in% promoter_terms)
      pos_genes <- unique(pos_differential_analysis_results_file$SYMBOL)
      print(paste0("Number of genes is ", length(pos_genes)))
      pos_results <- run_fmd_on_snATAC(pos_genes)
      pos_fmd_promoter_list[[snATAC_cell_type_for_file_name]][[as.character(fc)]] <- pos_results
    }
  }
}

# saveRDS(pos_fmd_promoter_list, file = paste0(snATAC_peak_annotated_dir, "pos_fmd_promoter.RDS"))
# pos_fmd_promoter_list <- readRDS(paste0(snATAC_peak_annotated_dir, "pos_fmd_promoter.RDS"))

neg_fmd_promoter_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  neg_fmd_promoter_list[[snATAC_cell_type_for_file_name]] <- list()
  for(fc in c(0.1, 0.2, 0.3, 0.585, 1, 2)) {
    neg_differential_analysis_results_file <- paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_0.01_fc_", fc, "_downregulated_annotated.tsv")
    if(file.exists(neg_differential_analysis_results_file) && file.size(neg_differential_analysis_results_file) != 1 && file.size(neg_differential_analysis_results_file) != 75) {
      neg_differential_analysis_results_file <- read.table(neg_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      neg_differential_analysis_results_file <- subset(neg_differential_analysis_results_file, annotation %in% promoter_terms)
      neg_genes <- unique(neg_differential_analysis_results_file$SYMBOL)
      print(paste0("Number of genes is ", length(neg_genes)))
      neg_results <- run_fmd_on_snATAC(neg_genes)
      neg_fmd_promoter_list[[snATAC_cell_type_for_file_name]][[as.character(fc)]] <- neg_results
    }
  }
}

# saveRDS(neg_fmd_promoter_list, file = paste0(snATAC_peak_annotated_dir, "neg_fmd_promoter.RDS"))
# neg_fmd_promoter_list <- readRDS(paste0(snATAC_peak_annotated_dir, "neg_fmd_promoter.RDS"))
