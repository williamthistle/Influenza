# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

snATAC_cell_types <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "CD4 Naive", "CD8 Naive", "cDC", "MAIT", "Proliferating", "pDC")
promoter_terms <- c("Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)")

# Step 1: Add pct.1 and pct.2 to final list of peaks
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  for(pct in c(0.01, 0.05, 0.1)) {
    sc_differential_analysis_results_file <- paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc_pct_", pct, ".tsv")
    sc_differential_analysis_results <- read.table(sc_differential_analysis_results_file, sep = "\t", header = TRUE)
    final_differential_analysis_results_file <- paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_", pct, ".tsv")
    final_differential_analysis_results <- read.table(final_differential_analysis_results_file, sep = "\t", header = TRUE)
    pct_1_vec <- c()
    pct_2_vec <- c()
    for(current_row in 1:nrow(final_differential_analysis_results)) {
      current_sc_peak <- sc_differential_analysis_results[rownames(sc_differential_analysis_results) == final_differential_analysis_results[current_row,]$Peak_Name,]
      pct_1_vec <- c(pct_1_vec, current_sc_peak$pct.1)
      pct_2_vec <- c(pct_2_vec, current_sc_peak$pct.2)
    }
    final_differential_analysis_results$pct.1 <- pct_1_vec
    final_differential_analysis_results$pct.2 <- pct_2_vec
    write.table(final_differential_analysis_results, file = final_differential_analysis_results_file, quote = FALSE, sep = "\t")
  }
}

# Step 2: Remove any peaks that don't have 1% of cells in BOTH pct.1 and pct.2
# This eliminates weird outliers (e.g., crazy high FC due to 0 min.pct value for one group)
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  for(analysis_type in c("sc", "final")) {
    for(pct in c(0.01, 0.05, 0.1)) {
      differential_analysis_results_file <- paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", pct, ".tsv")
      differential_analysis_results <- read.table(differential_analysis_results_file, sep = "\t", header = TRUE)
      differential_analysis_results <- differential_analysis_results[differential_analysis_results$pct.1 > 0 & differential_analysis_results$pct.2 > 0,]
      write.table(differential_analysis_results, file = differential_analysis_results_file, quote = FALSE, sep = "\t")
    }
  }
}

# Step 3: Annotate peaks
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
snATAC_peak_annotated_dir <- paste0(scATAC_hvl_placebo_das_dir, "annotated/")
if (!dir.exists(snATAC_peak_annotated_dir)) {dir.create(snATAC_peak_annotated_dir)}

for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  for(analysis_type in c("sc", "final")) {
    for(pct in c(0.01, 0.05, 0.1)) {
      differential_analysis_results_file <- paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", pct, ".tsv")
      differential_analysis_results <- read.table(differential_analysis_results_file, sep = "\t", header = TRUE)
      if(analysis_type == "final") {
        current_peaks <- differential_analysis_results$Peak_Name
      } else {
        current_peaks <- rownames(differential_analysis_results)
      }
      chromosomes <- sapply(strsplit(current_peaks, "-"), `[`, 1)
      start_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 2))
      end_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 3))
      differential_analysis_results$seqnames <- chromosomes
      differential_analysis_results$start <- start_coords
      differential_analysis_results$end <- end_coords
      if(analysis_type == "final") {
        upregulated_differential_analysis_results <- differential_analysis_results[differential_analysis_results$sc_log2FC > 0,]
        downregulated_differential_analysis_results <- differential_analysis_results[differential_analysis_results$sc_log2FC < 0,]
      } else {
        upregulated_differential_analysis_results <- differential_analysis_results[differential_analysis_results$avg_log2FC > 0,]
        downregulated_differential_analysis_results <- differential_analysis_results[differential_analysis_results$avg_log2FC < 0,]
      }
      for(fc_threshold in c(0.1, 0.2, 0.3, 0.585, 1, 2)) {
        if(analysis_type == "final") {
          upregulated_differential_analysis_results_fc <- upregulated_differential_analysis_results[upregulated_differential_analysis_results$sc_log2FC >= fc_threshold,]
          downregulated_differential_analysis_results_fc <- downregulated_differential_analysis_results[downregulated_differential_analysis_results$sc_log2FC <= -fc_threshold,]
        } else {
          upregulated_differential_analysis_results_fc <- upregulated_differential_analysis_results[upregulated_differential_analysis_results$avg_log2FC >= fc_threshold,]
          downregulated_differential_analysis_results_fc <- downregulated_differential_analysis_results[downregulated_differential_analysis_results$avg_log2FC <= -fc_threshold,]
        }
       if(nrow(upregulated_differential_analysis_results_fc) > 0) {
          upregulated_differential_analysis_results_fc_annotated <- annotatePeak(makeGRangesFromDataFrame(upregulated_differential_analysis_results_fc), TxDb = txdb, annoDb = "org.Hs.eg.db")
          upregulated_differential_analysis_results_fc_annotated <- as.data.frame(upregulated_differential_analysis_results_fc_annotated)
          if(analysis_type == "final") {
            upregulated_differential_analysis_results_fc_annotated$sc_pval <- upregulated_differential_analysis_results_fc$sc_pval
            upregulated_differential_analysis_results_fc_annotated$sc_log2FC <- upregulated_differential_analysis_results_fc$sc_log2FC
            upregulated_differential_analysis_results_fc_annotated$pseudo_bulk_pval <- upregulated_differential_analysis_results_fc$pseudo_bulk_pval
            upregulated_differential_analysis_results_fc_annotated$pseudo_bulk_robust_pval <- upregulated_differential_analysis_results_fc$pseudo_bulk_robust_pval
            upregulated_differential_analysis_results_fc_annotated$pseudo_bulk_log2FC <- upregulated_differential_analysis_results_fc$pseudo_bulk_log2FC
            upregulated_differential_analysis_results_fc_annotated$pct.1 <- upregulated_differential_analysis_results_fc$pct.1
            upregulated_differential_analysis_results_fc_annotated$pct.2 <- upregulated_differential_analysis_results_fc$pct.2
          } else {
            upregulated_differential_analysis_results_fc_annotated$avg_log2FC <- upregulated_differential_analysis_results_fc$avg_log2FC
            upregulated_differential_analysis_results_fc_annotated$pct.1 <- upregulated_differential_analysis_results_fc$pct.1
            upregulated_differential_analysis_results_fc_annotated$pct.2 <- upregulated_differential_analysis_results_fc$pct.2
            upregulated_differential_analysis_results_fc_annotated$p_val_adj <- upregulated_differential_analysis_results_fc$p_val_adj
          }
          write.table(upregulated_differential_analysis_results_fc_annotated, file = paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", pct, "_fc_", fc_threshold, "_upregulated_annotated.tsv"), 
                      sep = "\t", quote = FALSE, row.names = FALSE)
        }
        if(nrow(downregulated_differential_analysis_results_fc) > 0) {
          downregulated_differential_analysis_results_fc_annotated <- annotatePeak(makeGRangesFromDataFrame(downregulated_differential_analysis_results_fc), TxDb = txdb, annoDb = "org.Hs.eg.db")
          downregulated_differential_analysis_results_fc_annotated <- as.data.frame(downregulated_differential_analysis_results_fc_annotated)
          if(analysis_type == "final") {
            downregulated_differential_analysis_results_fc_annotated$sc_pval <- downregulated_differential_analysis_results_fc$sc_pval
            downregulated_differential_analysis_results_fc_annotated$sc_log2FC <- downregulated_differential_analysis_results_fc$sc_log2FC
            downregulated_differential_analysis_results_fc_annotated$pseudo_bulk_pval <- downregulated_differential_analysis_results_fc$pseudo_bulk_pval
            downregulated_differential_analysis_results_fc_annotated$pseudo_bulk_robust_pval <- downregulated_differential_analysis_results_fc$pseudo_bulk_robust_pval
            downregulated_differential_analysis_results_fc_annotated$pseudo_bulk_log2FC <- downregulated_differential_analysis_results_fc$pseudo_bulk_log2FC
            downregulated_differential_analysis_results_fc_annotated$pct.1 <- downregulated_differential_analysis_results_fc$pct.1
            downregulated_differential_analysis_results_fc_annotated$pct.2 <- downregulated_differential_analysis_results_fc$pct.2
          } else {
            downregulated_differential_analysis_results_fc_annotated$avg_log2FC <- downregulated_differential_analysis_results_fc$avg_log2FC
            downregulated_differential_analysis_results_fc_annotated$pct.1 <- downregulated_differential_analysis_results_fc$pct.1
            downregulated_differential_analysis_results_fc_annotated$pct.2 <- downregulated_differential_analysis_results_fc$pct.2
            downregulated_differential_analysis_results_fc_annotated$p_val_adj <- downregulated_differential_analysis_results_fc$p_val_adj
          }
          write.table(downregulated_differential_analysis_results_fc_annotated, file = paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", pct, "_fc_", fc_threshold, "_downregulated_annotated.tsv"), 
                      sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }
    }
  }
}

# Step 4: Print peaks that fall in promoter regions specifically
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  for(analysis_type in c("sc", "final")) {
    for(pct in c(0.01, 0.05, 0.1)) {
      for(fc_threshold in c(0.1, 0.2, 0.3, 0.585, 1, 2)) {
        upregulated_annotated_file <- paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", pct, "_fc_", fc_threshold, "_upregulated_annotated.tsv")
        downregulated_annotated_file <- paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", pct, "_fc_", fc_threshold, "_downregulated_annotated.tsv")
      
        if(file.exists(upregulated_annotated_file)) {
          upregulated_annotated_results <- read.table(upregulated_annotated_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
          upregulated_annotated_results <- upregulated_annotated_results[upregulated_annotated_results$annotation %in% promoter_terms,]
          write.table(upregulated_annotated_results, file = paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", pct, "_fc_", fc_threshold, "_upregulated_promoter_subset_annotated.tsv"),
                      sep = "\t", quote = FALSE, row.names = FALSE)
        }
      
        if(file.exists(downregulated_annotated_file)) {
          downregulated_annotated_results <- read.table(downregulated_annotated_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
          downregulated_annotated_results <- downregulated_annotated_results[downregulated_annotated_results$annotation %in% promoter_terms,]
          write.table(downregulated_annotated_results, file = paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", pct, "_fc_", fc_threshold, "_downregulated_promoter_subset_annotated.tsv"),
                      sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }
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
pos_enrichment_promoter_list <- list()
for(snATAC_cell_type in snATAC_cell_types) {
  snATAC_cell_type_for_file_name <- sub(" ", "_", snATAC_cell_type)
  pos_enrichment_promoter_list[[snATAC_cell_type_for_file_name]] <- list()
  for(fc in c(0.1, 0.2, 0.3, 0.585, 1, 2)) {
    pos_differential_analysis_results_file <- paste0(snATAC_peak_annotated_dir, "D28-vs-D_minus_1-degs-", snATAC_cell_type_for_file_name, "-time_point-controlling_for_subject_id_final_pct_0.01_fc_", fc, "_upregulated_annotated.tsv")
    if(file.exists(pos_differential_analysis_results_file) && file.size(pos_differential_analysis_results_file) != 1 && file.size(pos_differential_analysis_results_file) != 75) {
      pos_differential_analysis_results_file <- read.table(pos_differential_analysis_results_file, sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      pos_differential_analysis_results_file <- subset(pos_differential_analysis_results_file, annotation %in% promoter_terms)
      pos_genes <- unique(pos_differential_analysis_results_file$SYMBOL)
      print(paste0("Number of genes is ", length(pos_genes)))
      pos_results <- run_fmd_on_snATAC(pos_genes)
      pos_enrichment_promoter_list[[snATAC_cell_type_for_file_name]][[as.character(fc)]] <- pos_results
    }
  }
}

# saveRDS(pos_enrichment_promoter_list, file = paste0(snATAC_peak_annotated_dir, "pos_fmd_promoter.RDS"))
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

saveRDS(neg_fmd_promoter_list, file = paste0(snATAC_peak_annotated_dir, "neg_fmd_promoter.RDS"))
# neg_fmd_promoter_list <- readRDS(paste0(snATAC_peak_annotated_dir, "neg_fmd_promoter.RDS"))
