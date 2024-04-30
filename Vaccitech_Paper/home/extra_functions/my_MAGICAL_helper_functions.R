
combineRows <- function(df1, df2, matches) {
  combined <- data.frame()
  for (i in 1:nrow(matches)) {
    queryRow <- matches[i, "queryHits"]
    subjectRow <- matches[i, "subjectHits"]
    combinedRow <- cbind(df1[queryRow,], df2[subjectRow,])
    combined <- rbind(combined, combinedRow)
  }
  return(combined)
}

# Fix hg19 coordinates to be hg38 for a subset of marker peaks
add_hg38_coordinates_to_marker_peaks <- function(marker_table, hg38_table) {
  # create GRanges objects for differentially accessible marker peaks and ALL marker peaks (with hg38 start/end attached as extra columns)
  marker_peaks_granges <- makeGRangesFromDataFrame(df = marker_table, keep.extra.columns = TRUE)
  marker_peaks_with_hg38_coordinates_granges <- makeGRangesFromDataFrame(df = hg38_table, keep.extra.columns = TRUE)
  # Find overlap based on the hg19 coordinates
  marker_overlap <- as.data.frame(findOverlaps(marker_peaks_granges, marker_peaks_with_hg38_coordinates_granges))
  marker_overlap_df <- combineRows(marker_table, hg38_table, marker_overlap)
  # Now we can replace the hg19 coordinates with hg38 coordinates
  marker_overlap_df$start <- marker_overlap_df$hg38_start
  marker_overlap_df$end <- marker_overlap_df$hg38_end
  marker_overlap_df <- marker_overlap_df[,-c(12:17)]
  return(marker_overlap_df)
}

add_info_stage_1_MAGICAL <- function(magical_table) {
  current_sc_fc <- c()
  current_pseudobulk_fc <- c()
  current_sc_p_val <- c()
  current_pseudobulk_p_val <- c()
  current_pseudobulk_robust_p_val <- c()
  current_min_pct_1 <- c()
  current_min_pct_2 <- c()
  for(current_index in 1:nrow(magical_table)) {
    current_row <- magical_table[current_index,]
    current_cell_type <- current_row$Cell_Type
    current_sc_file <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", sub(" ", "_", current_cell_type), "-time_point-controlling_for_subject_id_sc_pct_0.01.tsv"), sep = "\t",
                              header = TRUE)
    current_pseudobulk_file <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", sub(" ", "_", current_cell_type), "-time_point-controlling_for_subject_id_pseudobulk_unfiltered.tsv"), sep = "\t",
                                  header = TRUE)
    current_peak <- paste0(current_row$Peak_chr, "-", current_row$Peak_start, "-", current_row$Peak_end)
    current_sc_row <- current_sc_file[rownames(current_sc_file) == current_peak,]
    current_pseudobulk_row <- current_pseudobulk_file[rownames(current_pseudobulk_file) == current_peak,]
    current_sc_fc <- c(current_sc_fc, current_sc_row$avg_log2FC)
    current_min_pct_1 <- c(current_min_pct_1, current_sc_row$pct.1)
    current_min_pct_2 <- c(current_min_pct_2, current_sc_row$pct.2)
    current_sc_p_val <- c(current_sc_p_val, current_sc_row$p_val)
    if(nrow(current_pseudobulk_row) == 1) {
      current_pseudobulk_fc <- c(current_pseudobulk_fc, current_pseudobulk_row$logFC_pseudobulk)
      current_pseudobulk_p_val <- c(current_pseudobulk_p_val, current_pseudobulk_row$p_value_pseudobulk)
      current_pseudobulk_robust_p_val <- c(current_pseudobulk_robust_p_val, current_pseudobulk_row$robust_p_value_pseudobulk)
    } else {
      current_pseudobulk_fc <- c(current_pseudobulk_fc, 0)
      current_pseudobulk_p_val <- c(current_pseudobulk_p_val, 1)
      current_pseudobulk_robust_p_val <- c(current_pseudobulk_robust_p_val, 1)
    }
  }
  magical_table$sc_p_val <- current_sc_p_val
  magical_table$pseudobulk_p_val <- current_pseudobulk_p_val
  magical_table$pseudobulk_robust_p_val <- current_pseudobulk_robust_p_val
  magical_table$pct.1 <- current_min_pct_1
  magical_table$pct.2 <- current_min_pct_2
  magical_table$sc_fc <- current_sc_fc
  magical_table$pseudobulk_fc <- current_pseudobulk_fc
  return(magical_table)
}

create_magical_gene_overlap_df <- function(magical_table, bulk_LRT) {
  # Grab mintchip data
  pos_mintchip_das <- list()
  neg_mintchip_das <- list()
  for(mintchip_marker in mintchip_markers) {
    current_pos_mintchip_das <- read.table(paste0(mintchip_das_dir, mintchip_marker, "/upregulated/", mintchip_marker, 
                                                  "_DESeq2_FC_0.1_upregulated_annotated.tsv"),
                                           sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    current_neg_mintchip_das <- read.table(paste0(mintchip_das_dir, mintchip_marker, "/downregulated/", mintchip_marker, 
                                                  "_DESeq2_FC_0.1_downregulated_annotated.tsv"),
                                           sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    pos_mintchip_das[[mintchip_marker]] <- unique(current_pos_mintchip_das$SYMBOL)
    neg_mintchip_das[[mintchip_marker]] <- unique(current_neg_mintchip_das$SYMBOL)
  }
  # Fill in vectors for different columns
  pos_H3K4me1 <- c()
  neg_H3K4me1 <- c()
  pos_H3K4me3 <- c()
  neg_H3K4me3 <- c()
  pos_H3K9me3 <- c()
  neg_H3K9me3 <- c()
  pos_H3K27Ac <- c()
  neg_H3K27Ac <- c()
  pos_H3K27me3 <- c()
  neg_H3K27me3 <- c()
  pos_H3K36me3 <- c()
  neg_H3K36me3 <- c()
  pos_snME <- c()
  neg_snME <- c()
  found_in_bulk_LRT <- c()
  # Parse MAGICAL table
  for(current_index in 1:nrow(magical_table)) {
    # Grab current gene
    current_row <- magical_table[current_index,]
    current_cell_type <- current_row$Cell_Type
    current_gene <- current_row$Gene_symbol
    # Is this gene found as closest gene to any mintchip site for any marker?
    pos_H3K4me1 <- c(pos_H3K4me1, ifelse(current_gene %in% pos_mintchip_das[["H3K4me1"]], TRUE, FALSE))
    neg_H3K4me1 <- c(neg_H3K4me1, ifelse(current_gene %in% neg_mintchip_das[["H3K4me1"]], TRUE, FALSE))
    pos_H3K4me3 <- c(pos_H3K4me3, ifelse(current_gene %in% pos_mintchip_das[["H3K4me3"]], TRUE, FALSE))
    neg_H3K4me3 <- c(neg_H3K4me3, ifelse(current_gene %in% neg_mintchip_das[["H3K4me3"]], TRUE, FALSE))
    pos_H3K9me3 <- c(pos_H3K9me3, ifelse(current_gene %in% pos_mintchip_das[["H3K9me3"]], TRUE, FALSE))
    neg_H3K9me3 <- c(neg_H3K9me3, ifelse(current_gene %in% neg_mintchip_das[["H3K9me3"]], TRUE, FALSE))
    pos_H3K27Ac <- c(pos_H3K27Ac, ifelse(current_gene %in% pos_mintchip_das[["H3K27Ac"]], TRUE, FALSE))
    neg_H3K27Ac <- c(neg_H3K27Ac, ifelse(current_gene %in% neg_mintchip_das[["H3K27Ac"]], TRUE, FALSE))
    pos_H3K27me3 <- c(pos_H3K27me3, ifelse(current_gene %in% pos_mintchip_das[["H3K27me3"]], TRUE, FALSE))
    neg_H3K27me3 <- c(neg_H3K27me3, ifelse(current_gene %in% neg_mintchip_das[["H3K27me3"]], TRUE, FALSE))
    pos_H3K36me3 <- c(pos_H3K36me3, ifelse(current_gene %in% pos_mintchip_das[["H3K36me3"]], TRUE, FALSE))
    neg_H3K36me3 <- c(neg_H3K36me3, ifelse(current_gene %in% neg_mintchip_das[["H3K36me3"]], TRUE, FALSE))
    # Is this gene found in snME for the associated cell type(s)?
    if(current_cell_type == "B") {
      snME_cell_types <- c("B-Mem", "B-Naive")
    } else if(current_cell_type == "CD14_Mono" | current_cell_type == "CD16_Mono") {
      snME_cell_types <- "Monocyte"
    } else if(current_cell_type == "NK") {
      snME_cell_types <- "NK-cell2"
    } else if(current_cell_type == "CD4_Naive") {
      snME_cell_types <- "Th-Naive"
    } else if(current_cell_type == "CD8_Naive") {
      snME_cell_types <- "Tc-Naive"
    } else if(current_cell_type == "CD4_Memory") {
      snME_cell_types <- c("Th-Mem")
    } else if(current_cell_type == "CD8_Memory") {
      snME_cell_types <- c("Tc-Mem")
    } else {
      snME_cell_types <- c()
    }
    current_pos_snME <- FALSE
    current_neg_snME <- FALSE
    for(snME_cell_type in snME_cell_types) {
      current_pos_snME_dmrs <- read.table(paste0(snME_results_dir, "Closest_Gene_Analysis_No_Promoter_Subset/", snME_cell_type, "_D28_hypermethylated_annotated_genes.tsv"),
                                          sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      current_neg_snME_dmrs <- read.table(paste0(snME_results_dir, "Closest_Gene_Analysis_No_Promoter_Subset/", snME_cell_type, "_D_minus_1_hypermethylated_annotated_genes.tsv"),
                                          sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      if(current_gene %in% current_pos_snME_dmrs$SYMBOL) {
        current_pos_snME <- TRUE
      }
      
      if(current_gene %in% current_neg_snME_dmrs$SYMBOL) {
        current_neg_snME <- TRUE
      }
    }
    pos_snME <- c(pos_snME, current_pos_snME)
    neg_snME <- c(neg_snME, current_neg_snME)
    found_in_bulk_LRT <- c(found_in_bulk_LRT, ifelse(current_gene %in% rownames(hvl_placebo_LRT_analysis_results_filtered), TRUE, FALSE))
  }
  magical_gene_overlap_df <- data.frame(Cell_Type = magical_table$Cell_Type, Gene_Name = magical_table$Gene_symbol,
                                        H3K4me1_pos = pos_H3K4me1, H3K4me1_neg = neg_H3K4me1, 
                                        H3K4me3_pos = pos_H3K4me3, H3K4me3_neg = neg_H3K4me3,
                                        H3K9me3_pos = pos_H3K9me3, H3K9me3_neg = neg_H3K9me3,
                                        H3K27Ac_pos = pos_H3K27Ac, H3K27Ac_neg = neg_H3K27Ac,
                                        H3K27me3_pos = pos_H3K27me3, H3K27me3_neg = neg_H3K27me3,
                                        H3K36me3_pos = pos_H3K36me3, H3K36me3_neg = neg_H3K36me3,
                                        pos_snME = pos_snME, neg_snME = neg_snME, bulk_LRT = found_in_bulk_LRT)
  return(magical_gene_overlap_df)
}

create_magical_site_overlap_df <- function(magical_table) {
  # Grab mintchip data
  pos_mintchip_das <- list()
  neg_mintchip_das <- list()
  for(mintchip_marker in mintchip_markers) {
    current_pos_mintchip_das <- read.table(paste0(mintchip_das_dir, mintchip_marker, "/upregulated/", mintchip_marker, 
                                                  "_DESeq2_FC_0.1_upregulated.tsv"),
                                           sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    current_neg_mintchip_das <- read.table(paste0(mintchip_das_dir, mintchip_marker, "/downregulated/", mintchip_marker, 
                                                  "_DESeq2_FC_0.1_downregulated.tsv"),
                                           sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    current_hg38_table <- read.table(paste0(mintchip_das_dir, mintchip_marker, "/", mintchip_marker, "_all_peaks_with_hg38_coordinates.tsv"),
               sep = "\t", header = TRUE)
    current_pos_mintchip_das_with_hg38 <- add_hg38_coordinates_to_marker_peaks(current_pos_mintchip_das,
                                                                               current_hg38_table)
    current_neg_mintchip_das_with_hg38 <- add_hg38_coordinates_to_marker_peaks(current_neg_mintchip_das,
                                                                               current_hg38_table)
    
    pos_mintchip_das[[mintchip_marker]] <- current_pos_mintchip_das_with_hg38
    neg_mintchip_das[[mintchip_marker]] <- current_neg_mintchip_das_with_hg38
  }
  # Fill in vectors for different columns
  overlap_pos_H3K4me1 <- c()
  overlap_neg_H3K4me1 <- c()
  overlap_pos_H3K4me3 <- c()
  overlap_neg_H3K4me3 <- c()
  overlap_pos_H3K9me3 <- c()
  overlap_neg_H3K9me3 <- c()
  overlap_pos_H3K27Ac <- c()
  overlap_neg_H3K27Ac <- c()
  overlap_pos_H3K27me3 <- c()
  overlap_neg_H3K27me3 <- c()
  overlap_pos_H3K36me3 <- c()
  overlap_neg_H3K36me3 <- c()
  overlap_pos_snME <- c()
  overlap_neg_snME <- c()
  # Parse MAGICAL table
  for(current_index in 1:nrow(magical_table)) {
    # Grab current gene
    current_row <- magical_table[current_index,]
    current_cell_type <- current_row$Cell_Type
    current_row$chr <- current_row$Peak_chr
    current_row$start <- current_row$Peak_start
    current_row$end <- current_row$Peak_end
    current_peak_granges <- makeGRangesFromDataFrame(df = current_row, keep.extra.columns = TRUE)
    overlap_pos_H3K4me1 <- c(overlap_pos_H3K4me1, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = pos_mintchip_das[["H3K4me1"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    overlap_neg_H3K4me1 <- c(overlap_neg_H3K4me1, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = neg_mintchip_das[["H3K4me1"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    
    overlap_pos_H3K4me3 <- c(overlap_pos_H3K4me3, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = pos_mintchip_das[["H3K4me3"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    overlap_neg_H3K4me3 <- c(overlap_neg_H3K4me3, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = neg_mintchip_das[["H3K4me3"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    
    overlap_pos_H3K9me3 <- c(overlap_pos_H3K9me3, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = pos_mintchip_das[["H3K9me3"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    overlap_neg_H3K9me3 <- c(overlap_neg_H3K9me3, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = neg_mintchip_das[["H3K9me3"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    
    overlap_pos_H3K27Ac <- c(overlap_pos_H3K27Ac, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = pos_mintchip_das[["H3K27Ac"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    overlap_neg_H3K27Ac <- c(overlap_neg_H3K27Ac, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = neg_mintchip_das[["H3K27Ac"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    
    overlap_pos_H3K27me3 <- c(overlap_pos_H3K27me3, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = pos_mintchip_das[["H3K27me3"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    overlap_neg_H3K27me3 <- c(overlap_neg_H3K27me3, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = neg_mintchip_das[["H3K27me3"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    
    overlap_pos_H3K36me3 <- c(overlap_pos_H3K36me3, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = pos_mintchip_das[["H3K36me3"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    overlap_neg_H3K36me3 <- c(overlap_neg_H3K36me3, ifelse(nrow(as.data.frame(findOverlaps(current_peak_granges, makeGRangesFromDataFrame(df = neg_mintchip_das[["H3K36me3"]], keep.extra.columns = TRUE)))) > 0, TRUE, FALSE))
    
    # Is this gene found in snME for the associated cell type(s)?
    if(current_cell_type == "B") {
      snME_cell_types <- c("B-Mem", "B-Naive")
    } else if(current_cell_type == "CD14_Mono" | current_cell_type == "CD16_Mono") {
      snME_cell_types <- "Monocyte"
    } else if(current_cell_type == "NK") {
      snME_cell_types <- "NK-cell2"
    } else if(current_cell_type == "CD4_Naive") {
      snME_cell_types <- "Th-Naive"
    } else if(current_cell_type == "CD8_Naive") {
      snME_cell_types <- "Tc-Naive"
    } else if(current_cell_type == "CD4_Memory") {
      snME_cell_types <- c("Th-Mem")
    } else if(current_cell_type == "CD8_Memory") {
      snME_cell_types <- c("Tc-Mem")
    } else {
      snME_cell_types <- c()
    }
    current_pos_snME <- FALSE
    current_neg_snME <- FALSE
    for(snME_cell_type in snME_cell_types) {
      current_pos_snME_dmrs <- read.table(paste0(snME_results_dir, "Closest_Gene_Analysis_No_Promoter_Subset/", snME_cell_type, "_D28_hypermethylated_annotated_genes.tsv"),
                                          sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      current_neg_snME_dmrs <- read.table(paste0(snME_results_dir, "Closest_Gene_Analysis_No_Promoter_Subset/", snME_cell_type, "_D_minus_1_hypermethylated_annotated_genes.tsv"),
                                          sep = "\t", header = TRUE, comment.char = "", quote = "\"")
      current_pos_snME_dmrs$start <- current_pos_snME_dmrs$start - 100
      current_pos_snME_dmrs$end <- current_pos_snME_dmrs$end + 100
      current_neg_snME_dmrs$start <- current_neg_snME_dmrs$start - 100
      current_neg_snME_dmrs$end <- current_neg_snME_dmrs$end + 100
      
      current_pos_snME_dmrs <- makeGRangesFromDataFrame(df = current_pos_snME_dmrs, keep.extra.columns = TRUE)
      current_neg_snME_dmrs <- makeGRangesFromDataFrame(df = current_neg_snME_dmrs, keep.extra.columns = TRUE)
      
      if(nrow(as.data.frame(findOverlaps(current_peak_granges, current_pos_snME_dmrs))) > 0) {
        current_pos_snME <- TRUE
      }
      
      if(nrow(as.data.frame(findOverlaps(current_peak_granges, current_neg_snME_dmrs))) > 0) {
        current_neg_snME <- TRUE
      }
    }
    overlap_pos_snME <- c(overlap_pos_snME, current_pos_snME)
    overlap_neg_snME <- c(overlap_neg_snME, current_neg_snME)
  }
  
  magical_site_overlap_df <- data.frame(Cell_Type = magical_table$Cell_Type, MAGICAL_peak_chr = magical_table$Peak_chr,
                                        MAGICAL_peak_start <- magical_table$Peak_start, MAGICAL_peak_end <- magical_table$Peak_end,
                                        overlap_pos_H3K4me1 = overlap_pos_H3K4me1, overlap_neg_H3K4me1 = overlap_neg_H3K4me1, 
                                        overlap_pos_H3K4me3 = overlap_pos_H3K4me3, overlap_neg_H3K4me3 = overlap_neg_H3K4me3, 
                                        overlap_pos_H3K9me3 = overlap_pos_H3K9me3, overlap_neg_H3K9me3 = overlap_neg_H3K9me3, 
                                        overlap_pos_H3K27Ac = overlap_pos_H3K27Ac, overlap_neg_H3K27Ac = overlap_neg_H3K27Ac, 
                                        overlap_pos_H3K27me3 = overlap_pos_H3K27me3, overlap_neg_H3K27me3 = overlap_neg_H3K27me3, 
                                        overlap_pos_H3K36me3 = overlap_pos_H3K36me3, overlap_neg_H3K36me3 = overlap_neg_H3K36me3, 
                                        overlap_pos_snME = overlap_pos_snME, overlap_neg_snME = overlap_neg_snME)
  return(magical_site_overlap_df)
}
