
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

add_info_stage_1_MAGICAL <- function(magical_table, deg_table) {
  # Add RNA info
  rna_sc_fc <- c()
  rna_pseudobulk_fc <- c()
  rna_sc_p_val <- c()
  rna_pseudobulk_p_val <- c()
  for(rna_index in 1:nrow(magical_table)) {
    rna_row <- magical_table[rna_index,]
    rna_cell_type <- rna_row$Cell_Type
    rna_cell_type <- gsub("_", " ", rna_cell_type)
    cell_type_subset_deg_table <- deg_table[deg_table$Cell_Type == rna_cell_type,]
    cell_type_subset_deg_table <- cell_type_subset_deg_table[cell_type_subset_deg_table$Gene_Name == rna_row$Gene_symbol,]
    rna_sc_fc <- c(rna_sc_fc, cell_type_subset_deg_table$sc_log2FC)
    rna_pseudobulk_fc <- c(rna_pseudobulk_fc, cell_type_subset_deg_table$pseudo_bulk_log2FC)
    rna_sc_p_val <- c(rna_sc_p_val, cell_type_subset_deg_table$sc_pval_adj)
    rna_pseudobulk_p_val <- c(rna_pseudobulk_p_val, cell_type_subset_deg_table$pseudo_bulk_pval)
  }
  
  magical_table$rna_sc_p_val <- rna_sc_p_val
  magical_table$rna_pseudobulk_p_val <- rna_pseudobulk_p_val
  magical_table$rna_sc_fc <- rna_sc_fc
  magical_table$rna_pseudobulk_fc <- rna_pseudobulk_fc
  
  # Add ATAC info
  atac_sc_fc <- c()
  atac_pseudobulk_fc <- c()
  atac_sc_p_val <- c()
  atac_pseudobulk_p_val <- c()
  atac_pseudobulk_robust_p_val <- c()
  atac_min_pct_1 <- c()
  atac_min_pct_2 <- c()
  for(atac_index in 1:nrow(magical_table)) {
    atac_row <- magical_table[atac_index,]
    atac_cell_type <- atac_row$Cell_Type
    atac_sc_file <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", sub(" ", "_", atac_cell_type), "-time_point-controlling_for_subject_id_sc_pct_0.01.tsv"), sep = "\t",
                              header = TRUE)
    atac_pseudobulk_file <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", sub(" ", "_", atac_cell_type), "-time_point-controlling_for_subject_id_pseudobulk_unfiltered.tsv"), sep = "\t",
                                  header = TRUE)
    atac_peak <- paste0(atac_row$Peak_chr, "-", atac_row$Peak_start, "-", atac_row$Peak_end)
    atac_sc_row <- atac_sc_file[rownames(atac_sc_file) == atac_peak,]
    atac_pseudobulk_row <- atac_pseudobulk_file[rownames(atac_pseudobulk_file) == atac_peak,]
    atac_sc_fc <- c(atac_sc_fc, atac_sc_row$avg_log2FC)
    atac_min_pct_1 <- c(atac_min_pct_1, atac_sc_row$pct.1)
    atac_min_pct_2 <- c(atac_min_pct_2, atac_sc_row$pct.2)
    atac_sc_p_val <- c(atac_sc_p_val, atac_sc_row$p_val)
    if(nrow(atac_pseudobulk_row) == 1) {
      atac_pseudobulk_fc <- c(atac_pseudobulk_fc, atac_pseudobulk_row$logFC_pseudobulk)
      atac_pseudobulk_p_val <- c(atac_pseudobulk_p_val, atac_pseudobulk_row$p_value_pseudobulk)
      atac_pseudobulk_robust_p_val <- c(atac_pseudobulk_robust_p_val, atac_pseudobulk_row$robust_p_value_pseudobulk)
    } else {
      atac_pseudobulk_fc <- c(atac_pseudobulk_fc, 0)
      atac_pseudobulk_p_val <- c(atac_pseudobulk_p_val, 1)
      atac_pseudobulk_robust_p_val <- c(atac_pseudobulk_robust_p_val, 1)
    }
  }
  magical_table$atac_sc_p_val <- atac_sc_p_val
  magical_table$atac_pseudobulk_p_val <- atac_pseudobulk_p_val
  magical_table$atac_pseudobulk_robust_p_val <- atac_pseudobulk_robust_p_val
  magical_table$atac_pct.1 <- atac_min_pct_1
  magical_table$atac_pct.2 <- atac_min_pct_2
  magical_table$atac_sc_fc <- atac_sc_fc
  magical_table$atac_pseudobulk_fc <- atac_pseudobulk_fc
  return(magical_table)
}

create_magical_gene_overlap_df <- function(magical_table, bulk_D5 = NULL, bulk_D8 = NULL) {
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
  bulk_D5_adj_p_value <- c()
  bulk_D5_fc <- c()
  bulk_D8_adj_p_value <- c()
  bulk_D8_fc <- c()
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
    # Day 5 bulk
    if(current_gene %in% rownames(bulk_D5[[2]])) {
      current_gene_bulk_info <- bulk_D5[[2]][rownames(bulk_D5[[2]]) %in% current_gene,]
      bulk_D5_adj_p_value <- c(bulk_D5_adj_p_value, current_gene_bulk_info$padj)
      bulk_D5_fc <- c(bulk_D5_fc, current_gene_bulk_info$log2FoldChange)
    } else if(current_gene %in% rownames(bulk_D5[[8]])) {
      current_gene_bulk_info <- bulk_D5[[8]][rownames(bulk_D5[[8]]) %in% current_gene,]
      bulk_D5_adj_p_value <- c(bulk_D5_adj_p_value, current_gene_bulk_info$padj)
      bulk_D5_fc <- c(bulk_D5_fc, current_gene_bulk_info$log2FoldChange)
    } else {
      bulk_D5_adj_p_value <- c(bulk_D5_adj_p_value, 1)
      bulk_D5_fc <- c(bulk_D5_fc, 0)
    }
    # Day 8 bulk
    if(current_gene %in% rownames(bulk_D8[[2]])) {
      current_gene_bulk_info <- bulk_D8[[2]][rownames(bulk_D8[[2]]) %in% current_gene,]
      bulk_D8_adj_p_value <- c(bulk_D8_adj_p_value, current_gene_bulk_info$padj)
      bulk_D8_fc <- c(bulk_D8_fc, current_gene_bulk_info$log2FoldChange)
    } else if(current_gene %in% rownames(bulk_D8[[8]])) {
      current_gene_bulk_info <- bulk_D8[[8]][rownames(bulk_D8[[8]]) %in% current_gene,]
      bulk_D8_adj_p_value <- c(bulk_D8_adj_p_value, current_gene_bulk_info$padj)
      bulk_D8_fc <- c(bulk_D8_fc, current_gene_bulk_info$log2FoldChange)
    } else {
      bulk_D8_adj_p_value <- c(bulk_D8_adj_p_value, 1)
      bulk_D8_fc <- c(bulk_D8_fc, 0)
    }
  }
  magical_gene_overlap_df <- data.frame(Cell_Type = magical_table$Cell_Type, Gene_Name = magical_table$Gene_symbol,
                                        H3K4me1_pos = pos_H3K4me1, H3K4me1_neg = neg_H3K4me1, 
                                        H3K4me3_pos = pos_H3K4me3, H3K4me3_neg = neg_H3K4me3,
                                        H3K9me3_pos = pos_H3K9me3, H3K9me3_neg = neg_H3K9me3,
                                        H3K27Ac_pos = pos_H3K27Ac, H3K27Ac_neg = neg_H3K27Ac,
                                        H3K27me3_pos = pos_H3K27me3, H3K27me3_neg = neg_H3K27me3,
                                        H3K36me3_pos = pos_H3K36me3, H3K36me3_neg = neg_H3K36me3,
                                        pos_snME = pos_snME, neg_snME = neg_snME, bulk_D5_adj_p_value = bulk_D5_adj_p_value,
                                        bulk_D5_fc = bulk_D5_fc, bulk_D8_adj_p_value = bulk_D8_adj_p_value, bulk_D8_fc = bulk_D8_fc)
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
                                        MAGICAL_peak_start = magical_table$Peak_start, MAGICAL_peak_end = magical_table$Peak_end,
                                        overlap_pos_H3K4me1 = overlap_pos_H3K4me1, overlap_neg_H3K4me1 = overlap_neg_H3K4me1, 
                                        overlap_pos_H3K4me3 = overlap_pos_H3K4me3, overlap_neg_H3K4me3 = overlap_neg_H3K4me3, 
                                        overlap_pos_H3K9me3 = overlap_pos_H3K9me3, overlap_neg_H3K9me3 = overlap_neg_H3K9me3, 
                                        overlap_pos_H3K27Ac = overlap_pos_H3K27Ac, overlap_neg_H3K27Ac = overlap_neg_H3K27Ac, 
                                        overlap_pos_H3K27me3 = overlap_pos_H3K27me3, overlap_neg_H3K27me3 = overlap_neg_H3K27me3, 
                                        overlap_pos_H3K36me3 = overlap_pos_H3K36me3, overlap_neg_H3K36me3 = overlap_neg_H3K36me3, 
                                        overlap_pos_snME = overlap_pos_snME, overlap_neg_snME = overlap_neg_snME)
  return(magical_site_overlap_df)
}

create_tf_targets_df <- function(magical_table) {
  split_df <- magical_table %>%
    mutate(TFs.binding.prob. = strsplit(as.character(TFs.binding.prob.), ",")) %>%
    unnest(TFs.binding.prob.) %>%
    mutate(TFs.binding.prob. = trimws(TFs.binding.prob.)) %>%
    filter(!is.na(TFs.binding.prob.))
  
  split_df <- subset(split_df, !grepl("^\\(", TFs.binding.prob.))
  TF_fc <- c()
  TF_sc_pval <- c()
  TF_pseudobulk_pval <- c()
  cell_types <- unique(split_df$Cell_Type)
  unfiltered_sc_results <- list()
  unfiltered_pseudobulk_results <- list()
  for(cell_type in cell_types) {
    current_unfiltered_sc_result <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", 
                                                      cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                               sep = "\t", header = TRUE)
    current_unfiltered_pseudobulk_result <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", 
                                                      cell_type, "-time_point-controlling_for_subject_id_pseudobulk_unfiltered.tsv"),
                                               sep = "\t", header = TRUE)
    unfiltered_sc_results[[cell_type]] <- current_unfiltered_sc_result
    unfiltered_pseudobulk_results[[cell_type]] <- current_unfiltered_pseudobulk_result
  }
  for(current_row_index in 1:nrow(split_df)) {
    current_row <- split_df[current_row_index,]
    current_cell_type <- current_row$Cell_Type
    current_sc_results <- unfiltered_sc_results[[current_cell_type]]
    current_sc_results <- current_sc_results[rownames(current_sc_results) == current_row$TFs.binding.prob.,]
    current_pseudobulk_results <- unfiltered_pseudobulk_results[[current_cell_type]]
    current_pseudobulk_results <- current_pseudobulk_results[rownames(current_pseudobulk_results) == current_row$TFs.binding.prob.,]
    TF_fc <- c(TF_fc, current_sc_results$avg_log2FC)
    TF_sc_pval <- c(TF_sc_pval, current_sc_results$p_val_adj)
    if(is.na(current_pseudobulk_results$pvalue)) {
      current_pseudobulk_results$pvalue <- 1
    }
    TF_pseudobulk_pval <- c(TF_pseudobulk_pval, current_pseudobulk_results$pvalue)
  }
  split_df$TF_fc <- TF_fc
  split_df$TF_sc_pval <- TF_sc_pval
  split_df$TF_pseudobulk_pval <- TF_pseudobulk_pval
  return(split_df)
}

create_tf_vs_cell_type_df <- function(magical_table) {
  split_df <- magical_table %>%
    mutate(TFs.binding.prob. = strsplit(as.character(TFs.binding.prob.), ",")) %>%
    unnest(TFs.binding.prob.) %>%
    mutate(TFs.binding.prob. = trimws(TFs.binding.prob.)) %>%
    filter(!is.na(TFs.binding.prob.))
  
  # Count the occurrences of each TF for each Cell_Type
  count_df <- split_df %>%
    group_by(Cell_Type, TFs.binding.prob.) %>%
    summarise(Count = n()) %>%
    arrange(Cell_Type, desc(Count))
  
  colnames(count_df) <- c("Cell_Type", "TF", "Count")
  count_df <- subset(count_df, !grepl("^\\(", TF))
  return(count_df)
}

find_tfs_in_scRNA_data <- function(tf_table, deg_table) {
  found_in_scRNA <- c()
  for(row_index in 1:nrow(tf_table)) {
    current_tf_row <- tf_table[row_index,]
    current_tf_cell_type <- gsub("_", " ", current_tf_row$Cell_Type)
    deg_table_subset <- deg_table[deg_table$Cell_Type == current_tf_cell_type,]
    found_in_scRNA <- c(found_in_scRNA, ifelse(current_tf_row$TF %in% deg_table_subset$Gene_Name, TRUE, FALSE))
  }
  tf_table$found_in_scRNA <- found_in_scRNA
  return(tf_table)
}

create_tf_heatmap_plot <- function(tf_table) {
  
  # Subset to top 5 TF for each cell type
  #tf_table <- tf_table %>%
  #  group_by(Cell_Type) %>%
  #  slice_head(n = 5)
  
  interesting_tfs <- c("BACH1", "BACH2", "CTCF", "ETS2", "FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "PATZ1", "PLAG1", "ZFX", "ZNF143", "ZNF281", "ZNF589")
  
  tf_table <- tf_table[tf_table$TF %in% interesting_tfs,]
  
  tf_table$Cell_Type <- gsub("_", " ",  tf_table$Cell_Type)
  tf_table$Cell_Type <- factor(tf_table$Cell_Type, levels = c("CD14 Mono", "CD16 Mono", "cDC", "pDC", "NK"))
  #ggplot(tf_table, aes(Cell_Type, TF, fill = Count)) + 
  #  geom_tile() + scale_fill_gradient(low="white", high="red") + 
  #  theme(axis.title.y=element_blank(),
  #        axis.ticks.y=element_blank())
  
  ggplot(tf_table, aes(Cell_Type, TF, fill = Count)) + 
    geom_tile() + scale_fill_gradientn(
      colors = c("ivory1", "lightcoral", "red"),
      values = scales::rescale(c(0, 20, 100))) +
    theme_minimal(base_size = 24) +
    theme(
      panel.background = element_rect(fill = "white", color = NA), # Set background to white
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1)
    ) +
    labs(title = NULL,
         x = NULL,
         y = NULL,
         fill = "Target Count")
}

find_overlapping_circuits <- function(magical_table) {
  current_circuits <- c()
  current_cell_types <- c()
  for(current_row_idx in 1:nrow(magical_table)) {
    current_row <- magical_table[current_row_idx,]
    current_circuits <- c(current_circuits, paste0(current_row$Gene_symbol, "-", current_row$Peak_chr, "-", current_row$Peak_start, "-", 
                              current_row$Peak_end))
    magical_table_subset <- magical_table[magical_table$Gene_symbol == current_row$Gene_symbol & 
                                            magical_table$Peak_chr == current_row$Peak_chr & 
                                            magical_table$Peak_start == current_row$Peak_start &
                                            magical_table$Peak_end == current_row$Peak_end,]
    current_cell_types <- c(current_cell_types, paste(magical_table_subset$Cell_Type, collapse = ","))
  }
  overlapping_circuit_df <- data.frame(circuits = current_circuits, cell_types = current_cell_types)
  return(overlapping_circuit_df)
}

find_overlapping_circuit_cell_type_matrix <- function(overlapping_circuits_df) {
  # Extract unique cell types
  vec <- overlapping_circuits_df$cell_types
  cell_types <- unique(unlist(strsplit(vec, ",")))
  
  # Create an empty matrix with rows and columns representing cell types
  mat <- matrix(0, nrow = length(cell_types), ncol = length(cell_types), dimnames = list(cell_types, cell_types))
  
  # Populate the matrix
  for (i in 1:length(vec)) {
    cells <- unlist(strsplit(vec[i], ","))
    for (j in 1:length(cells)) {
      for (k in 1:length(cells)) {
        mat[cells[j], cells[k]] <- mat[cells[j], cells[k]] + 1
      }
    }
  }
  
  diagonal_counts <- diag(mat)
  
  # Sort the cell types based on the counts on the diagonal
  sorted_cell_types <- names(sort(diagonal_counts, decreasing = TRUE))
  
  # Rearrange the rows and columns of the matrix
  mat <- mat[sorted_cell_types, sorted_cell_types]
  
  return(mat)
}
