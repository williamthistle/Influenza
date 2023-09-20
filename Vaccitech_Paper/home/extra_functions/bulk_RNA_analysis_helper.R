# Method to set up bulk analysis
setup_bulk_analysis=function(base_dir, data_dir) {
  # Read in count and metadata files
  gene_counts <<- fread(paste0(data_dir, "rsem_genes_count.processed.txt"), header = T, sep = "\t")
  all_metadata_file <<- paste0(base_dir, "metadata/all_metadata_sheet.tsv")
  all_metadata <<- read.table(all_metadata_file, header = TRUE, sep = "\t")
  # Currently, only capturing viral load for placebo subjects
  viral_load_file <<- paste0(base_dir, "metadata/bulk_RNA_viral_load.tsv")
  viral_load <<- read.table(viral_load_file, header = TRUE, sep = "\t")
  viral_load_primary <<- viral_load[viral_load$PARAMCD == "QPCRAUC",]
  viral_load_primary <<- viral_load_primary[viral_load_primary$TRT01A == "PLACEBO",]
  viral_load_primary$AVAL <<- as.numeric(viral_load_primary$AVAL)
  # Organize by viral load (high to low)
  viral_load_primary <<- viral_load_primary[order(viral_load_primary$AVAL, decreasing = TRUE),]
  # Take gene_id column from gene_counts and use contents as the rownames of gene_counts
  gene.row.names <<- as.character(gene_counts$gene_id)
  gene_counts <<- gene_counts[,2:ncol(gene_counts)]
  gene_counts <<- as.data.frame(gene_counts)
  rownames(gene_counts) <<- gene.row.names
  # Add period into time point (by itself, time point isn't unique - we need period information
  # to distinguish between D-1 in period 1 vs D-1 in period 2, for example)
  all_metadata$time_point <<- paste0(all_metadata$period, "_", all_metadata$time_point)
  # Make time point names safe for DESeq2
  all_metadata$time_point[all_metadata$time_point == '1_D1 predose'] <<- '1_D_minus_1'
  all_metadata$time_point[all_metadata$time_point == '2_D-2'] <<- '2_D_minus_2'
  all_metadata$time_point[all_metadata$time_point == '2_D-1'] <<- '2_D_minus_1'
  # Order factor levels for period 1, period 2, and all time points
  period_1_factors <<- c("1_D_minus_1", "1_D2", "1_D8", "1_D28")
  period_2_factors <<- c("2_D_minus_2", "2_D_minus_1", "2_D2", "2_D5", "2_D8", "2_D28")
  all_factors <<- c(period_1_factors, period_2_factors)
  # Set factor levels
  all_metadata$subject_id <<- as.factor(all_metadata$subject_id)
  all_metadata$time_point <<- factor(all_metadata$time_point, levels = all_factors)
  all_metadata$sex <<- as.factor(all_metadata$sex)
  all_metadata$age <<- as.factor(all_metadata$age)
  # Only keep metadata for bulk RNA-seq aliquots
  bulk_metadata <<- all_metadata[all_metadata$bulkRNA_seq == TRUE,]
  # Divide metadata into placebo and vaccinated
  placebo_metadata <<- bulk_metadata[bulk_metadata$treatment == "PLACEBO",]
  vaccinated_metadata <<- bulk_metadata[bulk_metadata$treatment == "MVA-NP+M1",]
  # Find placebo-associated and vaccinated-associated gene_counts
  kept_aliquots <<- placebo_metadata$aliquot_id
  placebo_counts <<- gene_counts[kept_aliquots]
  kept_aliquots <<- vaccinated_metadata$aliquot_id
  vaccinated_counts <<- gene_counts[kept_aliquots]
  # Sort columns in gene_counts and rows for each so they're in same order (for DESeq2)
  # ALL 
  sorted_col_names <<- sort(colnames(gene_counts))
  gene_counts <<- gene_counts[, sorted_col_names] 
  rownames(bulk_metadata) <<- bulk_metadata$aliquot_id
  sorted_row_names <<- sort(rownames(bulk_metadata))
  bulk_metadata <<- bulk_metadata[sorted_row_names,]
  # PLACEBO
  sorted_col_names <<- sort(colnames(placebo_counts))
  placebo_counts <<- placebo_counts[, sorted_col_names]
  rownames(placebo_metadata) <<- placebo_metadata$aliquot_id
  sorted_row_names <<- sort(rownames(placebo_metadata))
  placebo_metadata <<- placebo_metadata[sorted_row_names,]
  # VACCINATED
  sorted_col_names <<- sort(colnames(vaccinated_counts))
  vaccinated_counts <<- vaccinated_counts[, sorted_col_names]
  rownames(vaccinated_metadata) <<- vaccinated_metadata$aliquot_id
  sorted_row_names <<- sort(rownames(vaccinated_metadata))
  vaccinated_metadata <<- vaccinated_metadata[sorted_row_names,]
  # Drop aliquot ID column (it's stored in rownames)
  bulk_metadata <<- subset(bulk_metadata, select = -c(aliquot_id))
  placebo_metadata <<- subset(placebo_metadata, select = -c(aliquot_id))
  vaccinated_metadata <<- subset(vaccinated_metadata, select = -c(aliquot_id))
  # Probably OK to round expected counts from RSEM data. DESeq2 expects integers
  gene_counts <<- round(gene_counts)
  placebo_counts <<- round(placebo_counts)
  vaccinated_counts <<- round(vaccinated_counts)
  # Grab the main 23 subjects (13 high viral load, 10 low viral load) that have all 10 time points
  placebo_full_time_series_metadata <<- placebo_metadata[placebo_metadata$subject_id 
                                                        %in% names(table(placebo_metadata$subject_id)
                                                                   [table(placebo_metadata$subject_id) == 10]),]
  # Grab subject IDs for main 23 subjects
  placebo_full_time_series_subjects <<- unique(placebo_full_time_series_metadata$subject_id)
  # Reorder subject IDs according to viral load (high to low)
  placebo_full_time_series_subjects <<- placebo_full_time_series_subjects[order(match(placebo_full_time_series_subjects,viral_load_primary$SUBJID))]
  # Top 13 will be high viral load and bottom 10 will be low viral load 
  high_viral_load_subjects <<- placebo_full_time_series_subjects[1:13]
  low_viral_load_subjects <<- tail(placebo_full_time_series_subjects, n = 10)
  # Grab high viral load placebo counts and metadata
  high_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$subject_id %in% high_viral_load_subjects,])
  high_placebo_counts <<- placebo_counts[,high_placebo_aliquots]
  high_placebo_metadata <<- placebo_metadata[high_placebo_aliquots,]
  # Grab low viral load placebo counts and metadata
  low_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$subject_id %in% low_viral_load_subjects,])
  low_placebo_counts <<- placebo_counts[,low_placebo_aliquots]
  low_placebo_metadata <<- placebo_metadata[low_placebo_aliquots,]
  # Grab combination of high and low viral load placebo counts and metadata
  both_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$subject_id %in% high_viral_load_subjects | placebo_metadata$subject_id %in% low_viral_load_subjects,])
  both_placebo_counts <<- placebo_counts[,both_placebo_aliquots]
  both_placebo_metadata <<- placebo_metadata[both_placebo_aliquots,]
  viral_load_for_metadata <<- both_placebo_metadata$subject_id %in% high_viral_load_subjects
  viral_load_for_metadata <<- replace(viral_load_for_metadata, viral_load_for_metadata == TRUE, "HIGH")
  viral_load_for_metadata <<- replace(viral_load_for_metadata, viral_load_for_metadata == "FALSE", "LOW")
  both_placebo_metadata$viral_load <<- viral_load_for_metadata
  # Remove questionable low viral load individual
  removed_low_viral_aliquots <- rownames(placebo_metadata[placebo_metadata$subject_id == "f18c54d93cef4a4e",])
  placebo_metadata <<- placebo_metadata[!(placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
  placebo_counts <<- placebo_counts[,!(colnames(placebo_counts) %in% removed_low_viral_aliquots)]
  low_placebo_metadata <<- low_placebo_metadata[!(low_placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
  low_placebo_counts <<- low_placebo_counts[,!(colnames(low_placebo_counts) %in% removed_low_viral_aliquots)]
}

run_deseq_bulk_analysis_time_series=function(sample_type, counts, metadata, test_time, baseline_time, output_dir, output_name_prefix=NA) {
  # Select the two relevant time points from our metadata
  metadata_subset <- metadata[metadata$time_point == test_time | metadata$time_point == baseline_time,]
  # Remove subjects that only have one time point (not both)
  metadata_subset <- metadata_subset[metadata_subset$subject_id  %in% names(table(metadata_subset$subject_id)[table(metadata_subset$subject_id) == 2]),]
  # Select subset of counts associated with subjects
  counts_subset <- counts[rownames(metadata_subset)]
  # Run DESeq2
  current_analysis <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + time_point)
  current_analysis <- DESeq(current_analysis)
  # Grab results with alpha = 0.05 and lfcThreshold = 0.1
  current_analysis_results <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 0.1)
  current_analysis_results <- current_analysis_results[order(current_analysis_results$padj),]
  current_analysis_results <- subset(current_analysis_results, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 0.585 (1.5 fold increase)
  current_analysis_results_1.5 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 0.585)
  current_analysis_results_1.5 <- current_analysis_results_1.5[order(current_analysis_results_1.5$padj),]
  current_analysis_results_1.5 <- subset(current_analysis_results_1.5, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)   
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 1
  current_analysis_results_2 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 1)
  current_analysis_results_2 <- current_analysis_results_2[order(current_analysis_results_2$padj),]
  current_analysis_results_2 <- subset(current_analysis_results_2, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 2
  current_analysis_results_4 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 2)
  current_analysis_results_4 <- current_analysis_results_4[order(current_analysis_results_4$padj),]
  current_analysis_results_4 <- subset(current_analysis_results_4, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_4), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_4), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with alpha = 0.05 and no lfcThreshold set
  current_analysis_results_none <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.05)
  current_analysis_results_none <- current_analysis_results_none[order(current_analysis_results_none$padj),]
  current_analysis_results_none <- subset(current_analysis_results_none, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  return(list(current_analysis_results, current_analysis_results_1.5, current_analysis_results_2, current_analysis_results_4, current_analysis_results_none))
}

run_deseq_bulk_analysis_viral_load=function(sample_type, counts, metadata, test_time, test_cond, baseline_cond, output_dir, output_name_prefix=NA) {
  # Select the relevant time point from our metadata
  metadata_subset <- metadata[metadata$time_point == test_time,]  
  # Select subset of counts associated with subjects
  counts_subset <- counts[rownames(metadata_subset)]
  # Run DESeq2
  current_analysis <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ sex + age + viral_load)
  current_analysis <- DESeq(current_analysis)
  # Grab results with alpha = 0.05 and lfcThreshold = 0.1
  current_analysis_results <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05, lfcThreshold = 0.1)
  current_analysis_results <- current_analysis_results[order(current_analysis_results$padj),]
  current_analysis_results <- subset(current_analysis_results, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 0.585 (1.5 fold increase)
  current_analysis_results_1.5 <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05, lfcThreshold = 0.585)
  current_analysis_results_1.5 <- current_analysis_results_1.5[order(current_analysis_results_1.5$padj),]
  current_analysis_results_1.5 <- subset(current_analysis_results_1.5, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)   
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 1
  current_analysis_results_2 <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05, lfcThreshold = 1)
  current_analysis_results_2 <- current_analysis_results_2[order(current_analysis_results_2$padj),]
  current_analysis_results_2 <- subset(current_analysis_results_2, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 2
  current_analysis_results_4 <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05, lfcThreshold = 2)
  current_analysis_results_4 <- current_analysis_results_4[order(current_analysis_results_4$padj),]
  current_analysis_results_4 <- subset(current_analysis_results_4, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_4), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_4), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  return(list(current_analysis_results, current_analysis_results_1.5, current_analysis_results_2, current_analysis_results_4))
}

run_deseq2_LRT <- function(counts, metadata) {
  LRT_metadata <- metadata[metadata$time_point == "2_D28" | metadata$time_point == "2_D8" | 
                             metadata$time_point == "2_D5" | metadata$time_point == "2_D2" |
                             metadata$time_point == "2_D_minus_1",]
  LRT_metadata <- LRT_metadata[LRT_metadata$subject_id  %in% names(table(LRT_metadata$subject_id)[table(LRT_metadata$subject_id) == 5]),]
  LRT_counts <- counts[rownames(LRT_metadata)]
  LRT_analysis <- DESeqDataSetFromMatrix(countData = LRT_counts,
                                                               colData = LRT_metadata,
                                                               design = ~ subject_id + time_point)
  LRT_analysis <- DESeq(LRT_analysis, test = "LRT", reduced = ~ subject_id)
  LRT_analysis_results <- results(LRT_analysis, alpha = 0.05)
  LRT_analysis_results <- LRT_analysis_results[order(LRT_analysis_results$padj),]
  LRT_analysis_results <- subset(LRT_analysis_results, padj < 0.05)
  return(list(LRT_analysis, LRT_analysis_results))
}

plot_lrt_heatmap <- function(gene_list, LRT_analysis, heatmap_threshold, output_path) {
  passing_genes <- c()
  for(gene in gene_list) {
    if(gene %in% rownames(LRT_analysis[[2]])) {
      passing_genes <- c(passing_genes, gene)
    }
  }
  
  betas <- coef(LRT_analysis[[1]])
  betas <- betas[, -c(1:13)] # TODO: Remove hard-coding of column numbers
  betas <- betas[rownames(betas) %in% passing_genes,]
  colnames(betas) <- c("Day 2 vs Day -1", "Day 5 vs Day -1", "Day 8 vs Day -1", "Day 28 vs Day -1")
  pheatmap(betas, breaks=seq(from=-heatmap_threshold, to=heatmap_threshold, length=101),
           cluster_col=FALSE, fontsize = 30, width = 16, height = 12, filename = output_path, cex = 0.7)
}

find_degs_across_time_points_for_gene_list <- function(D2_df, D5_df, D8_df, D28_df, overall_df, gene_df) {
  for(gene in unique(gene_df$Gene_Name)) {
    gene_vector <- c(gene)
    specific_gene_df <- gene_df[gene_df$Gene_Name == gene,]
    cell_types <- specific_gene_df$Cell_Type
    cell_types <- paste0(cell_types, collapse = ", ")
    gene_vector <- c(gene_vector, cell_types)
    # D2
    if(gene %in% rownames(D2_df)) {
      gene_vector <- c(gene_vector, as.numeric(D2_df[rownames(D2_df) == gene,]$log2FoldChange))
    } else {
      gene_vector <- c(gene_vector, 0)
    }
    # D5
    if(gene %in% rownames(D5_df)) {
      gene_vector <- c(gene_vector, as.numeric(D5_df[rownames(D5_df) == gene,]$log2FoldChange))
    } else {
      gene_vector <- c(gene_vector, 0)
    }
    # D8
    if(gene %in% rownames(D8_df)) {
      gene_vector <- c(gene_vector, as.numeric(D8_df[rownames(D8_df) == gene,]$log2FoldChange))
    } else {
      gene_vector <- c(gene_vector, 0)
    }
    # D28
    if(gene %in% rownames(D28_df)) {
      gene_vector <- c(gene_vector, as.numeric(D28_df[rownames(D28_df) == gene,]$log2FoldChange))
    } else {
      gene_vector <- c(gene_vector, 0)
    }
    gene_vector <- as.data.frame(t(gene_vector))
    names(gene_vector) <- c("gene", "cell_types", "D2_fc", "D5_fc", "D8_fc", "D28_fc")
    overall_df <- rbind(overall_df, gene_vector)
  }
  overall_df$D2_fc <- as.numeric(overall_df$D2_fc) 
  overall_df$D5_fc <- as.numeric(overall_df$D5_fc) 
  overall_df$D8_fc <- as.numeric(overall_df$D8_fc) 
  overall_df$D28_fc <- as.numeric(overall_df$D28_fc) 
  return(overall_df)
}

fill_in_special_notes <- function(gene_df, viral_load = "HVL") {
  special_notes_vec <- c()
  for(gene in gene_df$gene) {
    specific_gene_df <- gene_df[gene_df$gene == gene,]
    specific_special_notes <- ""
    specific_special_notes <- add_day_fc_info(specific_special_notes, specific_gene_df, "D2", viral_load)
    specific_special_notes <- add_day_fc_info(specific_special_notes, specific_gene_df, "D5", viral_load)
    specific_special_notes <- add_day_fc_info(specific_special_notes, specific_gene_df, "D8", viral_load)
    specific_special_notes <- add_day_fc_info(specific_special_notes, specific_gene_df, "D28", viral_load)
    special_notes_vec <- c(special_notes_vec, specific_special_notes)
  }
  gene_df$special_notes <- special_notes_vec
  return(gene_df)
}

add_day_fc_info <- function(special_notes, gene_df, day, viral_load = "HVL") {
  # Positive FC
  if(gene_df[[paste0(day, "_fc")]] > 2) {
    special_notes <- paste0(special_notes, "FC > 2 for ", day, " for ", viral_load, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] > 1) {
    special_notes <- paste0(special_notes, "FC > 1 for ", day, " for ", viral_load, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] > 0.585) {
    special_notes <- paste0(special_notes, "FC > 0.585 for ", day, " for ", viral_load, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] > 0) {
    special_notes <- paste0(special_notes, "FC > 0 for ", day, " for ", viral_load, ". ")
  }
  # Negative FC
  if(gene_df[[paste0(day, "_fc")]] < -2) {
    special_notes <- paste0(special_notes, "FC < -2 for ", day, " for ", viral_load, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] < -1) {
    special_notes <- paste0(special_notes, "FC < -1 for ", day, " for ", viral_load, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] < -0.585) {
    special_notes <- paste0(special_notes, "FC < -0.585 for ", day, " for ", viral_load, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] < 0) {
    special_notes <- paste0(special_notes, "FC < 0 for ", day, " for ", viral_load, ". ")
  }
  return(special_notes)
}