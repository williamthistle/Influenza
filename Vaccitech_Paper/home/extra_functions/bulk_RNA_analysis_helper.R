# Method to set up bulk analysis
setup_bulk_analysis=function(metadata_dir, data_dir) {
  # Read in count and metadata files
  gene_counts <<- fread(paste0(data_dir, "rsem_genes_count.processed.txt"), header = T, sep = "\t")
  all_metadata_file <<- paste0(metadata_dir, "all_metadata_sheet.tsv")
  all_metadata <<- read.table(all_metadata_file, header = TRUE, sep = "\t")
  # Currently, only capturing viral load for placebo subjects
  viral_load_file <<- paste0(metadata_dir, "bulk_RNA_viral_load.tsv")
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
  # Grab the 14 vaccinated subjects (7 high T cell response, 7 low T cell response) that have all 10 time points
  vaccinated_full_time_series_metadata <<- vaccinated_metadata[vaccinated_metadata$subject_id 
                                                               %in% names(table(vaccinated_metadata$subject_id)
                                                                          [table(vaccinated_metadata$subject_id) == 10]),]
  ### VACCINATED ###
  immunogenicity_data <<- read.table(paste0(metadata_dir, "Immunogenicity_Data.tsv"), sep = "\t", header = TRUE)
  immunogenicity_data <<- immunogenicity_data[immunogenicity_data$SUBJID %in% vaccinated_full_time_series_metadata$subject_id,]
  immunogenicity_data <<- immunogenicity_data[,c(1,105,106,115,116)]
  immunogenicity_data <<- immunogenicity_data[order(immunogenicity_data$Vaccination.Day8_IFNg_NP.Background_SFC.10.6.cells, decreasing = TRUE),]
  high_t_cell_response_subjects <<- immunogenicity_data$SUBJID[1:7]
  low_t_cell_response_subjects <<- immunogenicity_data$SUBJID[8:14]
  # Grab high t cell response vaccinated counts and metadata
  high_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$subject_id %in% high_t_cell_response_subjects,])
  high_vaccinated_counts <<- vaccinated_counts[,high_vaccinated_aliquots]
  high_vaccinated_metadata <<- vaccinated_metadata[high_vaccinated_aliquots,]
  # Grab low viral load placebo counts and metadata
  low_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$subject_id %in% low_t_cell_response_subjects,])
  low_vaccinated_counts <<- vaccinated_counts[,low_vaccinated_aliquots]
  low_vaccinated_metadata <<- vaccinated_metadata[low_vaccinated_aliquots,]
  # Grab combination of high and low viral load placebo counts and metadata
  both_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$subject_id %in% high_t_cell_response_subjects | vaccinated_metadata$subject_id %in% low_t_cell_response_subjects,])
  both_vaccinated_counts <<- vaccinated_counts[,both_vaccinated_aliquots]
  both_vaccinated_metadata <<- vaccinated_metadata[both_vaccinated_aliquots,]
  t_cell_for_metadata <<- both_vaccinated_metadata$subject_id %in% high_t_cell_response_subjects
  t_cell_for_metadata <<- replace(t_cell_for_metadata, t_cell_for_metadata == TRUE, "HIGH")
  t_cell_for_metadata <<- replace(t_cell_for_metadata, t_cell_for_metadata == "FALSE", "LOW")
  both_vaccinated_metadata$t_cell_response <<- t_cell_for_metadata
  ### PLACEBO ###
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
  both_placebo_metadata <<-  both_placebo_metadata[!(both_placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
  both_placebo_counts <<- both_placebo_counts[,!(colnames(both_placebo_counts) %in% removed_low_viral_aliquots)]
  low_placebo_metadata <<- low_placebo_metadata[!(low_placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
  low_placebo_counts <<- low_placebo_counts[,!(colnames(low_placebo_counts) %in% removed_low_viral_aliquots)]
}

run_deseq_bulk_analysis_time_series=function(sample_type, counts, metadata, test_time, baseline_time, output_dir, output_name_prefix=NA, alpha = 0.05) {
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  # Select the two relevant time points from our metadata
  metadata_subset <- metadata[metadata$time_point == test_time | metadata$time_point == baseline_time,]
  # Remove subjects that only have one time point (not both)
  metadata_subset <- metadata_subset[metadata_subset$subject_id  %in% names(table(metadata_subset$subject_id)[table(metadata_subset$subject_id) == 2]),]
  # Select subset of counts associated with subjects
  counts_subset <- counts[rownames(metadata_subset)]
  # Run DESeq2
  current_analysis <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + time_point)
  current_analysis <- DESeq(current_analysis)
  # Grab results with lfcThreshold = 0.1
  current_analysis_results <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 0.1)
  current_analysis_results <- current_analysis_results[order(current_analysis_results$padj),]
  current_analysis_results <- subset(current_analysis_results, padj < alpha)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with lfcThreshold = 0.585 (1.5 fold increase)
  current_analysis_results_1.5 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 0.585)
  current_analysis_results_1.5 <- current_analysis_results_1.5[order(current_analysis_results_1.5$padj),]
  current_analysis_results_1.5 <- subset(current_analysis_results_1.5, padj < alpha)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)   
  }
  # Grab results with lfcThreshold = 1
  current_analysis_results_2 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 1)
  current_analysis_results_2 <- current_analysis_results_2[order(current_analysis_results_2$padj),]
  current_analysis_results_2 <- subset(current_analysis_results_2, padj < alpha)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with lfcThreshold = 2
  current_analysis_results_4 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 2)
  current_analysis_results_4 <- current_analysis_results_4[order(current_analysis_results_4$padj),]
  current_analysis_results_4 <- subset(current_analysis_results_4, padj < alpha)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_4), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_4), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with no lfcThreshold set
  current_analysis_results_none <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha)
  current_analysis_results_none <- current_analysis_results_none[order(current_analysis_results_none$padj),]
  current_analysis_results_none <- subset(current_analysis_results_none, padj < alpha)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  return(list(current_analysis_results, current_analysis_results_1.5, current_analysis_results_2, current_analysis_results_4, current_analysis_results_none))
}

run_deseq_bulk_analysis_viral_load=function(sample_type, counts, metadata, test_time, test_cond, baseline_cond, output_dir, output_name_prefix=NA) {
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
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
  # Grab results with alpha = 0.05 and no lfcThreshold
  current_analysis_results_none <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05)
  current_analysis_results_none <- current_analysis_results_none[order(current_analysis_results_none$padj),]
  current_analysis_results_none <- subset(current_analysis_results_none, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  return(list(current_analysis_results, current_analysis_results_1.5, current_analysis_results_2, current_analysis_results_4, current_analysis_results_none))
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
    cell_types <- paste0(cell_types, collapse = ",")
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

fill_in_special_notes <- function(gene_df) {
  special_notes_vec <- c()
  for(gene in gene_df$gene) {
    specific_gene_df <- gene_df[gene_df$gene == gene,]
    specific_special_notes <- ""
    specific_special_notes <- add_day_fc_info(specific_special_notes, specific_gene_df, "D2")
    specific_special_notes <- add_day_fc_info(specific_special_notes, specific_gene_df, "D5")
    specific_special_notes <- add_day_fc_info(specific_special_notes, specific_gene_df, "D8")
    specific_special_notes <- add_day_fc_info(specific_special_notes, specific_gene_df, "D28")
    special_notes_vec <- c(special_notes_vec, specific_special_notes)
  }
  gene_df$special_notes <- special_notes_vec
  return(gene_df)
}

add_day_fc_info <- function(special_notes, gene_df, day) {
  # Positive FC
  if(gene_df[[paste0(day, "_fc")]] > 2) {
    special_notes <- paste0(special_notes, "FC > 2 for ", day, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] > 1) {
    special_notes <- paste0(special_notes, "FC > 1 for ", day, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] > 0.585) {
    special_notes <- paste0(special_notes, "FC > 0.585 for ", day, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] > 0) {
    special_notes <- paste0(special_notes, "FC > 0 for ", day, ". ")
  }
  # Negative FC
  if(gene_df[[paste0(day, "_fc")]] < -2) {
    special_notes <- paste0(special_notes, "FC < -2 for ", day, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] < -1) {
    special_notes <- paste0(special_notes, "FC < -1 for ", day, ". ")
  } else if(gene_df[[paste0(day, "_fc")]] < -0.585) {
    special_notes <- paste0(special_notes, "FC < -0.585 for ", ". ")
  } else if(gene_df[[paste0(day, "_fc")]] < 0) {
    special_notes <- paste0(special_notes, "FC < 0 for ", day, ". ")
  }
  return(special_notes)
}

find_overlapping_motifs_between_atac_and_rna <- function(peak_dir, sc_gene_table, cell_types, pos_bulk_genes = NULL, neg_bulk_genes = NULL) {
  overlapping_motif_df <- data.frame(tf = character(), cell_types = character(), found_in_bulk = logical())
  for(cell_type in cell_types) {
    print(cell_type)
    if(file.exists(paste0(peak_dir, cell_type, "_D28_D1_motif_up_pseudobulk_corrected.tsv"))) {
      # Up-regulated enriched TFs
      current_motifs_up <- read.table(paste0(peak_dir, cell_type, "_D28_D1_motif_up_pseudobulk_corrected.tsv"), sep = "\t", header = TRUE)
      current_motifs_up <- current_motifs_up[current_motifs_up$p_adj < 0.05,]
      current_tfs_up <- current_motifs_up$TF
      # Cut off _ID so we can match with gene names from RNA
      current_tfs_up <- sub("_.*", "", current_tfs_up)
      # Down-regulated enriched TFs
      current_motifs_down <- read.table(paste0(peak_dir, cell_type, "_D28_D1_motif_down_pseudobulk_corrected.tsv"), sep = "\t", header = TRUE)
      current_motifs_down <- current_motifs_down[current_motifs_down$p_adj < 0.05,]
      current_tfs_down <- current_motifs_down$TF
      # Cut off _ID so we can match with gene names from RNA
      current_tfs_down <- sub("_.*", "", current_tfs_down)
      # Grab genes from SC table for current cell type
      cell_type_in_df <- sub("_", " ", cell_type)
      cell_type_sc_gene_table <- sc_gene_table[sc_gene_table$Cell_Type == cell_type_in_df,]
      # Grab overlap between pos SC genes and upregulated TF motifs
      overlapping_pos_tfs <- intersect(current_tfs_up, cell_type_sc_gene_table[cell_type_sc_gene_table$sc_log2FC > 0,]$Gene_Name)
      # Grab overlap between neg SC genes and downregulated TF motifs
      overlapping_neg_tfs <- intersect(current_tfs_down, cell_type_sc_gene_table[cell_type_sc_gene_table$sc_log2FC < 0,]$Gene_Name)
      for(pos_tf in overlapping_pos_tfs) {
        # If the TF is already present in our DF, then just add the new cell type to cell_types
        if(pos_tf %in% overlapping_motif_df$tf) {
          current_cell_types <- overlapping_motif_df[overlapping_motif_df$tf == pos_tf,]$cell_types
          current_cell_types <- paste0(current_cell_types, ", ", cell_type_in_df, " (Up)")
          overlapping_motif_df[overlapping_motif_df$tf == pos_tf,]$cell_types <- current_cell_types
        } else {
          # If the tf was found in the bulk data, note that
          found_in_bulk <- FALSE
          if(pos_tf %in% pos_bulk_genes) {
            found_in_bulk <- TRUE
          }
          # Add tf to df (and mark it as up)
          motif_vector <- c(pos_tf, paste0(cell_type_in_df, " (Up)"), found_in_bulk)
          motif_vector <- as.data.frame(t(motif_vector))
          names(motif_vector) <- c("tf", "cell_types", "found_in_bulk")
          overlapping_motif_df <- rbind(overlapping_motif_df, motif_vector)
        }
      }
      for(neg_tf in overlapping_neg_tfs) {
        # If the TF is already present in our DF, then just add the new cell type to cell_types
        if(neg_tf %in% overlapping_motif_df$tf) {
          current_cell_types <- overlapping_motif_df[overlapping_motif_df$tf == neg_tf,]$cell_types
          current_cell_types <- paste0(current_cell_types, ",", cell_type_in_df, " (Down)")
          overlapping_motif_df[overlapping_motif_df$tf == neg_tf,]$cell_types <- current_cell_types
        } else {
          # If the tf was found in the bulk data, note that
          found_in_bulk <- FALSE
          if(neg_tf %in% neg_bulk_genes) {
            found_in_bulk <- TRUE
          }
          # Add tf to df (and mark it as down)
          motif_vector <- c(neg_tf, paste0(cell_type_in_df, " (Down)"), found_in_bulk)
          motif_vector <- as.data.frame(t(motif_vector))
          names(motif_vector) <- c("tf", "cell_types", "found_in_bulk")
          overlapping_motif_df <- rbind(overlapping_motif_df, motif_vector)
        }
      }
    }
  }
  return(overlapping_motif_df)
}

fill_in_info_for_magical_output <- function(overall_magical_df, das_dir, sc_pseudobulk_deg_combined_cell_types_table, sc_pseudobulk_deg_table,
                                            pos_bulk_genes, neg_bulk_genes, sc_peaks, overall_pseudobulk_motif_enrichment_df) {
  # Distance between relevant gene TSS and peak 
  dist_between_gene_tss_and_peak <- c()
  # FC for circuit gene from SC data (positive or negative?)
  gene_fcs <- c()
  # FC for circuit site from SC (pseudobulk) data (positive or negative?)
  site_fcs <- c()
  # Cell types where a circuit was found with the same gene
  magical_cell_types <- c()
  # Cell types where the gene was found to be a DEG (combined cell types used as input for MAGICAL)
  combined_single_cell_types <- c()
  # Cell types where the gene was found to be a DEG (more granular cell types from the original DEG analysis)
  original_single_cell_types <- c()
  # Was the gene found to be a DEG in bulk data?
  # If yes, note that its FC direction always matches the single cell FC direction 
  bulk_boolean <- c()
  # What kind of peak is the site?
  # What is the closest gene to the site?
  site_type <- c()
  nearest_gene_to_site <- c()
  #dist_to_nearest_gene_tss <- c()
  #dist_to_nearest_gene_start <- c()
  pseudobulk_motif_enrichment <- c()
  for(current_circuit_index in 1:nrow(overall_magical_df)) {
    current_circuit <- overall_magical_df[current_circuit_index,]
    cell_type <- current_circuit$Cell_Type
    cell_type_in_df <- sub("_", " ", cell_type)
    current_gene <- current_circuit$Gene_symbol
    current_tfs <- current_circuit$TFs.binding.prob.
    current_tfs <- unlist(strsplit(current_tfs, ","))
    current_tfs <- current_tfs[c(TRUE, FALSE)]
    all_pseudobulk_info <- ""
    for(current_tf in current_tfs) {
      pseudobulk_enrichment_subset_df <- overall_pseudobulk_motif_enrichment_df[overall_pseudobulk_motif_enrichment_df$TF == current_tf,]
      pseudobulk_enrichment_subset_df <- pseudobulk_enrichment_subset_df[pseudobulk_enrichment_subset_df$Cell_Type == cell_type,]
      if(nrow(pseudobulk_enrichment_subset_df) > 0) {
        current_pseudobulk_info <- paste0(current_tf, ": ")
        current_pseudobulk_tf_info <- paste(pseudobulk_enrichment_subset_df$Direction, " (Rank ", pseudobulk_enrichment_subset_df$rank, ")", sep = "")
        current_pseudobulk_tf_info <- paste(current_pseudobulk_tf_info, collapse = ",")
        current_pseudobulk_info <- paste0(current_pseudobulk_info, current_pseudobulk_tf_info, ".")
      } else {
        current_pseudobulk_info <- paste0(current_tf, ": N/A.")
      }
      all_pseudobulk_info <- paste0(all_pseudobulk_info, current_pseudobulk_info)
    }
    pseudobulk_motif_enrichment <- c(pseudobulk_motif_enrichment, all_pseudobulk_info)
    # Capture distance between gene TSS and peak
    gene_tss <- current_circuit$Gene_TSS
    peak_start <- current_circuit$Peak_start
    peak_end <- current_circuit$Peak_end
    peak_mid <- (peak_start + peak_end) / 2
    dist_between_gene_tss_and_peak <- c(dist_between_gene_tss_and_peak, abs(gene_tss - peak_mid))
    # Capture gene FC
    gene_fc <- sc_pseudobulk_deg_combined_cell_types_table[sc_pseudobulk_deg_combined_cell_types_table$Gene_Name == current_circuit$Gene_symbol,]
    gene_fc <- gene_fc[gene_fc$Cell_Type == cell_type_in_df,]$sc_log2FC
    gene_fcs <- c(gene_fcs, gene_fc)
    # Capture site FC
    site_fc <- read.table(paste0(sc_das_dir, "diff_peaks/", sub(" ", "_", cell_type), "_D28_D1_diff_pseudo_filtered.tsv"), sep = "\t",
                                                      header = TRUE)
    site_fc <- site_fc[site_fc$chr == current_circuit$Peak_chr & site_fc$start == current_circuit$Peak_start & site_fc$end == current_circuit$Peak_end,]
    site_fc <- site_fc$log2FoldChange
    site_fcs <- c(site_fcs, site_fc)
    # Capture cell types from MAGICAL circuits for current gene
    magical_gene_cell_types <- overall_magical_df[overall_magical_df$Gene_symbol == current_gene,]$Cell_Type
    magical_gene_cell_types <- sort(unique(magical_gene_cell_types))
    magical_gene_cell_types <- paste(magical_gene_cell_types, collapse = ",")
    magical_cell_types <- c(magical_cell_types, magical_gene_cell_types)
    # Capture cell types from DEGs for combined SC cell types
    combined_gene_cell_types <- sc_pseudobulk_deg_combined_cell_types_table[sc_pseudobulk_deg_combined_cell_types_table$Gene_Name == current_gene,]$Cell_Type
    combined_gene_cell_types <- sort(combined_gene_cell_types)
    combined_gene_cell_types <- paste(combined_gene_cell_types, collapse = ",")
    combined_single_cell_types <- c(combined_single_cell_types, combined_gene_cell_types)
    # Capture cell types from DEGS for original more granular SC cell types
    original_gene_cell_types <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Gene_Name == current_gene,]$Cell_Type
    if(length(original_gene_cell_types) > 0) {
      original_gene_cell_types <- sort(original_gene_cell_types)
      original_gene_cell_types <- paste(original_gene_cell_types, collapse = ",")
      original_single_cell_types <- c(original_single_cell_types, original_gene_cell_types)
    } else {
      original_single_cell_types <- c(original_single_cell_types, "N/A")
    }
    # Capture whether gene was found in bulk RNA-seq data
    if(current_gene %in% pos_bulk_genes | current_gene %in% neg_bulk_genes) {
      bulk_boolean <- c(bulk_boolean, TRUE)
    } else {
      bulk_boolean <- c(bulk_boolean, FALSE)
    }
    current_peak_info <- sc_peaks[sc_peaks$value == current_circuit$Peak_chr & sc_peaks$start == current_circuit$Peak_start & sc_peaks$end == current_circuit$Peak_end,]
    site_type <- c(site_type, current_peak_info$peakType)
    nearest_gene_to_site <- c(nearest_gene_to_site, current_peak_info$nearestGene)
    #dist_to_nearest_gene_tss <- c(dist_to_nearest_gene_tss, current_peak_info$distToTSS)
    #dist_to_nearest_gene_start <- c(dist_to_nearest_gene_start, current_peak_info$distToGeneStart)
  }
  overall_magical_df$dist_between_gene_tss_and_peak <- dist_between_gene_tss_and_peak
  overall_magical_df$gene_fc <- gene_fcs
  overall_magical_df$site_fc <- site_fcs
  overall_magical_df$nearest_gene_to_site <- nearest_gene_to_site
  #overall_magical_df$dist_to_nearest_gene_start <- dist_to_nearest_gene_start
  #overall_magical_df$dist_to_nearest_gene_tss <- dist_to_nearest_gene_tss
  overall_magical_df$site_type_for_nearest_gene <- site_type
  overall_magical_df$bulk <- bulk_boolean
  overall_magical_df$magical_cell_types <- magical_cell_types
  overall_magical_df$combined_single_cell_types <- combined_single_cell_types
  overall_magical_df$original_single_cell_types <- original_single_cell_types
  overall_magical_df$pseudobulk_motif_enrichment_for_tf <- pseudobulk_motif_enrichment
  return(overall_magical_df)
}

fill_in_info_for_magical_tf_output <- function(overall_magical_df, overall_pseudobulk_motif_enrichment_df) {
  # TF 
  tf <- c()
  # Total Circuit Count
  total_circuit_count <- c()
  # Circuit Cell Types (Count)
  circuit_cell_types <- c()
  # Bound_Genes (Cell Type)
  bound_genes <- c()
  # Found_as_Circuit_Genes (Cell Type)
  found_as_circuit_genes <- c()
  # Found_as_DEGs_in_Combined_Cell_Types (Cell Type and Direction)
  found_as_DEGs_in_combined_cell_types <- c()
  # Found_as_DEGs_in_Original_SC_Cell_Types (Cell Type and Direction)
  found_as_DEGs_in_original_cell_types <- c()
  # Found_as_DEGs_in_Bulk (Direction)
  found_as_DEGs_in_bulk <- c()
  # Pseudobulk_Motif_Enrichment (Cell Type and Direction)
  pseudobulk_motif_enrichment <- c()
  # Grab all TFs found in MAGICAL circuits
  for(current_circuit_index in 1:nrow(overall_magical_df)) {
    current_circuit <- overall_magical_df[current_circuit_index,]
    current_tfs <- current_circuit$TFs.binding.prob.
    current_tfs <- unlist(strsplit(current_tfs, ","))
    current_tfs <- current_tfs[c(TRUE, FALSE)]
    tf <- c(tf, current_tfs)
  }
  # Sort TFs by name and find total count
  tf_summary <- table(tf)
  tf <- names(tf_summary)
  total_circuit_count <- as.numeric(unname(tf_summary))
  # All bulk passing genes
  passing_bulk_genes <- c(high_passing_pos_genes, high_passing_neg_genes)
  for(current_tf in tf) {
    # Find cell types for circuits that have current TF
    tf_rows <- overall_magical_df[grepl(current_tf, overall_magical_df$TFs.binding.prob.),]
    current_circuit_cell_types <- sort(unique(tf_rows$Cell_Type))
    current_circuit_cell_types <- paste0(current_circuit_cell_types, collapse = ",")
    circuit_cell_types <- c(circuit_cell_types, current_circuit_cell_types)
    # Find genes that are bound by current TF in MAGICAL circuits (and record cell types for each gene)
    current_bound_genes <- sort(unique(tf_rows$Gene_symbol))
    current_bound_genes <- paste0(current_bound_genes, collapse = ",")
    bound_genes <- c(bound_genes, current_bound_genes)
    # If the current TF is also found as a circuit gene, record that info (gene name and cell type)
    # If the current TF is not found as a circuit gene, just write N/A
    found_as_circuit_genes_tf <- overall_magical_df[overall_magical_df$Gene_symbol == current_tf,]
    if(nrow(found_as_circuit_genes_tf) > 0) {
      current_found_as_circuit_gene_cell_types <- found_as_circuit_genes_tf$Cell_Type
      current_found_as_circuit_gene_cell_types <- unique(current_found_as_circuit_gene_cell_types)
      current_found_as_circuit_gene_cell_types <- paste(current_found_as_circuit_gene_cell_types, collapse = ",")
      found_as_circuit_genes <- c(found_as_circuit_genes, current_found_as_circuit_gene_cell_types)
    } else {
      found_as_circuit_genes <- c(found_as_circuit_genes, "N/A")
    }
    # If the current TF is found as a DEG in the combined cell type data (input for MAGICAL), record that info (cell type and FC direction)
    # If the current TF is not found, just write N/A
    found_as_combined_deg_tf <- sc_pseudobulk_deg_combined_cell_types_table[sc_pseudobulk_deg_combined_cell_types_table$Gene_Name == current_tf,]
    if(nrow(found_as_combined_deg_tf) > 0) {
      current_combined_deg_cell_types <- found_as_combined_deg_tf$Cell_Type
      current_combined_deg_fcs <- found_as_combined_deg_tf$sc_log2FC
      current_combined_deg_fcs <- ifelse(current_combined_deg_fcs > 0, "Pos", "Neg")
      current_deg_gene_info <- paste(current_combined_deg_cell_types, " (", current_combined_deg_fcs, ")", sep = "")
      current_deg_gene_info <- paste(current_deg_gene_info, collapse = ",")
      found_as_DEGs_in_combined_cell_types <- c(found_as_DEGs_in_combined_cell_types, current_deg_gene_info)
    } else {
      found_as_DEGs_in_combined_cell_types <- c(found_as_DEGs_in_combined_cell_types, "N/A")
    }
    # If the current TF is found as a DEG in the original sc data, record that info (cell type and FC direction)
    # If the current TF is not found, just write N/A
    found_as_original_deg_tf <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Gene_Name == current_tf,]
    if(nrow(found_as_original_deg_tf) > 0) {
      current_original_deg_cell_types <- found_as_original_deg_tf$Cell_Type
      current_original_deg_fcs <- found_as_original_deg_tf$sc_log2FC
      current_original_deg_fcs <- ifelse(current_original_deg_fcs > 0, "Pos", "Neg")
      current_deg_gene_info <- paste(current_original_deg_cell_types, " (", current_original_deg_fcs, ")", sep = "")
      current_deg_gene_info <- paste(current_deg_gene_info, collapse = ",")
      found_as_DEGs_in_original_cell_types <- c(found_as_DEGs_in_original_cell_types, current_deg_gene_info)
    } else {
      found_as_DEGs_in_original_cell_types <- c(found_as_DEGs_in_original_cell_types, "N/A")
    }
    # If the current TF is found as a DEG in the bulk data, record that info (FC direction)
    # If the current TF is not found, just write N/A
    if(current_tf %in% high_passing_pos_genes) {
      found_as_DEGs_in_bulk <- c(found_as_DEGs_in_bulk, "Pos")
    } else if(current_tf %in% high_passing_neg_genes) {
      found_as_DEGs_in_bulk <- c(found_as_DEGs_in_bulk, "Neg")
    } else {
      found_as_DEGs_in_bulk <- c(found_as_DEGs_in_bulk, "N/A")
    }
    # If the current TF is found as enriched in the pseudobulk motif analysis, record that info (cell type and direction)
    # If the current TF is not found, just write N/A
    current_pseudobulk_motif_enrichment_tf <- overall_pseudobulk_motif_enrichment_df[overall_pseudobulk_motif_enrichment_df$TF == current_tf,]
    if(nrow(current_pseudobulk_motif_enrichment_tf) > 0) {
      current_pseudobulk_motif_enrichment_cell_types <- current_pseudobulk_motif_enrichment_tf$Cell_Type
      current_pseudobulk_motif_enrichment_directions <- current_pseudobulk_motif_enrichment_tf$Direction
      current_pseudobulk_motif_enrichment_ranks <- current_pseudobulk_motif_enrichment_tf$rank
      current_pseudobulk_motif_enrichment_info <- paste(current_pseudobulk_motif_enrichment_cell_types, " (", current_pseudobulk_motif_enrichment_directions, "; Rank ", current_pseudobulk_motif_enrichment_ranks, ")", sep = "")
      current_pseudobulk_motif_enrichment_info <- paste(current_pseudobulk_motif_enrichment_info, collapse = ",")
      pseudobulk_motif_enrichment <- c(pseudobulk_motif_enrichment, current_pseudobulk_motif_enrichment_info)
    } else {
      pseudobulk_motif_enrichment <- c(pseudobulk_motif_enrichment, "N/A")
    }
  }
  overall_magical_tf_df <- data.frame(tf = tf, total_circuit_count = total_circuit_count, circuit_cell_types = circuit_cell_types,
                                      bound_genes = bound_genes, found_as_circuit_genes = found_as_circuit_genes,
                                      found_as_DEGs_in_combined_cell_types = found_as_DEGs_in_combined_cell_types,
                                      found_as_DEGs_in_original_cell_types = found_as_DEGs_in_original_cell_types,
                                      found_as_DEGs_in_bulk = found_as_DEGs_in_bulk,
                                      pseudobulk_motif_enrichment = pseudobulk_motif_enrichment)
  return(overall_magical_tf_df)
}






fill_in_sc_deg_info_for_time_series <- function(sc_gene_df, high_placebo_counts, high_placebo_metadata, output_dir, sc_fc_direction, alpha = 0.05) {
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  final_df <- data.frame(Gene = character(), Day = character(), Fold.Change.Abs = numeric(), 
                                         Fold.Change.Direction.Raw = character(), Fold.Change.Direction = character(), 
                                         Adjusted.P.Value = numeric())
  # Grab possible genes
  if(sc_fc_direction == "up") {
    possible_genes <- unique(sc_gene_df[sc_gene_df$sc_log2FC > 0,]$Gene_Name)
  } else {
    possible_genes <- unique(sc_gene_df[sc_gene_df$sc_log2FC < 0,]$Gene_Name)
  }
  
  # Find bulk RNA-seq validated genes
  high_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                             "2_D28", "2_D_minus_1", data_dir, "high", alpha = alpha)
  raw_high_placebo_period_2_D28_vs_D_minus_1_results <- high_placebo_period_2_D28_vs_D_minus_1_results[[5]]
  filtered_sc_trained_immunity_genes <- intersect(possible_genes,
                                                                   rownames(raw_high_placebo_period_2_D28_vs_D_minus_1_results))
  
  # Verify that FC in bulk agrees with FC in SC (if not, discard)
  final_filtered_sc_trained_immunity_genes <- c()
  
  for(gene in filtered_sc_trained_immunity_genes) {
    sc_fc <- sc_gene_df[sc_gene_df$Gene_Name == gene,]$sc_log2FC
    bulk_fc <- raw_high_placebo_period_2_D28_vs_D_minus_1_results[rownames(raw_high_placebo_period_2_D28_vs_D_minus_1_results) == gene,]$log2FoldChange
    if(all(sign(sc_fc) == sign(bulk_fc))) {
      final_filtered_sc_trained_immunity_genes <- c(final_filtered_sc_trained_immunity_genes, gene)
    }
  }
  
  # Compute results for D2 / D5 / D8 (with alpha parameter and without alpha parameter)
  # WITH ALPHA SET
  high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                       "2_D2", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D2_vs_D_minus_1_", alpha, "/"), "high", alpha = alpha)
  high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered[[5]]
  
  high_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                       "2_D5", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D5_vs_D_minus_1_", alpha, "/"), "high", alpha = alpha)
  high_placebo_period_2_D5_vs_D_minus_1_results <- high_placebo_period_2_D5_vs_D_minus_1_results[[5]]
  
  high_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                       "2_D8", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D8_vs_D_minus_1_", alpha, "/"), "high", alpha = alpha)
  high_placebo_period_2_D8_vs_D_minus_1_results <- high_placebo_period_2_D8_vs_D_minus_1_results[[5]]
  # WITHOUT ALPHA SET
  high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                       "2_D2", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D2_vs_D_minus_1_", alpha, "/"), "high", alpha = 0.99999)
  high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered[[5]]
  
  high_placebo_period_2_D5_vs_D_minus_1_results_unfiltered <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                       "2_D5", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D5_vs_D_minus_1_", alpha, "/"), "high", alpha = 0.99999)
  high_placebo_period_2_D5_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D5_vs_D_minus_1_results_unfiltered[[5]]
  
  high_placebo_period_2_D8_vs_D_minus_1_results_unfiltered <- run_deseq_bulk_analysis_time_series("placebo", high_placebo_counts, high_placebo_metadata,
                                                                                       "2_D8", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D8_vs_D_minus_1_", alpha, "/"), "high", alpha = 0.99999)
  high_placebo_period_2_D8_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D8_vs_D_minus_1_results_unfiltered[[5]]
  
  

}