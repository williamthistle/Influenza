# Method to set up bulk analysis
setup_bulk_rna_analysis=function(metadata_dir, data_dir) {
  # Read in count and metadata files
  gene_counts <<- fread(paste0(data_dir, "rsem_genes_count.processed.txt"), header = T, sep = "\t")
  gene_counts_normalized <<- fread(paste0(data_dir, "rsem_genes_count.processed.normalized.txt"), header = T, sep = "\t")
  gene_counts_normalized_without_scale <<- fread(paste0(data_dir, "rsem_genes_count.processed.normalized.without.scale.txt"), header = T, sep = "\t")
  
  gene_counts_normalized <<- as.data.frame(gene_counts_normalized)
  rownames(gene_counts_normalized) <<- gene_counts_normalized$gene_id
  gene_counts_normalized <<- gene_counts_normalized[,-c(1)]
  
  gene_counts_normalized_without_scale <<- as.data.frame(gene_counts_normalized_without_scale)
  rownames(gene_counts_normalized_without_scale) <<- gene_counts_normalized_without_scale$gene_id
  gene_counts_normalized_without_scale <<- gene_counts_normalized_without_scale[,-c(1)]
  
  all_metadata_file <<- paste0(metadata_dir, "all_metadata_sheet.tsv")
  all_metadata <<- read.table(all_metadata_file, header = TRUE, sep = "\t")
  # Read in viral load
  viral_load_file <<- paste0(metadata_dir, "bulk_RNA_viral_load.tsv")
  viral_load <<- read.table(viral_load_file, header = TRUE, sep = "\t")
  viral_load_primary <<- viral_load[viral_load$PARAMCD == "QPCRAUC",]
  viral_load_primary$AVAL <<- as.numeric(viral_load_primary$AVAL)
  # Read in immunogenicity data
  immunogenicity_data <<- read.table(paste0(metadata_dir, "Immunogenicity_Data.tsv"), sep = "\t", header = TRUE)
  immunogenicity_data$subject_id <<- immunogenicity_data$SUBJID
  # Read in NAI data
  nai_data <<- read.table(paste0(metadata_dir, "placebo_nai_recoded_processed.tsv"), sep = "\t", header = TRUE)
  # Organize by viral load (high to low)
  viral_load_primary <<- viral_load_primary[order(viral_load_primary$AVAL, decreasing = TRUE),]
  # Separate viral loads into vaccinated / placebo
  vaccinated_viral_load_primary <<- viral_load_primary[viral_load_primary$TRT01A == "MVA-NP+M1",]
  placebo_viral_load_primary <<- viral_load_primary[viral_load_primary$TRT01A == "PLACEBO",]
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
  
  # Add cell type proportion content
  rows_to_merge <- bulk_metadata$aliquot_id
  bulk_metadata <<- merge(bulk_metadata[bulk_metadata$aliquot_id %in% rows_to_merge, ], 
                          cibersort_cell_type_proportions[cibersort_cell_type_proportions$aliquot_id %in% rows_to_merge, ], 
                          by = "aliquot_id", 
                          all.x = TRUE, 
                          all.y = TRUE)
  
  # Divide metadata into placebo and vaccinated (and add qPCRAUC viral load values and t cell info for vaccinated)
  placebo_metadata <<- bulk_metadata[bulk_metadata$treatment == "PLACEBO",]
  names(placebo_viral_load_primary)[names(placebo_viral_load_primary) == 'SUBJID'] <<- 'subject_id'
  kept_columns <- colnames(placebo_metadata)
  kept_columns <- c(kept_columns, "AVAL")
  placebo_metadata <<- merge(placebo_metadata,  placebo_viral_load_primary, by = "subject_id")
  placebo_metadata <<- placebo_metadata[,kept_columns]
  vaccinated_metadata <<- bulk_metadata[bulk_metadata$treatment == "MVA-NP+M1",]
  names(vaccinated_viral_load_primary)[names(vaccinated_viral_load_primary) == 'SUBJID'] <<- 'subject_id'
  kept_columns <- colnames(vaccinated_metadata)
  kept_columns <- c(kept_columns, "AVAL")
  vaccinated_metadata <<- merge(vaccinated_metadata, vaccinated_viral_load_primary, by = "subject_id")
  vaccinated_metadata <<- vaccinated_metadata[,kept_columns]
  
  immunogenicity_data_t_cell_subset <<- immunogenicity_data[,c(1,105,106,115,116)]
  immunogenicity_data_t_cell_subset <<- immunogenicity_data_t_cell_subset[immunogenicity_data_t_cell_subset$SUBJID %in% vaccinated_metadata$subject_id,]
  immunogenicity_data_t_cell_subset$subject_id <<- immunogenicity_data_t_cell_subset$SUBJID
  immunogenicity_data_t_cell_subset <<- immunogenicity_data_t_cell_subset[,-c(1)]
  
  vaccinated_metadata <<- merge(vaccinated_metadata, 
                               immunogenicity_data_t_cell_subset, 
                               by = "subject_id")
  
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
  # Get high T cell and low T cell response subjects that were profiled with multiomics
  immunogenicity_data_full_time_series <<- immunogenicity_data[immunogenicity_data$SUBJID %in% vaccinated_full_time_series_metadata$subject_id,]
  immunogenicity_data_full_time_series <<- immunogenicity_data_full_time_series[,c(1,105,106,115,116)]
  immunogenicity_data_full_time_series <<- immunogenicity_data_full_time_series[order(immunogenicity_data_full_time_series$Vaccination.Day8_IFNg_NP.Background_SFC.10.6.cells, decreasing = TRUE),]
  high_t_cell_response_subjects <<- immunogenicity_data_full_time_series$SUBJID[1:7]
  low_t_cell_response_subjects <<- immunogenicity_data_full_time_series$SUBJID[8:14]
  # Grab high t cell response vaccinated counts and metadata
  high_t_cell_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$subject_id %in% high_t_cell_response_subjects,])
  high_t_cell_vaccinated_counts <<- vaccinated_counts[,high_t_cell_vaccinated_aliquots]
  high_t_cell_vaccinated_metadata <<- vaccinated_metadata[high_t_cell_vaccinated_aliquots,]
  # Grab low t cell response vaccinated counts and metadata
  low_t_cell_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$subject_id %in% low_t_cell_response_subjects,])
  low_t_cell_vaccinated_counts <<- vaccinated_counts[,low_t_cell_vaccinated_aliquots]
  low_t_cell_vaccinated_metadata <<- vaccinated_metadata[low_t_cell_vaccinated_aliquots,]
  # Grab combination of high and low t cell response vaccinated counts and metadata
  both_t_cell_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$subject_id %in% high_t_cell_response_subjects | vaccinated_metadata$subject_id %in% low_t_cell_response_subjects,])
  both_t_cell_vaccinated_counts <<- vaccinated_counts[,both_t_cell_vaccinated_aliquots]
  both_t_cell_vaccinated_metadata <<- vaccinated_metadata[both_t_cell_vaccinated_aliquots,]
  t_cell_for_metadata <<- both_t_cell_vaccinated_metadata$subject_id %in% high_t_cell_response_subjects
  t_cell_for_metadata <<- replace(t_cell_for_metadata, t_cell_for_metadata == TRUE, "HIGH")
  t_cell_for_metadata <<- replace(t_cell_for_metadata, t_cell_for_metadata == FALSE, "LOW")
  both_t_cell_vaccinated_metadata$t_cell_response <<- t_cell_for_metadata
  ### PLACEBO ###
  # Grab subject IDs for main 23 subjects
  placebo_full_time_series_subjects <<- unique(placebo_full_time_series_metadata$subject_id)
  # Reorder subject IDs according to viral load (high to low)
  placebo_full_time_series_subjects <<- placebo_full_time_series_subjects[order(match(placebo_full_time_series_subjects,placebo_viral_load_primary$subject_id))]
  # Top 13 will be high viral load and bottom 10 will be low viral load 
  high_viral_load_subjects <<- placebo_full_time_series_subjects[1:13]
  low_viral_load_subjects <<- tail(placebo_full_time_series_subjects, n = 10)
  # Grab high viral load placebo counts and metadata
  hvl_full_time_series_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$subject_id %in% high_viral_load_subjects,])
  hvl_full_time_series_placebo_counts <<- placebo_counts[,hvl_full_time_series_placebo_aliquots]
  hvl_full_time_series_placebo_metadata <<- placebo_metadata[hvl_full_time_series_placebo_aliquots,]
  # Grab low viral load placebo counts and metadata
  lvl_full_time_series_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$subject_id %in% low_viral_load_subjects,])
  lvl_full_time_series_placebo_counts <<- placebo_counts[,lvl_full_time_series_placebo_aliquots]
  lvl_full_time_series_placebo_metadata <<- placebo_metadata[lvl_full_time_series_placebo_aliquots,]
  # Grab combination of high and low viral load placebo counts and metadata
  both_full_time_series_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$subject_id %in% high_viral_load_subjects | placebo_metadata$subject_id %in% low_viral_load_subjects,])
  both_full_time_series_placebo_counts <<- placebo_counts[,both_full_time_series_placebo_aliquots]
  both_full_time_series_placebo_metadata <<- placebo_metadata[both_full_time_series_placebo_aliquots,]
  # Remove questionable low viral load individual (0 qPCRAUC) from placebo
  removed_low_viral_aliquots <- rownames(placebo_metadata[placebo_metadata$subject_id == "f18c54d93cef4a4e",])
  placebo_metadata <<- placebo_metadata[!(placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
  placebo_counts <<- placebo_counts[,!(colnames(placebo_counts) %in% removed_low_viral_aliquots)]
  both_full_time_series_placebo_metadata <<-  both_full_time_series_placebo_metadata[!(both_full_time_series_placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
  both_full_time_series_placebo_counts <<- both_full_time_series_placebo_counts[,!(colnames(both_full_time_series_placebo_counts) %in% removed_low_viral_aliquots)]
  lvl_full_time_series_placebo_metadata <<- lvl_full_time_series_placebo_metadata[!(lvl_full_time_series_placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
  lvl_full_time_series_placebo_counts <<- lvl_full_time_series_placebo_counts[,!(colnames(lvl_full_time_series_placebo_counts) %in% removed_low_viral_aliquots)]
  # Remove questionable low viral load individuals (NA or 0 qPCRAUC) from vaccinated
  vaccinated_metadata <<- vaccinated_metadata[!is.na(vaccinated_metadata$AVAL),]
  removed_low_viral_aliquots <- rownames(vaccinated_metadata[vaccinated_metadata$AVAL == 0,])
  removed_low_viral_subjects <- unique(vaccinated_metadata[vaccinated_metadata$AVAL == 0,]$subject_id)
  vaccinated_metadata <<- vaccinated_metadata[!(vaccinated_metadata$subject_id %in% removed_low_viral_subjects),]
  vaccinated_counts <<- vaccinated_counts[,!(colnames(vaccinated_counts) %in% removed_low_viral_aliquots)]
  both_t_cell_vaccinated_metadata <<-  both_t_cell_vaccinated_metadata[!(both_t_cell_vaccinated_metadata$subject_id %in% removed_low_viral_subjects),]
  both_t_cell_vaccinated_counts <<- both_t_cell_vaccinated_counts[,!(colnames(both_t_cell_vaccinated_counts) %in% removed_low_viral_aliquots)]
  low_t_cell_vaccinated_metadata <<- low_t_cell_vaccinated_metadata[!(low_t_cell_vaccinated_metadata$subject_id %in% removed_low_viral_subjects),]
  low_t_cell_vaccinated_counts <<- low_t_cell_vaccinated_counts[,!(colnames(low_t_cell_vaccinated_counts) %in% removed_low_viral_aliquots)]
  # Subset to HVL and LVL for vaccinated (based on thresholds in placebo)
  # Grab high viral load vaccinated counts and metadata
  hvl_threshold <- min(hvl_full_time_series_placebo_metadata$AVAL)
  lvl_threshold <- max(lvl_full_time_series_placebo_metadata$AVAL)
  hvl_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$AVAL >= hvl_threshold,])
  hvl_vaccinated_counts <<- vaccinated_counts[,hvl_vaccinated_aliquots]
  hvl_vaccinated_metadata <<- vaccinated_metadata[hvl_vaccinated_aliquots,]
  # Grab mvl viral load vaccinated counts and metadata
  mvl_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$AVAL > lvl_threshold & vaccinated_metadata$AVAL < hvl_threshold,])
  mvl_vaccinated_counts <<- vaccinated_counts[,mvl_vaccinated_aliquots]
  mvl_vaccinated_metadata <<- vaccinated_metadata[mvl_vaccinated_aliquots,]
  # Grab lvl viral load vaccinated counts and metadata
  lvl_vaccinated_aliquots <<- rownames(vaccinated_metadata[vaccinated_metadata$AVAL <= lvl_threshold,])
  lvl_vaccinated_counts <<- vaccinated_counts[,lvl_vaccinated_aliquots]
  lvl_vaccinated_metadata <<- vaccinated_metadata[lvl_vaccinated_aliquots,]
  # Get full list of HVL, MVL, and LVL for placebo
  hvl_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$AVAL >= hvl_threshold,])
  hvl_placebo_counts <<- placebo_counts[,hvl_placebo_aliquots]
  hvl_placebo_metadata <<- placebo_metadata[hvl_placebo_aliquots,]
  mvl_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$AVAL > lvl_threshold & placebo_metadata$AVAL < hvl_threshold,])
  mvl_placebo_counts <<- placebo_counts[,mvl_placebo_aliquots]
  mvl_placebo_metadata <<- placebo_metadata[mvl_placebo_aliquots,]
  lvl_placebo_aliquots <<- rownames(placebo_metadata[placebo_metadata$AVAL <= lvl_threshold,])
  lvl_placebo_counts <<- placebo_counts[,lvl_placebo_aliquots]
  lvl_placebo_metadata <<- placebo_metadata[lvl_placebo_aliquots,]
  # Get full list of all subjects for full time series
  all_full_time_series_metadata <<- rbind(both_t_cell_vaccinated_metadata[,-c(51,52,53,54,55)], both_full_time_series_placebo_metadata)
  all_full_time_series_counts <<- cbind(both_t_cell_vaccinated_counts, both_full_time_series_placebo_counts)
}

run_deseq_bulk_analysis_time_series=function(sample_type, counts, metadata, test_time, baseline_time, output_dir, output_name_prefix=NA, alpha = 0.05, apply_correction = "sv") {
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  print(sample_type)
  print(test_time)
  print(baseline_time)
  print(output_dir)
  print(output_name_prefix)
  print(alpha)
  if(alpha > 0.9) {
    alpha_for_subsetting_df <- 1.5
  } else {
    alpha_for_subsetting_df <- alpha
  }
  # Select the two relevant time points from our metadata
  metadata_subset <- metadata[metadata$time_point == test_time | metadata$time_point == baseline_time,]
  # Remove subjects that only have one time point (not both)
  metadata_subset <- metadata_subset[metadata_subset$subject_id  %in% names(table(metadata_subset$subject_id)[table(metadata_subset$subject_id) == 2]),]
  # Print distribution of time points to make sure we're doing
  print(table(metadata_subset$time_point))
  # Select subset of counts associated with subjects
  counts_subset <- counts[rownames(metadata_subset)]
  # Run DESeq2
  median_value <- median(metadata_subset$Absolute.score..sig.score.)
  
  # Replace values below median with "LOW" and above median with "HIGH"
  #metadata_subset$Absolute.score..sig.score. <- ifelse(metadata_subset$Absolute.score..sig.score. < median_value, "LOW", ifelse(metadata_subset$Absolute.score..sig.score. > median_value, "HIGH", metadata_subset$Absolute.score..sig.score. ))
  #metadata_subset$Absolute.score..sig.score. <- factor(metadata_subset$Absolute.score..sig.score., levels = c("LOW", "HIGH"))
  #metadata_subset$Absolute.score..sig.score. <- scale(metadata_subset$Absolute.score..sig.score.)
  #current_analysis <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + Absolute.score..sig.score. + time_point)
  if(apply_correction == "sv") {
    print("Applying SV correction")
    base_model <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + time_point)
    base_model <- estimateSizeFactors(base_model)
    dat <- counts(base_model, normalized = TRUE)
    idx <- rowMeans(dat) > 1
    dat <- dat[idx, ]
    mod  <- model.matrix(~ subject_id + time_point, colData(base_model))
    mod0 <- model.matrix(~ subject_id, colData(base_model))
    print(paste0("Number of surrogate variables recommended by be: ", num.sv(dat, mod, method = "be")))
    print(paste0("Number of surrogate variables recommended by leek: ", num.sv(dat, mod, method = "leek")))
    svseq <- svaseq(dat, mod, mod0, n.sv = 1)
    metadata_subset$SV1 <- svseq$sv[,1]
    metadata_subset$SV1_discrete <- svseq$sv[,1]
    median_value <- median(metadata_subset$SV1_discrete)
    metadata_subset$SV1_discrete <- ifelse(metadata_subset$SV1_discrete < median_value, "LOW", ifelse(metadata_subset$SV1_discrete > median_value, "HIGH", metadata_subset$SV1_discrete ))
    metadata_subset$SV1_discrete <- factor(metadata_subset$SV1_discrete, levels = c("LOW", "HIGH"))
    current_analysis <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + SV1 + time_point)
    sva_model_discrete <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + SV1_discrete + time_point)
    vsd <- vst(sva_model_discrete, blind = FALSE)
    pcaData <- plotPCA(vsd, intgroup = c( "subject_id", "SV1_discrete", "time_point"), pcsToUse = c(1,2), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    sv_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = SV1_discrete, shape = time_point)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      ggtitle("PCA with VST data (SV)")
    subject_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = subject_id, shape = time_point)) +
      geom_point(size = 3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      ggtitle("PCA with VST data (Subject ID)")
  } else if(apply_correction == "absolute_score") {
    metadata_subset$Absolute.score..sig.score. <- scale(metadata_subset$Absolute.score..sig.score.)
    current_analysis <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + Absolute.score..sig.score. + time_point)
    sv_plot <- NULL
    subject_plot <- NULL
    t_cell_plot <- NULL
  } else if(apply_correction == "none") {
    print("Applying no correction")
    current_analysis <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + time_point)
    sv_plot <- NULL
    subject_plot <- NULL
    t_cell_plot <- NULL
  }
  
  current_analysis <- DESeq(current_analysis)
  save(current_analysis, file = paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, ".rds"))
  # Grab results with no lfcThreshold set
  current_analysis_results_none <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha)
  current_analysis_results_none <- current_analysis_results_none[order(current_analysis_results_none$padj),]
  current_analysis_results_none <- subset(current_analysis_results_none, padj < alpha_for_subsetting_df)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_fc_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_fc_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with lfcThreshold = 0.1
  current_analysis_results_0.1 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 0.1)
  current_analysis_results_0.1 <- current_analysis_results_0.1[order(current_analysis_results_0.1$padj),]
  current_analysis_results_0.1 <- subset(current_analysis_results_0.1, padj < alpha_for_subsetting_df)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_0.1), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_fc_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_0.1), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_fc_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with lfcThreshold = 0.2
  current_analysis_results_0.2 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 0.2)
  current_analysis_results_0.2 <- current_analysis_results_0.2[order(current_analysis_results_0.2$padj),]
  current_analysis_results_0.2 <- subset(current_analysis_results_0.2, padj < alpha_for_subsetting_df)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_0.2), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_fc_0.2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_0.2), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_fc_0.2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with lfcThreshold = 0.3
  current_analysis_results_0.3 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 0.3)
  current_analysis_results_0.3 <- current_analysis_results_0.3[order(current_analysis_results_0.3$padj),]
  current_analysis_results_0.3 <- subset(current_analysis_results_0.3, padj < alpha_for_subsetting_df)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_0.3), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_fc_0.3.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_0.3), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_fc_0.3.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with lfcThreshold = 0.585 (1.5 fold increase)
  current_analysis_results_0.585 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 0.585)
  current_analysis_results_0.585 <- current_analysis_results_0.585[order(current_analysis_results_0.585$padj),]
  current_analysis_results_0.585 <- subset(current_analysis_results_0.585, padj < alpha_for_subsetting_df)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_0.585), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_fc_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_0.585), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_fc_585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)   
  }
  # Grab results with lfcThreshold = 1
  current_analysis_results_1 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 1)
  current_analysis_results_1 <- current_analysis_results_1[order(current_analysis_results_1$padj),]
  current_analysis_results_1 <- subset(current_analysis_results_1, padj < alpha_for_subsetting_df)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_1), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_fc_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_1), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_fc_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with lfcThreshold = 2
  current_analysis_results_2 <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = alpha, lfcThreshold = 2)
  current_analysis_results_2 <- current_analysis_results_2[order(current_analysis_results_2$padj),]
  current_analysis_results_2 <- subset(current_analysis_results_2, padj < alpha_for_subsetting_df)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_fc_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_fc_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with no p-value thresholding (grab every gene)
  current_analysis_results_unfiltered <- results(current_analysis, contrast = c("time_point", test_time, baseline_time), alpha = 0.99999)
  current_analysis_results_unfiltered <- current_analysis_results_unfiltered[order(current_analysis_results_unfiltered$padj),]
  current_analysis_results_unfiltered <- subset(current_analysis_results_unfiltered, padj < 1.5)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_unfiltered), paste0(output_dir, test_time, "_vs_", baseline_time, "_", sample_type, "_fc_unfiltered.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_unfiltered), paste0(output_dir, test_time, "_vs_", baseline_time, "_", output_name_prefix, "_", sample_type, "_fc_unfiltered.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  return(list(current_analysis_results_none, current_analysis_results_0.1, current_analysis_results_0.2, current_analysis_results_0.3, current_analysis_results_0.585, current_analysis_results_1, current_analysis_results_2, current_analysis_results_unfiltered, sv_plot, subject_plot))
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
    write.table(rownames(current_analysis_results), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_fc_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_0.1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 0.585 (1.5 fold increase)
  current_analysis_results_1.5 <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05, lfcThreshold = 0.585)
  current_analysis_results_1.5 <- current_analysis_results_1.5[order(current_analysis_results_1.5$padj),]
  current_analysis_results_1.5 <- subset(current_analysis_results_1.5, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_fc_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_1.5), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_0.585.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)   
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 1
  current_analysis_results_2 <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05, lfcThreshold = 1)
  current_analysis_results_2 <- current_analysis_results_2[order(current_analysis_results_2$padj),]
  current_analysis_results_2 <- subset(current_analysis_results_2, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_fc_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_2), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with alpha = 0.05 and lfcThreshold = 2
  current_analysis_results_4 <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05, lfcThreshold = 2)
  current_analysis_results_4 <- current_analysis_results_4[order(current_analysis_results_4$padj),]
  current_analysis_results_4 <- subset(current_analysis_results_4, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_4), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_fc_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    write.table(rownames(current_analysis_results_4), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", output_name_prefix, "_", sample_type, "_2.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  # Grab results with alpha = 0.05 and no lfcThreshold
  current_analysis_results_none <- results(current_analysis, contrast = c("viral_load", test_cond, baseline_cond), alpha = 0.05)
  current_analysis_results_none <- current_analysis_results_none[order(current_analysis_results_none$padj),]
  current_analysis_results_none <- subset(current_analysis_results_none, padj < 0.05)
  if(is.na(output_name_prefix)) {
    write.table(rownames(current_analysis_results_none), paste0(output_dir, test_cond, "_vs_", baseline_cond, "_", sample_type, "_fc_none.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
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
  
  print("Applying SV correction")
  base_model <- DESeqDataSetFromMatrix(countData = LRT_counts, colData = LRT_metadata, design = ~ subject_id + time_point)
  base_model <- estimateSizeFactors(base_model)
  dat <- counts(base_model, normalized = TRUE)
  idx <- rowMeans(dat) > 1
  dat <- dat[idx, ]
  mod  <- model.matrix(~ subject_id + time_point, colData(base_model))
  mod0 <- model.matrix(~ subject_id, colData(base_model))
  print(paste0("Number of surrogate variables recommended by be: ", num.sv(dat, mod, method = "be")))
  print(paste0("Number of surrogate variables recommended by leek: ", num.sv(dat, mod, method = "leek")))
  svseq <- svaseq(dat, mod, mod0, n.sv = 1)
  LRT_metadata$SV1 <- svseq$sv[,1]
  
  LRT_analysis <- DESeqDataSetFromMatrix(countData = LRT_counts,
                                         colData = LRT_metadata,
                                         design = ~ subject_id + SV1 + time_point)
  LRT_analysis <- DESeq(LRT_analysis, test = "LRT", reduced = ~ subject_id + SV1)
  LRT_analysis_results <- results(LRT_analysis, alpha = 0.05)
  LRT_analysis_results <- LRT_analysis_results[order(LRT_analysis_results$padj),]
  LRT_analysis_results_filtered <- subset(LRT_analysis_results, padj < 0.05)
  return(list(LRT_analysis, LRT_analysis_results_filtered))
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

# This method finds overlap between enriched transcription factors in pseudobulk ATAC-seq and those genes expressed in the scRNA-seq data
find_overlapping_motifs_between_atac_and_rna <- function(peak_dir, sc_gene_table, cell_types, pos_bulk_genes = NULL, neg_bulk_genes = NULL) {
  overlapping_motif_df <- data.frame(tf = character(), cell_types = character(), found_in_bulk = logical())
  for(cell_type in cell_types) {
    if(file.exists(paste0(peak_dir, cell_type, "_D28_D1_motif_up_pseudobulk_only.tsv"))) {
      # Up-regulated enriched TFs
      current_motifs_up <- read.table(paste0(peak_dir, cell_type, "_D28_D1_motif_up_pseudobulk_only.tsv"), sep = "\t", header = TRUE)
      current_motifs_up <- current_motifs_up[current_motifs_up$p_adj < 0.05,]
      current_tfs_up <- current_motifs_up$TF
      # Cut off _ID so we can match with gene names from RNA
      current_tfs_up <- sub("_.*", "", current_tfs_up)
      # Down-regulated enriched TFs
      current_motifs_down <- read.table(paste0(peak_dir, cell_type, "_D28_D1_motif_down_pseudobulk_only.tsv"), sep = "\t", header = TRUE)
      current_motifs_down <- current_motifs_down[current_motifs_down$p_adj < 0.05,]
      current_tfs_down <- current_motifs_down$TF
      # Cut off _ID so we can match with gene names from RNA
      current_tfs_down <- sub("_.*", "", current_tfs_down)
      # Grab genes from SC table for current cell type
      cell_type_in_df <- sub("_", " ", cell_type)
      cell_type_sc_gene_table <- sc_gene_table[sc_gene_table$Cell_Type == cell_type_in_df,]
      print(current_tfs_up)
      print(current_tfs_down)
      print(cell_type_sc_gene_table[cell_type_sc_gene_table$sc_log2FC > 0,]$Gene_Name)
      print(cell_type_sc_gene_table[cell_type_sc_gene_table$sc_log2FC < 0,]$Gene_Name)
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

# This method finds overlap between enriched transcription factors in mintchip and those genes expressed in the scRNA-seq data
find_overlapping_motifs_between_mintchip_and_rna <- function(peak_dir, sc_gene_table, markers, pos_bulk_genes = NULL, neg_bulk_genes = NULL) {
  overlapping_motif_df <- data.frame(marker = character(), tf = character(), cell_types = character(), found_in_bulk = logical())
  for(marker in markers) {
    specific_marker_dir_pos <- paste0(marker, "_D28_D1_diff_filtered_homer_pos/")
    specific_marker_dir_neg <- paste0(marker, "_D28_D1_diff_filtered_homer_neg/")
    if(file.exists(paste0(peak_dir, specific_marker_dir_pos, "knownResults.txt"))) {
      # Up-regulated enriched TFs
      current_motifs_up <- read.table(paste0(peak_dir, specific_marker_dir_pos, "knownResults.txt"), sep = "\t", header = TRUE, comment.char = "")
      current_motifs_up <- current_motifs_up[current_motifs_up$q.value..Benjamini. < 0.05,]
      current_tfs_up <- current_motifs_up$Motif.Name
      # Cut off _ID so we can match with gene names from RNA
      current_tfs_up <- sub("\\(.*", "", current_tfs_up)
      # Down-regulated enriched TFs
      current_motifs_down <- read.table(paste0(peak_dir, specific_marker_dir_neg, "knownResults.txt"), sep = "\t", header = TRUE, comment.char = "")
      current_motifs_down <- current_motifs_down[current_motifs_down$q.value..Benjamini. < 0.05,]
      current_tfs_down <- current_motifs_down$Motif.Name
      # Cut off _ID so we can match with gene names from RNA
      current_tfs_down <- sub("\\(.*", "", current_tfs_down)
      # Grab overlap between pos SC genes and upregulated TF motifs
      overlapping_pos_tfs <- intersect(current_tfs_up, sc_gene_table[sc_gene_table$sc_log2FC > 0,]$Gene_Name)
      overlapping_pos_tfs_upper <- intersect(toupper(current_tfs_up), sc_gene_table[sc_gene_table$sc_log2FC > 0,]$Gene_Name)
      overlapping_pos_tfs <- unique(union(overlapping_pos_tfs, overlapping_pos_tfs_upper))
      # Grab overlap between neg SC genes and downregulated TF motifs
      overlapping_neg_tfs <- intersect(current_tfs_down, sc_gene_table[sc_gene_table$sc_log2FC < 0,]$Gene_Name)
      overlapping_neg_tfs_upper <- intersect(toupper(current_tfs_down), sc_gene_table[sc_gene_table$sc_log2FC < 0,]$Gene_Name)
      overlapping_neg_tfs <- unique(union(overlapping_neg_tfs, overlapping_neg_tfs_upper))
      for(pos_tf in overlapping_pos_tfs) {
        cell_types <- sc_gene_table[sc_gene_table$Gene_Name == pos_tf,]$Cell_Type
        # If the tf was found in the bulk data, note that
        found_in_bulk <- FALSE
        if(pos_tf %in% pos_bulk_genes) {
          found_in_bulk <- TRUE
        }
        # Add tf to df (and mark it as up)
        motif_vector <- c(marker, paste0(pos_tf, " (Up)"), paste(cell_types, collapse = ","), found_in_bulk)
        motif_vector <- as.data.frame(t(motif_vector))
        names(motif_vector) <- c("marker", "tf", "cell_types", "found_in_bulk")
        overlapping_motif_df <- rbind(overlapping_motif_df, motif_vector)
      }
      for(neg_tf in overlapping_neg_tfs) {
        cell_types <- sc_gene_table[sc_gene_table$Gene_Name == neg_tf,]$Cell_Type
        # If the tf was found in the bulk data, note that
        found_in_bulk <- FALSE
        if(neg_tf %in% neg_bulk_genes) {
          found_in_bulk <- TRUE
        }
        # Add tf to df (and mark it as down)
        motif_vector <- c(marker, paste0(neg_tf, " (Down)"), paste(cell_types, collapse = ","), found_in_bulk)
        motif_vector <- as.data.frame(t(motif_vector))
        names(motif_vector) <- c("marker", "tf", "cell_types", "found_in_bulk")
        overlapping_motif_df <- rbind(overlapping_motif_df, motif_vector)
      }
    }
  }
  return(overlapping_motif_df)
}



fill_in_info_for_magical_output <- function(overall_magical_df, das_dir, sc_pseudobulk_deg_combined_cell_types_table, sc_pseudobulk_deg_table,
                                            pos_bulk_genes, neg_bulk_genes, sc_peaks, overall_motif_enrichment_df) {
  # Distance between relevant gene TSS and peak 
  dist_between_gene_tss_and_peak <- c()
  # FC for circuit gene from SC data (positive or negative?)
  gene_fcs <- c()
  # FC for circuit site from SC data (positive or negative?)
  site_fcs <- c()
  # PCT 1 for circuit site from SC data
  site_pct_1s <- c()
  # PCT 2 for circuit site from SC data
  site_pct_2s <- c()
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
  motif_enrichment <- c()
  for(current_circuit_index in 1:nrow(overall_magical_df)) {
    current_circuit <- overall_magical_df[current_circuit_index,]
    cell_type <- current_circuit$Cell_Type
    cell_type_in_df <- sub("_", " ", cell_type)
    current_gene <- current_circuit$Gene_symbol
    current_tfs <- current_circuit$TFs.binding.prob.
    current_tfs <- unlist(strsplit(current_tfs, ","))
    current_tfs <- current_tfs[c(TRUE, FALSE)]
    motif_info <- ""
    for(current_tf in current_tfs) {
      motif_enrichment_subset_df <- overall_motif_enrichment_df[overall_motif_enrichment_df$motif.name == current_tf,]
      motif_enrichment_subset_df <- motif_enrichment_subset_df[motif_enrichment_subset_df$Cell_Type == cell_type_in_df,]
      if(nrow(motif_enrichment_subset_df) > 0) {
        current_motif_info <- paste0(current_tf, ": ")
        current_motif_tf_info <- paste(motif_enrichment_subset_df$Direction, " (Rank ", motif_enrichment_subset_df$rank, ")", sep = "")
        current_motif_tf_info <- paste(current_motif_tf_info, collapse = ",")
        current_motif_info <- paste0(current_motif_info, current_motif_tf_info, ".")
      } else {
        current_motif_info <- paste0(current_tf, ": N/A. ")
      }
      motif_info <- paste0(motif_info, current_motif_info)
    }
    motif_enrichment <- c(motif_enrichment, motif_info)
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
    site_info <- read.table(paste0(sc_das_dir, "diff_peaks/D28-vs-D_minus_1-degs-", sub(" ", "_", cell_type), "-time_point-controlling_for_subject_id_sc_filtered_pct_0.01.tsv"), sep = "\t",
                            header = TRUE)
    current_peak_name <- paste0(current_circuit$Peak_chr, "-", current_circuit$Peak_start, "-", current_circuit$Peak_end)
    site_info <- site_info[rownames(site_info) == current_peak_name,]
    site_fc <- site_info$avg_log2FC
    site_pct_1 <- site_info$pct.1
    site_pct_2 <- site_info$pct.2
    site_fcs <- c(site_fcs, site_fc)
    site_pct_1s <- c(site_pct_1s, site_pct_1)
    site_pct_2s <- c(site_pct_2s, site_pct_2)
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
  overall_magical_df$site_pct_1 <- site_pct_1s
  overall_magical_df$site_pct_2 <- site_pct_2s
  overall_magical_df$nearest_gene_to_site <- nearest_gene_to_site
  #overall_magical_df$dist_to_nearest_gene_start <- dist_to_nearest_gene_start
  #overall_magical_df$dist_to_nearest_gene_tss <- dist_to_nearest_gene_tss
  overall_magical_df$site_type_for_nearest_gene <- site_type
  overall_magical_df$bulk <- bulk_boolean
  overall_magical_df$magical_cell_types <- magical_cell_types
  overall_magical_df$combined_single_cell_types <- combined_single_cell_types
  overall_magical_df$original_single_cell_types <- original_single_cell_types
  overall_magical_df$motif_enrichment_for_tf <- motif_enrichment
  return(overall_magical_df)
}

fill_in_info_for_magical_tf_output <- function(overall_magical_df, overall_motif_enrichment_df) {
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
  # Motif_Enrichment (Cell Type and Direction)
  motif_enrichment <- c()
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
    current_motif_enrichment_tf <- overall_motif_enrichment_df[overall_motif_enrichment_df$motif.name == current_tf,]
    if(nrow(current_motif_enrichment_tf) > 0) {
      current_motif_enrichment_cell_types <- current_motif_enrichment_tf$Cell_Type
      current_motif_enrichment_directions <- current_motif_enrichment_tf$Direction
      current_motif_enrichment_ranks <- current_motif_enrichment_tf$rank
      current_motif_enrichment_info <- paste(current_motif_enrichment_cell_types, " (", current_motif_enrichment_directions, "; Rank ", current_motif_enrichment_ranks, ")", sep = "")
      current_motif_enrichment_info <- paste(current_motif_enrichment_info, collapse = ",")
      motif_enrichment <- c(motif_enrichment, current_motif_enrichment_info)
    } else {
      motif_enrichment <- c(motif_enrichment, "N/A")
    }
  }
  overall_magical_tf_df <- data.frame(tf = tf, total_circuit_count = total_circuit_count, circuit_cell_types = circuit_cell_types,
                                      bound_genes = bound_genes, found_as_circuit_genes = found_as_circuit_genes,
                                      found_as_DEGs_in_combined_cell_types = found_as_DEGs_in_combined_cell_types,
                                      found_as_DEGs_in_original_cell_types = found_as_DEGs_in_original_cell_types,
                                      found_as_DEGs_in_bulk = found_as_DEGs_in_bulk,
                                      motif_enrichment = motif_enrichment)
  return(overall_magical_tf_df)
}






fill_in_sc_deg_info_for_time_series <- function(sc_gene_df, counts, metadata, output_dir, sc_fc_direction, alpha = 0.05, filter_D28 = TRUE) {
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
  if(filter_D28) {
    high_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                          "2_D28", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D28_vs_D_minus_1_alpha_", alpha, "/"), "high", alpha = alpha)
    high_placebo_period_2_D28_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D28_vs_D_minus_1_results
  } else {
    high_placebo_period_2_D28_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                          "2_D28", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D28_vs_D_minus_1_alpha_", alpha, "/"), "high", alpha = alpha)
    high_placebo_period_2_D28_vs_D_minus_1_results_unfiltered <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                                     "2_D28", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D28_vs_D_minus_1_alpha_", alpha, "/"), "high", alpha = 0.99999)
  }
  high_placebo_period_2_D28_vs_D_minus_1_results <- high_placebo_period_2_D28_vs_D_minus_1_results[[1]]
  high_placebo_period_2_D28_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D28_vs_D_minus_1_results_unfiltered[[1]]
  
  if(filter_D28) {
    filtered_sc_trained_immunity_genes <- intersect(possible_genes,
                                                    rownames(high_placebo_period_2_D28_vs_D_minus_1_results))
  } else {
    filtered_sc_trained_immunity_genes <- intersect(possible_genes,
                                                    rownames(high_placebo_period_2_D28_vs_D_minus_1_results_unfiltered))
  }
  
  # Verify that FC in bulk agrees with FC in SC (if not, discard)
  final_filtered_sc_trained_immunity_genes <- c()
  
  for(gene in filtered_sc_trained_immunity_genes) {
    sc_fc <- sc_gene_df[sc_gene_df$Gene_Name == gene,]$sc_log2FC
    if(filter_D28) {
      bulk_fc <- high_placebo_period_2_D28_vs_D_minus_1_results[rownames(high_placebo_period_2_D28_vs_D_minus_1_results) == gene,]$log2FoldChange
    } else {
      bulk_fc <- high_placebo_period_2_D28_vs_D_minus_1_results_unfiltered[rownames(high_placebo_period_2_D28_vs_D_minus_1_results_unfiltered) == gene,]$log2FoldChange
    }
    
    if(all(sign(sc_fc) == sign(bulk_fc))) {
      final_filtered_sc_trained_immunity_genes <- c(final_filtered_sc_trained_immunity_genes, gene)
    } else {
      print(paste0("Gene ", gene, " has mismatched FC (positive FC in SC and negative FC in bulk or vice versa)"))
    }
  }
  
  # Compute results for D2 / D5 / D8 (with alpha parameter and without alpha parameter)
  # WITH ALPHA SET
  high_placebo_period_2_D2_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                       "2_D2", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D2_vs_D_minus_1_alpha_", alpha, "/"), "high", alpha = alpha)
  high_placebo_period_2_D2_vs_D_minus_1_results <- high_placebo_period_2_D2_vs_D_minus_1_results[[1]]
  
  high_placebo_period_2_D5_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                       "2_D5", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D5_vs_D_minus_1_alpha_", alpha, "/"), "high", alpha = alpha)
  high_placebo_period_2_D5_vs_D_minus_1_results <- high_placebo_period_2_D5_vs_D_minus_1_results[[1]]
  
  high_placebo_period_2_D8_vs_D_minus_1_results <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                       "2_D8", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D8_vs_D_minus_1_alpha_", alpha, "/"), "high", alpha = alpha)
  high_placebo_period_2_D8_vs_D_minus_1_results <- high_placebo_period_2_D8_vs_D_minus_1_results[[1]]
  # WITHOUT ALPHA SET
  high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                                  "2_D2", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D2_vs_D_minus_1_alpha_0.99999/"), "high", alpha = 0.99999)
  high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered[[1]]
  
  high_placebo_period_2_D5_vs_D_minus_1_results_unfiltered <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                                  "2_D5", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D5_vs_D_minus_1_alpha_0.99999/"), "high", alpha = 0.99999)
  high_placebo_period_2_D5_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D5_vs_D_minus_1_results_unfiltered[[1]]
  
  high_placebo_period_2_D8_vs_D_minus_1_results_unfiltered <- run_deseq_bulk_analysis_time_series("placebo", counts, metadata,
                                                                                                  "2_D8", "2_D_minus_1", paste0(output_dir, "high_placebo_period_2_D8_vs_D_minus_1_alpha_0.99999/"), "high", alpha = 0.99999)
  high_placebo_period_2_D8_vs_D_minus_1_results_unfiltered <- high_placebo_period_2_D8_vs_D_minus_1_results_unfiltered[[1]]
  
  # Parse our list of genes
  for(current_gene in final_filtered_sc_trained_immunity_genes) {
    cell_types <- sc_gene_df[sc_gene_df$Gene_Name == current_gene,]$Cell_Type
    cell_types <- paste0("(", paste(cell_types, collapse = ", "), ")")
    final_df <- rbind(final_df, grab_deg_info_for_sc_gene(current_gene, cell_types, high_placebo_period_2_D2_vs_D_minus_1_results, 
                                                          high_placebo_period_2_D2_vs_D_minus_1_results_unfiltered, "D2"))
    final_df <- rbind(final_df, grab_deg_info_for_sc_gene(current_gene, cell_types, high_placebo_period_2_D5_vs_D_minus_1_results, 
                                                          high_placebo_period_2_D5_vs_D_minus_1_results_unfiltered, "D5"))
    final_df <- rbind(final_df, grab_deg_info_for_sc_gene(current_gene, cell_types, high_placebo_period_2_D8_vs_D_minus_1_results, 
                                                          high_placebo_period_2_D8_vs_D_minus_1_results_unfiltered, "D8"))
    final_df <- rbind(final_df, grab_deg_info_for_sc_gene(current_gene, cell_types, high_placebo_period_2_D28_vs_D_minus_1_results, 
                                                          high_placebo_period_2_D28_vs_D_minus_1_results_unfiltered, "D28"))
  }
  final_df$Day <- replace(final_df$Day, final_df$Day == "D2", "Day.2")
  final_df$Day <- replace(final_df$Day, final_df$Day == "D5", "Day.5")
  final_df$Day <- replace(final_df$Day, final_df$Day == "D8", "Day.8")
  final_df$Day <- replace(final_df$Day, final_df$Day == "D28", "Day.28")
  final_df$Day <- factor(final_df$Day, levels = c("Day.2","Day.5","Day.8","Day.28"))
  return(final_df)
}

grab_deg_info_for_sc_gene <- function(current_gene, cell_types, current_deseq2_results, current_deseq2_results_unfiltered, current_day) {
  if(current_gene %in% rownames(current_deseq2_results)) {
    gene_result <- current_deseq2_results[rownames(current_deseq2_results) == current_gene,]
    fold.change.abs <- abs(gene_result$log2FoldChange)
    fold.change.direction.raw <- ifelse(gene_result$log2FoldChange > 0, "Positive", "Negative")
    fold.change.direction <- ifelse(gene_result$log2FoldChange > 0, "Positive", "Negative")
    adjusted.p.value <- gene_result$padj
  } else {
    gene_result <- current_deseq2_results_unfiltered[rownames(current_deseq2_results_unfiltered) == current_gene,]
    fold.change.abs <- abs(gene_result$log2FoldChange)
    fold.change.direction.raw <- ifelse(gene_result$log2FoldChange > 0, "Positive", "Negative")
    fold.change.direction <- "Not Significant"
    adjusted.p.value <- gene_result$padj
  }
  return(data.frame(Gene = paste0(cell_types, " ", current_gene), Day = current_day, Fold.Change.Abs = fold.change.abs, 
                    Fold.Change.Direction.Raw = fold.change.direction.raw, Fold.Change.Direction = fold.change.direction, 
                    Adjusted.P.Value = adjusted.p.value))
}

# Reads motif table (necessary now because file names include total peak counts)
read_motif_table <- function(current_dir, analysis_type, pct, fc_threshold, direction) {
  # List files in the directory
  files <- list.files(current_dir)
  
  # Define the tokens you want to match
  pct <- paste0("_", pct, "_")
  if(direction == "pos") {
    fc_threshold <- paste0("_", fc_threshold, "_")
  } else {
    fc_threshold <- paste0("_-", fc_threshold, "_")
  }
  
  direction <- paste0("_", direction, "_")
  tokens <- c(analysis_type, pct, fc_threshold, direction)
  
  # Function to check if a file contains all tokens
  check_all_tokens <- function(file, tokens) {
    all(sapply(tokens, function(token) grepl(token, file)))
  }
  
  
  # Filter files that contain all tokens
  matching_file <- files[ sapply(files, check_all_tokens, tokens) ]
  
  # Read in table
  matching_file <- read.table(paste0(current_dir, matching_file), sep = "\t", header = TRUE)
  return(matching_file)
}

# Evaluate bulk cell type proportion changes from CIBERSORTx
evaluate_bulk_cell_type_proportion_changes <- function(metadata, bulk_cell_types, test_time, baseline_time, absolute_score = FALSE) {
  metadata_subset <- metadata[metadata$time_point == test_time | metadata$time_point == baseline_time,]
  # Remove subjects that only have one time point (not both)
  metadata_subset <- metadata_subset[metadata_subset$subject_id %in% names(table(metadata_subset$subject_id)[table(metadata_subset$subject_id) == 2]),]
  #if(absolute_score) {
  #  cols_to_divide <- 23:44
    # Iterate over each row and divide values in specified columns by 'Absolute.score..sig.score.'
  #  for (i in 1:nrow(metadata_subset)) {
  #    metadata_subset[i, cols_to_divide] <- metadata_subset[i, cols_to_divide] / metadata_subset[i, "Absolute.score..sig.score."]
  #  }
  # }
  # Find unadjusted p-values for each cell type of interest
  # We use coin::wilcox_test because it handles ties (which happen when there are 0s in the cell type proportions)
  # Unfortunately, we can't use loops because of the way the function call works (I think)
  cell_type_proportion_p_values <- c()
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(B.cells.naive ~ time_point, data = metadata_subset)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(T.cells.CD8 ~ time_point, data = metadata_subset)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(T.cells.CD4.naive ~ time_point, data = metadata_subset)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(T.cells.CD4.memory.resting ~ time_point, data = metadata_subset)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(T.cells.regulatory..Tregs. ~ time_point, data = metadata_subset)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(NK.cells.resting ~ time_point, data = metadata_subset)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(Monocytes ~ time_point, data = metadata_subset)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(Mast.cells.resting ~ time_point, data = metadata_subset)))
  cell_type_proportion_p_values <- c(cell_type_proportion_p_values, coin::pvalue(coin::wilcox_test(Neutrophils ~ time_point, data = metadata_subset)))
  # Adjust for multiple hypothesis testing
  cell_type_proportion_p_values_adjusted <- p.adjust(cell_type_proportion_p_values, method = "BH")
  names(cell_type_proportion_p_values_adjusted) <- bulk_cell_types
  return(cell_type_proportion_p_values_adjusted)
}

remove_text_in_parentheses <- function(input_string) {
  output_string <- gsub("\\([^\\)]+\\)", "", input_string)
  return(trimws(output_string))
}

compare_absolute_score_and_sv <- function(counts, metadata, test_time, baseline_time) {
  metadata_subset <- metadata[metadata$time_point == test_time | metadata$time_point == baseline_time,]
  # Remove subjects that only have one time point (not both)
  metadata_subset <- metadata_subset[metadata_subset$subject_id  %in% names(table(metadata_subset$subject_id)[table(metadata_subset$subject_id) == 2]),]
  # Select subset of counts associated with subjects
  counts_subset <- counts[rownames(metadata_subset)]
  # Run base analysis with no correction
  base_model <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + time_point)
  base_run <- DESeq(base_model)
  base_res <- run_fc_thresholding(base_run, test_time, baseline_time)
  # Run analysis with log transformed absolute score
  metadata_subset$Absolute.score..sig.score. <- scale(metadata_subset$Absolute.score..sig.score.)
  absolute_score_model <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + Absolute.score..sig.score. + time_point)
  absolute_score_run <- DESeq(absolute_score_model)
  absolute_score_res <- run_fc_thresholding(absolute_score_run, test_time, baseline_time)
  # Run analysis with SV1
  base_model <- estimateSizeFactors(base_model)
  dat  <- counts(base_model, normalized = TRUE)
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  mod  <- model.matrix(~ subject_id + time_point, colData(base_model))
  mod0 <- model.matrix(~ subject_id, colData(base_model))
  svseq <- svaseq(dat, mod, mod0, n.sv = 1)
  metadata_subset$SV1 <- svseq$sv[,1]
  sv_model <- DESeqDataSetFromMatrix(countData = counts_subset, colData = metadata_subset, design = ~ subject_id + time_point + SV1)
  sv_run <- DESeq(sv_model)
  sv_res <- run_fc_thresholding(sv_run, test_time, baseline_time)
  # Check correlation between SV1 and other numeric covariates
  metadata_subset_for_cor <- metadata_subset[,23:51]
  SV1_cors <- metadata_subset_for_cor %>% 
    correlate() %>% 
    focus(SV1) %>%
    filter(!is.na(SV1)) %>%
    arrange(SV1)
  # Return everything
  return(list(base_res, absolute_score_res, sv_res, SV1_cors))
}

run_fc_thresholding <- function(run, test_time, baseline_time) {
  # No FC threshold
  fc_0 <- results(run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05)
  fc_0 <- fc_0[order(fc_0$padj),]
  fc_0 <- subset(fc_0, padj < 0.05)
  # lfcThreshold = 0.1
  fc_0.1 <- results(run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 0.1)
  fc_0.1 <- fc_0.1[order(fc_0.1$padj),]
  fc_0.1 <- subset(fc_0.1, padj < 0.05)
  # lfcThreshold = 0.2
  fc_0.2 <- results(run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 0.2)
  fc_0.2 <- fc_0.2[order(fc_0.2$padj),]
  fc_0.2 <- subset(fc_0.2, padj < 0.05)
  # lfcThreshold = 0.3
  fc_0.3 <- results(run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 0.3)
  fc_0.3 <- fc_0.3[order(fc_0.3$padj),]
  fc_0.3 <- subset(fc_0.3, padj < 0.05)
  # lfcThreshold = 0.585
  fc_0.585 <- results(run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 0.585)
  fc_0.585 <- fc_0.585[order(fc_0.585$padj),]
  fc_0.585 <- subset(fc_0.585, padj < 0.05)
  # lfcThreshold = 1
  fc_1 <- results(run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 1)
  fc_1 <- fc_1[order(fc_1$padj),]
  fc_1 <- subset(fc_1, padj < 0.05)
  # lfcThreshold = 2
  fc_2 <- results(run, contrast = c("time_point", test_time, baseline_time), alpha = 0.05, lfcThreshold = 2)
  fc_2 <- fc_2[order(fc_2$padj),]
  fc_2 <- subset(fc_2, padj < 0.05)
  # Raw (no p-value filtering)
  fc_unfiltered <- results(run, contrast = c("time_point", test_time, baseline_time), alpha = 0.99999)
  fc_unfiltered <- fc_unfiltered[order(fc_unfiltered$padj),]
  fc_unfiltered <- subset(fc_unfiltered, padj < 1.5)
  return(list(fc_0, fc_0.1, fc_0.2, fc_0.3, fc_0.585, fc_1, fc_2, fc_unfiltered))
}

# Wayne approach
apply_wayne_classifier <- function(counts, metadata, contrast) {
  all_results <- list()
  # Subset to relevant counts for samples
  counts <- counts[,rownames(metadata)]
  rownames(counts) <- gsub("-", "_", rownames(counts))
  counts <- t(counts)
  # Record correct labels
  label <- as.vector(metadata$time_point)
  
  # Set up task
  if(length(contrast) == 2) {
    task <- "classification"
  } else {
    task <- "classification_multi"
  }
  
  # Perform classification task
  current_cv_report <- cv_report(counts, label, task = task, contrast = contrast)
  all_results[[1]] <- current_cv_report
  # Grab predictions and add metadata
  predictions <- as.data.frame(current_cv_report[[1]])
  predictions$label <- label
  predictions$AVAL <- metadata$AVAL
  predictions$viral_load_category <- metadata$viral_load_category
  predictions$MNT_Baseline <- metadata$MNT_Baseline
  predictions$MNT_Day_Minus_1_Pre_Challenge <- metadata$MNT_Day_Minus_1_Pre_Challenge
  predictions$MNT_Day_28_Post_Challenge <- metadata$MNT_Day_28_Post_Challenge
  predictions$MNT_Day_28_Change <- metadata$MNT_Day_28_Change
  predictions$HAI_Day_Minus_1_Pre_Challenge <- metadata$HAI_Day_Minus_1_Pre_Challenge
  predictions$HAI_Day_28_Post_Challenge <- metadata$HAI_Day_28_Post_Challenge
  predictions$HAI_Day_28_Change <- metadata$HAI_Day_28_Change
  predictions$NAI_Day_Minus_1_Pre_Challenge <- metadata$NAI_Day_Minus_1_Pre_Challenge
  predictions$NAI_Day_28_Post_Challenge <- metadata$NAI_Day_28_Post_Challenge
  predictions$NAI_Day_28_Change <- metadata$NAI_Day_28_Change
  predictions$sample <- rownames(metadata)
  
  # Rename columns appropriately
  prediction_colnames <- paste("time_", contrast, sep="")
  if(task == "classification_multi") {
    colnames(predictions) <- c(prediction_colnames, "correct_label", "AVAL", "viral_load_category", "MNT_Baseline",
                               "MNT_Day_Minus_1_Pre_Challenge", "MNT_Day_28_Post_Challenge",
                               "MNT_Day_28_Change", "HAI_Day_Minus_1_Pre_Challenge", "HAI_Day_28_Post_Challenge",
                               "HAI_Day_28_Change", "NAI_Day_Minus_1_Pre_Challenge", "NAI_Day_28_Post_Challenge",
                               "NAI_Day_28_Change", "sample")
  } else {
    colnames(predictions) <- c("prediction_prob", "correct_label", "AVAL", "viral_load_category", "MNT_Baseline",
                               "MNT_Day_Minus_1_Pre_Challenge", "MNT_Day_28_Post_Challenge", 
                               "MNT_Day_28_Change", "HAI_Day_Minus_1_Pre_Challenge", "HAI_Day_28_Post_Challenge", 
                               "HAI_Day_28_Change", "NAI_Day_Minus_1_Pre_Challenge", "NAI_Day_28_Post_Challenge",
                               "NAI_Day_28_Change", "sample")
    predictions[[prediction_colnames[1]]] <- 1 - predictions$prediction_prob
    predictions[[prediction_colnames[2]]] <- predictions$prediction_prob
  }
  
  # Subset to controls
  predictions_control_subset <- predictions[predictions$correct_label == contrast[1],]
  
  # Get info for barplot
  predictions_control_subset_barplot <- predictions_control_subset
  
  # Put in correct format for ggplot2
  predictions_control_subset_barplot <- predictions_control_subset_barplot %>% 
    pivot_longer(
      cols = prediction_colnames,
      names_to = "time_point",
      values_to = "probability"
    )
  
  # Set up ordering of time points
  barplot_time_point_order <- rev(prediction_colnames)
  predictions_control_subset_barplot$time_point <- factor(predictions_control_subset_barplot$time_point,
                                                          levels = barplot_time_point_order)
  
  # Set up ordering of samples
  sample_order <- predictions_control_subset_barplot[predictions_control_subset_barplot$time_point == prediction_colnames[1],]
  sample_order <- sample_order[order(sample_order$probability, decreasing = TRUE), ]$sample
  predictions_control_subset_barplot$sample <- factor(predictions_control_subset_barplot$sample,
                                                                                   levels = sample_order)
  
  # Set up legend
  barplot_time_point_order_for_legend <- gsub("time_2_D28|time_1_D28", "28 Days Post", barplot_time_point_order)
  barplot_time_point_order_for_legend <- gsub("time_2_D8|time_1_D8", "8 Days Post", barplot_time_point_order_for_legend)
  barplot_time_point_order_for_legend <- gsub("time_2_D5|time_1_D5", "5 Days Post", barplot_time_point_order_for_legend)
  barplot_time_point_order_for_legend <- gsub("time_2_D2|time_1_D2", "2 Days Post", barplot_time_point_order_for_legend)
  barplot_time_point_order_for_legend <- gsub("time_2_D_minus_1|time_1_D_minus_1", "Control", barplot_time_point_order_for_legend)
  
  # Create barplot
  predictions_control_subset_barplot <- ggplot(data = predictions_control_subset_barplot, aes(fill=time_point, y=probability, x=sample)) + 
    geom_bar(position="fill", stat="identity") +  theme(axis.text.x=element_blank(),
                                                        axis.ticks.x=element_blank()) + 
    xlab("All Control Samples") + ylab("Probability") + guides(fill=guide_legend(title="Time Points")) + 
    scale_fill_discrete(labels=barplot_time_point_order_for_legend)
  
  all_results[[2]] <- predictions_control_subset_barplot
  
  # Create correlation plots
  end_coord <- length(barplot_time_point_order) - 1
  j <- 3
  for(i in 1:end_coord) {
    # Viral load
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "AVAL", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("Viral Load")
    all_results[[j]] <- cor_plot
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$AVAL, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # MNT_Baseline
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "MNT_Baseline", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("MNT Baseline")
    all_results[[j]] <- cor_plot
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$MNT_Baseline, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # MNT_Day_Minus_1_Pre_Challenge
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "MNT_Day_Minus_1_Pre_Challenge", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("MNT (Day -1, Pre-Challenge)")
    all_results[[j]] <- cor_plot
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$MNT_Day_Minus_1_Pre_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # MNT_Day_28_Post_Challenge
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "MNT_Day_28_Post_Challenge", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("MNT (Day 28, Post-Challenge)")
    all_results[[j]] <- cor_plot
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$MNT_Day_28_Post_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # MNT_Day_28_Change
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "MNT_Day_28_Change", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("MNT (Day 28, Change from Baseline)")
    all_results[[j]] <- cor_plot
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$MNT_Day_28_Post_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # HAI_Day_Minus_1_Pre_Challenge
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "HAI_Day_Minus_1_Pre_Challenge", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("HAI (Day -1, Pre-Challenge)")
    all_results[[j]] <- cor_plot
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$HAI_Day_Minus_1_Pre_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # HAI_Day_28_Post_Challenge
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "HAI_Day_28_Post_Challenge", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("HAI (Day 28, Post-Challenge)")
    all_results[[j]] <- cor_plot
    
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$HAI_Day_28_Post_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # HAI_Day_28_Change
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "HAI_Day_28_Change", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("HAI (Day 28, Change from Baseline)")
    all_results[[j]] <- cor_plot
    
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$HAI_Day_28_Post_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # NAI_Day_Minus_1_Pre_Challenge
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "NAI_Day_Minus_1_Pre_Challenge", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("NAI (Day -1, Pre-Challenge)")
    all_results[[j]] <- cor_plot
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$NAI_Day_Minus_1_Pre_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # NAI_Day_28_Post_Challenge
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "NAI_Day_28_Post_Challenge", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("NAI (Day 28, Post-Challenge)")
    all_results[[j]] <- cor_plot
    
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$NAI_Day_28_Post_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
    # NAI_Day_28_Change
    cor_plot <- ggplot(data = predictions_control_subset, mapping = aes_string(x = barplot_time_point_order[[i]], y = "NAI_Day_28_Change", color = "viral_load_category", group = 1)) +
      geom_point(size = 2) +
      sm_statCorr(corr_method = "spearman") + xlab(paste0("Probability of Misclassification as ", barplot_time_point_order_for_legend[[i]])) + ylab("NAI (Day 28, Change from Baseline)")
    all_results[[j]] <- cor_plot
    
    
    correlation_val <- cor.test(predictions_control_subset[[barplot_time_point_order[[i]]]], predictions_control_subset$NAI_Day_28_Post_Challenge, method = "spearman")
    print(correlation_val$estimate)
    print(correlation_val$p.value)
    
    j <- j + 1
    
  }
  
  return(all_results)
}

# Create classifier based on multiclassPairs method
create_multiclassPairs_classifier_for_bulk_data <- function(counts, metadata) {
  counts <- counts[,rownames(metadata)]
  rownames(counts) <- gsub("-", "_", rownames(counts))
  time_points <- metadata$time_point
  classes <- as.vector(unique(time_points))
  
  # Create classifier object
  classifier_object <- ReadData(Data = counts,
                                Labels = time_points,
                                verbose = TRUE)
  
  # Get list of filtered genes
  filtered_genes <- filter_genes_TSP(data_object = classifier_object,
                                     filter = "one_vs_one",
                                     platform_wise = FALSE,
                                     featureNo = 1000,
                                     UpDown = TRUE,
                                     verbose = TRUE)
  
  # Train classifier
  classifier <- train_one_vs_rest_TSP(data_object = classifier_object,
                                      filtered_genes = filtered_genes,
                                      k_range = 5:50,
                                      include_pivot = FALSE,
                                      one_vs_one_scores = TRUE,
                                      platform_wise_scores = FALSE,
                                      seed = get_speedi_seed(),
                                      verbose = TRUE)
  
  # Test classifier on training data
  results_train <- predict_one_vs_rest_TSP(classifier = classifier,
                                           Data = classifier_object,
                                           tolerate_missed_genes = TRUE,
                                           weighted_votes = TRUE,
                                           classes = classes,
                                           verbose = TRUE)
  
  # Check that predictions are accurate
  results_train$max_score == time_points
  
  caret::confusionMatrix(data = factor(results_train$max_score, 
                                       levels = unique(classifier_object$data$Labels)),
                         reference = factor(classifier_object$data$Labels, 
                                            levels = unique(classifier_object$data$Labels)),
                         mode="everything")
  return(list(classifier, results_train))
}
