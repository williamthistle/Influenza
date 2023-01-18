# Method to set up bulk analysis
setup_bulk_analysis=function() {
  # Read in count and metadata files
  base_dir <<- "C:/Users/wat2/Documents/GitHub/Influenza/"
  data_dir <<- "C:/Users/wat2/Documents/local_data_files/"
  gene_counts <<- fread(paste0(data_dir, "rsem_genes_count.processed.txt"), header = T, sep = "\t")
  all_metadata_file <<- paste0(base_dir, "all_metadata_sheet.tsv")
  all_metadata <<- read.table(all_metadata_file, header = TRUE, sep = "\t")
  # Currently, only capturing viral load for placebo subjects
  viral_load_file <<- paste0(base_dir, "bulk_RNA_viral_load.tsv")
  viral_load <<- read.table(viral_load_file, header = TRUE, sep = "\t")
  viral_load_primary <<- viral_load[viral_load$PARAMCD == "QPCRAUC",]
  viral_load_primary <<- viral_load_primary[viral_load_primary$TRT01A == "PLACEBO",]
  viral_load_primary$AVAL <<- as.numeric(viral_load_primary$AVAL)
  # Organize by viral load (high to low).
  viral_load_primary <<- viral_load_primary[order(viral_load_primary$AVAL, decreasing = TRUE),]
  # Take gene_id column from gene_counts and use contents as the rownames of gene_counts
  gene.row.names <<- as.character(gene_counts$gene_id)
  gene_counts <<- gene_counts[,2:ncol(gene_counts)]
  gene_counts <<- as.data.frame(gene_counts)
  rownames(gene_counts) <<- gene.row.names
  # Only keep metadata for bulk RNA-seq aliquots
  all_metadata <<- all_metadata[all_metadata$bulkRNA_seq == TRUE,]
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
  # Divide metadata into placebo and vaccinated
  placebo_metadata <<- all_metadata[all_metadata$treatment == "PLACEBO",]
  vaccinated_metadata <<- all_metadata[all_metadata$treatment == "MVA-NP+M1",]
  # Find placebo-associated and vaccinated-associated gene_counts
  kept_aliquots <<- placebo_metadata$aliquot_id
  placebo_counts <<- gene_counts[kept_aliquots]
  kept_aliquots <<- vaccinated_metadata$aliquot_id
  vaccinated_counts <<- gene_counts[kept_aliquots]
  # Sort columns in gene_counts and rows for each so they're in same order (for DESeq2)
  colnames(gene_counts) <<- sort(colnames(gene_counts))
  rownames(all_metadata) <<- all_metadata$aliquot_id
  rownames(all_metadata) <<- sort(rownames(all_metadata))
  colnames(placebo_counts) <<- sort(colnames(placebo_counts))
  rownames(placebo_metadata) <<- placebo_metadata$aliquot_id
  rownames(placebo_metadata) <<- sort(rownames(placebo_metadata))
  colnames(vaccinated_counts) <<- sort(colnames(vaccinated_counts))
  rownames(vaccinated_metadata) <<- vaccinated_metadata$aliquot_id
  rownames(vaccinated_metadata) <<- sort(rownames(vaccinated_metadata))
  # Drop aliquot ID column (it's stored in rownames)
  all_metadata <<- subset(all_metadata, select = -c(aliquot_id))
  placebo_metadata <<- subset(placebo_metadata, select = -c(aliquot_id))
  vaccinated_metadata <<- subset(vaccinated_metadata, select = -c(aliquot_id))
  # Probably OK to round expected counts from RSEM data. DESeq2 expects integers
  gene_counts <<- round(gene_counts)
  placebo_counts <<- round(placebo_counts)
  vaccinated_counts <<- round(vaccinated_counts)
  # Currently not filtering sex associated genes
  #sex_associated_genes <<- find_sex_associated_genes(paste0(data_dir, "sex_associated_genes/"))
}

run_deseq_bulk_analysis=function(sample_type, counts, metadata, test_time, baseline_time, output_dir, output_name_prefix=NA) {
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
  return(list(current_analysis_results, current_analysis_results_1.5, current_analysis_results_2, current_analysis_results_4))
}

# Method to find sex associated genes - not currently used
find_sex_associated_genes=function(sex_associated_dir, padj_threshold = 0.05, log2fc_threshold = 0.1) {
  sex_associated_gene_files <- list.files(sex_associated_dir, pattern = ".csv")
  sex_associated_genes <- c()
  number_of_studies <- length(sex_associated_gene_files)
  for (sex_associated_gene_file in sex_associated_gene_files) {
    print(sex_associated_gene_file)
    sex_associated_gene_file_contents <- read.csv(paste0(sex_associated_dir, sex_associated_gene_file))
    row.names <- as.character(sex_associated_gene_file_contents$X)
    sex_associated_gene_file_contents <- sex_associated_gene_file_contents[,2:ncol(sex_associated_gene_file_contents)]
    sex_associated_gene_file_contents <- as.data.frame(sex_associated_gene_file_contents)
    rownames(sex_associated_gene_file_contents) <- row.names
    sex_associated_gene_file_contents <- subset(sex_associated_gene_file_contents, padj < padj_threshold & 
                                                  abs(log2FoldChange) > log2fc_threshold)
    for (entry in rownames(sex_associated_gene_file_contents)) {
      if(entry %in% genemap$ensembl_gene_id) {
        sex_associated_genes <- c(sex_associated_genes, unique(genemap[genemap$ensembl_gene_id == entry,]$hgnc_symbol))
      } else{
        print(entry)
      }
    }
  }
  # Remove repeated genes and the "" entries (no HGNC symbol)
  #sex_associated_genes <- unique(sex_associated_genes[sex_associated_genes != ""])
  return(sex_associated_genes)
}
