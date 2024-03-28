run_differential_expression_cluster <- function(sc_obj, marker_dir) {
  print(paste0("Performing differential expression for each cluster"))
  # Make reading data parallel
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }
  
  n.cores <- 16
  
  message(paste0("Number of cores: ", n.cores))
  
  registerDoMC(n.cores)
  message("Begin parallelizing...")
  cluster_ids <- unique(sc_obj$seurat_clusters)
  
  compiled_output <- foreach(
    i = 1:length(cluster_ids),
    .packages = c("Seurat", "base")
  ) %dopar% {
    cluster_id <- cluster_ids[[i]]
    cluster.markers <- FindMarkers(sc_obj, test.use="LR", latent.vars = 'subject_id', ident.1 = cluster_id, assay = "SCT", recorrect_umi = FALSE)
    write.table(cluster.markers, paste0(marker_dir, "markers_", cluster_id, ".txt"), quote = FALSE, sep = "\t")
    return(i)
  }
  message("All done!")
}

run_differential_expression_controlling_for_subject_id <- function(sc_obj, analysis_dir, sample_metadata_for_SPEEDI_df, group, cell_types) {
  print(paste0("Performing differential expression for group ", group, " for each cell type (controlling for subject ID)"))
  expected_num_samples <- length(unique(sc_obj$sample))
  
  if(group == "viral_load_category") {
    first_group <- "high"
    second_group <- "low"
  } else if(group == "time_point") {
    first_group <- "D28"
    second_group <- "D_minus_1"
  } else if(group == "sex") {
    first_group <- "M"
    second_group <- "F"
  }
  
  # Make reading data parallel
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }
  
  n.cores <- 12
  
  message(paste0("Number of cores: ", n.cores))
  
  registerDoMC(n.cores)
  message("Begin parallelizing...")
  
  compiled_output <- foreach(
    i = 1:length(cell_types),
    .packages = c("Seurat", "base")
  ) %dopar% {
    current_cell_type <- cell_types[[i]]
    cell_type_for_file_name <- sub(" ", "_", current_cell_type)
    print(current_cell_type)
    if(current_cell_type %in% unique(sc_obj$predicted_celltype_majority_vote)) {
      idxPass <- which(sc_obj$predicted_celltype_majority_vote %in% current_cell_type)
    } else {
      idxPass <- which(sc_obj$magical_cell_types %in% current_cell_type)
      print("This cell type is for MAGICAL processing")
    }
    if(length(idxPass) == 0) {
      print("No cells found, so skipping cell type")
    } else {
      cellsPass <- names(sc_obj$orig.ident[idxPass])
      cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
      if(length(unique(cells_subset$sample)) != expected_num_samples) {
        print("Missing at least one sample for analysis, so skipping cell type")
      } else {
        DefaultAssay(cells_subset) <- "SCT"
        Idents(cells_subset) <- group
        # Run SC DE with no thresholding and write to file
        current_de <- FindMarkers(cells_subset, test.use="LR", latent.vars = 'subject_id', ident.1 = first_group, ident.2 = second_group, logfc.threshold = 0, min.pct = 0, assay = "SCT", recorrect_umi = FALSE)
        write.table(current_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_sc_unfiltered.tsv"), quote = FALSE, sep = "\t")
        # Filter results and write to file
        current_de <- current_de[current_de$p_val_adj < 0.05,]
        current_de <- current_de[abs(current_de$avg_log2FC) >= 0.1,]
        current_de <- current_de[current_de$pct.1 >= 0.1 | current_de$pct.2 >= 0.1, ]
        if(nrow(current_de) > 0) {
          write.table(current_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_sc.tsv"), quote = FALSE, sep = "\t")
        }
        # Run pseudobulk
        Seurat::DefaultAssay(cells_subset) <- "RNA"
        pseudobulk_counts <- SPEEDI::create_pseudobulk_counts(cells_subset, log_flag = FALSE)
        pseudobulk_metadata <- sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$subject_id %in% cells_subset$subject_id,]
        pseudobulk_metadata$aliquots <- rownames(pseudobulk_metadata)
        pseudobulk_metadata <- pseudobulk_metadata[match(colnames(pseudobulk_counts), pseudobulk_metadata$aliquots),]
        pseudobulk_analysis <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_counts, colData = pseudobulk_metadata, design = stats::formula(paste("~ subject_id + ",group)))
        pseudobulk_analysis <- DESeq2::DESeq(pseudobulk_analysis)
        pseudobulk_analysis_results_contrast <- utils::tail(DESeq2::resultsNames(pseudobulk_analysis), n=1)
        print(pseudobulk_analysis_results_contrast)
        pseudobulk_analysis_results <- DESeq2::results(pseudobulk_analysis, name=pseudobulk_analysis_results_contrast)
        write.table(pseudobulk_analysis_results, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_pseudobulk_unfiltered.tsv"), quote = FALSE, sep = "\t")
        pseudobulk_analysis_results <- pseudobulk_analysis_results[rowSums(is.na(pseudobulk_analysis_results)) == 0, ] # Remove NAs
        pseudobulk_analysis_results <- pseudobulk_analysis_results[pseudobulk_analysis_results$pvalue < 0.05,]
        if(nrow(pseudobulk_analysis_results) > 0) {
          write.table(pseudobulk_analysis_results, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_pseudobulk.tsv"), quote = FALSE, sep = "\t")
        }
        # Find SC genes that are pseudobulk filtered
        final_genes <- intersect(rownames(current_de), rownames(pseudobulk_analysis_results))
        # Record information about remaining genes in final_current_de
        final_current_de <- data.frame(Cell_Type = character(), Gene_Name = character(), sc_pval_adj = character(), sc_log2FC = character(), pseudo_bulk_pval = character(),
                                       pseudo_bulk_log2FC = character())
        for(current_gene in final_genes) {
          current_sc_pval_adj <- current_de[rownames(current_de) == current_gene,]$p_val_adj
          current_sc_log2FC <- current_de[rownames(current_de) == current_gene,]$avg_log2FC
          current_pseudo_bulk_pval <- pseudobulk_analysis_results[rownames(pseudobulk_analysis_results) == current_gene,]$pvalue
          current_pseudo_bulk_log2FC <- pseudobulk_analysis_results[rownames(pseudobulk_analysis_results) == current_gene,]$log2FoldChange
          current_row <- data.frame(current_cell_type, current_gene, current_sc_pval_adj, current_sc_log2FC, current_pseudo_bulk_pval, current_pseudo_bulk_log2FC)
          names(current_row) <- c("Cell_Type", "Gene_Name", "sc_pval_adj", "sc_log2FC", "pseudo_bulk_pval", "pseudo_bulk_log2FC")
          final_current_de <- rbind(final_current_de, current_row)
        }
        if(nrow(final_current_de) > 0) {
          write.table(final_current_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_final.tsv"), quote = FALSE, sep = "\t")
        }
      }
    }
    return(current_cell_type)
  }
  
  # Combine final results into one table
  final_pseudobulk_corrected_de <- list()
  for(current_cell_type in cell_types) {
    cell_type_for_file_name <- sub(" ", "_", current_cell_type)
    if(file.exists(paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_final.tsv"))) {
      current_cell_type_file <- read.table(paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_final.tsv"), sep = "\t", header = TRUE)
      final_pseudobulk_corrected_de[[current_cell_type]] <- current_cell_type_file
    }
  }
  final_pseudobulk_corrected_de <- do.call(rbind, final_pseudobulk_corrected_de)
  
  write.table(final_pseudobulk_corrected_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", group, ".final.list.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  pos_final_pseudobulk_corrected_de <- final_pseudobulk_corrected_de[final_pseudobulk_corrected_de$sc_log2FC > 0,]
  write.table(pos_final_pseudobulk_corrected_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", group, ".final.pos.list.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  neg_final_pseudobulk_corrected_de <- final_pseudobulk_corrected_de[final_pseudobulk_corrected_de$sc_log2FC < 0,]
  write.table(neg_final_pseudobulk_corrected_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", group, ".final.neg.list.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  print(paste0("Done performing differential expression for group ", group, " for each cell type"))
}

