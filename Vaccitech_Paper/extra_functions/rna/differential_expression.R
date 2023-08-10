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
    cluster.markers <- FindMarkers(sc_obj, ident.1 = cluster_id, assay = "SCT", recorrect_umi = FALSE)
    write.table(cluster.markers, paste0(marker_dir, "markers_", cluster_id, ".txt"), quote = FALSE, sep = "\t")
    return(i)
  }
  message("All done!")
}

run_differential_expression_group <- function(sc_obj, analysis_dir, group) {
  print(paste0("Performing differential expression for group ", group, " for each cell type"))
  
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
  
  all_cell_types <- union(unique(sc_obj$predicted_celltype_majority_vote), unique(sc_obj$magical_cell_types))
  
  compiled_output <- foreach(
    i = 1:length(all_cell_types),
    .packages = c("Seurat", "base")
  ) %dopar% {
    cell_type <- all_cell_types[[i]]
    print(cell_type)
    if(cell_type %in% unique(sc_obj$predicted_celltype_majority_vote)) {
      idxPass <- which(sc_obj$predicted_celltype_majority_vote %in% cell_type)
    } else {
      idxPass <- which(sc_obj$magical_cell_types %in% cell_type)
      print("This cell type is for MAGICAL processing")
    }
    print(length(idxPass))
    cellsPass <- names(sc_obj$orig.ident[idxPass])
    cells_subset <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    # TODO: Make this print table for relevant group
    #print(table(cells_subset$viral_load))
    DefaultAssay(cells_subset) <- "SCT"
    Idents(cells_subset) <- group
    if(group == "viral_load") {
      first_group <- "high"
      second_group <- "low"
    } else if(group == "time_point") {
      first_group <- "D28"
      second_group <- "D_minus_1"
    } else if(group == "sex") {
      first_group <- "M"
      second_group <- "F"
    }
    #diff_markers = FindMarkers(cells_subset, test.use="LR", latent.vars = 'subject_id', ident.1 = first_group, ident.2 = second_group, logfc.threshold = 0, min.pct = 0, assay = "SCT", recorrect_umi = FALSE)
    diff_markers <- FindMarkers(cells_subset, ident.1 = first_group, ident.2 = second_group, logfc.threshold = 0.1, min.pct = 0.1, assay = "SCT", recorrect_umi = FALSE)
    cell_type <- sub(" ", "_", cell_type)
    write.csv(diff_markers, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type, "-", group, ".csv"), quote = FALSE)
    return(i)
  }
  print(paste0("Done performing differential expression for group ", group, " for each cell type"))
}

run_differential_expression_controlling_for_subject_id <- function(sc_obj, analysis_dir, sample_metadata_for_SPEEDI_df, group, cell_types) {
  print(paste0("Performing differential expression for group ", group, " for each cell type (controlling for subject ID)"))
  final_current_de <- data.frame(Cell_Type = character(), Gene_Name = character(), sc_pval_adj = character(), sc_log2FC = character(), pseudo_bulk_pval = character(),
                                 pseudo_bulk_log2FC = character())  
  for(current_cell_type in cell_types) {
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
      DefaultAssay(cells_subset) <- "SCT"
      Idents(cells_subset) <- group
      if(group == "viral_load") {
        first_group <- "high"
        second_group <- "low"
      } else if(group == "time_point") {
        first_group <- "D28"
        second_group <- "D_minus_1"
      } else if(group == "sex") {
        first_group <- "M"
        second_group <- "F"
      }
      current_de = FindMarkers(cells_subset, test.use="LR", latent.vars = 'subject_id', ident.1 = first_group, ident.2 = second_group, logfc.threshold = 0.1, min.pct = 0.1, assay = "SCT", recorrect_umi = FALSE)
      current_de <- current_de[current_de$p_val_adj < 0.05,]
      cell_type_for_file_name <- sub(" ", "_", current_cell_type)
      write.table(current_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_sc.tsv"), quote = FALSE, sep = "\t")
      # Run pseudobulk
      Seurat::DefaultAssay(cells_subset) <- "RNA"
      pseudobulk_counts <- SPEEDI::create_pseudobulk_counts(cells_subset, log_flag = FALSE)
      pseudobulk_metadata <- sample_metadata_for_SPEEDI_df[sample_metadata_for_SPEEDI_df$subject_id %in% cells_subset$subject_id,]
      pseudobulk_metadata$aliquots <- rownames(pseudobulk_metadata)
      pseudobulk_metadata <- pseudobulk_metadata[match(colnames(pseudobulk_counts), pseudobulk_metadata$aliquots),]
      pseudobulk_analysis <- DESeq2::DESeqDataSetFromMatrix(countData = pseudobulk_counts, colData = pseudobulk_metadata, design = stats::formula(paste("~ subject_id + ",group)))
      pseudobulk_analysis <- DESeq2::DESeq(pseudobulk_analysis)
      pseudobulk_analysis_results_contrast <- utils::tail(DESeq2::resultsNames(pseudobulk_analysis), n=1)
      pseudobulk_analysis_results <- DESeq2::results(pseudobulk_analysis, name=pseudobulk_analysis_results_contrast)
      pseudobulk_analysis_results <- pseudobulk_analysis_results[rowSums(is.na(pseudobulk_analysis_results)) == 0, ] # Remove NAs
      pseudobulk_analysis_results <- pseudobulk_analysis_results[pseudobulk_analysis_results$pvalue < 0.05,]
      pseudobulk_analysis_results <- pseudobulk_analysis_results[pseudobulk_analysis_results$log2FoldChange < -0.3 | pseudobulk_analysis_results$log2FoldChange > 0.3,]
      write.table(pseudobulk_analysis_results, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", cell_type_for_file_name, "-", group, "-controlling_for_subject_id_pseudobulk.tsv"), quote = FALSE, sep = "\t")
      final_genes <- intersect(rownames(current_de), rownames(pseudobulk_analysis_results))
      # Record information about remaining genes in final_current_de
      for(current_gene in final_genes) {
        current_sc_pval_adj <- current_de[rownames(current_de) == current_gene,]$p_val_adj
        current_sc_log2FC <- current_de[rownames(current_de) == current_gene,]$avg_log2FC
        current_pseudo_bulk_pval <- pseudobulk_analysis_results[rownames(pseudobulk_analysis_results) == current_gene,]$pvalue
        current_pseudo_bulk_log2FC <- pseudobulk_analysis_results[rownames(pseudobulk_analysis_results) == current_gene,]$log2FoldChange
        current_row <- data.frame(current_cell_type, current_gene, current_sc_pval_adj, current_sc_log2FC, current_pseudo_bulk_pval, current_pseudo_bulk_log2FC)
        names(current_row) <- c("Cell_Type", "Gene_Name", "sc_pval_adj", "sc_log2FC", "pseudo_bulk_pval", "pseudo_bulk_log2FC")
        final_current_de <- rbind(final_current_de, current_row)
      }
    }
  }
  write.table(final_current_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", group, ".final.list.tsv"), quote = FALSE, sep = "\t")
  pos_final_current_de <- final_current_de[final_current_de$sc_log2FC > 0,]
  write.table(pos_final_current_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", group, ".final.pos.list.tsv"), quote = FALSE, sep = "\t")
  neg_final_current_de <- final_current_de[final_current_de$sc_log2FC < 0,]
  write.table(neg_final_current_de, paste0(analysis_dir, first_group, "-vs-", second_group, "-degs-", group, ".final.neg.list.tsv"), quote = FALSE, sep = "\t")
  print(paste0("Done performing differential expression for group ", group, " for each cell type"))
}

