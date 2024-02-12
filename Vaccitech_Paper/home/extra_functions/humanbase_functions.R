assign_cell_types_to_humanbase_results <- function(humanbase_file_path, gene_table) {
  # Table parsed into parts (using \t\t\t\t\t as delimiter)
  parsed_table <- list()
  
  # Read the file as table_lines
  table_lines <- readLines(humanbase_file_path)
  
  # Initialize an empty string to store the current component
  current_component <- ""
  
  # Loop through each line
  for (line in table_lines) {
    if (line == "\t\t\t\t\t") {
      # If the line is the separator, add the current component to the list
      if (nchar(current_component) > 0) {
        parsed_table <- append(parsed_table, list(current_component))
      }
      # Reset the current component
      current_component <- ""
    } else {
      # Append the line to the current component
      current_component <- paste(current_component, line, sep = "\n")
    }
  }
  
  # Add the last component to the list
  if (nchar(current_component) > 0) {
    parsed_table <- append(parsed_table, list(current_component))
  }
  # Info about overall genes for each module
  overall_module_data <- read.table(text = parsed_table[[2]], header = TRUE, sep = "\t")
  overall_module_data <- overall_module_data[,1:2]
  # Info about individual GO terms for each module and associated genes
  individual_go_term_data <- read.table(text = parsed_table[[3]], header = TRUE, sep = "\t")
  
  # Overall modules
  overall_module_cell_types <- c()
  for(current_genes in overall_module_data$CLUSTER_GENES) {
    current_genes <- strsplit(current_genes, ",")[[1]]
    final_gene_list <- ""
    for(current_gene in current_genes) {
      current_cell_types <- gene_table[gene_table$Gene_Name %in% current_gene,]$Cell_Type
      current_gene <- paste0(current_gene, " (", paste0(current_cell_types, collapse = ", "), ")")
      if(nchar(final_gene_list) == 0) {
        final_gene_list <- current_gene
      } else {
        final_gene_list <- paste0(final_gene_list, ", ", current_gene)
      }
    }
    overall_module_cell_types <- c(overall_module_cell_types, final_gene_list)
  }
  overall_module_data$CLUSTER_GENES <- overall_module_cell_types
  
  # Individual GO terms
  individual_go_term_cell_types <- c()
  for(current_genes in individual_go_term_data$TERM_GENES) {
    current_genes <- strsplit(current_genes, ",")[[1]]
    final_gene_list <- ""
    for(current_gene in current_genes) {
      current_cell_types <- gene_table[gene_table$Gene_Name %in% current_gene,]$Cell_Type
      current_gene <- paste0(current_gene, " (", paste0(current_cell_types, collapse = ", "), ")")
      if(nchar(final_gene_list) == 0) {
        final_gene_list <- current_gene
      } else {
        final_gene_list <- paste0(final_gene_list, ", ", current_gene)
      }
    }
    individual_go_term_cell_types <- c(individual_go_term_cell_types, final_gene_list)
  }
  individual_go_term_data$TERM_GENES <- individual_go_term_cell_types
  
  # Return parsed info with cell types
  return(list(overall_module_data, individual_go_term_data))
}

assign_cell_types_to_pathway_results <- function(pathway_info, gene_table) {
  individual_pathway_term_cell_types <- c()
  for(current_genes in pathway_info$Genes) {
    current_genes <- strsplit(current_genes, ";")[[1]]
    current_cell_types <- gene_table[gene_table$Gene_Name %in% current_genes,]$Cell_Type
    current_cell_type_summary <- sort(table(current_cell_types), decreasing = TRUE)
    formatted_elements <- character(length(current_cell_type_summary))
    
    # Loop through the named vector and format each element
    for (i in 1:length(current_cell_type_summary)) {
      formatted_elements[i] <- paste0(names(current_cell_type_summary)[i], " (", current_cell_type_summary[i], ")")
    }
    
    # Combine the formatted elements into a single string
    final_cell_type_str <- paste(formatted_elements, collapse = ", ")
    individual_pathway_term_cell_types <- c(individual_pathway_term_cell_types, final_cell_type_str)
  }
  pathway_info$Cell_Types <- individual_pathway_term_cell_types
  
  # Return parsed info with cell types
  return(pathway_info)
}

run_fmd_on_flu_data <- function(list_of_deseq_results_at_different_fc_thresholds) {
  fmd_results <- list()
  log_fc_thresholds <- c(0,0.1,0.2,0.3,0.585,1,2)
  index <- 1
  for(current_result in list_of_deseq_results_at_different_fc_thresholds) {
    log_fc_threshold <- log_fc_thresholds[index]
    current_gene_list <- rownames(current_result[current_result$log2FoldChange > 0,])
    if(length(current_gene_list) > 0 && length(current_gene_list) < 2000) {
      current_fmd_result <- SPEEDI::RunFMD_RNA(current_gene_list, "blood")
      fmd_result_name <- paste0("logfc_", log_fc_threshold, "_upregulated")
      fmd_results[[fmd_result_name]] <- current_fmd_result
    }
    current_gene_list <- rownames(current_result[current_result$log2FoldChange < 0,])
    if(length(current_gene_list) > 0 && length(current_gene_list) < 2000) {
      current_fmd_result <- SPEEDI::RunFMD_RNA(current_gene_list, "blood")
      fmd_result_name <- paste0("logfc_", log_fc_threshold, "_downregulated")
      fmd_results[[fmd_result_name]] <- current_fmd_result
    }
    index <- index + 1
  }
  fmd_links <- c()
  for(entry in fmd_results) { 
    fmd_links <- c(fmd_links, entry[[1]])
  }
  fmd_result_df <- data.frame(fc_threshold_name = names(fmd_results), link = fmd_links)
  return(fmd_result_df)
}

get_modules_for_downstream_analysis <- function(humanbase_file_path) {
  final_modules <- c()
  # Table parsed into parts (using \t\t\t\t\t as delimiter)
  parsed_table <- list()
  
  # Read the file as table_lines
  table_lines <- readLines(humanbase_file_path)
  
  # Initialize an empty string to store the current component
  current_component <- ""
  
  # Loop through each line
  for (line in table_lines) {
    if (line == "\t\t\t\t\t") {
      # If the line is the separator, add the current component to the list
      if (nchar(current_component) > 0) {
        parsed_table <- append(parsed_table, list(current_component))
      }
      # Reset the current component
      current_component <- ""
    } else {
      # Append the line to the current component
      current_component <- paste(current_component, line, sep = "\n")
    }
  }
  
  # Add the last component to the list
  if (nchar(current_component) > 0) {
    parsed_table <- append(parsed_table, list(current_component))
  }
  # Info about overall genes for each module
  overall_module_data <- read.table(text = parsed_table[[2]], header = TRUE, sep = "\t")
  overall_module_data <- overall_module_data[,1:2]
  
  # Get list of modules with >= 15 genes
  for(current_row in 1:nrow(overall_module_data)) {
    current_module <- overall_module_data[current_row,]$CLUSTER_NAME
    current_genes <- overall_module_data[current_row,]$CLUSTER_GENES
    current_genes <- strsplit(current_genes, ",")[[1]]
    if(length(current_genes) >= 15) {
      final_modules <- c(final_modules, current_module)
    }
  }
  return(final_modules)
}


get_lists_of_commands <- function(analysis_dir, direction = "upregulated") {
  output_dir <- paste0(analysis_dir, "output/")
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  compare_module_dir <- paste0(output_dir, "module_comparisons/")
  if (!dir.exists(compare_module_dir)) {dir.create(compare_module_dir)}
  pathway_analyses_dir <- paste0(output_dir, "pathway_analyses/")
  if (!dir.exists(pathway_analyses_dir)) {dir.create(pathway_analyses_dir)}
  analyzed_files <- list.files(analysis_dir)
  analyzed_files <- analyzed_files[!analyzed_files %in% "output"]
  common_suffix <- unique(sub(paste0(".*_([^_]+_", direction, ".tsv)$"), "\\1", analyzed_files))
  common_suffix_without_file_ext <- tools::file_path_sans_ext(common_suffix)
  terms_for_upregulated_scRNA_comparisons <- sub(paste0("_([^_]+_", direction, ".tsv)$"), "", analyzed_files)
  # Generate all possible pairwise combinations
  combinations <- combn(terms_for_upregulated_scRNA_comparisons, 2, simplify = TRUE)
  total_prefixes <- unique(c(combinations[1, ], combinations[2, ]))
  # We only want to run module summary on modules that have at least 15 genes
  hb_module_summaries <- list()
  for(prefix in total_prefixes) {
    hb_module_summaries[[prefix]] <- get_modules_for_downstream_analysis(paste0(analysis_dir, prefix, "_", common_suffix))
  }
  # Create vector of commands for compareHBfiles
  commands <- sprintf("python compareHBfiles.py -n1 %s%s_%s -n2 %s%s_%s -tsim 0.2 > %s%s_vs_%s_%s",
                      analysis_dir, combinations[1, ], common_suffix, analysis_dir, combinations[2, ], common_suffix, 
                      compare_module_dir, combinations[1, ], combinations[2, ], common_suffix)
  
  # Print the commands without quotes and leading numbers
  cat(commands, sep = "\n")
  
  # Create vector of commands for HBmodulesummary
  module_summary_commands <- c()
  for(prefix in names(hb_module_summaries)) {
    current_modules <- hb_module_summaries[[prefix]]
    for(current_module in current_modules) {
      module_summary_commands <- c(module_summary_commands, paste0("python HBmodulesummary.py -n ", analysis_dir, prefix, "_", common_suffix, " -gl Reactome_2022 -m ",
                                                                   current_module, " -o ", pathway_analyses_dir, prefix, "_", common_suffix_without_file_ext, "_", 
                                                                   current_module, "_gsea_info_output.tsv"))
    }
  }
  
  # Print the commands without quotes and leading numbers
  cat(module_summary_commands, sep = "\n")
}



