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
    current_cell_types <- gene_table[gene_table$Gene_Name %in% current_genes,]$Cell_Type
    current_cell_type_summary <- sort(table(current_cell_types), decreasing = TRUE)
    formatted_elements <- character(length(current_cell_type_summary))
    
    # Loop through the named vector and format each element
    for (i in 1:length(current_cell_type_summary)) {
      formatted_elements[i] <- paste0(names(current_cell_type_summary)[i], " (", current_cell_type_summary[i], ")")
    }
    
    # Combine the formatted elements into a single string
    final_cell_type_str <- paste(formatted_elements, collapse = ", ")
    overall_module_cell_types <- c(overall_module_cell_types, final_cell_type_str)
  }
  overall_module_data$Cell_Types <- overall_module_cell_types
  
  # Individual GO terms
  individual_go_term_cell_types <- c()
  for(current_genes in individual_go_term_data$TERM_GENES) {
    current_genes <- strsplit(current_genes, ",")[[1]]
    current_cell_types <- gene_table[gene_table$Gene_Name %in% current_genes,]$Cell_Type
    current_cell_type_summary <- sort(table(current_cell_types), decreasing = TRUE)
    formatted_elements <- character(length(current_cell_type_summary))
    
    # Loop through the named vector and format each element
    for (i in 1:length(current_cell_type_summary)) {
      formatted_elements[i] <- paste0(names(current_cell_type_summary)[i], " (", current_cell_type_summary[i], ")")
    }
    
    # Combine the formatted elements into a single string
    final_cell_type_str <- paste(formatted_elements, collapse = ", ")
    individual_go_term_cell_types <- c(individual_go_term_cell_types, final_cell_type_str)
  }
  individual_go_term_data$Cell_Types <- individual_go_term_cell_types
  
  # Return parsed info with cell types
  return(list(overall_module_data, individual_go_term_data))
}