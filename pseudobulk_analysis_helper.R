library(DESeq2)
library(org.Hs.eg.db)

grab_transformed_pseudobulk_counts <- function(pseudobulk_dir, cell_types) {
  total_pseudobulk_df <- NULL
  if(grepl("Single Cell", pseudobulk_dir)) {
    # We start with B and add the rest - create sample-level pseudobulk counts
    total_pseudobulk_df <- read.table(paste0(pseudobulk_dir, "pseudo_bulk_RNA_count_B.txt"), header = TRUE, sep = "\t")
    for(cell_type in cell_types) {
      if(file.exists(paste0("pseudo_bulk_RNA_count_", cell_type, ".txt"))) {
        cell_type_pseudobulk_df <- read.table(paste0(pseudobulk_dir, "pseudo_bulk_RNA_count_", cell_type, ".txt"), header = TRUE, sep = "\t")
        for(current_col_index in 2:ncol(cell_type_pseudobulk_df)) {
          total_pseudobulk_df[,current_col_index] <- total_pseudobulk_df[,current_col_index] + cell_type_pseudobulk_df[,current_col_index]
        }
      }
    }
  } else {
    # We start with B and add the rest - create sample-level pseudobulk counts
    total_pseudobulk_df <- read.table(paste0(pseudobulk_dir, "pseudo_bulk_RNA_count_B.csv"), header = TRUE, sep = ",")
    for(cell_type in cell_types) {
      if(file.exists(paste0("pseudo_bulk_RNA_count_", cell_type, ".csv"))) {
        cell_type_pseudobulk_df <- read.table(paste0(pseudobulk_dir, "pseudo_bulk_RNA_count_", cell_type, ".csv"), header = TRUE, sep = ",")
        for(current_col_index in 2:ncol(cell_type_pseudobulk_df)) {
          total_pseudobulk_df[,current_col_index] <- total_pseudobulk_df[,current_col_index] + cell_type_pseudobulk_df[,current_col_index]
        }
      }
    }
  }
  
  
  # Grab gene names and set them to row names
  rownames(total_pseudobulk_df) <- total_pseudobulk_df$X
  total_pseudobulk_df <- total_pseudobulk_df[,-1]
  
  # Use DESeq2 varianceStabilizingTransformation function to log-normalize gene counts (as requested by MetaIntegrator)
  # I picked varianceStabilizingTransformation instead of rlog because rlog gives negative values (MetaIntegrator doesn't like those)
  # "Also, negative gene expression values are problematic for geometric mean calculation."
  # IMPORTANT NOTE: For old SC data, we have SCT counts that have ALREADY been normalized in some way.
  # It may make sense to start with RNA data for single cell and multiome.
  total_pseudobulk_df_log <- varianceStabilizingTransformation(as.matrix(total_pseudobulk_df))
  
  return(total_pseudobulk_df_log)
}

create_metaintegrator_obj <- function(token, counts, metadata = NULL, case_time_point = NULL, control_time_point = NULL) {
  # Take appropriate subset of metadata (based on timepoints) and log-normalize counts for use in MetaIntegrator
  if(token == "bulk") {
    metadata <- metadata[metadata$time_point == case_time_point | metadata$time_point == control_time_point,]
    counts <- counts[rownames(metadata)]
    counts <- varianceStabilizingTransformation(as.matrix(counts))
  }
  # MetaIntegrator (or at least the Data Compendium) expects ENTREZ IDs so we map all our gene names to their respective ENTREZ IDs
  # Do we actually have to do this?
  entrez_gene_list <- mapIds(org.Hs.eg.db, rownames(counts), "ENTREZID", "SYMBOL")
  entrez_gene_list[is.na(entrez_gene_list)] <- names(entrez_gene_list[is.na(entrez_gene_list)]) # Replace NAs with gene names
  counts_entrez <- counts
  rownames(counts_entrez) <- unname(entrez_gene_list)
  
  # Create ENTREZ -> ENTREZ keys (including some gene names because they don't have corresponding ENTREZ IDs)
  metaintegrator_keys <- entrez_gene_list
  names(metaintegrator_keys) <- metaintegrator_keys
  
  if(token == "multiome") {
    #Remove "Sample_" from each column name
    for (col in 1:ncol(counts_entrez)){
      colnames(counts_entrez)[col] <-  sub("Sample_", "", colnames(counts_entrez)[col])
    }
    sample_names <- colnames(counts_entrez)
    sample_class <- c()
    for(sample_name in sample_names) {
      sample_class <- c(sample_class, sample_metadata[sample_metadata$aliquot == sample_name,]$time_point)
    }
    sample_class <- as.numeric(grepl("D28", sample_class))
  } else if(token == "single cell") {
    # Create class vector (1 is disease, 0 is control)
    sample_names <- colnames(counts_entrez)
    sample_class <- as.numeric(grepl("D28", sample_names))
  } else if (token == "bulk") {
    sample_names <- colnames(counts_entrez)
    sample_class <- as.character(metadata$time_point)
    sample_class[sample_class == control_time_point] <- 0
    sample_class[sample_class == case_time_point] <- 1
  }
  
  # Create pheno dataframe
  pheno_class <- sample_class
  pheno_class[pheno_class == 0] <- "Healthy"
  pheno_class[pheno_class == 1] <- "Virus"
  
  pathogen_class <- pheno_class
  pathogen_class[pathogen_class == "Virus"] <- "Influenza virus"
  
  pheno_df <- data.frame(Class = pheno_class, Pathogen = pathogen_class)
  rownames(pheno_df) <- sample_names
  
  # Create MetaIntegrator object
  metaintegrator_obj <- list()
  metaintegrator_obj$expr <- counts_entrez
  metaintegrator_obj$keys <- metaintegrator_keys
  metaintegrator_obj$formattedName <- "FLU"
  metaintegrator_obj$class <- sample_class
  names(metaintegrator_obj$class) <- colnames(metaintegrator_obj$expr)
  metaintegrator_obj$pheno <- pheno_df
  
  # Method expects a list containing datasets, so we'll just stick our obj in a list
  internal_data_list <- list()
  internal_data_list[[1]] <- metaintegrator_obj
  return(internal_data_list)
}
auc_names <- c("Filtering_Assay", "Filtering_Method", "Discovery_Assay", "Discovery_Dataset", "Pos_Genes", "Neg_Genes", "Total_Passing_Genes", "Total_Genes", "Percentage_of_Passing_Genes")

add_auc_row <- function(auc_df, auc_names, filtering_assay, filtering_method, discovery_assay, discovery_dataset, discovery_auc_df, auc_col_name) {
  current_row <- data.frame(filtering_assay, filtering_method, discovery_assay, discovery_dataset, nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] > 0.7,]),
                            nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] < 0.3,]), 
                            nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] > 0.7,]) + nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] < 0.3,]),
                            nrow(discovery_auc_df), (nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] > 0.7,]) + 
                                                              nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] < 0.3,])) / nrow(discovery_auc_df))
  names(current_row) <- auc_names
  auc_df <- rbind(auc_df, current_row)
}
