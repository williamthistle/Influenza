library(DESeq2)
library(org.Hs.eg.db)

grab_transformed_pseudobulk_counts <- function(pseudobulk_dir, cell_types) {
  total_pseudobulk_df <- NULL
  if(grepl("Aliza", pseudobulk_dir)) {
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
  
  
  if(token == "mine") {
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
  } else if(token == "Aliza") {
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

add_auc_row <- function(auc_df, auc_names, filtering_assay, filtering_method, discovery_assay, discovery_dataset, discovery_auc_df, auc_col_name) {
  current_row <- data.frame(filtering_assay, filtering_method, discovery_assay, discovery_dataset, nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] > 0.7,]),
                            nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] < 0.3,]), 
                            nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] > 0.7,]) + nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] < 0.3,]),
                            nrow(discovery_auc_df), (nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] > 0.7,]) + 
                                                              nrow(discovery_auc_df[discovery_auc_df[[auc_col_name]] < 0.3,])) / nrow(discovery_auc_df))
  names(current_row) <- auc_names
  auc_df <- rbind(auc_df, current_row)
}

find_aucs_of_interest <- function(gene_table, metaintegrator_obj, source) {
  gene_list <- NULL
  if("Gene_symbol" %in% colnames(gene_table)) {
    colnames(gene_table)[1] <- "Cell_Type"
    colnames(gene_table)[2] <- "Gene_Name"
    colnames(gene_table)[7] <- "sc_log2FC"
  }
  gene_list <- unique(gene_table$Gene_Name)
  all_aucs <- na.omit(test_individual_genes_on_datasets(gene_list, metaintegrator_obj, source))
  cell_type_df <- data.frame(Gene_Name = character(), Cell_Type = character(), stringsAsFactors = FALSE)
  gene_table_subset <- gene_table[gene_table$Gene_Name %in% all_aucs$gene_name,]
  for(gene_name in unique(gene_table_subset$Gene_Name)) {
    gene_subset <- gene_table_subset[gene_table_subset$Gene_Name %in% gene_name,]
    cell_types <- paste(gene_subset$Cell_Type, collapse = ",")
    current_row <- data.frame(unique(gene_subset$Gene_Name), cell_types)
    names(current_row) <- c("Gene_Name", "Cell_Type")
    cell_type_df <- rbind(cell_type_df, current_row)
  }
  all_aucs$gene_name <- all_aucs$gene_name[order(match(all_aucs$gene_name,cell_type_df$Gene_Name))]
  all_aucs$cell_type <- cell_type_df$Cell_Type
  # Only use genes that have AUC that matches fold change direction (pos FC = pos AUC, neg FC = neg AUC)
  # for further analysis.
  curated_gene_list <- c()
  for(current_gene in all_aucs$gene_name) {
    associated_fc <- gene_table[gene_table$Gene_Name == current_gene,]$sc_log2FC
    associated_auc <- all_aucs[all_aucs$gene_name == current_gene,]$gene_auc
    if((all(associated_fc < 0) & all(associated_auc < 0.5)) | (all(associated_fc > 0) & all(associated_auc > 0.5))) {
      curated_gene_list <- c(curated_gene_list, current_gene)
    }
  }
  curated_all_aucs <- all_aucs[all_aucs$gene_name %in% curated_gene_list,]
  curated_gene_list <- curated_all_aucs[curated_all_aucs$gene_auc < 0.3 | curated_all_aucs$gene_auc > 0.7,]$gene_name
  high_genes <- curated_all_aucs[curated_all_aucs$gene_auc > 0.7,]$gene_name
  low_genes <- curated_all_aucs[curated_all_aucs$gene_auc < 0.3,]$gene_name
  return(list(all_aucs, curated_gene_list, high_genes, low_genes))
}