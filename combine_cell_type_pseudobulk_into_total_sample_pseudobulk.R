library(DESeq2)
library(data.table)
library(org.Hs.eg.db)

base_dir <- "~/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
setup_bulk_analysis()
sample_metadata <- read.table(paste0(base_dir, "sample_metadata.tsv"), sep = "\t", header = TRUE)

pseudobulk_dirs <- c("../OneDrive - Princeton University/Influenza Analysis/Single Cell RNA-Seq/MAGICAL Analyses/Placebo 6 Sample (Run by Aliza)/scRNA/pseudo_bulk",
                     "../OneDrive - Princeton University/Influenza Analysis/True Multiome/MAGICAL Analyses/14 Placebo Sample (Final)/scRNA_pseudobulk")

for(pseudobulk_dir in pseudobulk_dirs) {
  setwd("~/")
  setwd(pseudobulk_dir)
  cell_types <- c("CD4_Memory", "CD8_Memory", "cDC", "HSPC", "pDC", "Plasmablast", "Proliferating", "NK", "T_Naive", "CD14_Mono", "CD16_Mono", "MAIT")
  
  if(grepl("Single Cell", pseudobulk_dir)) {
    # We start with B and add the rest - create sample-level pseudobulk counts
    total_pseudobulk_df <- read.table(paste0("pseudo_bulk_RNA_count_B.txt"), header = TRUE, sep = "\t")
    
    for(cell_type in cell_types) {
      print(cell_type)
      if(file.exists(paste0("pseudo_bulk_RNA_count_", cell_type, ".txt"))) {
        cell_type_pseudobulk_df <- read.table(paste0("pseudo_bulk_RNA_count_", cell_type, ".txt"), header = TRUE, sep = "\t")
        for(current_col_index in 2:ncol(cell_type_pseudobulk_df)) {
          total_pseudobulk_df[,current_col_index] <- total_pseudobulk_df[,current_col_index] + cell_type_pseudobulk_df[,current_col_index]
        }
      }
    }
  } else {
    # We start with B and add the rest - create sample-level pseudobulk counts
    total_pseudobulk_df <- read.table(paste0("pseudo_bulk_RNA_count_B.csv"), header = TRUE, sep = ",")
    
    for(cell_type in cell_types) {
      print(cell_type)
      if(file.exists(paste0("pseudo_bulk_RNA_count_", cell_type, ".csv"))) {
        cell_type_pseudobulk_df <- read.table(paste0("pseudo_bulk_RNA_count_", cell_type, ".csv"), header = TRUE, sep = ",")
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
  total_pseudobulk_df_log <- varianceStabilizingTransformation(as.matrix(total_pseudobulk_df))
  if(grepl("Single Cell", pseudobulk_dir)) {
    write.table(total_pseudobulk_df_log, file = "total_pseudobulk.txt", sep = "\t", quote = FALSE)
  } else {
    write.table(total_pseudobulk_df_log, file = "total_pseudobulk.csv", sep = ",", quote = FALSE)
  }
  
  # MetaIntegrator (or at least the Data Compendium) expects ENTREZ IDs so we map all our gene names to their respective ENTREZ IDs
  # Do we actually have to do this?
  entrez_gene_list <- mapIds(org.Hs.eg.db, rownames(total_pseudobulk_df_log), "ENTREZID", "SYMBOL")
  entrez_gene_list[is.na(entrez_gene_list)] <- names(entrez_gene_list[is.na(entrez_gene_list)]) # Replace NAs with gene names
  total_pseudobulk_df_log_entrez <- total_pseudobulk_df_log
  rownames(total_pseudobulk_df_log_entrez) <- unname(entrez_gene_list)
  
  # Create ENTREZ -> ENTREZ keys (including some gene names because they don't have corresponding ENTREZ IDs)
  metaintegrator_keys <- entrez_gene_list
  names(metaintegrator_keys) <- metaintegrator_keys
  
  if(!grepl("Single Cell", pseudobulk_dir)) {
    #Remove "Sample_" from each column name
    for (col in 1:ncol(total_pseudobulk_df_log_entrez)){
      colnames(total_pseudobulk_df_log_entrez)[col] <-  sub("Sample_", "", colnames(total_pseudobulk_df_log_entrez)[col])
    }
    sample_names <- colnames(total_pseudobulk_df_log_entrez)
    sample_class <- c()
    for(sample_name in sample_names) {
      sample_class <- c(sample_class, sample_metadata[sample_metadata$aliquot == sample_name,]$time_point)
    }
    sample_class <- as.numeric(grepl("D28", sample_class))
  } else {
    # Create class vector (1 is disease, 0 is control)
    sample_names <- colnames(total_pseudobulk_df_log_entrez)
    sample_class <- as.numeric(grepl("D28", sample_names))
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
  metaintegrator_obj$expr <- total_pseudobulk_df_log_entrez
  metaintegrator_obj$keys <- metaintegrator_keys
  metaintegrator_obj$formattedName <- "FLU"
  metaintegrator_obj$class <- sample_class
  names(metaintegrator_obj$class) <- colnames(metaintegrator_obj$expr)
  metaintegrator_obj$pheno <- pheno_df
  
  # Method expects a list containing datasets, so we'll just stick our obj in a list
  internal_data_list <- list()
  internal_data_list[[1]] <- metaintegrator_obj
  
  # Calculate individual AUCs for our genes from MAGICAL (single cell)
  # IDEA - separate into HVL and LVL and see whether it has higher AUC for higher viral load ppl
  sc_pseudobulk_aucs <- test_individual_genes_on_datasets(pseudobulk_multiome_14_genes, internal_data_list, "Single_Cell_Paired", "sc_pseudobulk")
  # We're mainly interested in ones that have AUC > 0.7 or AUC < 0.3
  nrow(sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc > 0.7,])
  sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc > 0.7,]
  sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc > 0.7,]$gene_name
  nrow(sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc < 0.3,])
  sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc < 0.3,]
  sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc < 0.3,]$gene_name
}


# Interesting idea - compare AUCs in our gene list(s) to other datasets (e.g., Data Compendium discovery or our own bulk RNA-seq)
comparison_aucs <- sc_pseudobulk_aucs
comparison_aucs$discovery_aucs <- sc_discovery_flu_aucs$flu_discovery_gene_auc
comparison_aucs <- comparison_aucs[,c(1,2,4,3)]