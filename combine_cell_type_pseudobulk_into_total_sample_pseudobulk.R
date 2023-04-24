library(DESeq2)
library(org.Hs.eg.db)

setwd("~/")
setwd("../OneDrive - Princeton University/Influenza Analysis/Single Cell RNA-Seq/MAGICAL Analyses/Placebo 6 Sample (Run by Aliza)/scRNA/pseudo_bulk")
cell_types <- c("CD4_Memory", "CD8_Memory", "cDC", "HSPC", "pDC", "Plasmablast", "Proliferating", "NK", "T_Naive", "CD14_Mono", "CD16_Mono", "MAIT")

# We start with B and add the rest - create sample-level pseudobulk counts
total_pseudobulk_df <- read.table(paste0("pseudo_bulk_RNA_count_B.txt"), header = TRUE, sep = "\t")

for(cell_type in cell_types) {
  print(cell_type)
  cell_type_pseudobulk_df <- read.table(paste0("pseudo_bulk_RNA_count_", cell_type, ".txt"), header = TRUE, sep = "\t")
  for(current_col_index in 2:ncol(cell_type_pseudobulk_df)) {
    total_pseudobulk_df[,current_col_index] <- total_pseudobulk_df[,current_col_index] + cell_type_pseudobulk_df[,current_col_index]
  }
}

# Grab gene names and set them to row names
rownames(total_pseudobulk_df) <- total_pseudobulk_df$X
total_pseudobulk_df <- total_pseudobulk_df[,-1]

# Use DESeq2 rlog function to log-normalize gene counts (as requested by MetaIntegrator)
total_pseudobulk_df_log <- rlog(as.matrix(total_pseudobulk_df))
write.table(total_pseudobulk_df_log, file = "total_pseudobulk.tsv", sep = "\t", quote = FALSE)

# MetaIntegrator (or at least the Data Compendium) expects ENTREZ IDs so we map all our gene names to their respective ENTREZ IDs
# Do we actually have to do this?
entrez_gene_list <- mapIds(org.Hs.eg.db, rownames(total_pseudobulk_df_log), "ENTREZID", "SYMBOL")
entrez_gene_list[is.na(entrez_gene_list)] <- names(entrez_gene_list[is.na(entrez_gene_list)]) # Replace NAs with gene names
total_pseudobulk_df_log_entrez <- total_pseudobulk_df_log
rownames(total_pseudobulk_df_log_entrez) <- unname(entrez_gene_list)

# Create ENTREZ -> ENTREZ keys (including some gene names because they don't have corresponding ENTREZ IDs)
metaintegrator_keys <- entrez_gene_list
names(metaintegrator_keys) <- metaintegrator_keys

# Create class vector (1 is disease, 0 is control)
sample_names <- colnames(total_pseudobulk_df_log_entrez)
sample_class <- as.numeric(grepl("D28", sample_names))

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
metaintegrator_obj$formattedName <- "FLU_SC"
metaintegrator_obj$class <- sample_class
names(metaintegrator_obj$class) <- colnames(metaintegrator_obj$expr)
metaintegrator_obj$pheno <- pheno_df

# Method expects a list containing datasets, so we'll just stick our obj in a list
internal_data_list <- list()
internal_data_list[[1]] <- metaintegrator_obj

# Calculate individual AUCs for our genes from MAGICAL (single cell)
sc_pseudobulk_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, internal_data_list, "Single_Cell_Paired", "sc_pseudobulk")
# We're mainly interested in ones that have AUC > 0.7 or AUC < 0.3
nrow(sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc > 0.7,])
sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc > 0.7,]
sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc > 0.7,]$gene_name
nrow(sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc < 0.3,])
sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc < 0.3,]
sc_pseudobulk_aucs[sc_pseudobulk_aucs$sc_pseudobulk_gene_auc < 0.3,]$gene_name

# Interesting idea - compare AUCs in our gene list(s) to other datasets (e.g., Data Compendium discovery or our own bulk RNA-seq)
comparison_aucs <- sc_pseudobulk_aucs
comparison_aucs$discovery_aucs <- sc_discovery_flu_aucs$flu_discovery_gene_auc
comparison_aucs <- comparison_aucs[,c(1,2,4,3)]