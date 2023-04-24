library(DESeq2)
library(org.Hs.eg.db)

setwd("~/")
setwd("../OneDrive - Princeton University/Influenza Analysis/Single Cell RNA-Seq/MAGICAL Analyses/Placebo 6 Sample (Run by Aliza)/scRNA/pseudo_bulk")

high_metadata_subset <- high_placebo_metadata[high_placebo_metadata$time_point == "2_D5" | high_placebo_metadata$time_point == "2_D_minus_1",]
low_metadata_subset <- low_placebo_metadata[low_placebo_metadata$time_point == "2_D5" | low_placebo_metadata$time_point == "2_D_minus_1",]

high_counts_subset <- high_placebo_counts[rownames(high_metadata_subset)]
low_counts_subset <- low_placebo_counts[rownames(low_metadata_subset)]

high_counts_subset <- rlog(as.matrix(high_counts_subset))
low_counts_subset <- rlog(as.matrix(low_counts_subset))

entrez_gene_list <- mapIds(org.Hs.eg.db, rownames(high_counts_subset), "ENTREZID", "SYMBOL")
entrez_gene_list[is.na(entrez_gene_list)] <- names(entrez_gene_list[is.na(entrez_gene_list)]) # Replace NAs with gene names
high_counts_subset_entrez <- high_counts_subset
rownames(high_counts_subset_entrez) <- unname(entrez_gene_list)
low_counts_subset_entrez <- low_counts_subset
rownames(low_counts_subset_entrez) <- unname(entrez_gene_list)

# HIGH VIRAL LOAD AUC TESTING
# Create entrez -> entrez keys (some gene names because they don't have corresponding entrez IDs)
metaintegrator_keys <- entrez_gene_list
names(metaintegrator_keys) <- metaintegrator_keys

# Create class vector (1 is disease, 0 is control)
sample_names <- colnames(high_counts_subset_entrez)
sample_class <- as.character(high_metadata_subset$time_point)
sample_class[sample_class == "2_D_minus_1"] <- 0
sample_class[sample_class == "2_D5"] <- 1

# Create pheno dataframe
pheno_class <- sample_class
pheno_class[pheno_class == 0] <- "Healthy"
pheno_class[pheno_class == 1] <- "Virus"

pathogen_class <- pheno_class
pathogen_class[pathogen_class == "Virus"] <- "Influenza virus"

pheno_df <- data.frame(Class = pheno_class, Pathogen = pathogen_class)
rownames(pheno_df) <- sample_names

# Create metaintegrator object
bulk_high_metaintegrator_obj <- list()
bulk_high_metaintegrator_obj$expr <- high_counts_subset_entrez
bulk_high_metaintegrator_obj$keys <- metaintegrator_keys
bulk_high_metaintegrator_obj$formattedName <- "FLU_BULK_HIGH"
bulk_high_metaintegrator_obj$class <- sample_class
names(bulk_high_metaintegrator_obj$class) <- colnames(bulk_high_metaintegrator_obj$expr)
bulk_high_metaintegrator_obj$pheno <- pheno_df

bulk_high_internal_data_list <- list()
bulk_high_internal_data_list[[1]] <- bulk_high_metaintegrator_obj

sc_bulk_high_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, bulk_high_internal_data_list, "Single_Cell_Paired", "bulk_high")
nrow(sc_bulk_high_aucs[sc_bulk_high_aucs$bulk_high_gene_auc > 0.7,])
sc_bulk_high_aucs[sc_bulk_high_aucs$bulk_high_gene_auc > 0.7,]
sc_bulk_high_aucs[sc_bulk_high_aucs$bulk_high_gene_auc > 0.7,]$gene_name
nrow(sc_bulk_high_aucs[sc_bulk_high_aucs$bulk_high_gene_auc < 0.3,])
sc_bulk_high_aucs[sc_bulk_high_aucs$bulk_high_gene_auc < 0.3,]
sc_bulk_high_aucs[sc_bulk_high_aucs$bulk_high_gene_auc < 0.3,]$gene_name

# LOW VIRAL LOAD AUC TESTING
# Create entrez -> entrez keys (some gene names because they don't have corresponding entrez IDs)
metaintegrator_keys <- entrez_gene_list
names(metaintegrator_keys) <- metaintegrator_keys

# Create class vector (1 is disease, 0 is control)
sample_names <- colnames(low_counts_subset_entrez)
sample_class <- as.character(low_metadata_subset$time_point)
sample_class[sample_class == "2_D_minus_1"] <- 0
sample_class[sample_class == "2_D5"] <- 1

# Create pheno dataframe
pheno_class <- sample_class
pheno_class[pheno_class == 0] <- "Healthy"
pheno_class[pheno_class == 1] <- "Virus"

pathogen_class <- pheno_class
pathogen_class[pathogen_class == "Virus"] <- "Influenza virus"

pheno_df <- data.frame(Class = pheno_class, Pathogen = pathogen_class)
rownames(pheno_df) <- sample_names

# Create metaintegrator object
bulk_low_metaintegrator_obj <- list()
bulk_low_metaintegrator_obj$expr <- low_counts_subset_entrez
bulk_low_metaintegrator_obj$keys <- metaintegrator_keys
bulk_low_metaintegrator_obj$formattedName <- "FLU_BULK_LOW"
bulk_low_metaintegrator_obj$class <- sample_class
names(bulk_low_metaintegrator_obj$class) <- colnames(bulk_low_metaintegrator_obj$expr)
bulk_low_metaintegrator_obj$pheno <- pheno_df

bulk_low_internal_data_list <- list()
bulk_low_internal_data_list[[1]] <- bulk_low_metaintegrator_obj

sc_bulk_low_aucs <- test_individual_genes_on_datasets(MAGICAL_single_cell_genes, bulk_low_internal_data_list, "Single_Cell_Paired", "bulk_low")
nrow(sc_bulk_low_aucs[sc_bulk_low_aucs$bulk_low_gene_auc > 0.7,])
sc_bulk_low_aucs[sc_bulk_low_aucs$bulk_low_gene_auc > 0.7,]
sc_bulk_low_aucs[sc_bulk_low_aucs$bulk_low_gene_auc > 0.7,]$gene_name
nrow(sc_bulk_low_aucs[sc_bulk_low_aucs$bulk_low_gene_auc < 0.3,])
sc_bulk_low_aucs[sc_bulk_low_aucs$bulk_low_gene_auc < 0.3,]
sc_bulk_low_aucs[sc_bulk_low_aucs$bulk_low_gene_auc < 0.3,]$gene_name