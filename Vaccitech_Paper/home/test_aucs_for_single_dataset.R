current_dir <- "C:/Users/wat2/OneDrive - Princeton University/Influenza Analysis/Single Cell RNA-Seq/MAGICAL Analyses/Placebo 4 Subject HVL (SPEEDI) - SCT/"

current_sc_pseudobulk_gene_table <- read.table(paste0(current_dir, "D28_D1_DESeq2_pseudobulk_genes.tsv"), sep = "\t", header = TRUE)

#cd14_monocytes_gene_table <- read.table(paste0(current_dir, "scRNA_DEGs/D28-vs-D_minus_1-degs-CD14_Mono-time_point.csv"), sep = ",", header = TRUE)

#current_sc_magical_gene_table <- read.table(paste0(current_dir, "D28_D1_MAGICAL.txt"), sep = "\t", header = TRUE)

#current_strict_sc_pseudobulk_gene_table <- read.table(paste0(current_dir, "D28_D1_MAGICAL_pseudobulk_genes_strict.txt"), sep = "\t", header = TRUE)
#current_strict_sc_magical_gene_table <- read.table(paste0(current_dir, "D28_D1_MAGICAL_strict.txt"), sep = "\t", header = TRUE)

current_sc_pseudobulk_genes <- unique(current_sc_pseudobulk_gene_table$Gene_Name)
print(paste0("Number of genes that pass pseudobulk (scRNA): ", length(current_sc_pseudobulk_genes)))
#current_sc_magical_genes <- unique(current_sc_magical_gene_table$Gene_symbol)
#print(paste0("Number of genes that pass MAGICAL (scRNA): ", length(current_sc_magical_genes)))

#current_strict_sc_pseudobulk_genes <- unique(current_strict_sc_pseudobulk_gene_table$gene)
#print(paste0("Number of genes that pass pseudobulk (scRNA, stricter): ", length(current_strict_sc_pseudobulk_genes)))
#current_strict_sc_magical_genes <- unique(current_strict_sc_magical_gene_table$Gene_symbol)
#print(paste0("Number of genes that pass MAGICAL (scRNA, stricter): ", length(current_strict_sc_magical_genes)))

temp_auc_df <- data.frame(Filtering_Assay = character(), Filtering_Method = character(), Discovery_Assay = character(), 
                     Discovery_Dataset = character(), Pos_Genes = integer(), Neg_Genes = integer(), Total_Passing_Genes = integer(), 
                     Total_Genes = integer(), Percentage_of_Passing_Genes = double(), stringsAsFactors = FALSE)
temp_auc_names <- c("Filtering_Assay", "Filtering_Method", "Discovery_Assay", "Discovery_Dataset", "Passing_Pos_Genes", "Passing_Neg_Genes", "Total_Passing_Genes", "Total_Genes", "Percentage_of_Passing_Genes")

# Calculate individual AUCs for our gene lists on their respective pseudobulk data
# Why is it 439 instead of 452 (length of pseudobulk genes)?
current_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_sc_pseudobulk_genes, sc_pseudobulk_metaintegrator_obj, "sc_paired", "sc_pseudobulk"))
temp_auc_df <- add_auc_row(temp_auc_df, temp_auc_names, "Single Cell", "Cell Type Pseudobulk", "Single Cell", "Total Pseudobulk", current_sc_pseudobulk_gene_aucs, "sc_pseudobulk_gene_auc")

#current_strict_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_strict_sc_pseudobulk_genes, sc_pseudobulk_metaintegrator_obj, "sc_paired", "sc_pseudobulk"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk (Strict)", "Single Cell", "Total Pseudobulk", current_strict_sc_pseudobulk_gene_aucs, "sc_pseudobulk_gene_auc")

#current_sc_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_sc_magical_genes, sc_pseudobulk_metaintegrator_obj, "sc_paired", "sc_magical"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "MAGICAL", "Single Cell", "Total Pseudobulk", current_sc_magical_gene_aucs, "sc_magical_gene_auc")

#current_strict_sc_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_strict_sc_magical_genes, sc_pseudobulk_metaintegrator_obj, "sc_paired", "sc_magical"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "MAGICAL (Strict)", "Single Cell", "Total Pseudobulk", current_strict_sc_magical_gene_aucs, "sc_magical_gene_auc")

### CALCULATE AUCS ON BULK (SAME SUBJECTS AS SINGLE CELL) - basically an alternative for pseudobulk
current_sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_sc_pseudobulk_genes, sc_D28_bulk_metaintegrator_obj, "sc_paired", "sc_bulk_D28"))
temp_auc_df <- add_auc_row(temp_auc_df, temp_auc_names, "Single Cell", "Cell Type Pseudobulk", "Single Cell", "D28 Bulk for Single Cell Subjects", current_sc_bulk_D28_sc_pseudobulk_gene_aucs, "sc_bulk_D28_gene_auc")

#current_strict_sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_strict_sc_pseudobulk_genes, sc_D28_bulk_metaintegrator_obj, "sc_paired", "sc_bulk_D28"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk (Strict)", "Single Cell", "D28 Bulk for Single Cell Subjects", current_strict_sc_bulk_D28_sc_pseudobulk_gene_aucs, "sc_bulk_D28_gene_auc")

#current_sc_bulk_D28_sc_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_sc_magical_genes, sc_D28_bulk_metaintegrator_obj, "sc_paired", "sc_bulk_D28_magical"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "MAGICAL", "Single Cell", "D28 Bulk for Single Cell Subjects", current_sc_bulk_D28_sc_magical_gene_aucs, "sc_bulk_D28_magical_gene_auc")

#current_strict_sc_bulk_D28_sc_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_strict_sc_magical_genes, sc_D28_bulk_metaintegrator_obj, "sc_paired", "sc_bulk_D28_magical"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "MAGICAL (Strict)", "Single Cell", "D28 Bulk for Single Cell Subjects", current_strict_sc_bulk_D28_sc_magical_gene_aucs, "sc_bulk_D28_magical_gene_auc")

### CALCULATE AUCS ON D28 HVL BULK
current_high_sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_sc_pseudobulk_genes, high_D28_bulk_metaintegrator_obj, "sc_paired", "high_bulk_D28"))
temp_auc_df <- add_auc_row(temp_auc_df, temp_auc_names, "Single Cell", "Cell Type Pseudobulk", "Single Cell", "D28 Bulk for HVL Subjects", current_high_sc_bulk_D28_sc_pseudobulk_gene_aucs, "high_bulk_D28_gene_auc")

#current_high_strict_sc_bulk_D28_sc_pseudobulk_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_strict_sc_pseudobulk_genes, high_D28_bulk_metaintegrator_obj, "sc_paired", "high_bulk_D28"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "Cell Type Pseudobulk (Strict)", "Single Cell", "D28 Bulk for HVL Subjects", current_high_strict_sc_bulk_D28_sc_pseudobulk_gene_aucs, "high_bulk_D28_gene_auc")

#current_high_sc_bulk_D28_sc_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_sc_magical_genes, high_D28_bulk_metaintegrator_obj, "sc_paired", "high_bulk_D28_magical"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "MAGICAL", "Single Cell", "D28 Bulk for HVL Subjects", current_high_sc_bulk_D28_sc_magical_gene_aucs, "high_bulk_D28_magical_gene_auc")

#current_high_strict_sc_bulk_D28_sc_magical_gene_aucs <- na.omit(test_individual_genes_on_datasets(current_strict_sc_magical_genes, high_D28_bulk_metaintegrator_obj, "sc_paired", "high_bulk_D28_magical"))
#temp_auc_df <- add_auc_row(temp_auc_df, auc_names, "Single Cell", "MAGICAL (Strict)", "Single Cell", "D28 Bulk for HVL Subjects", current_high_strict_sc_bulk_D28_sc_magical_gene_aucs, "high_bulk_D28_magical_gene_auc")
