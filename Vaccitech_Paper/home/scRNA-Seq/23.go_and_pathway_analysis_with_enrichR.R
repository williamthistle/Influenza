# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Function to count genes
count_genes <- function(gene_string) {
  if (is.na(gene_string)) {
    return(0)
  } else {
    return(length(unlist(strsplit(gene_string, ";"))))
  }
}

# Function to remove the token from terms
remove_token <- function(term) {
  return(sub(" \\(GO:[0-9]+\\)$", "", term))
}

# CD14 Mono
cd14_mono_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "CD14 Mono",]
cd14_mono_upregulated_genes <- cd14_mono_genes[cd14_mono_genes$sc_log2FC > 0,]$Gene_Name
cd14_mono_downregulated_genes <- cd14_mono_genes[cd14_mono_genes$sc_log2FC < 0,]$Gene_Name
# CD16 Mono
cd16_mono_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD16_Mono_Upregulated.tsv"))
cd16_mono_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD16_Mono_Downregulated.tsv"))
# NK
nk_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "NK_Upregulated.tsv"))
nk_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "NK_Downregulated.tsv"))
# cDC
cDC_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "cDC_Upregulated.tsv"))
cDC_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "cDC_Downregulated.tsv"))
# CD4 Memory
cd4_memory_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD4_Memory_Upregulated.tsv"))
cd4_memory_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD4_Memory_Downregulated.tsv"))
# CD8 Memory
cd8_memory_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD8_Memory_Upregulated.tsv"))
cd8_memory_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD8_Memory_Downregulated.tsv"))
# CD8 Naive
cd8_naive_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD8_Naive_Upregulated.tsv"))
cd8_naive_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "CD8_Naive_Downregulated.tsv"))
# MAIT
mait_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "MAIT_Upregulated.tsv"))
mait_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "MAIT_Downregulated.tsv"))
# B Memory
b_memory_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "B_Memory_Upregulated.tsv"))
b_memory_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "B_Memory_Downregulated.tsv"))
# B Naive
b_naive_upregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "B_Naive_Upregulated.tsv"))
b_naive_downregulated_modules <- get_module_genes_for_downstream_analysis(paste0(sc_humanbase_dir, "B_Naive_Downregulated.tsv"))

# 
# dbs <- listEnrichrDbs()
go_dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023")
pathway_dbs <- c("Reactome_2022")

go_results <- enrichr(cd14_mono_upregulated_genes, go_dbs)
pathway_results <- enrichr(cd14_mono_upregulated_genes, pathway_dbs)
  
# Sleep for 1 second so we don't overload enrichR server with requests
Sys.sleep(1)


# Adjusted.P.value cutoff
cutoff <- 0.05


# Subset the dataframe based on the cutoff
subset_df <- go_results[[1]][go_results[[1]]$Adjusted.P.value <= cutoff, ]

# Count the number of genes for each row
subset_df$Num_Genes <- sapply(subset_df$Genes, count_genes)

# Remove the token from terms
subset_df$Term <- remove_token(subset_df$Term)

subset_df$Term <- factor(subset_df$Term, levels = subset_df$Term)

# Create the barplot
ggplot(subset_df, aes(x = Num_Genes, y = fct_rev(Term), fill = Adjusted.P.value)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") +  # Reversed color scale
  labs(x = "Number of Genes", y = "Term") +
  guides(fill = guide_colorbar(reverse = TRUE))  # Reverse the legend
  