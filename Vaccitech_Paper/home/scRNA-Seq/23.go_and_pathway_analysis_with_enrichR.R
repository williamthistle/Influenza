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

create_enrichment_plot <- function(results, cutoff = 0.05) {
  # Subset the dataframe based on the cutoff
  subset_df <- results[results$Adjusted.P.value < cutoff, ]
  
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
}

get_enrichr_results <- function(gene_list) {
  # Set up databases
  # Use dbs <- listEnrichrDbs() to get full list of databases
  dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "Reactome_2022")
  
  # Get results
  results <- enrichr(gene_list, dbs)
  results[[1]] <- results[[1]][results[[1]]$Adjusted.P.value < 0.05,]
  results[[2]] <- results[[2]][results[[2]]$Adjusted.P.value < 0.05,]
  results[[3]] <- results[[3]][results[[3]]$Adjusted.P.value < 0.05,]
  results[[4]] <- results[[4]][results[[4]]$Adjusted.P.value < 0.05,]
  return(results)
}

# CD14 Mono
cd14_mono_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "CD14 Mono",]
cd14_mono_upregulated_genes <- cd14_mono_genes[cd14_mono_genes$sc_log2FC > 0,]$Gene_Name
cd14_mono_downregulated_genes <- cd14_mono_genes[cd14_mono_genes$sc_log2FC < 0,]$Gene_Name
# CD16 Mono
cd16_mono_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "CD16 Mono",]
cd16_mono_upregulated_genes <- cd16_mono_genes[cd16_mono_genes$sc_log2FC > 0,]$Gene_Name
cd16_mono_downregulated_genes <- cd16_mono_genes[cd16_mono_genes$sc_log2FC < 0,]$Gene_Name
# NK
nk_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "NK",]
nk_upregulated_genes <- nk_genes[nk_genes$sc_log2FC > 0,]$Gene_Name
nk_downregulated_genes <- nk_genes[nk_genes$sc_log2FC < 0,]$Gene_Name
# cDC
cDC_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "cDC",]
cDC_upregulated_genes <- cDC_genes[cDC_genes$sc_log2FC > 0,]$Gene_Name
cDC_downregulated_genes <- cDC_genes[cDC_genes$sc_log2FC < 0,]$Gene_Name
# CD4 Memory
cd4_memory_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "CD4 Memory",]
cd4_memory_upregulated_genes <- cd4_memory_genes[cd4_memory_genes$sc_log2FC > 0,]$Gene_Name
cd4_memory_downregulated_genes <- cd4_memory_genes[cd4_memory_genes$sc_log2FC < 0,]$Gene_Name
# CD8 Memory
cd8_memory_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "CD8 Memory",]
cd8_memory_upregulated_genes <- cd8_memory_genes[cd8_memory_genes$sc_log2FC > 0,]$Gene_Name
cd8_memory_downregulated_genes <- cd8_memory_genes[cd8_memory_genes$sc_log2FC < 0,]$Gene_Name
# CD8 Naive
cd8_naive_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "CD8 Naive",]
cd8_naive_upregulated_genes <- cd8_naive_genes[cd8_naive_genes$sc_log2FC > 0,]$Gene_Name
cd8_naive_downregulated_genes <- cd8_naive_genes[cd8_naive_genes$sc_log2FC < 0,]$Gene_Name
# MAIT
mait_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "MAIT",]
mait_upregulated_genes <- mait_genes[mait_genes$sc_log2FC > 0,]$Gene_Name
mait_downregulated_genes <- mait_genes[mait_genes$sc_log2FC < 0,]$Gene_Name
# B Memory
b_memory_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "B memory",]
b_memory_upregulated_genes <- b_memory_genes[b_memory_genes$sc_log2FC > 0,]$Gene_Name
b_memory_downregulated_genes <- b_memory_genes[b_memory_genes$sc_log2FC < 0,]$Gene_Name
# B Naive
b_naive_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "B naive",]
b_naive_upregulated_genes <- b_naive_genes[b_naive_genes$sc_log2FC > 0,]$Gene_Name
b_naive_downregulated_genes <- b_naive_genes[b_naive_genes$sc_log2FC < 0,]$Gene_Name

cd14_mono_enrichr_upregulated_results <- get_enrichr_results(cd14_mono_upregulated_genes)
cd14_mono_enrichr_downregulated_results <- get_enrichr_results(cd14_mono_downregulated_genes)

cd16_mono_enrichr_upregulated_results <- get_enrichr_results(cd16_mono_upregulated_genes)
cd16_mono_enrichr_downregulated_results <- get_enrichr_results(cd16_mono_downregulated_genes)

nk_enrichr_upregulated_results <- get_enrichr_results(nk_upregulated_genes)
nk_enrichr_downregulated_results <- get_enrichr_results(nk_downregulated_genes)

cDC_enrichr_upregulated_results <- get_enrichr_results(cDC_upregulated_genes)
cDC_enrichr_downregulated_results <- get_enrichr_results(cDC_downregulated_genes)

cd4_memory_enrichr_upregulated_results <- get_enrichr_results(cd4_memory_upregulated_genes)
cd4_memory_enrichr_downregulated_results <- get_enrichr_results(cd4_memory_downregulated_genes)

cd8_memory_enrichr_upregulated_results <- get_enrichr_results(cd8_memory_upregulated_genes)
cd8_memory_enrichr_downregulated_results <- get_enrichr_results(cd8_memory_downregulated_genes)

cd8_naive_enrichr_upregulated_results <- get_enrichr_results(cd8_naive_upregulated_genes)
cd8_naive_enrichr_downregulated_results <- get_enrichr_results(cd8_naive_downregulated_genes)

mait_enrichr_upregulated_results <- get_enrichr_results(mait_upregulated_genes)
mait_enrichr_downregulated_results <- get_enrichr_results(mait_downregulated_genes)

b_memory_enrichr_upregulated_results <- get_enrichr_results(b_memory_upregulated_genes)
b_memory_enrichr_downregulated_results <- get_enrichr_results(b_memory_downregulated_genes)

b_naive_enrichr_upregulated_results <- get_enrichr_results(b_naive_upregulated_genes)
b_naive_enrichr_downregulated_results <- get_enrichr_results(b_naive_downregulated_genes)





create_enrichment_plot(pathway_results[[1]])
  
# Sleep for 1 second so we don't overload enrichR server with requests
Sys.sleep(1)



  