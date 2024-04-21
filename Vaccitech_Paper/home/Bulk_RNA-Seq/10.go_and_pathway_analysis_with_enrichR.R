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
  term <- sub(" \\(GO:[0-9]+\\)$", "", term)
  term <- gsub(" R-HSA-\\d+", "", term)
  return(term)
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
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "red", high = "blue") +  # Reversed color scale
    labs(x = "Number of Genes", y = "Term") +
    guides(fill = guide_colorbar(reverse = TRUE)) +  # Reverse the legend  
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1)))))
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

count_genes <- function(x) {
  matches <- gregexpr(";", x)
  count <- sapply(matches, function(match) sum(match != -1))
  return(count + 1)
}

# HVL Naive D5 
hvl_naive_d5 <- hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results[[5]]
hvl_naive_d5_upregulated_genes <- rownames(hvl_naive_d5[hvl_naive_d5$log2FoldChange > 0,])
hvl_naive_d5_downregulated_genes <- rownames(hvl_naive_d5[hvl_naive_d5$log2FoldChange < 0,])

bulk_D5_enrichr_upregulated_results <- get_enrichr_results(hvl_naive_d5_upregulated_genes)
create_enrichment_plot(bulk_D5_enrichr_upregulated_results[[4]])
bulk_D5_enrichr_downregulated_results <- get_enrichr_results(hvl_naive_d5_downregulated_genes)

# HVL Naive D8
hvl_naive_d8 <- hvl_placebo_period_2_D8_vs_D_minus_1_results[[5]]
hvl_naive_d8_upregulated_genes <- rownames(hvl_naive_d8[hvl_naive_d8$log2FoldChange > 0,])
hvl_naive_d8_downregulated_genes <- rownames(hvl_naive_d8[hvl_naive_d8$log2FoldChange < 0,])

bulk_D8_enrichr_upregulated_results <- get_enrichr_results(hvl_naive_d8_upregulated_genes)
create_enrichment_plot(bulk_D8_enrichr_upregulated_results[[4]][1:20,])


cd14_mono_enrichr_upregulated_results <- get_enrichr_results(cd14_mono_upregulated_genes)
cd14_mono_enrichr_downregulated_results <- get_enrichr_results(cd14_mono_downregulated_genes)

cd16_mono_enrichr_upregulated_results <- get_enrichr_results(cd16_mono_upregulated_genes)
cd16_mono_enrichr_downregulated_results <- get_enrichr_results(cd16_mono_downregulated_genes)

nk_enrichr_upregulated_results <- get_enrichr_results(nk_upregulated_genes)
nk_enrichr_downregulated_results <- get_enrichr_results(nk_downregulated_genes)

cDC_enrichr_upregulated_results <- get_enrichr_results(cDC_upregulated_genes)
cDC_enrichr_downregulated_results <- get_enrichr_results(cDC_downregulated_genes)

cd4_memory_enrichr_upregulated_results <- get_enrichr_results(cd4_memory_upregulated_genes)
# Not uninteresting pathways, but maybe not enough space
cd4_memory_enrichr_downregulated_results <- get_enrichr_results(cd4_memory_downregulated_genes)

cd8_memory_enrichr_upregulated_results <- get_enrichr_results(cd8_memory_upregulated_genes)
# Not uninteresting pathways, but maybe not enough space - maybe discuss downregulated pathways in both CD4 and CD8? Shared?
cd8_memory_enrichr_downregulated_results <- get_enrichr_results(cd8_memory_downregulated_genes)

cd8_naive_enrichr_upregulated_results <- get_enrichr_results(cd8_naive_upregulated_genes)
cd8_naive_enrichr_downregulated_results <- get_enrichr_results(cd8_naive_downregulated_genes)

mait_enrichr_upregulated_results <- get_enrichr_results(mait_upregulated_genes)
mait_enrichr_downregulated_results <- get_enrichr_results(mait_downregulated_genes)

b_memory_enrichr_upregulated_results <- get_enrichr_results(b_memory_upregulated_genes)
b_memory_enrichr_downregulated_results <- get_enrichr_results(b_memory_downregulated_genes)

b_naive_enrichr_upregulated_results <- get_enrichr_results(b_naive_upregulated_genes)
b_naive_enrichr_downregulated_results <- get_enrichr_results(b_naive_downregulated_genes)

cd14_test <- cd14_mono_enrichr_upregulated_results[[1]]
cd14_test$Gene_Count <- sapply(cd14_test$Genes, count_genes)
cd14_test <- cd14_test[cd14_test$Gene_Count >= 5,]
create_enrichment_plot(cd14_test)

cd16_test <- cd16_mono_enrichr_downregulated_results[[1]]
cd16_test$Gene_Count <- sapply(cd16_test$Genes, count_genes)
cd16_test <- cd16_test[cd16_test$Gene_Count >= 5,]
create_enrichment_plot(cd16_test)

create_enrichment_plot(cd14_mono_enrichr_upregulated_results[[1]])
 
create_enrichment_plot(cd16_mono_enrichr_upregulated_results[[1]])

create_enrichment_plot(cd14_mono_enrichr_upregulated_results[[4]])

# Sleep for 1 second so we don't overload enrichR server with requests
Sys.sleep(1)



  