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
bulk_D5_enrichr_downregulated_results <- get_enrichr_results(hvl_naive_d5_downregulated_genes)

# HVL Naive D8
hvl_naive_d8 <- hvl_placebo_period_2_D8_vs_D_minus_1_results[[5]]
hvl_naive_d8_upregulated_genes <- rownames(hvl_naive_d8[hvl_naive_d8$log2FoldChange > 0,])
hvl_naive_d8_downregulated_genes <- rownames(hvl_naive_d8[hvl_naive_d8$log2FoldChange < 0,])

bulk_D8_enrichr_upregulated_results <- get_enrichr_results(hvl_naive_d8_upregulated_genes)

# Plot upregulated pathways for D5 / D8

bulk_D5_plotted_processes <- bulk_D5_enrichr_upregulated_results[[4]]
overlapping_process_indices <- which(bulk_D8_enrichr_upregulated_results[[4]]$Term %in% bulk_D5_plotted_processes$Term)
overlapping_process_indices <- c(overlapping_process_indices, 1, 3, 5, 12)
bulk_D8_plotted_processes <- bulk_D8_enrichr_upregulated_results[[4]][overlapping_process_indices,]  

bulk_D5_plotted_processes$Day <- "Day 5"
bulk_D8_plotted_processes$Day <- "Day 8"

all_plotted_processes <- rbind(bulk_D5_plotted_processes, bulk_D8_plotted_processes)
all_plotted_processes$Gene.Count <- sapply(all_plotted_processes$Genes, count_genes)
all_plotted_processes <- all_plotted_processes[,c(1,4,10,11)]
all_plotted_processes$Term <- remove_token(all_plotted_processes$Term)
colnames(all_plotted_processes) <- c("Reactome.Pathway", "Adjusted.P.Value", "Day", "Gene.Count")

all_plotted_processes$Day <- factor(all_plotted_processes$Day, levels = c("Day 5", "Day 8"))
all_plotted_processes$Reactome.Pathway <- factor(all_plotted_processes$Reactome.Pathway, levels = unique(all_plotted_processes$Reactome.Pathway))

up_pathway_plot <- ggplot(data = all_plotted_processes, aes(x = Day, y = Reactome.Pathway, size = Gene.Count, color = Adjusted.P.Value)) +
  geom_point() +
  theme_minimal() + guides(fill = guide_colorbar(reverse = TRUE)) +
  labs(
    title = "Pathway Analaysis of Acute Phase of Influenza",
    x = "Day",
    y = "Reactome Pathway",
    size = "Gene Count",
    color = "Adjusted P Value"
  ) +
  theme(plot.title = element_text(hjust = 1))

ggsave(filename = paste0(monocyte_dir, "HB_monocyte_upregulated_genes_monocyte_gsea_top_upregulated_pathways.tiff"), plot = up_pathway_plot_file, device='tiff', dpi=300)
