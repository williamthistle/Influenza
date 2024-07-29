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
    theme_classic(base_size = 24) +
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "red", high = "blue") +  # Reversed color scale
    labs(x = "Number of Genes", y = "Term") +
    guides(fill = guide_colorbar(reverse = TRUE)) +  # Reverse the legend  
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1))))) #+
    #theme(legend.position="none")
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
# pDC
pDC_genes <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == "pDC",]
pDC_upregulated_genes <- pDC_genes[pDC_genes$sc_log2FC > 0,]$Gene_Name
pDC_downregulated_genes <- pDC_genes[pDC_genes$sc_log2FC < 0,]$Gene_Name
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

cd14_mono_enrichr_all_results <- get_enrichr_results(unique(cd14_mono_genes$Gene_Name))
cd14_mono_enrichr_upregulated_results <- get_enrichr_results(cd14_mono_upregulated_genes)
cd14_mono_enrichr_downregulated_results <- get_enrichr_results(cd14_mono_downregulated_genes)

cd16_mono_enrichr_all_results <- get_enrichr_results(unique(cd16_mono_genes$Gene_Name))
cd16_mono_enrichr_upregulated_results <- get_enrichr_results(cd16_mono_upregulated_genes)
cd16_mono_enrichr_downregulated_results <- get_enrichr_results(cd16_mono_downregulated_genes)

nk_enrichr_all_results <- get_enrichr_results(unique(nk_genes$Gene_Name))
nk_enrichr_upregulated_results <- get_enrichr_results(nk_upregulated_genes)
nk_enrichr_downregulated_results <- get_enrichr_results(nk_downregulated_genes)

cDC_enrichr_upregulated_results <- get_enrichr_results(cDC_upregulated_genes)
cDC_enrichr_downregulated_results <- get_enrichr_results(cDC_downregulated_genes)

pDC_enrichr_upregulated_results <- get_enrichr_results(pDC_upregulated_genes)
pDC_enrichr_downregulated_results <- get_enrichr_results(pDC_downregulated_genes)

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

nk_enrichr_upregulated_results <- get_enrichr_results(nk_upregulated_genes)
nk_enrichr_downregulated_results <- get_enrichr_results(nk_downregulated_genes)

cd14_mono_plotted_processes <- cd14_mono_enrichr_upregulated_results[[1]][c(5,15,17,22,26),]
cd14_mono_plotted_processes$Cell_Type <- "CD14 Mono"
cd16_mono_plotted_processes <- cd16_mono_enrichr_upregulated_results[[1]][c(1,2,3,5,6,8),]
cd16_mono_plotted_processes$Cell_Type <- "CD16 Mono"
cd4_memory_plotted_processes <- cd4_memory_enrichr_upregulated_results[[1]][c(1,2,3),]
cd4_memory_plotted_processes$Cell_Type <- "CD4 Memory"
cd8_memory_plotted_processes <- cd8_memory_enrichr_upregulated_results[[1]][c(1,2,4),]
cd8_memory_plotted_processes$Cell_Type <- "CD8 Memory"

all_plotted_processes <- rbind(cd14_mono_plotted_processes, cd16_mono_plotted_processes, cd4_memory_plotted_processes,
                               cd8_memory_plotted_processes)
all_plotted_processes$Gene.Count <- sapply(all_plotted_processes$Genes, count_genes)
all_plotted_processes <- all_plotted_processes[,c(1,4,10,11)]
all_plotted_processes$Term <- remove_token(all_plotted_processes$Term)
colnames(all_plotted_processes) <- c("Biological.Process", "Adjusted.P.Value", "Cell.Type", "Gene.Count")

all_plotted_processes$Cell.Type <- factor(all_plotted_processes$Cell.Type, levels = c("CD14 Mono", "CD16 Mono", "CD4 Memory", "CD8 Memory"))
all_plotted_processes$Biological.Process <- factor(all_plotted_processes$Biological.Process, levels = unique(all_plotted_processes$Biological.Process))

up_pathway_plot <- ggplot(data = all_plotted_processes, aes(x = Cell.Type, y = Biological.Process, size = Gene.Count, color = Adjusted.P.Value)) +
  geom_point() +
  theme_minimal() + guides(fill = guide_colorbar(reverse = TRUE)) +
  labs(
    title = "Biological Processes Associated with Upregulated Genes from Immune Cell Types",
    x = "Cell Type",
    y = "Biological Process",
    size = "Gene Count",
    color = "Adjusted P Value"
  ) +
  theme(plot.title = element_text(hjust = 1))

ggsave(filename = paste0("C:/Users/willi/Desktop/", "HB_monocyte_upregulated_genes_monocyte_gsea_top_upregulated_pathways.png"), plot = up_pathway_plot, device='png', dpi=300, width = 8)


cd14_test <- cd14_mono_enrichr_upregulated_results[[1]]
cd14_test$Gene_Count <- sapply(cd14_test$Genes, count_genes)
cd14_test <- cd14_test[cd14_test$Gene_Count >= 5,]
create_enrichment_plot(cd14_test)

cd16_test <- cd16_mono_enrichr_upregulated_results[[1]]
cd16_test$Gene_Count <- sapply(cd16_test$Genes, count_genes)
cd16_test <- cd16_test[cd16_test$Gene_Count >= 4,]
create_enrichment_plot(cd16_test)


cd16_test <- cd16_mono_enrichr_downregulated_results[[1]]
cd16_test$Gene_Count <- sapply(cd16_test$Genes, count_genes)
cd16_test <- cd16_test[cd16_test$Gene_Count >= 5,]
create_enrichment_plot(cd16_test)


cd14_mono_enrichr_upregulated_results_go <- cd14_mono_enrichr_upregulated_results[[1]]
cd14_mono_enrichr_upregulated_results_go$Pathway_Type <- "GO"
cd14_mono_enrichr_upregulated_results_go
cd14_mono_enrichr_upregulated_results_reactome <- cd14_mono_enrichr_upregulated_results[[4]]
cd14_mono_enrichr_upregulated_results_reactome$Pathway_Type <- "Reactome"




cd14_mono_enrichr_upregulated_results_for_barplot <- cd14_mono_enrichr_upregulated_results




cd14_mono_enrichr_upregulated_results_for_barplot <- cd14_mono_enrichr_upregulated_results_for_barplot[c(1,5,10,15,17,22,26),]


create_enrichment_plot(cd14_mono_enrichr_upregulated_results[[1]])
 
create_enrichment_plot(cd16_mono_enrichr_upregulated_results[[1]])

cd14_reactome_plot <- create_enrichment_plot(cd14_mono_enrichr_upregulated_results[[4]])
ggsave("C:/Users/willi/Desktop/CD14_reactome_plot.png", plot = cd14_reactome_plot, device = "png", width = 10, height = 5, units = "in")

cd14_go_plot <- create_enrichment_plot(cd14_mono_enrichr_upregulated_results[[1]])
ggsave("C:/Users/wat2/Desktop/CD14_go_plot_with_legend.png", plot = cd14_go_plot, device = "png", width = 18, height = 9, units = "in")

cd16_go_plot <- create_enrichment_plot(NK_enrichr_upregulated_results[[4]])
ggsave("C:/Users/wat2/Desktop/CD16_go_plot.png", plot = cd16_go_plot, device = "png", width = 14, height = 9, units = "in")

# Sleep for 1 second so we don't overload enrichR server with requests
Sys.sleep(1)



  