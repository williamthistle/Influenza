### THIS FUNCTION:
### RUNS MAGICAL AND CREATES 1) OVERALL MAGICAL DF, AND 2) TF-CENTRIC MAGICAL DF

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))
source(paste0(base_dir, "extra_functions/my_MAGICAL_helper_functions.R"))

# Grab starting DFs
hvl_magical_df <- read.table(paste0(MAGICAL_hvl_placebo_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"), sep = "\t", header = TRUE)
vaccinated_magical_df <- read.table(paste0(MAGICAL_hvl_vaccinated_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"), sep = "\t", header = TRUE)

correlation_cell_types <- c("CD16_Mono", "cDC", "NK", "CD14_Mono")

sc_correlations <- list()
sc_correlation_plots <- list()

for(cell_type in correlation_cell_types) {
  print(cell_type)
  sc_correlations[[cell_type]] <- list()
  sc_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  if(cell_type_no_underscore == "NK CD56bright") {
    cell_type_no_underscore <- "NK_CD56bright"
  }
  # Unfiltered SC
  unfiltered_cell_type_sc_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                             sep = "\t", header = TRUE)
  unfiltered_cell_type_vaccinated_sc_degs <- read.table(paste0(scRNA_hvl_vaccinated_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                 sep = "\t", header = TRUE)
  # Filtered SC
  cell_type_hvl_magical_df <- hvl_magical_df[hvl_magical_df$Cell_Type == cell_type,]
  cell_type_hvl_magical_genes <- unique(cell_type_hvl_magical_df$Gene_symbol)
  
  cell_type_vaccinated_magical_df <- vaccinated_magical_df[vaccinated_magical_df$Cell_Type == cell_type,]
  cell_type_vaccinated_magical_genes <- unique(cell_type_vaccinated_magical_df$Gene_symbol)
  
  # Check correlation of union of HVL and vaccinated pseudobulk corrected genes
  union_sc_pseudobulk_corrected_degs <- union(cell_type_hvl_magical_genes, cell_type_vaccinated_magical_genes)
  print(paste0("Number of SC DEGs (pseudobulk corrected) in primary data is: ", length(union_sc_pseudobulk_corrected_degs)))
  
  unfiltered_cell_type_sc_degs <- unfiltered_cell_type_sc_degs[rownames(unfiltered_cell_type_sc_degs) %in% union_sc_pseudobulk_corrected_degs,]
  unfiltered_cell_type_vaccinated_sc_degs <- unfiltered_cell_type_vaccinated_sc_degs[rownames(unfiltered_cell_type_vaccinated_sc_degs) %in% union_sc_pseudobulk_corrected_degs,]
  
  unfiltered_cell_type_sc_degs <- unfiltered_cell_type_sc_degs[order(rownames(unfiltered_cell_type_sc_degs)),]
  unfiltered_cell_type_vaccinated_sc_degs <- unfiltered_cell_type_vaccinated_sc_degs[order(rownames(unfiltered_cell_type_vaccinated_sc_degs)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(unfiltered_cell_type_sc_degs), first_fc = unfiltered_cell_type_vaccinated_sc_degs$avg_log2FC,
                                             second_fc = unfiltered_cell_type_sc_degs$avg_log2FC)
  
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_pseudobulk_corrected_vs_sc_vaccinated_union"]] <- correlation_val
  
  # Plot correlation
  # Idea! Use data to set limits. 
  sc_correlation_plots[[cell_type]][["sc_primary_pseudobulk_corrected_vs_sc_vaccinated_union"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman", text_size = 6) + xlab("Vaccinated HVL FC") + ylab("Naive HVL FC") + labs(title = cell_type_no_underscore) + 
    xlim(-3.5, 3.5) + ylim(-3.5, 3.5) + theme(aspect.ratio = 1, text=element_text(size=15)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
}

hvl_vs_vaccinated_pseudobulk_corrected_plots <- lapply(sc_correlation_plots, function(x) x[[1]])
ggsave("C:/Users/willi/Desktop/hvl_vs_vaccinated_fc_correlation.png", plot = patchwork::wrap_plots(hvl_vs_vaccinated_pseudobulk_corrected_plots, ncol = 3, nrow = 2), height = 10, width = 10)
