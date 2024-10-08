# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

find_atac_and_rna_correlation <- function(atac_df, rna_df) {
  rna_df <- rna_df[rownames(rna_df) %in% atac_df$SYMBOL,]
  atac_df <- atac_df[atac_df$SYMBOL %in% rownames(rna_df),]
  
  # Grab FC associated with each DAS symbol
  corresponding_rna_fc <- c()
  
  for(current_row_index in 1:nrow(atac_df)) {
    current_symbol <- atac_df[current_row_index,]$SYMBOL
    current_rna_fc <- rna_df[rownames(rna_df) == current_symbol,]$avg_log2FC
    corresponding_rna_fc <- c(corresponding_rna_fc, current_rna_fc)
  }
  
  comparing_first_vs_second_df <- data.frame(gene_name = atac_df$SYMBOL, first_fc = atac_df$sc_log2FC,
                                             second_fc = corresponding_rna_fc)
  return(comparing_first_vs_second_df)
}

correlation_cell_types <- c("CD14_Mono", "CD16_Mono", "cDC", "pDC", "NK")

sc_correlations <- list()
sc_correlation_plots <- list()

for(cell_type in correlation_cell_types) {
  print(cell_type)
  sc_correlations[[cell_type]] <- list()
  sc_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  # Unfiltered scRNA
  unfiltered_cell_type_sc_deg <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                            sep = "\t", header = TRUE)
  # Annotated scATAC
  upregulated_cell_type_sc_das <- read.table(paste0("C:/Users/wat2/Desktop/final_curated_list_for_me/D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final_pct_0.01_fc_0.1_upregulated_promoter_subset_annotated.tsv"),
                                               sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  downregulated_cell_type_sc_das <- read.table(paste0("C:/Users/wat2/Desktop/final_curated_list_for_me/D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final_pct_0.01_fc_0.1_downregulated_promoter_subset_annotated.tsv"),
                                 sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  cell_type_sc_das <- rbind(upregulated_cell_type_sc_das, downregulated_cell_type_sc_das)
  # Check correlation for primary significant SC genes in validation set
  comparing_first_vs_second_df <- find_atac_and_rna_correlation(cell_type_sc_das, unfiltered_cell_type_sc_deg)
  
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc,
                              method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_atac_vs_rna"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_atac_vs_rna"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("ATAC FC") + ylab("RNA FC") + labs(title = cell_type_no_underscore) + xlim(-6, 6) + ylim(-6, 6)
}

pseudobulk_corrected_plots <- lapply(sc_correlation_plots, function(x) x[[1]])
ggsave("C:/Users/wat2/Desktop/atac_vs_rna_correlations.png", plot = patchwork::wrap_plots(pseudobulk_corrected_plots, ncol = 2, nrow = 3), height = 10, width = 10)
