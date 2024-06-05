# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

find_sc_correlation_final <- function(first_gene_df, second_gene_df) {
  compare_first_df <- first_gene_df
  compare_second_df <- second_gene_df[rownames(second_gene_df) %in% compare_first_df$Peak_Name,]
  compare_first_df <- compare_first_df[compare_first_df$Peak_Name %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% compare_first_df$Peak_Name,]
  compare_first_df <- compare_first_df[order(compare_first_df$Peak_Name),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$sc_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  return(comparing_first_vs_second_df)
}

correlation_cell_types <- c("CD14_Mono", "CD16_Mono", "cDC", "pDC", "NK")

sc_correlations <- list()
sc_correlation_plots <- list()

total_overlap <- 0
total_found <- 0

for(cell_type in correlation_cell_types) {
  print(cell_type)
  sc_correlations[[cell_type]] <- list()
  sc_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  # Unfiltered SC
  unfiltered_cell_type_sc_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                             sep = "\t", header = TRUE)
  unfiltered_cell_type_validation_sc_das <- read.table(paste0(scATAC_lvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                       sep = "\t", header = TRUE)
  # Filtered SC
  cell_type_sc_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final_pct_0.1.tsv"),
                                          sep = "\t", header = TRUE)
  # Set fold change threshold (only keep high FC peaks)
  # cell_type_sc_das <- cell_type_sc_das[abs(cell_type_sc_das$sc_log2FC) > 1,]
  # cell_type_sc_das <- cell_type_sc_das[cell_type_sc_das$sc_pval < 0.01,]
  cell_type_validation_sc_das <- read.table(paste0(scATAC_hvl_vaccinated_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final_pct_0.1.tsv"),
                                                     sep = "\t", header = TRUE)
  # Check correlation for primary significant SC genes in validation set
  comparing_first_vs_second_df <- find_sc_correlation_final(cell_type_sc_das, unfiltered_cell_type_validation_sc_das)
  
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc,
                              method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_hvl_vs_lvl"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_hvl_vs_lvl"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("Naive FC") + ylab("Vaccinated FC") + labs(title = cell_type_no_underscore) + xlim(-5, 5) + ylim(-5, 5)
}

pseudobulk_corrected_plots <- lapply(sc_correlation_plots, function(x) x[[1]])
ggsave("C:/Users/willi/Desktop/hvl_vs_lvl_das_correlations.png", plot = patchwork::wrap_plots(pseudobulk_corrected_plots, ncol = 2, nrow = 3), height = 10, width = 10)
