# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

find_sc_correlation_final <- function(sc_das, validation_das, unfiltered_sc_das, unfiltered_validation_das) {
  union_sc_pseudobulk_corrected_das <- union(sc_das$Peak_Name, validation_das$Peak_Name)
  print(paste0("Number of SC DAS (pseudobulk corrected) is: ", length(union_sc_pseudobulk_corrected_das)))
  
  unfiltered_sc_das <- unfiltered_sc_das[rownames(unfiltered_sc_das) %in% union_sc_pseudobulk_corrected_das,]
  unfiltered_validation_das <- unfiltered_validation_das[rownames(unfiltered_validation_das) %in% union_sc_pseudobulk_corrected_das,]
  
  unfiltered_sc_das <- unfiltered_sc_das[rownames(unfiltered_sc_das) %in% rownames(unfiltered_validation_das),]
  unfiltered_validation_das <- unfiltered_validation_das[rownames(unfiltered_validation_das) %in% rownames(unfiltered_sc_das),]
  
  unfiltered_sc_das <- unfiltered_sc_das[order(rownames(unfiltered_sc_das)),]
  unfiltered_validation_das <- unfiltered_validation_das[order(rownames(unfiltered_validation_das)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(unfiltered_sc_das), first_fc = unfiltered_validation_das$avg_log2FC,
                                             second_fc = unfiltered_sc_das$avg_log2FC)
  return(comparing_first_vs_second_df)
}

correlation_cell_types <- c("pDC", "CD16_Mono", "cDC", "NK", "CD14_Mono")

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
  cell_type_sc_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final_pct_0.01.tsv"),
                                          sep = "\t", header = TRUE)
  cell_type_validation_sc_das <- read.table(paste0(scATAC_lvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final_pct_0.01.tsv"),
                                                     sep = "\t", header = TRUE)
  # Remove any entries that are min.pct of 0 that will mess up fold change comparisons
  unfiltered_cell_type_sc_das <- unfiltered_cell_type_sc_das[unfiltered_cell_type_sc_das$pct.1 > 0 & unfiltered_cell_type_sc_das$pct.2 > 0,]
  unfiltered_cell_type_validation_sc_das <- unfiltered_cell_type_validation_sc_das[unfiltered_cell_type_validation_sc_das$pct.1 > 0 & unfiltered_cell_type_validation_sc_das$pct.2 > 0,]
  # cell_type_sc_das <- cell_type_sc_das[cell_type_sc_das$Peak_Name %in% rownames(unfiltered_cell_type_validation_sc_das),]
  # Check correlation for primary significant SC genes in validation set
  comparing_first_vs_second_df <- find_sc_correlation_final(cell_type_sc_das, cell_type_validation_sc_das, unfiltered_cell_type_sc_das, unfiltered_cell_type_validation_sc_das)
  
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc,
                              method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_hvl_vs_lvl"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_hvl_vs_lvl"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman", text_size = 6) + xlab("Naive LVL FC") + ylab("Naive HVL FC") + labs(title = cell_type_no_underscore) +
    xlim(-8, 8) + ylim(-8, 8) + theme(aspect.ratio = 1, text=element_text(size=15)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
}

pseudobulk_corrected_plots <- lapply(sc_correlation_plots, function(x) x[[1]])
ggsave("C:/Users/willi/Desktop/hvl_vs_lvl_das_correlations.png", plot = patchwork::wrap_plots(pseudobulk_corrected_plots, ncol = 3, nrow = 2), height = 10, width = 10)
