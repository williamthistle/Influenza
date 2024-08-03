### THIS FUNCTION:
### RUNS MAGICAL AND CREATES 1) OVERALL MAGICAL DF, AND 2) TF-CENTRIC MAGICAL DF

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))
source(paste0(base_dir, "extra_functions/my_MAGICAL_helper_functions.R"))

find_sc_correlation_final <- function(sc_das, validation_das, unfiltered_sc_das, unfiltered_validation_das) {
  union_sc_pseudobulk_corrected_das <- union(rownames(sc_das), rownames(validation_das))
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
  # Unfiltered SC
  unfiltered_cell_type_sc_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                            sep = "\t", header = TRUE)
  unfiltered_cell_type_validation_sc_das <- read.table(paste0(scATAC_hvl_vaccinated_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                       sep = "\t", header = TRUE)
  # Filtered SC
  hvl_magical_peak_set <- hvl_magical_df[hvl_magical_df$Cell_Type == cell_type,]
  hvl_magical_peak_set <- paste0(hvl_magical_peak_set$Peak_chr, "-", hvl_magical_peak_set$Peak_start, "-", hvl_magical_peak_set$Peak_end)
  cell_type_sc_das <- unfiltered_cell_type_sc_das[rownames(unfiltered_cell_type_sc_das) %in% hvl_magical_peak_set,]
  
  vaccinated_magical_peak_set <- vaccinated_magical_df[vaccinated_magical_df$Cell_Type == cell_type,]
  vaccinated_magical_peak_set <- paste0(vaccinated_magical_peak_set$Peak_chr, "-", vaccinated_magical_peak_set$Peak_start, "-", vaccinated_magical_peak_set$Peak_end)
  cell_type_validation_sc_das <- unfiltered_cell_type_validation_sc_das[rownames(unfiltered_cell_type_validation_sc_das) %in% vaccinated_magical_peak_set,]
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
  sc_correlations[[cell_type]][["sc_primary_hvl_vs_vaccinated"]] <- correlation_val
  
  # How does fold change compare between HVL and vaccinated?
  print("Number of HVL > Vaccinated points")
  print(sum(abs(comparing_first_vs_second_df$second_fc) > abs(comparing_first_vs_second_df$first_fc)))
  print("Number of Vaccinated > HVL points")
  print(sum(abs(comparing_first_vs_second_df$first_fc) > abs(comparing_first_vs_second_df$second_fc)))
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_hvl_vs_vaccinated"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman", text_size = 6) + xlab("Naive vaccinated FC") + ylab("Naive HVL FC") + labs(title = cell_type_no_underscore) +
    xlim(-8, 8) + ylim(-8, 8) + theme(aspect.ratio = 1, text=element_text(size=15)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
}

hvl_vs_vaccinated_pseudobulk_corrected_plots <- lapply(sc_correlation_plots, function(x) x[[1]])
ggsave("C:/Users/willi/Desktop/hvl_vs_vaccinated_fc_magical_das_correlation.png", plot = patchwork::wrap_plots(hvl_vs_vaccinated_pseudobulk_corrected_plots, ncol = 2, nrow = 2), height = 10, width = 10)
