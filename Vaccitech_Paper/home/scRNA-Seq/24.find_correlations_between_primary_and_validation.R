# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

correlation_cell_types <- c("CD14_Mono", "CD16_Mono", "cDC", "pDC", "NK", "NK_CD56bright")

sc_correlations <- list()
sc_correlation_plots <- list()

for(cell_type in correlation_cell_types) {
  print(cell_type)
  sc_correlations[[cell_type]] <- list()
  sc_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  if(cell_type_no_underscore == "B memory") {
    cell_type_no_underscore <- "B Memory"
  } else if(cell_type_no_underscore == "B naive") {
    cell_type_no_underscore <- "B Naive"
  }
  # Unfiltered SC
  unfiltered_cell_type_sc_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                             sep = "\t", header = TRUE)
  unfiltered_cell_type_validation_sc_degs <- read.table(paste0(scRNA_hvl_vaccinated_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                        sep = "\t", header = TRUE)
  # Filtered SC
  cell_type_sc_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc.tsv"),
                                          sep = "\t", header = TRUE)
  cell_type_validation_sc_degs <- read.table(paste0(scRNA_hvl_vaccinated_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc.tsv"),
                                                     sep = "\t", header = TRUE)
  # Pseudobulk corrected
  cell_type_sc_pseudobulk_corrected_degs <- read.table(paste0(scRNA_hvl_placebo_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final.tsv"),
                                  sep = "\t", header = TRUE)
  rownames(cell_type_sc_pseudobulk_corrected_degs) <- cell_type_sc_pseudobulk_corrected_degs$Gene_Name
  cell_type_validation_sc_pseudobulk_corrected_degs <- read.table(paste0(scRNA_hvl_vaccinated_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final.tsv"),
                                             sep = "\t", header = TRUE)
  rownames(cell_type_validation_sc_pseudobulk_corrected_degs) <- cell_type_validation_sc_pseudobulk_corrected_degs$Gene_Name
  
  # Check correlation for primary significant SC genes in validation set
  primary_sc_degs <- rownames(cell_type_sc_degs)
  print(paste0("Number of SC DEGs in primary data is: ", length(primary_sc_degs)))
  compare_first_df <- cell_type_sc_degs
  compare_second_df <- unfiltered_cell_type_validation_sc_degs[rownames(unfiltered_cell_type_validation_sc_degs) %in% primary_sc_degs,]
  compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$avg_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  #print(correlation_val$estimate)
  #print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_vs_sc_validation"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_vs_sc_validation"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr() + xlab("Naive FC") + ylab("Vaccinated FC") + labs(title = cell_type_no_underscore)
  
  # Check correlation for primary significant SC genes (pseudobulk_corrected) in validation set
  primary_sc_pseudobulk_corrected_degs <- rownames(cell_type_sc_pseudobulk_corrected_degs)
  print(paste0("Number of SC DEGs (pseudobulk corrected) in primary data is: ", length(primary_sc_pseudobulk_corrected_degs)))
  compare_first_df <- cell_type_sc_pseudobulk_corrected_degs
  compare_second_df <- unfiltered_cell_type_validation_sc_degs[rownames(unfiltered_cell_type_validation_sc_degs) %in% primary_sc_pseudobulk_corrected_degs,]
  compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$sc_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc, method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_pseudobulk_corrected_vs_sc_validation"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_pseudobulk_corrected_vs_sc_validation"]] <- ggplot(data = comparing_first_vs_second_df, mapping = aes(x = first_fc, y = second_fc)) +
    geom_point(size = 2) +
    sm_statCorr(corr_method = "spearman") + xlab("Naive FC") + ylab("Vaccinated FC") + labs(title = cell_type_no_underscore) +  xlim(-2.5, 2.5) + ylim(-2.5, 2.5)
  
  
  
  # Check correlation for union of significant SC genes in primary and validation sets
  union_of_significant_degs <- union(rownames(cell_type_sc_degs), rownames(cell_type_validation_sc_degs))
  print(paste0("Number of union SC DEGs is: ", length(union_of_significant_degs)))
  compare_first_df <- unfiltered_cell_type_sc_degs[rownames(unfiltered_cell_type_sc_degs) %in% union_of_significant_degs,]
  compare_second_df <- unfiltered_cell_type_validation_sc_degs[rownames(unfiltered_cell_type_validation_sc_degs) %in% union_of_significant_degs,]
  compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$avg_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  #print(correlation_val$estimate)
  #print(correlation_val$p.value)
  sc_correlations[[cell_type]][["union_sc_primary_vs_sc_validation"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["union_sc_primary_vs_sc_validation"]] <- ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
    geom_point() +
    geom_smooth(method=lm) + xlab("Naive FC") + ylab("Vaccinated FC")
  
  
  # Check correlation for union of significant SC genes (pseudobulk corrected) in primary and validation sets
  union_of_significant_degs <- union(rownames(cell_type_sc_pseudobulk_corrected_degs), rownames(cell_type_validation_sc_pseudobulk_corrected_degs))
  print(paste0("Number of union SC (pseudobulk corrected) DEGs is: ", length(union_of_significant_degs)))
  compare_first_df <- unfiltered_cell_type_sc_degs[rownames(unfiltered_cell_type_sc_degs) %in% union_of_significant_degs,]
  compare_second_df <- unfiltered_cell_type_validation_sc_degs[rownames(unfiltered_cell_type_validation_sc_degs) %in% union_of_significant_degs,]
  compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$avg_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  #print(correlation_val$estimate)
  #print(correlation_val$p.value)
  sc_correlations[[cell_type]][["union_pseudobulk_corrected_sc_primary_vs_sc_validation"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["union_pseudobulk_corrected_sc_primary_vs_sc_validation"]] <- ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
    geom_point() +
    geom_smooth(method=lm) + xlab("Naive FC") + ylab("Vaccinated FC")
}

pseudobulk_corrected_plots <- lapply(sc_correlation_plots, function(x) x[[2]])
ggsave("C:/Users/willi/Desktop/hvl_vaccinated_vs_naive_deg_correlation.png", plot = patchwork::wrap_plots(pseudobulk_corrected_plots, ncol = 2, nrow = 3), height = 10, width = 10)
