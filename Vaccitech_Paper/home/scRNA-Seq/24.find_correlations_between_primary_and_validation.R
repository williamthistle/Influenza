# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# I picked my 9 favorite cell types so I can have a 3x3 grid of correlation plots
correlation_cell_types <- c("CD4_Memory", "CD8_Memory", "cDC", "NK", "CD14_Mono", "CD16_Mono", "MAIT", "B_naive", "B_memory")

sc_correlations <- list()
sc_correlation_plots <- list()

for(cell_type in correlation_cell_types) {
  print(cell_type)
  sc_correlations[[cell_type]] <- list()
  sc_correlation_plots[[cell_type]] <- list()
  cell_type_no_underscore <- gsub("_", " ", cell_type)
  # Unfiltered SC
  unfiltered_cell_type_sc_degs <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                             sep = "\t", header = TRUE)
  unfiltered_cell_type_validation_sc_degs <- read.table(paste0(sc_deg_validation_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                        sep = "\t", header = TRUE)
  # Filtered SC
  cell_type_sc_degs <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc.tsv"),
                                          sep = "\t", header = TRUE)
  cell_type_validation_sc_degs <- read.table(paste0(sc_deg_validation_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc.tsv"),
                                                     sep = "\t", header = TRUE)
  # Pseudobulk corrected
  cell_type_sc_pseudobulk_corrected_degs <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final.tsv"),
                                  sep = "\t", header = TRUE)
  rownames(cell_type_sc_pseudobulk_corrected_degs) <- cell_type_sc_pseudobulk_corrected_degs$Gene_Name
  cell_type_validation_sc_pseudobulk_corrected_degs <- read.table(paste0(sc_deg_validation_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_final.tsv"),
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
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_vs_sc_validation"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_vs_sc_validation"]] <- ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
    geom_point() +
    geom_smooth(method=lm) + xlab("Primary FC") + ylab("Validation FC")
  
  
  
  # Check correlation for primary significant SC genes (pseudobulk_corrected) in validation set
  primary_sc_pseudobulk_corrected_degs <- rownames(cell_type_sc_pseudobulk_corrected_degs)
  print(paste0("Number of SC DEGs (pseudobulk corrected) in primary data is: ", length(primary_sc_pseudobulk_corrected_degs)))
  compare_first_df <- primary_sc_pseudobulk_corrected_degs
  compare_second_df <- unfiltered_cell_type_validation_sc_degs[rownames(unfiltered_cell_type_validation_sc_degs) %in% primary_sc_pseudobulk_corrected_degs,]
  compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$sc_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  # Calculate correlation
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc)
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["sc_primary_pseudobulk_corrected_vs_sc_validation"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["sc_primary_pseudobulk_corrected_vs_sc_validation"]] <- ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
    geom_point() +
    geom_smooth(method=lm) + xlab("Primary FC") + ylab("Validation FC")
  
  
  
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
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["union_sc_primary_vs_sc_validation"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["union_sc_primary_vs_sc_validation"]] <- ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
    geom_point() +
    geom_smooth(method=lm) + xlab("Primary FC") + ylab("Validation FC")
  
  
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
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  sc_correlations[[cell_type]][["union_pseudobulk_corrected_sc_primary_vs_sc_validation"]] <- correlation_val
  
  # Plot correlation
  sc_correlation_plots[[cell_type]][["union_pseudobulk_corrected_sc_primary_vs_sc_validation"]] <- ggplot(comparing_first_vs_second_df, aes(x=first_fc, y=second_fc)) + 
    geom_point() +
    geom_smooth(method=lm) + xlab("Primary FC") + ylab("Validation FC")
}


