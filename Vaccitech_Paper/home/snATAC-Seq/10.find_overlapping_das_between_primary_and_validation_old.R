# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# First argument is primary (placebo), second argument is validation (vaccinated)
find_sc_correlation_final <- function(first_gene_df, second_gene_df) {
  compare_first_df <- first_gene_df
  compare_second_df <- second_gene_df[rownames(second_gene_df) %in% compare_first_df$Peak_Name,]
  compare_first_df <- compare_first_df[compare_first_df$Peak_Name %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% compare_first_df$Peak_Name,]
  compare_first_df <- compare_first_df[order(compare_first_df$Peak_Name),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$sc_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  
  print("Calculating Spearman correlation")
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc, method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  
  return(as.vector(correlation_val$estimate))
}

# First argument is primary (placebo), second argument is validation (vaccinated)
find_sc_correlation_sc <- function(first_gene_df, second_gene_df) {
  compare_first_df <- first_gene_df
  compare_second_df <- second_gene_df[rownames(second_gene_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[rownames(compare_first_df) %in% rownames(compare_second_df),]
  compare_second_df <- compare_second_df[rownames(compare_second_df) %in% rownames(compare_first_df),]
  compare_first_df <- compare_first_df[order(rownames(compare_first_df)),]
  compare_second_df <- compare_second_df[order(rownames(compare_second_df)),]
  
  comparing_first_vs_second_df <- data.frame(gene_name = rownames(compare_first_df), first_fc = compare_first_df$avg_log2FC,
                                             second_fc = compare_second_df$avg_log2FC)
  
  print("Calculating Spearman correlation")
  correlation_val <- cor.test(comparing_first_vs_second_df$first_fc, comparing_first_vs_second_df$second_fc, method = "spearman")
  print(correlation_val$estimate)
  print(correlation_val$p.value)
  
  return(correlation_val)
}


overlapping_das_cell_types <- c("CD14_Mono", "CD16_Mono", "cDC","NK", "CD4_Memory", "CD8_Memory", "CD8_Naive", "MAIT", "B", "Proliferating", "CD4_Naive")

for(cell_type in overlapping_das_cell_types) {
  # Unfiltered
  unfiltered_cell_type_sc_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "Old/D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                            sep = "\t", header = TRUE)
  unfiltered_cell_type_validation_sc_das <- read.table(paste0(scATAC_hvl_vaccinated_das_dir, "Old/D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_sc_unfiltered.tsv"),
                                                       sep = "\t", header = TRUE)
  for(analysis_type in c("sc", "final")) {
    for(min_pct in c(0.01, 0.05, 0.1)) {
      print(paste0("Overlap and correlations for ", cell_type, " (", analysis_type, ", ", min_pct, ")"))
      # Filtered
      placebo_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", min_pct, ".tsv"),
                                     sep = "\t", header = TRUE)
      vaccination_das <- read.table(paste0(scATAC_hvl_vaccinated_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", min_pct, ".tsv"),
                                    sep = "\t", header = TRUE)
      if(analysis_type == "sc") {
        overlapping_das <- intersect(rownames(placebo_das), rownames(vaccination_das))
        print(paste0("Total overlapping DAS: ", length(overlapping_das)))
        print(paste0("Total placebo DAS: ", nrow(placebo_das)))
        print(paste0("Percent of placebo DAS that overlap: ", (length(overlapping_das) / nrow(placebo_das) * 100)))
        #find_sc_correlation_sc(placebo_das, unfiltered_cell_type_validation_sc_das)
      } else {
        overlapping_das_with_robust <- intersect(placebo_das$Peak_Name, vaccination_das$Peak_Name)
        print(paste0("Total overlapping DAS (with DESeq2 pseudobulk): ", length(overlapping_das_with_robust)))
        print(paste0("Total placebo DAS (with DESeq2 pvalue): ", nrow(placebo_das)))
        print(paste0("Percent of placebo DAS that overlap (with DESeq2 pvalue): ", (length(overlapping_das_with_robust) / nrow(placebo_das) * 100)))
        #find_sc_correlation_final(placebo_das, unfiltered_cell_type_validation_sc_das)
      }
    }
  }
}
