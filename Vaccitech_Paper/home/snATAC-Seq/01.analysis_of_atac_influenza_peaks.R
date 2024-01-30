# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

snATAC_get_list_of_closest_genes <- function(file_path, sc_peaks) {
  current_peaks <- read.table(file_path, sep = "\t", header = TRUE)
  current_peaks <- current_peaks$Peak_Name
  chromosomes <- sapply(strsplit(current_peaks, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(current_peaks, "-"), `[`, 3))
  current_peaks_coords <- data.frame(seqnames = chromosomes, start = start_coords, end = end_coords)
  closest_genes <- dplyr::inner_join(sc_peaks, current_peaks_coords, by = c("seqnames", "start", "end"))
  return(closest_genes$nearestGene)
}

run_fmd_on_snATAC <- function(snATAC_list) {
  fmd_results <- list()
  index <- 1
  for(current_nearest_gene_list in snATAC_list) {
    if(length(current_nearest_gene_list) > 1 && length(current_nearest_gene_list) < 2000) {
      current_fmd_result <- SPEEDI::RunFMD_RNA(current_nearest_gene_list, "blood")
    } else { 
      current_fmd_result <- "EMPTY OR OVER 2000 PEAKS (TOO MANY)"
    }
    fmd_results[[index]] <- current_fmd_result
    index <- index + 1
  }
  return(fmd_results)
}

# Read in relevant peaks
snATAC_peaks_for_hb_files <- list()

# Write annotated peaks down
cell_types <- c("B", "CD4_Memory", "CD8_Memory", "CD14_Mono", "CD16_Mono", "MAIT", "NK", "Proliferating", "T_Naive")
for(cell_type in cell_types) {
  input_neg_file_path <- paste0(sc_das_dir, 
                                "diff_peaks/D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_neg.tsv")
  annotated_peaks_neg <- snATAC_get_list_of_closest_genes(file_path = input_neg_file_path, sc_peaks = sc_peaks)
  output_neg_file_path <- paste0(sc_das_dir, 
                                 "diff_peaks/D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_neg_annotated.tsv")
  write.table(annotated_peaks_neg, file = output_neg_file_path, sep = "\t", quote = FALSE)
  
  input_pos_file_path <- paste0(sc_das_dir, 
                                "diff_peaks/D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_pos.tsv")
  annotated_peaks_pos <- snATAC_get_list_of_closest_genes(file_path = input_pos_file_path, sc_peaks = sc_peaks)
  output_pos_file_path <- paste0(sc_das_dir, 
                                 "diff_peaks/D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_pos_annotated.tsv")
  write.table(annotated_peaks_pos, file = output_pos_file_path, sep = "\t", quote = FALSE)
}

# NOTE: I may want to do logFC < -1 for negative to get better signal. 
snATAC_nearest_gene_FMD_results <- run_fmd_on_snATAC(snATAC_peaks_for_hb_files)

# saveRDS(snATAC_nearest_gene_FMD_results, file = paste0(sc_das_dir, "snATAC_nearest_gene_FMD_results.rds"))
# snATAC_nearest_gene_FMD_results <- readRDS(file = paste0(sc_das_dir, "snATAC_nearest_gene_FMD_results.rds"))

snATAC_nearest_gene_FMD_results



























# Alternative strategy: using same peaks as motif discovery
# I think it's worse
# CD4 Memory

names(snATAC_peaks_for_hb_files)[[1]] <- "CD4_Memory_Negative"
# CD8 Memory
snATAC_peaks_for_hb_files[[2]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-CD8_Memory-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_neg.tsv"),
                                                                   logFC = -0.585, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[2]] <- "CD8_Memory_Negative"
snATAC_peaks_for_hb_files[[3]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-CD8_Memory-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_pos.tsv"),
                                                                   logFC = 0.585, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[3]] <- "CD8_Memory_Positive"
# CD14 Mono
snATAC_peaks_for_hb_files[[4]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-CD14_Mono-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_neg.tsv"),
                                                                   logFC = -1, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[4]] <- "CD14_Mono_Negative"
snATAC_peaks_for_hb_files[[5]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-CD14_Mono-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_pos.tsv"),
                                                                   logFC = 1, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[5]] <- "CD14_Mono_Positive"
# CD16 Mono
snATAC_peaks_for_hb_files[[6]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-CD16_Mono-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_neg.tsv"),
                                                                   logFC = -1, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[6]] <- "CD16_Mono_Negative"
snATAC_peaks_for_hb_files[[7]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-CD16_Mono-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_pos.tsv"),
                                                                   logFC = 1, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[7]] <- "CD16_Mono_Positive"
# T Naive
snATAC_peaks_for_hb_files[[8]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-T_Naive-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_neg.tsv"),
                                                                   logFC = -1, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[8]] <- "T_Naive_Negative"
snATAC_peaks_for_hb_files[[9]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-T_Naive-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_pos.tsv"),
                                                                   logFC = 1, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[9]] <- "T_Naive_Positive"
# T Naive
snATAC_peaks_for_hb_files[[10]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-B-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_neg.tsv"),
                                                                   logFC = -0.1, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[10]] <- "B_Negative"
snATAC_peaks_for_hb_files[[11]] <- snATAC_get_list_of_closest_genes(file_path = paste0(sc_das_dir, 
                                                                                      "diff_peaks/D28-vs-D_minus_1-degs-B-time_point-controlling_for_subject_id_sc_filtered_pct_0.01_pos.tsv"),
                                                                   logFC = 0.1, sc_peaks = sc_peaks)
names(snATAC_peaks_for_hb_files)[[11]] <- "B_Positive"


# I can also do liftover for mintchip and then compare overlap between mintchip markers and my peaks
# I can also look at overlap between closest genes and genes associataed with specific markers
# For example, are we seeing a lot of one or two types of markers in CD14 Mono?
# I guess I could also do innate immune system and adaptive, but I'll delay that for now



















# Subset sc_peaks to only have the lenient DAS
colnames(sc_peaks)[colnames(sc_peaks) == "value"] <- "chr"
sc_peaks_lenient_subset <- merge(sc_peaks, sc_das_lenient, by = c("chr", "idx"))
colnames(sc_peaks_lenient_subset)[colnames(sc_peaks_lenient_subset) == "start.x"] <- "start"
colnames(sc_peaks_lenient_subset)[colnames(sc_peaks_lenient_subset) == "end.x"] <- "end"
sc_peaks_lenient_subset <- sc_peaks_lenient_subset[, -which(names(sc_peaks_lenient_subset) %in% c("width", "score", "replicateScoreQuantile", "Reproducibility", "names", "GroupReplicate", "groupScoreQuantile", "N", "start.y", "end.y"))]

sc_peaks_lenient_subset_extended <- sc_peaks_lenient_subset
sc_peaks_lenient_subset_extended$start <- sc_peaks_lenient_subset_extended$start - 250
sc_peaks_lenient_subset_extended$end <- sc_peaks_lenient_subset_extended$end + 250

mintchip_table_extended <- mintchip_table
mintchip_table_extended$start <- mintchip_table_extended$start - 200
mintchip_table_extended$end <- mintchip_table_extended$end + 200

# OVERLAPPING PEAKS
overlap_ranges <- data.frame()
for(current_chr in unique(sc_peaks_lenient_subset_extended$chr)) {
  sc_peaks_lenient_subset_extended_chr_subset <- sc_peaks_lenient_subset_extended[sc_peaks_lenient_subset_extended$chr == current_chr,]
  mintchip_table_extended_subset <- mintchip_table_extended[mintchip_table_extended$seqnames == current_chr,]
  
  # Iterate through each pair of ranges and check for overlaps
  for (i in 1:nrow(sc_peaks_lenient_subset_extended_chr_subset)) {
    for (j in 1:nrow(mintchip_table_extended_subset)) {
      if (sc_peaks_lenient_subset_extended_chr_subset$start[i] <= mintchip_table_extended_subset$end[j] && mintchip_table_extended_subset$start[j] <= sc_peaks_lenient_subset_extended_chr_subset$end[i]) {
        overlap_ranges <- rbind(overlap_ranges, c(current_chr, sc_peaks_lenient_subset_extended_chr_subset$Cell_Type[i], sc_peaks_lenient_subset_extended_chr_subset$start[i], sc_peaks_lenient_subset_extended_chr_subset$end[i], mintchip_table_extended_subset$start[j], mintchip_table_extended_subset$end[j]))
      }
    }
  }
  
  # Rename columns of the overlap_ranges data frame
  if(nrow(overlap_ranges) > 0) {
    colnames(overlap_ranges) <- c("chr", "Cell_Type", "ATAC_peaks_start", "ATAC_peaks_end", "mintchip_peaks_start", "mintchip_peaks_end")
  }
}

# OVERLAPPING GENES
mintchip_validated_genes <- sc_peaks_lenient_subset[sc_peaks_lenient_subset$nearestGene %in% mintchip_table$SYMBOL,]
# Run pos and neg genes through HumanBase - really good groupings
mintchip_validated_genes_pos <- mintchip_validated_genes[mintchip_validated_genes$sc_log2FC > 0,]
mintchip_validated_genes_neg <- mintchip_validated_genes[mintchip_validated_genes$sc_log2FC < 0,]
# What about HB within cell types?
# CD14 Mono
mintchip_validated_genes_cd14 <- mintchip_validated_genes[mintchip_validated_genes$Cell_Type == "CD14_Mono",]
