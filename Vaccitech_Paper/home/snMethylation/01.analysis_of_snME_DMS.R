# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

find_matching_snATAC_das_and_snME_dms <- function(snATAC_cell_type, neg_peak_file_path, pos_peak_file_path, overlapping_snME_and_all_peaks, lenient = FALSE) {
  # Find matching snME cell types
  if(snATAC_cell_type == "B") {
    matching_snME_cell_types <- c("B-Mem", "B-Naive")
  } else if(snATAC_cell_type == "CD14 Mono") {
    matching_snME_cell_types <- c("Monocyte")
  } else if(snATAC_cell_type == "CD16 Mono") {
    matching_snME_cell_types <- c("Monocyte")
  } else if(snATAC_cell_type == "NK") {
    matching_snME_cell_types <- c("NK-cell2")
  } else if(snATAC_cell_type == "T Naive") {
    matching_snME_cell_types <- c("Tc-Naive", "Th-Naive")
  } else if(snATAC_cell_type == "CD4 Memory") {
    matching_snME_cell_types <- c("Tc-Mem", "Th-Mem")
  } else if(snATAC_cell_type == "CD8 Memory") {
    matching_snME_cell_types <- c("Tc-Mem", "Th-Mem")
  } else {
    raise("ERROR!")
  }
  
  # We expanded the snATAC peaks to find matches between snME and snATAC peaks
  # To find whether those snATAC peaks are differentially accessible, we need to shrink them again so they match with our DAS tables
  if(lenient) {
    overlapping_snME_and_all_peaks$start <- overlapping_snME_and_all_peaks$start + 250
    overlapping_snME_and_all_peaks$end <- overlapping_snME_and_all_peaks$end - 250
  }
  
  # Grab positive DAS and negative DAS
  pos_peak_table <- read.table(pos_peak_file_path, sep = "\t", header = TRUE)
  neg_peak_table <- read.table(neg_peak_file_path, sep = "\t", header = TRUE)
  
  # Grab coordinate information for each positive and negative peak so we can find overlap between DAS and DMS
  pos_peaks <- pos_peak_table$Peak_Name
  pos_chromosomes <- sapply(strsplit(pos_peaks, "-"), `[`, 1)
  pos_start_coords <- as.numeric(sapply(strsplit(pos_peaks, "-"), `[`, 2))
  pos_end_coords <- as.numeric(sapply(strsplit(pos_peaks, "-"), `[`, 3))
  pos_peak_table$value <- pos_chromosomes
  pos_peak_table$start <- pos_start_coords
  pos_peak_table$end <- pos_end_coords
  
  neg_peaks <- neg_peak_table$Peak_Name
  neg_chromosomes <- sapply(strsplit(neg_peaks, "-"), `[`, 1)
  neg_start_coords <- as.numeric(sapply(strsplit(neg_peaks, "-"), `[`, 2))
  neg_end_coords <- as.numeric(sapply(strsplit(neg_peaks, "-"), `[`, 3))
  neg_peak_table$value <- neg_chromosomes
  neg_peak_table$start <- neg_start_coords
  neg_peak_table$end <- neg_end_coords
  
  final_list_of_sites <- list()
  index <- 1
  for(snME_cell_type in matching_snME_cell_types) {
    # Subset to current cell type in snME data
    overlapping_snME_and_all_peaks_cell_type_subset <- overlapping_snME_and_all_peaks[overlapping_snME_and_all_peaks$dms_cell_type == snME_cell_type,]
    # Find overlap between pos/neg peaks and snME DMS
    pos_peaks_overlapping_with_snME <- dplyr::inner_join(pos_peak_table, overlapping_snME_and_all_peaks_cell_type_subset, by = c("value", "start", "end"))
    if(nrow(pos_peaks_overlapping_with_snME) > 0) {
      for(i in 1:nrow(pos_peaks_overlapping_with_snME)) {
        final_list_of_sites[[index]] <- pos_peaks_overlapping_with_snME[i,]
        index <- index + 1
      }
    }

    neg_peaks_overlapping_with_snME <- dplyr::inner_join(neg_peak_table, overlapping_snME_and_all_peaks_cell_type_subset, by = c("value", "start", "end"))
    if(nrow(neg_peaks_overlapping_with_snME) > 0) {
      for(i in 1:nrow(neg_peaks_overlapping_with_snME)) {
        final_list_of_sites[[index]] <- neg_peaks_overlapping_with_snME[i,]
        index <- index + 1
      }
    }
  }
  
  final_list_of_sites <- do.call(rbind, final_list_of_sites)
  final_list_of_sites <- final_list_of_sites %>%
    select(-value, -start, -end)
  return(final_list_of_sites)
}

# Find all overlap between snATAC-seq (hg38) and snMethylation (hg38)

# Create lenient version of sc_peaks (adds 250 to start / end of coordinates)
sc_peaks_lenient <- sc_peaks
sc_peaks_lenient$start <- sc_peaks_lenient$start - 250 
sc_peaks_lenient$end <- sc_peaks_lenient$end + 250

# Initialize an empty list to store the filtered rows
# filtered_rows uses the original peaks
# filtered_rows_lenient is more lenient (adds 500 to start / end)
filtered_rows <- list()
filtered_rows_lenient <- list()

# Iterate through rows in snME_dms
for (i in 1:nrow(snME_dms)) {
  # Keep track of which row we're on
  if(i %% 10000 == 0) {
    print(i)
  }
  
  # Grab important info from snME_dms
  chr_val <- snME_dms[i, "chr"]
  coordinate_val <- snME_dms[i, "start"]
  cell_type <- snME_dms[i, "celltype"]
  methylation <- snME_dms[i, "methylation"]
  
  # Subset sc_peaks based on chr
  sc_peaks_chr_subset <- subset(sc_peaks, value == chr_val)
  sc_peaks_lenient_chr_subset <- subset(sc_peaks_lenient, value == chr_val)
  
  # Filter rows based on the condition (snME site within start and end)
  filtered_row <- subset(sc_peaks_chr_subset, coordinate_val >= start & coordinate_val <= end)
  filtered_row_lenient <- subset(sc_peaks_lenient_chr_subset, coordinate_val >= start & coordinate_val <= end)
  
  # If any matching rows are found, store them
  if (nrow(filtered_row) > 0) {
    filtered_row$dms_coordinate <- coordinate_val
    filtered_row$dms_cell_type <- cell_type
    filtered_row$dms_methylation <- methylation
    filtered_rows[[i]] <- filtered_row
  }
  if(nrow(filtered_row_lenient) > 0) {
    filtered_row_lenient$dms_coordinate <- coordinate_val
    filtered_row_lenient$dms_cell_type <- cell_type
    filtered_row_lenient$dms_methylation <- methylation
    filtered_rows_lenient[[i]] <- filtered_row_lenient
  }
}

# Combine the filtered rows into df
# 5494 peaks overlap with snME DMS
snATAC_snME_overlap <- do.call(rbind, filtered_rows)
write.table(snATAC_snME_overlap, file = paste0(snME_data_dir, "snME_dms_overlap_with_all_atac_peaks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 10420 peaks (lenient) overlap with snME DMS
snATAC_snME_lenient_overlap <- do.call(rbind, filtered_rows_lenient)
write.table(snATAC_snME_lenient_overlap, file = paste0(snME_data_dir, "snME_dms_overlap_with_all_atac_peaks_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The next question is: Are any of these peaks differentially accessible in the relevant cell type?
# If so, we have matching differential accessibility and differential methylation
atac_cell_types_for_snME_analysis <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "T Naive")

magical_output_dir <- paste0(sc_magical_dir, "Output/")
magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

for(cell_type in atac_cell_types_for_snME_analysis) {
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  cell_type_snATAC_snME_overlap_das_subset <- find_matching_snATAC_das_and_snME_dms(cell_type, paste0(sc_das_dir, 
                                                                                                                     "diff_peaks/D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_neg.tsv"),
                                                                          paste0(sc_das_dir, 
                                                                                 "diff_peaks/D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_pos.tsv"),
                                                                          snATAC_snME_overlap, lenient = FALSE)
  cell_type_snATAC_snME_overlap_das_subset_lenient <- find_matching_snATAC_das_and_snME_dms(cell_type, paste0(sc_das_dir, 
                                                                                                                             "diff_peaks/D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_neg.tsv"),
                                                                                  paste0(sc_das_dir, 
                                                                                         "diff_peaks/D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_pos.tsv"),
                                                                                  snATAC_snME_lenient_overlap, lenient = TRUE)
  write.table(cell_type_snATAC_snME_overlap_das_subset, file = paste0(snME_data_dir, cell_type_for_file_name, "_das_snME_dms_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cell_type_snATAC_snME_overlap_das_subset_lenient, file = paste0(snME_data_dir, cell_type_for_file_name, "_das_snME_dms_overlap_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Overlap between MAGICAL circuits and DAS peaks that contained DMS
  magical_results_cell_type_subset <- magical_results[magical_results$Cell_Type == cell_type_for_file_name,]
  
  # Normal
  peaks <- cell_type_snATAC_snME_overlap_das_subset$Peak_Name
  chromosomes <- sapply(strsplit(peaks, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 3))
  cell_type_snATAC_snME_overlap_das_subset$Peak_chr <- chromosomes
  cell_type_snATAC_snME_overlap_das_subset$Peak_start <- start_coords
  cell_type_snATAC_snME_overlap_das_subset$Peak_end <- end_coords
  magical_results_cell_type_subset_overlap <- dplyr::inner_join(magical_results_cell_type_subset, cell_type_snATAC_snME_overlap_das_subset, by = c("Peak_chr", "Peak_start", "Peak_end"))
  write.table(magical_results_cell_type_subset_overlap, file = paste0(snME_data_dir, cell_type_for_file_name, "_das_snME_dms_overlap_MAGICAL.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Lenient
  peaks <- cell_type_snATAC_snME_overlap_das_subset_lenient$Peak_Name
  chromosomes <- sapply(strsplit(peaks, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 3))
  cell_type_snATAC_snME_overlap_das_subset_lenient$Peak_chr <- chromosomes
  cell_type_snATAC_snME_overlap_das_subset_lenient$Peak_start <- start_coords
  cell_type_snATAC_snME_overlap_das_subset_lenient$Peak_end <- end_coords
  magical_results_cell_type_subset_overlap_lenient <- dplyr::inner_join(magical_results_cell_type_subset, cell_type_snATAC_snME_overlap_das_subset_lenient, by = c("Peak_chr", "Peak_start", "Peak_end"))
  write.table(magical_results_cell_type_subset_overlap_lenient, file = paste0(snME_data_dir, cell_type_for_file_name, "_das_snME_dms_overlap_MAGICAL_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
