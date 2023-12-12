# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Function to check for overlapping ranges
check_overlap <- function(start1, end1, start2, end2) {
  overlap <- (start1 <= end2) & (end1 >= start2)
  return(overlap)
}

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
    overlapping_snME_and_all_peaks_cell_type_subset <- overlapping_snME_and_all_peaks[overlapping_snME_and_all_peaks$cell_type_snME == snME_cell_type,]
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

for(chr in unique(sc_peaks$value)) {
  print(chr)
  sc_peaks_chr_subset <- sc_peaks[sc_peaks$value == chr,]
  sc_peaks_lenient_chr_subset <- sc_peaks_lenient[sc_peaks_lenient$value == chr,]
  snME_dms_chr_subset <- snME_dms[snME_dms$chr == chr,]
  
  overlap_indices_snME_dms <- list()
  overlap_indices_snME_dms_lenient <- list()
  
  # Iterate through rows of B and check for overlaps with A
  for (i in 1:nrow(snME_dms_chr_subset)) {
    overlap <- check_overlap(snME_dms_chr_subset$start[i], snME_dms_chr_subset$end[i], sc_peaks_chr_subset$start, sc_peaks_chr_subset$end)
    if (any(overlap)) {
      overlap_indices_snME_dms[[i]] <- which(overlap)
    } else {
      overlap_indices_snME_dms[[i]] <- NA
    }
    overlap <- check_overlap(snME_dms_chr_subset$start[i], snME_dms_chr_subset$end[i], sc_peaks_lenient_chr_subset$start, sc_peaks_lenient_chr_subset$end)
    if (any(overlap)) {
      overlap_indices_snME_dms_lenient[[i]] <- which(overlap)
    } else {
      overlap_indices_snME_dms_lenient[[i]] <- NA
    }
  }
  
  for(current_snME_row_index in 1:length(overlap_indices_snME_dms)) {
    current_peak_row_index <- overlap_indices_snME_dms[[current_snME_row_index]]
    if(!is.na(current_peak_row_index)) {
      current_snME_row <- snME_dms_chr_subset[current_snME_row_index,]
      current_peak_row <- sc_peaks_chr_subset[current_peak_row_index,]
      current_peak_row$start_methyl <- current_snME_row$start
      current_peak_row$end_methyl <- current_snME_row$end
      current_peak_row$methyl_status <- current_snME_row$methylation
      current_peak_row$number_of_dms <- current_snME_row$number_of_dms
      current_peak_row$cell_type_snME <- current_snME_row$celltype
      filtered_rows[[length(filtered_rows) + 1]] <- current_peak_row
    }
  }
  
  for(current_snME_row_index in 1:length(overlap_indices_snME_dms_lenient)) {
    current_peak_row_index <- overlap_indices_snME_dms_lenient[[current_snME_row_index]]
    if(sum(!is.na(current_peak_row_index)) >= 1) {
      for(entry in current_peak_row_index) {
        current_snME_row <- snME_dms_chr_subset[current_snME_row_index,]
        current_peak_row <- sc_peaks_lenient_chr_subset[entry,]
        current_peak_row$start_methyl <- current_snME_row$start
        current_peak_row$end_methyl <- current_snME_row$end
        current_peak_row$methyl_status <- current_snME_row$methylation
        current_peak_row$number_of_dms <- current_snME_row$number_of_dms
        current_peak_row$cell_type_snME <- current_snME_row$celltype
        filtered_rows_lenient[[length(filtered_rows_lenient) + 1]] <- current_peak_row
      } 
    }
  }
}

# 5498 peaks overlap with snME DMR
snATAC_snME_overlap <- do.call(rbind, filtered_rows)
write.table(snATAC_snME_overlap, file = paste0(snME_results_dir, "snME_dms_overlap_with_all_atac_peaks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 10424 peaks (lenient) overlap with snME DMR
snATAC_snME_overlap_lenient <- do.call(rbind, filtered_rows_lenient)
write.table(snATAC_snME_overlap_lenient, file = paste0(snME_results_dir, "snME_dms_overlap_with_all_atac_peaks_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

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
                                                                                  snATAC_snME_overlap_lenient, lenient = TRUE)
  write.table(cell_type_snATAC_snME_overlap_das_subset, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cell_type_snATAC_snME_overlap_das_subset_lenient, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
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
  if(nrow(magical_results_cell_type_subset_overlap) > 0) {
    write.table(magical_results_cell_type_subset_overlap, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_MAGICAL.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
 
  # Lenient
  peaks <- cell_type_snATAC_snME_overlap_das_subset_lenient$Peak_Name
  chromosomes <- sapply(strsplit(peaks, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 3))
  cell_type_snATAC_snME_overlap_das_subset_lenient$Peak_chr <- chromosomes
  cell_type_snATAC_snME_overlap_das_subset_lenient$Peak_start <- start_coords
  cell_type_snATAC_snME_overlap_das_subset_lenient$Peak_end <- end_coords
  magical_results_cell_type_subset_overlap_lenient <- dplyr::inner_join(magical_results_cell_type_subset, cell_type_snATAC_snME_overlap_das_subset_lenient, by = c("Peak_chr", "Peak_start", "Peak_end"))
  if(nrow(magical_results_cell_type_subset_overlap_lenient) > 0) {
    write.table(magical_results_cell_type_subset_overlap_lenient, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_MAGICAL_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
