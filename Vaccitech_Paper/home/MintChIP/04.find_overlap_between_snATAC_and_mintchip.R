# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Function to check for overlapping ranges
check_overlap <- function(start1, end1, start2, end2) {
  overlap <- (start1 <= end2) & (end1 >= start2)
  return(overlap)
}

find_matching_snATAC_das_and_mintchip_das <- function(snATAC_cell_type, neg_peak_file_path, pos_peak_file_path, marker, mintchip_file_path, mintchip_hg38_file_path, lenient = FALSE) {
  # Grab positive DAS and negative snaTAC DAS
  pos_peak_table <- read.table(pos_peak_file_path, sep = "\t", header = TRUE)
  neg_peak_table <- read.table(neg_peak_file_path, sep = "\t", header = TRUE)
  
  # Save coordinate information for each positive and negative peak so we can more easily find overlap
  pos_peaks <- pos_peak_table$Peak_Name
  pos_chromosomes <- sapply(strsplit(pos_peaks, "-"), `[`, 1)
  pos_start_coords <- as.numeric(sapply(strsplit(pos_peaks, "-"), `[`, 2))
  pos_end_coords <- as.numeric(sapply(strsplit(pos_peaks, "-"), `[`, 3))
  pos_peak_table$seqnames <- pos_chromosomes
  pos_peak_table$start <- pos_start_coords
  pos_peak_table$end <- pos_end_coords
  
  neg_peaks <- neg_peak_table$Peak_Name
  neg_chromosomes <- sapply(strsplit(neg_peaks, "-"), `[`, 1)
  neg_start_coords <- as.numeric(sapply(strsplit(neg_peaks, "-"), `[`, 2))
  neg_end_coords <- as.numeric(sapply(strsplit(neg_peaks, "-"), `[`, 3))
  neg_peak_table$seqnames <- neg_chromosomes
  neg_peak_table$start <- neg_start_coords
  neg_peak_table$end <- neg_end_coords
  
  # Grab mintchip DAS
  mintchip_table <- read.table(mintchip_file_path, sep = "\t", header = TRUE)
  mintchip_table_for_hg38_coordinates <- read.table(mintchip_hg38_file_path, sep = "\t", header = TRUE)
  hg38_start <- c()
  hg38_end <- c()
  for(current_row_num in 1:nrow(mintchip_table)) {
    current_mintchip_row <- mintchip_table[current_row_num,]
    matching_hg38_mintchip_row <- mintchip_table_for_hg38_coordinates[mintchip_table_for_hg38_coordinates$seqnames == current_mintchip_row$seqnames
                                                                      & mintchip_table_for_hg38_coordinates$start == current_mintchip_row$start
                                                                      & mintchip_table_for_hg38_coordinates$end == current_mintchip_row$end,]
    if(nrow(matching_hg38_mintchip_row) > 0) {
      components <- unlist(strsplit(matching_hg38_mintchip_row$hg38_coordinates, "[:-]"))
      hg38_start <- c(hg38_start, components[2])
      hg38_end <- c(hg38_end, components[3])
    } else {
      hg38_start <- c(hg38_start, "REMOVE")
      hg38_end <- c(hg38_end, "REMOVE")
    }
  }
  mintchip_table$hg38_start <- hg38_start
  mintchip_table$hg38_end <- hg38_end
  mintchip_table <- subset(mintchip_table, !grepl("REMOVE", hg38_start))
  mintchip_table$hg38_start <- as.numeric(mintchip_table$hg38_start)
  mintchip_table$hg38_end <- as.numeric(mintchip_table$hg38_end)
  
  
  #pos_mintchip_table <- mintchip_table[mintchip_table$Fold > 0,]
  #neg_mintchip_table <- mintchip_table[mintchip_table$Fold < 0,]
  
  # Expand peaks (double length) if lenient is TRUE
  if(lenient) {
    pos_peak_table$start <- pos_peak_table$start - 250
    pos_peak_table$end <- pos_peak_table$end + 250
    neg_peak_table$start <- neg_peak_table$start - 250
    neg_peak_table$end <- neg_peak_table$end + 250
    
    mintchip_table$hg38_start <- mintchip_table$hg38_start - 200
    mintchip_table$hg38_end <- mintchip_table$hg38_end + 200
    # pos_mintchip_table$hg38_start <- pos_mintchip_table$hg38_start - 200
    # pos_mintchip_table$hg38_end <- pos_mintchip_table$hg38_end + 200
    # neg_mintchip_table$hg38_start <- neg_mintchip_table$hg38_start - 200
    # neg_mintchip_table$hg38_end <- neg_mintchip_table$hg38_end + 200
  }
  
  filtered_rows <- list()
  print("Positive peaks")
  positive_overlap_indices <- list()
  for(chr in unique(pos_peak_table$seqnames)) {
    print(chr)
    pos_peak_table_chr_subset <- pos_peak_table[pos_peak_table$seqnames == chr,]
    mintchip_table_chr_subset <- mintchip_table[mintchip_table$seqnames == chr,]
    for (i in 1:nrow(mintchip_table_chr_subset)) {
      overlap <- check_overlap(mintchip_table_chr_subset$hg38_start[i], mintchip_table_chr_subset$hg38_end[i], pos_peak_table_chr_subset$start, pos_peak_table_chr_subset$end)
      if (any(overlap)) {
        print(overlap)
        positive_overlap_indices[[i]] <- which(overlap)
      } else {
        positive_overlap_indices[[i]] <- NA
      }
    }
    
    for(current_mintchip_row_index in 1:length(positive_overlap_indices)) {
      current_peak_row_index <- positive_overlap_indices[[current_mintchip_row_index]]
      if(!is.na(current_peak_row_index)) {
        current_mintchip_row <- mintchip_table[current_mintchip_row_index,]
        current_peak_row <- pos_peak_table_chr_subset[current_peak_row_index,]
        current_peak_row$mintchip_marker <- marker
        current_peak_row$mintchip_start <- current_mintchip_row$start
        current_peak_row$mintchip_end <- current_mintchip_row$end
        current_peak_row$mintchip_fc <- current_mintchip_row$Fold      
        current_peak_row$mintchip_pvalue <- current_mintchip_row$p.value   
        filtered_rows[[length(filtered_rows) + 1]] <- current_peak_row
      }
    }
  }
  
  snATAC_mintchip_overlap <- do.call(rbind, filtered_rows)
  

    
    sc_peaks_lenient_chr_subset <- sc_peaks_lenient[sc_peaks_lenient$value == chr,]
    snME_dms_chr_subset <- snME_dms[snME_dms$chr == chr,]
    
    overlap_indices_snME_dms <- list()
    overlap_indices_snME_dms_lenient <- list()
    
    
    
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


# We will use consensus peak FC 0 threshold for liftover (hg19 -> hg38)

mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")

for(marker in mintchip_markers) {
  current_hg19_sites <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0.tsv"), sep = "\t", header = TRUE)
  current_hg19_sites <- current_hg19_sites$coordinates
  write.table(current_hg19_sites, file = paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_coordinates.tsv"),
                                                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Go to UCSC site to perform liftover

# Add hg38 coordinates to each file
for(marker in mintchip_markers) {
  current_marker_sites <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0.tsv"), sep = "\t", header = TRUE)
  associated_hg38_coordinates <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_coordinates_hg38.tsv"), 
                                            sep = "\t")$V1
  if(file.size(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_coordinates_hg38_error.tsv")) > 2) {
    associated_hg38_errors <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_coordinates_hg38_error.tsv"), 
                                         sep = "\t")$V1
    # Remove any site that didn't have associated hg38 site in liftover
    current_marker_sites <- current_marker_sites[current_marker_sites$coordinates != associated_hg38_errors, ]
  }
  # Add hg38 info
  current_marker_sites$hg38_coordinates <- associated_hg38_coordinates
  write.table(current_marker_sites, file = paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_with_hg38_coordinates.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Next, we want to find overlap between each cell type and each marker
# Define FC threshold above method (maybe should be different FC threshold for each marker?)
atac_cell_types_for_mintchip_analysis <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "MAIT", "NK", "Proliferating", "T Naive")

for(atac_cell_type in atac_cell_types_for_mintchip_analysis) {
  print(atac_cell_type)
  for(marker in mintchip_markers) {
    print(marker)
    atac_cell_type_for_file_name <- sub(" ", "_", atac_cell_type)
    cell_type_snATAC_snME_overlap_das_subset <- find_matching_snATAC_das_and_mintchip_das(cell_type, paste0(sc_das_dir, 
                                                                                                        "diff_peaks/D28-vs-D_minus_1-degs-", atac_cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_neg.tsv"),
                                                                                      paste0(sc_das_dir, 
                                                                                             "diff_peaks/D28-vs-D_minus_1-degs-", atac_cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01_pos.tsv"),
                                                                                      marker, paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0.1.tsv"), paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_with_hg38_coordinates.tsv"), lenient = FALSE)
  }
}



