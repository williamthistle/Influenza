# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

combineRows <- function(df1, df2, matches) {
  combined <- data.frame()
  for (i in 1:nrow(matches)) {
    queryRow <- matches[i, "queryHits"]
    subjectRow <- matches[i, "subjectHits"]
    combinedRow <- cbind(df1[queryRow,], df2[subjectRow,])
    combined <- rbind(combined, combinedRow)
  }
  return(combined)
}

# Fix hg19 coordinates to be hg38
add_hg38_coordinates_to_marker_peaks <- function(marker_table, hg38_table) {
  marker_peaks_granges <- makeGRangesFromDataFrame(df = marker_table, keep.extra.columns = TRUE)
  marker_peaks_with_hg38_coordinates_granges <- makeGRangesFromDataFrame(df = hg38_table, keep.extra.columns = TRUE)
  marker_overlap <- as.data.frame(findOverlaps(marker_peaks_granges, marker_peaks_with_hg38_coordinates_granges))
  marker_overlap_df <- combineRows(marker_table, hg38_table, marker_overlap)
  components <- strsplit(marker_overlap_df$hg38_coordinates, "[:-]")
  hg38_start <- c()
  hg38_end <- c()
  for(component in components) {
    hg38_start <- c(hg38_start, component[2])
    hg38_end <- c(hg38_end, component[3])
  }
  marker_overlap_df$start <- hg38_start
  marker_overlap_df$end <- hg38_end
  marker_overlap_df <- marker_overlap_df[,-c(12:25)]
  return(marker_overlap_df)
}

find_matching_snATAC_das_and_mintchip_das <- function(snATAC_cell_type, neg_peak_file_path, pos_peak_file_path, mintchip_markers, fc_threshold, lenient = FALSE) {
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
  
  # Expand peaks (double length) if lenient is TRUE
  if(lenient) {
    pos_peak_table$start <- pos_peak_table$start - 250
    pos_peak_table$end <- pos_peak_table$end + 250
    neg_peak_table$start <- neg_peak_table$start - 250
    neg_peak_table$end <- neg_peak_table$end + 250
  }
  
  filtered_rows <- list()
  # Traverse markers one at a time
  for(marker in mintchip_markers) {
    print(marker)
    # Grab marker peaks and file containing hg38 coordinates (since markers originally were hg19)
    current_marker_peaks_file_path <- paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_", fc_threshold, ".tsv")
    current_marker_hg38_coordinate_file_path <- paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_with_hg38_coordinates.tsv")
    current_marker_peaks <- read.table(current_marker_peaks_file_path, sep = "\t", header = TRUE)
    current_marker_hg38_coordinates <- read.table(current_marker_hg38_coordinate_file_path, sep = "\t", header = TRUE)
    hg38_coords <- c()
    hg38_start <- c()
    hg38_end <- c()
    for(current_row_num in 1:nrow(current_marker_peaks)) {
      current_mintchip_row <- current_marker_peaks[current_row_num,]
      matching_hg38_mintchip_row <- current_marker_hg38_coordinates[current_marker_hg38_coordinates$seqnames == current_mintchip_row$seqnames
                                                                        & current_marker_hg38_coordinates$start == current_mintchip_row$start
                                                                        & current_marker_hg38_coordinates$end == current_mintchip_row$end,]
      if(nrow(matching_hg38_mintchip_row) > 0) {
        components <- unlist(strsplit(matching_hg38_mintchip_row$hg38_coordinates, "[:-]"))
        hg38_coords <- c(hg38_coords, matching_hg38_mintchip_row$hg38_coordinates)
        hg38_start <- c(hg38_start, components[2])
        hg38_end <- c(hg38_end, components[3])
      } else {
        hg38_coords <- c(hg38_coords, "REMOVE")
        hg38_start <- c(hg38_start, "REMOVE")
        hg38_end <- c(hg38_end, "REMOVE")
      }
    }
    current_marker_peaks$hg38_coordinates <- hg38_coords
    current_marker_peaks$hg38_start <- hg38_start
    current_marker_peaks$hg38_end <- hg38_end
    current_marker_peaks <- subset(current_marker_peaks, !grepl("REMOVE", hg38_start))
    current_marker_peaks$hg38_start <- as.numeric(current_marker_peaks$hg38_start)
    current_marker_peaks$hg38_end <- as.numeric(current_marker_peaks$hg38_end)
    
    # Is this OK?
    if(lenient) {
      current_marker_peaks$start <- current_marker_peaks$start - 200
      current_marker_peaks$end <- current_marker_peaks$end + 200
      current_marker_peaks$hg38_start <- current_marker_peaks$hg38_start - 200
      current_marker_peaks$hg38_end <- current_marker_peaks$hg38_end + 200
    }
    
    for(chr in unique(pos_peak_table$seqnames)) {
      positive_overlap_indices <- list()
      pos_peak_table_chr_subset <- pos_peak_table[pos_peak_table$seqnames == chr,]
      mintchip_table_chr_subset <- current_marker_peaks[current_marker_peaks$seqnames == chr,]
      if(nrow(mintchip_table_chr_subset) == 0) {
        next
      }
      for (i in 1:nrow(mintchip_table_chr_subset)) {
        overlap <- check_overlap(mintchip_table_chr_subset$hg38_start[i], mintchip_table_chr_subset$hg38_end[i], pos_peak_table_chr_subset$start, pos_peak_table_chr_subset$end)
        if (any(overlap)) {
          positive_overlap_indices[[i]] <- which(overlap)
        } else {
          positive_overlap_indices[[i]] <- NA
        }
      }
      
      for(current_mintchip_row_index in 1:length(positive_overlap_indices)) {
        current_peak_row_index <- positive_overlap_indices[[current_mintchip_row_index]]
        if(!is.na(current_peak_row_index)) {
          current_mintchip_row <- mintchip_table_chr_subset[current_mintchip_row_index,]
          current_peak_row <- pos_peak_table_chr_subset[current_peak_row_index,]
          current_peak_row$mintchip_marker <- marker
          current_peak_row$mintchip_start <- current_mintchip_row$hg38_start
          current_peak_row$mintchip_end <- current_mintchip_row$hg38_end
          current_peak_row$mintchip_fc <- current_mintchip_row$Fold      
          current_peak_row$mintchip_pvalue <- current_mintchip_row$p.value   
          filtered_rows[[length(filtered_rows) + 1]] <- current_peak_row
        }
      }
    }
    
    for(chr in unique(neg_peak_table$seqnames)) {
      negative_overlap_indices <- list()
      neg_peak_table_chr_subset <- neg_peak_table[neg_peak_table$seqnames == chr,]
      mintchip_table_chr_subset <- current_marker_peaks[current_marker_peaks$seqnames == chr,]
      if(nrow(mintchip_table_chr_subset) == 0) {
        next
      }
      for (i in 1:nrow(mintchip_table_chr_subset)) {
        overlap <- check_overlap(mintchip_table_chr_subset$hg38_start[i], mintchip_table_chr_subset$hg38_end[i], neg_peak_table_chr_subset$start, neg_peak_table_chr_subset$end)
        if (any(overlap)) {
          negative_overlap_indices[[i]] <- which(overlap)
        } else {
          negative_overlap_indices[[i]] <- NA
        }
      }
      
      for(current_mintchip_row_index in 1:length(negative_overlap_indices)) {
        current_peak_row_index <- negative_overlap_indices[[current_mintchip_row_index]]
        if(!is.na(current_peak_row_index)) {
          current_mintchip_row <- mintchip_table_chr_subset[current_mintchip_row_index,]
          current_peak_row <- neg_peak_table_chr_subset[current_peak_row_index,]
          current_peak_row$mintchip_marker <- marker
          current_peak_row$mintchip_start <- current_mintchip_row$hg38_start
          current_peak_row$mintchip_end <- current_mintchip_row$hg38_end
          current_peak_row$mintchip_fc <- current_mintchip_row$Fold      
          current_peak_row$mintchip_pvalue <- current_mintchip_row$p.value   
          filtered_rows[[length(filtered_rows) + 1]] <- current_peak_row
        }
      }
    }
  }
  snATAC_mintchip_overlap <- do.call(rbind, filtered_rows)
  return(snATAC_mintchip_overlap)
}

# MAGICAL info
magical_output_dir <- paste0(sc_magical_dir, "Output/")
magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

# Find all overlap between snATAC-seq (hg38) and MintChIP (hg38)

# Create lenient version of sc_peaks (adds 250 to start / end of coordinates)
colnames(sc_peaks)[1] <- "seqnames"
sc_peaks_lenient <- sc_peaks
sc_peaks_lenient$start <- sc_peaks_lenient$start - 250 
sc_peaks_lenient$end <- sc_peaks_lenient$end + 250
sc_peaks_more_lenient <- sc_peaks
sc_peaks_more_lenient$start <- sc_peaks_lenient$start - 250 
sc_peaks_more_lenient$end <- sc_peaks_lenient$end + 250

sc_peaks_granges <- makeGRangesFromDataFrame(df = sc_peaks, keep.extra.columns = TRUE)
sc_peaks_lenient_granges <- makeGRangesFromDataFrame(df = sc_peaks_lenient, keep.extra.columns = TRUE)
sc_peaks_more_lenient_granges <- makeGRangesFromDataFrame(df = sc_peaks_more_lenient, keep.extra.columns = TRUE)

# Mintchip markers
mintchip_markers <- c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27Ac", "H3K27me3", "H3K36me3")

# We will use consensus peak FC 0 threshold for liftover (hg19 -> hg38)
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

# Find overlap between ALL marker peaks and ALL snATAC-seq peaks
unfiltered_all_marker_overlap_nums <- c()
for(marker in mintchip_markers) {
  print(marker)
  all_marker_peaks <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_hg38.tsv"),
                               sep = "\t")
  colnames(all_marker_peaks) <- c("seqnames", "start", "end", "coordinates", "blah")
  all_marker_peaks_granges <- makeGRangesFromDataFrame(df = all_marker_peaks, keep.extra.columns = TRUE)
  all_marker_overlap <- as.data.frame(findOverlaps(all_marker_peaks_granges, sc_peaks_granges))
  percentage_peak_overlap <- nrow(all_marker_overlap) / nrow(all_marker_peaks)
  unfiltered_all_marker_overlap_nums <- c(unfiltered_all_marker_overlap_nums, paste0(marker, ",", nrow(all_marker_overlap),
                                                                                     ",", nrow(all_marker_peaks),
                                                                                     ",",
                                                                                     percentage_peak_overlap))
}

# Find overlap between DAS marker peaks and ALL snATAC-seq peaks
all_marker_overlap_df <- list()
for(marker in mintchip_markers) {
  print(marker)
  # Fix hg19 coordinates to be hg38
  marker_with_hg38_df <- add_hg38_coordinates_to_marker_peaks(read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0.1.tsv"),
                                                  sep = "\t", header = TRUE),
                                       read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_with_hg38_coordinates.tsv"),
                                                  sep = "\t", header = TRUE))
  # Find overlap between marker peaks and sc peaks
  marker_peaks_granges <- makeGRangesFromDataFrame(df = marker_with_hg38_df, keep.extra.columns = TRUE)
  marker_overlap <- as.data.frame(findOverlaps(marker_peaks_granges, sc_peaks_granges))
  marker_overlap_df <- combineRows(marker_with_hg38_df, sc_peaks, marker_overlap)
  marker_overlap_df$marker <- marker
  all_marker_overlap_df[[length(all_marker_overlap_df) + 1]] <- marker_overlap_df
}

# 1479 total marker peaks overlap with ALL snATAC-seq peaks
# H3K27Ac H3K27me3 H3K36me3  H3K4me1  H3K4me3  H3K9me3 
#   248       22       62      712      381       54 
all_marker_overlap_df <- do.call(rbind, all_marker_overlap_df)
  
# Next, we want to find overlap between each cell type and each marker
atac_cell_types_for_mintchip_analysis <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "MAIT", "NK", "Proliferating", "T Naive")
atac_mintchip_tables <- list()
for(atac_cell_type in atac_cell_types_for_mintchip_analysis) {
  print(atac_cell_type)
  atac_cell_type_for_file_name <- sub(" ", "_", atac_cell_type)
  atac_cell_type_peaks <- read.table(paste0(sc_das_dir, "diff_peaks/D28-vs-D_minus_1-degs-", atac_cell_type_for_file_name, 
                                            "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01.tsv"),
                                     sep = "\t", header = TRUE)
  components <- strsplit(atac_cell_type_peaks$Peak_Name, "[:-]")
  chr <- c()
  start <- c()
  end <- c()
  for(component in components) {
    chr <- c(chr, component[1])
    start <- c(start, component[2])
    end <- c(end, component[3])
  }
  atac_cell_type_peaks$seqnames <- chr
  atac_cell_type_peaks$start <- start
  atac_cell_type_peaks$end <- end
  atac_cell_type_peaks_granges <- makeGRangesFromDataFrame(df = atac_cell_type_peaks, keep.extra.columns = TRUE)
  cell_type_marker_overlap_df <- list()
  for(marker in mintchip_markers) {
    print(marker)
    marker_with_hg38_df <- add_hg38_coordinates_to_marker_peaks(read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0.1.tsv"),
                                                                           sep = "\t", header = TRUE),
                                                                read.table(paste0(mintchip_das_dir, marker, "/", marker, "_consensus_peak_set_FC_0_with_hg38_coordinates.tsv"),
                                                                           sep = "\t", header = TRUE))
    marker_peaks_granges <- makeGRangesFromDataFrame(df = marker_with_hg38_df, keep.extra.columns = TRUE)
    marker_overlap <- as.data.frame(findOverlaps(atac_cell_type_peaks_granges, marker_peaks_granges))
    if(nrow(marker_overlap) > 0) {
      marker_overlap_df <- combineRows(atac_cell_type_peaks, marker_with_hg38_df, marker_overlap)
      marker_overlap_df$marker <- marker
      cell_type_marker_overlap_df[[length(cell_type_marker_overlap_df) + 1]] <- marker_overlap_df
    }
  }
  cell_type_marker_overlap_df <- do.call(rbind, cell_type_marker_overlap_df)
  atac_mintchip_tables[[length(atac_mintchip_tables) + 1]] <- cell_type_marker_overlap_df
}

atac_mintchip_tables <- do.call(rbind, atac_mintchip_tables)

  
  
  
  
  

  
  # Overlap between MAGICAL circuits and DAS peaks that contained mintchip marker
  magical_results_cell_type_subset <- magical_results[magical_results$Cell_Type == atac_cell_type_for_file_name,]
  
  # Normal
  peaks <- cell_type_snATAC_mintchip_das_overlap$Peak_Name
  chromosomes <- sapply(strsplit(peaks, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 3))
  cell_type_snATAC_mintchip_das_overlap$Peak_chr <- chromosomes
  cell_type_snATAC_mintchip_das_overlap$Peak_start <- start_coords
  cell_type_snATAC_mintchip_das_overlap$Peak_end <- end_coords
  magical_results_cell_type_subset_overlap <- dplyr::inner_join(magical_results_cell_type_subset, cell_type_snATAC_mintchip_das_overlap, by = c("Peak_chr", "Peak_start", "Peak_end"))
  if(nrow(magical_results_cell_type_subset_overlap) > 0) {
    print(magical_results_cell_type_subset_overlap)
    #write.table(magical_results_cell_type_subset_overlap, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_MAGICAL.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}



