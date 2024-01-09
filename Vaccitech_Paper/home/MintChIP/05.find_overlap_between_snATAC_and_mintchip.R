# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))
magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

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

# Fix hg19 coordinates to be hg38 for a subset of marker peaks
add_hg38_coordinates_to_marker_peaks <- function(marker_table, hg38_table) {
  # create GRanges objects for differentially accessible marker peaks and ALL marker peaks (with hg38 start/end attached as extra columns)
  marker_peaks_granges <- makeGRangesFromDataFrame(df = marker_table, keep.extra.columns = TRUE)
  marker_peaks_with_hg38_coordinates_granges <- makeGRangesFromDataFrame(df = hg38_table, keep.extra.columns = TRUE)
  # Find overlap based on the hg19 coordinates
  marker_overlap <- as.data.frame(findOverlaps(marker_peaks_granges, marker_peaks_with_hg38_coordinates_granges))
  marker_overlap_df <- combineRows(marker_table, hg38_table, marker_overlap)
  # Now we can replace the hg19 coordinates with hg38 coordinates
  marker_overlap_df$start <- marker_overlap_df$hg38_start
  marker_overlap_df$end <- marker_overlap_df$hg38_end
  marker_overlap_df <- marker_overlap_df[,-c(12:17)]
  return(marker_overlap_df)
}

# Important note - in addition to being different assays, the mintchip was run on 11 subjects (versus 4 for our ATAC peaks).
# So that will definitely affect the overlap we find.

# Find overlap between ALL marker peaks and ALL snATAC-seq peaks
# Note that numbers may be slightly off because we don't do the combineRows step, but they are mostly accurate.
unfiltered_all_marker_overlap_nums <- data.frame(marker = character(), num_marker_peaks = numeric(),
                                                 num_atac_peaks = numeric(), num_of_marker_peaks_overlapping = numeric(),
                                                 percentage_of_marker_peaks_overlapping = numeric(),
                                                 num_of_atac_peaks_overlapping = numeric(),
                                                 percentage_of_atac_peaks_overlapping = numeric())
for(marker in mintchip_markers) {
  print(marker)
  all_marker_peaks <- read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_with_hg38_coordinates.tsv"),
                               sep = "\t", header = TRUE)
  all_marker_peaks$start <- all_marker_peaks$hg38_start
  all_marker_peaks$end <- all_marker_peaks$hg38_end
  all_marker_peaks_granges <- makeGRangesFromDataFrame(df = all_marker_peaks, keep.extra.columns = TRUE)
  all_marker_overlap <- as.data.frame(findOverlaps(all_marker_peaks_granges, sc_peaks_granges))
  marker_percentage_peak_overlap <- length(unique(all_marker_overlap$queryHits)) / nrow(all_marker_peaks)
  atac_percentage_peak_overlap <- length(unique(all_marker_overlap$subjectHits)) / nrow(sc_peaks)

  current_row <- c(marker, nrow(all_marker_peaks), nrow(sc_peaks), length(unique(all_marker_overlap$queryHits)),
                   marker_percentage_peak_overlap, length(unique(all_marker_overlap$subjectHits)), atac_percentage_peak_overlap)
  current_row <- data.frame(t(current_row))
  colnames(current_row) <- colnames(unfiltered_all_marker_overlap_nums)
  unfiltered_all_marker_overlap_nums <- rbind(unfiltered_all_marker_overlap_nums, current_row)
}

# The overlap is pretty good for all markers!
unfiltered_all_marker_overlap_nums

# Find overlap between DAS marker peaks and ALL snATAC-seq peaks
# Let's do DESeq2 with FC threshold 0 since that's the best case scenario I'd say
das_marker_overlap_nums <- data.frame(marker = character(), num_marker_peaks = numeric(),
                                      num_atac_peaks = numeric(), num_of_marker_peaks_overlapping = numeric(),
                                      percentage_of_marker_peaks_overlapping = numeric(),
                                      num_of_atac_peaks_overlapping = numeric(),
                                      percentage_of_atac_peaks_overlapping = numeric())

for(marker in mintchip_markers) {
  print(marker)
  # Fix hg19 coordinates to be hg38
  das_marker_with_hg38_df <- add_hg38_coordinates_to_marker_peaks(read.table(paste0(mintchip_das_dir, marker, "/", marker, "_DESeq2_FC_0.tsv"),
                                                  sep = "\t", header = TRUE),
                                       read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_with_hg38_coordinates.tsv"),
                                                  sep = "\t", header = TRUE))
  # Find overlap between marker peaks and sc peaks
  das_marker_peaks_granges <- makeGRangesFromDataFrame(df = das_marker_with_hg38_df, keep.extra.columns = TRUE)
  das_marker_overlap <- as.data.frame(findOverlaps(das_marker_peaks_granges, sc_peaks_granges))
  marker_percentage_peak_overlap <- length(unique(das_marker_overlap$queryHits)) / nrow(das_marker_with_hg38_df)
  atac_percentage_peak_overlap <- length(unique(das_marker_overlap$subjectHits)) / nrow(sc_peaks)
  
  current_row <- c(marker, nrow(das_marker_with_hg38_df), nrow(sc_peaks), length(unique(das_marker_overlap$queryHits)),
                   marker_percentage_peak_overlap, length(unique(das_marker_overlap$subjectHits)), atac_percentage_peak_overlap)
  current_row <- data.frame(t(current_row))
  colnames(current_row) <- colnames(das_marker_overlap_nums)
  das_marker_overlap_nums <- rbind(das_marker_overlap_nums, current_row)
}

# Not too bad!
das_marker_overlap_nums

# Finally, let's find overlap between the DAS in each cell type and each marker

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
    das_marker_with_hg38_df <- add_hg38_coordinates_to_marker_peaks(read.table(paste0(mintchip_das_dir, marker, "/", marker, "_DESeq2_FC_0.tsv"),
                                                                               sep = "\t", header = TRUE),
                                                                    read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_with_hg38_coordinates.tsv"),
                                                                               sep = "\t", header = TRUE))
    das_marker_peaks_granges <- makeGRangesFromDataFrame(df = das_marker_with_hg38_df, keep.extra.columns = TRUE)
    marker_overlap <- as.data.frame(findOverlaps(atac_cell_type_peaks_granges, das_marker_peaks_granges))
    if(nrow(marker_overlap) > 0) {
      marker_overlap_df <- combineRows(atac_cell_type_peaks, das_marker_with_hg38_df, marker_overlap)
      marker_overlap_df$marker <- marker
      cell_type_marker_overlap_df[[length(cell_type_marker_overlap_df) + 1]] <- marker_overlap_df
    }
  }
  cell_type_marker_overlap_df <- do.call(rbind, cell_type_marker_overlap_df)
  colnames(cell_type_marker_overlap_df)[10] <- "seqnames_marker"
  colnames(cell_type_marker_overlap_df)[11] <- "start_marker"
  colnames(cell_type_marker_overlap_df)[12] <- "end_marker"
  atac_mintchip_tables[[length(atac_mintchip_tables) + 1]] <- cell_type_marker_overlap_df
}
atac_mintchip_tables <- do.call(rbind, atac_mintchip_tables)

# Summarize overlap for different cell types and different markers
print(table(atac_mintchip_tables$marker))
for(cell_type in unique(atac_mintchip_tables$Cell_Type)) { 
  print(cell_type)
  cell_type_subset <- atac_mintchip_tables[atac_mintchip_tables$Cell_Type == cell_type,]
  print(table(cell_type_subset$marker))
}

# Find MAGICAL overlap
MAGICAL_overlap_with_mintchip <- list()
for(cell_type in unique(atac_mintchip_tables$Cell_Type)) {
  print(cell_type)
  cell_type_for_magical <- sub(" ", "_", cell_type)
  atac_mintchip_tables_cell_type_subset <- atac_mintchip_tables[atac_mintchip_tables$Cell_Type == cell_type,]
  magical_results_cell_type_subset <- magical_results[magical_results$Cell_Type == cell_type_for_magical,]
  magical_results_cell_type_subset$seqnames <- magical_results_cell_type_subset$Peak_chr
  magical_results_cell_type_subset$start <- magical_results_cell_type_subset$Peak_start
  magical_results_cell_type_subset$end <- magical_results_cell_type_subset$Peak_end
  result <- merge(magical_results_cell_type_subset, atac_mintchip_tables_cell_type_subset, by = c("seqnames", "start", "end"))
  MAGICAL_overlap_with_mintchip[[length(MAGICAL_overlap_with_mintchip) + 1]] <- result
}
MAGICAL_overlap_with_mintchip <- do.call(rbind, MAGICAL_overlap_with_mintchip)

