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

# Finally, let's find overlap between the DMR in each cell type and each marker
dmr_cell_types_for_mintchip_analysis <- unique(snME_dms$celltype)
dmr_mintchip_tables <- list()
for(dmr_cell_type in dmr_cell_types_for_mintchip_analysis) {
  print(dmr_cell_type)
  snME_dms_subset <- snME_dms[snME_dms$celltype == dmr_cell_type,]
  snME_dms_subset$start <- snME_dms_subset$start - 100
  snME_dms_subset$end <- snME_dms_subset$end + 100
  snME_dms_subset_granges <- makeGRangesFromDataFrame(df = snME_dms_subset, keep.extra.columns = TRUE)
  cell_type_marker_overlap_df <- list()
  for(marker in mintchip_markers) {
    print(marker)
    das_marker_with_hg38_df <- add_hg38_coordinates_to_marker_peaks(read.table(paste0(mintchip_das_dir, marker, "/", marker, "_DESeq2_FC_0.tsv"),
                                                                               sep = "\t", header = TRUE),
                                                                    read.table(paste0(mintchip_das_dir, marker, "/", marker, "_all_peaks_with_hg38_coordinates.tsv"),
                                                                               sep = "\t", header = TRUE))
    das_marker_peaks_granges <- makeGRangesFromDataFrame(df = das_marker_with_hg38_df, keep.extra.columns = TRUE)
    marker_overlap <- as.data.frame(findOverlaps(snME_dms_subset_granges, das_marker_peaks_granges))
    if(nrow(marker_overlap) > 0) {
      marker_overlap_df <- combineRows(snME_dms_subset, das_marker_with_hg38_df, marker_overlap)
      marker_overlap_df$marker <- marker
      cell_type_marker_overlap_df[[length(cell_type_marker_overlap_df) + 1]] <- marker_overlap_df
    }
  }
  cell_type_marker_overlap_df <- do.call(rbind, cell_type_marker_overlap_df)
  colnames(cell_type_marker_overlap_df)[35] <- "seqnames_marker"
  colnames(cell_type_marker_overlap_df)[36] <- "start_marker"
  colnames(cell_type_marker_overlap_df)[37] <- "end_marker"
  dmr_mintchip_tables[[length(dmr_mintchip_tables) + 1]] <- cell_type_marker_overlap_df
}
dmr_mintchip_tables <- do.call(rbind, dmr_mintchip_tables)

# Summarize overlap for different cell types and different markers
print(table(dmr_mintchip_tables$marker))
for(cell_type in unique(dmr_mintchip_tables$celltype)) { 
  print(cell_type)
  cell_type_subset <- dmr_mintchip_tables[dmr_mintchip_tables$celltype == cell_type,]
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

