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

find_matching_snATAC_das_and_snME_dms <- function(snATAC_cell_type, sc_peak_file_path, snME_dms, lenient = 0) {
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
  
  # Grab DAS
  peak_table <- read.table(sc_peak_file_path, sep = "\t", header = TRUE)
  
  # Grab coordinate information for each positive and negative peak so we can find overlap between DAS and DMS
  peaks <- peak_table$Peak_Name
  chromosomes <- sapply(strsplit(peaks, "-"), `[`, 1)
  start_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 2))
  end_coords <- as.numeric(sapply(strsplit(peaks, "-"), `[`, 3))
  peak_table$seqnames <- chromosomes
  peak_table$start <- start_coords
  peak_table$end <- end_coords
  
  peak_table$start <- peak_table$start - lenient
  peak_table$end <- peak_table$end + lenient
  
  peak_table_granges <- makeGRangesFromDataFrame(df = peak_table, keep.extra.columns = TRUE)
  
  final_df <- list()
  
  for(snME_cell_type in matching_snME_cell_types) {
    snME_dms_subset <- snME_dms[snME_dms$celltype == snME_cell_type,]
    snME_dms_subset_granges <- makeGRangesFromDataFrame(df = snME_dms_subset, keep.extra.columns = TRUE)
    overlap <- as.data.frame(findOverlaps(peak_table_granges, snME_dms_subset_granges))
    overlap_df <- combineRows(peak_table, snME_dms_subset,overlap)
    overlap_df <- overlap_df[,-c(27:52)]
    final_df[[length(final_df) + 1]] <- overlap_df
  }
  
  final_df <- do.call(rbind, final_df)
  return(final_df)
}

# Find all overlap between snATAC-seq (hg38) and snMethylation (hg38)

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

colnames(snME_dms)[2] <- "seqnames"
snME_dms_granges <- makeGRangesFromDataFrame(df = snME_dms, keep.extra.columns = TRUE)

# 5498
sc_peak_and_snME_dms_overlap <- as.data.frame(findOverlaps(sc_peaks_granges, snME_dms_granges))
sc_peak_and_snME_dms_overlap_df <- combineRows(sc_peaks, snME_dms, sc_peak_and_snME_dms_overlap)
sc_peak_and_snME_dms_overlap_df <- sc_peak_and_snME_dms_overlap_df[,-c(27:52)]

# 10424
sc_peak_lenient_and_snME_dms_overlap <- as.data.frame(findOverlaps(sc_peaks_lenient_granges, snME_dms_granges))
sc_peak_lenient_and_snME_dms_overlap_df <- combineRows(sc_peaks_lenient, snME_dms, sc_peak_lenient_and_snME_dms_overlap)
sc_peak_lenient_and_snME_dms_overlap_df <- sc_peak_lenient_and_snME_dms_overlap_df[,-c(27:52)]

# 14622
sc_peak_more_lenient_and_snME_dms_overlap <- as.data.frame(findOverlaps(sc_peaks_more_lenient_granges, snME_dms_granges))
sc_peak_more_lenient_and_snME_dms_overlap_df <- combineRows(sc_peaks_more_lenient, snME_dms, sc_peak_more_lenient_and_snME_dms_overlap)
sc_peak_more_lenient_and_snME_dms_overlap_df <- sc_peak_more_lenient_and_snME_dms_overlap_df[,-c(27:52)]

write.table(sc_peak_and_snME_dms_overlap_df, file = paste0(snME_results_dir, "snME_dms_overlap_with_all_atac_peaks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sc_peak_lenient_and_snME_dms_overlap_df, file = paste0(snME_results_dir, "snME_dms_overlap_with_all_atac_peaks_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sc_peak_more_lenient_and_snME_dms_overlap_df, file = paste0(snME_results_dir, "snME_dms_overlap_with_all_atac_peaks_more_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# The next question is: Are any of these peaks differentially accessible in the relevant cell type?
# If so, we have matching differential accessibility and differential methylation
atac_cell_types_for_snME_analysis <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "T Naive")

magical_output_dir <- paste0(sc_magical_dir, "Output/")
magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

for(cell_type in atac_cell_types_for_snME_analysis) {
  cell_type_for_file_name <- sub(" ", "_", cell_type)
  cell_type_snATAC_snME_overlap_das_subset <- find_matching_snATAC_das_and_snME_dms(cell_type, paste0(sc_das_dir, 
                                                                                                                     "diff_peaks/D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01.tsv"),
                                                                                    snME_dms, lenient = 0)
  cell_type_snATAC_snME_overlap_das_subset_lenient <- find_matching_snATAC_das_and_snME_dms(cell_type, paste0(sc_das_dir, 
                                                                                                              "diff_peaks/D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01.tsv"),
                                                                                            snME_dms, lenient = 250)
  cell_type_snATAC_snME_overlap_das_subset_more_lenient <- find_matching_snATAC_das_and_snME_dms(cell_type, paste0(sc_das_dir, 
                                                                                                              "diff_peaks/D28-vs-D_minus_1-degs-", cell_type_for_file_name, "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01.tsv"),
                                                                                            snME_dms, lenient = 500)
  
  write.table(cell_type_snATAC_snME_overlap_das_subset, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cell_type_snATAC_snME_overlap_das_subset_lenient, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cell_type_snATAC_snME_overlap_das_subset_more_lenient, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_more_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  sc_peaks_granges <- makeGRangesFromDataFrame(df = cell_type_snATAC_snME_overlap_das_subset, keep.extra.columns = TRUE)
  sc_peaks_lenient_granges <- makeGRangesFromDataFrame(df = cell_type_snATAC_snME_overlap_das_subset_lenient, keep.extra.columns = TRUE)
  sc_peaks_more_lenient_granges <- makeGRangesFromDataFrame(df = cell_type_snATAC_snME_overlap_das_subset_more_lenient, keep.extra.columns = TRUE)
  
  # Overlap between MAGICAL circuits and DAS peaks that contained DMS
  magical_results_cell_type_subset <- magical_results[magical_results$Cell_Type == cell_type_for_file_name,]
  colnames(magical_results_cell_type_subset)[5] <- "seqnames"
  colnames(magical_results_cell_type_subset)[6] <- "start"
  colnames(magical_results_cell_type_subset)[7] <- "end"
  
  magical_granges <- makeGRangesFromDataFrame(df = magical_results_cell_type_subset, keep.extra.columns = TRUE)
  
  # Normal
  magical_overlap <- as.data.frame(findOverlaps(sc_peaks_granges, magical_granges))
  magical_lenient_overlap <- as.data.frame(findOverlaps(sc_peaks_lenient_granges, magical_granges))
  magical_more_lenient_overlap <- as.data.frame(findOverlaps(sc_peaks_more_lenient_granges, magical_granges))
  
  if(nrow(magical_overlap) > 0) {
    magical_overlap_df <- combineRows(cell_type_snATAC_snME_overlap_das_subset, magical_results_cell_type_subset, magical_overlap)
    write.table(magical_overlap_df, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_MAGICAL.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  if(nrow(magical_lenient_overlap) > 0) {
    magical_overlap_lenient_df <- combineRows(cell_type_snATAC_snME_overlap_das_subset_lenient, magical_results_cell_type_subset, magical_lenient_overlap)
    write.table(magical_overlap_lenient_df, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_MAGICAL_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  if(nrow(magical_more_lenient_overlap) > 0) {
    magical_overlap_more_lenient_df <- combineRows(cell_type_snATAC_snME_overlap_das_subset_more_lenient, magical_results_cell_type_subset, magical_more_lenient_overlap)
    write.table(magical_more_lenient_overlap, file = paste0(snME_results_dir, "Cell_Type_Overlap/", cell_type_for_file_name, "_das_snME_dms_overlap_MAGICAL_more_lenient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
