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

find_matching_cell_types_snME <- function(snATAC_cell_type) {
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
  return(matching_snME_cell_types)
}

# Important note - in addition to being different assays, the snME was run on 13 subjects (versus 4 for our ATAC peaks).
# So that will definitely affect the overlap we find.

# Find overlap between DMRs and ALL snATAC-seq peaks
dmr_overlap_nums <- data.frame(cell_type = character(), num_dmr = numeric(),
                                      num_atac_peaks = numeric(), num_of_dmr_overlapping = numeric(),
                                      percentage_of_dmr_overlapping = numeric(),
                                      num_of_atac_peaks_overlapping = numeric(),
                                      percentage_of_atac_peaks_overlapping = numeric())

for(cell_type in snME_dms_cell_types) {
  print(cell_type)
  current_snME_dms_data <- snME_dms[snME_dms$celltype == cell_type,]
  # Extend DMR for overlap
  current_snME_dms_data$start <- current_snME_dms_data$start - 100
  current_snME_dms_data$end <- current_snME_dms_data$end + 100
  # Find overlap between DMR and sc peaks
  dmr_granges <- makeGRangesFromDataFrame(df = current_snME_dms_data, keep.extra.columns = TRUE)
  dmr_overlap <- as.data.frame(findOverlaps(dmr_granges, sc_peaks_granges))
  dmr_percentage_peak_overlap <- length(unique(dmr_overlap$queryHits)) / nrow(current_snME_dms_data)
  atac_percentage_peak_overlap <- length(unique(dmr_overlap$subjectHits)) / nrow(sc_peaks)
  
  current_row <- c(cell_type, nrow(current_snME_dms_data), nrow(sc_peaks), length(unique(dmr_overlap$queryHits)),
                   dmr_percentage_peak_overlap, length(unique(dmr_overlap$subjectHits)), atac_percentage_peak_overlap)
  current_row <- data.frame(t(current_row))
  colnames(current_row) <- colnames(dmr_overlap_nums)
  dmr_overlap_nums <- rbind(dmr_overlap_nums, current_row)
}

dmr_overlap_nums

# The next question is: Are any of these peaks differentially accessible in the relevant cell type?
# If so, we have matching differential accessibility and differential methylation
atac_cell_types_for_snME_analysis <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "T Naive")
atac_snME_tables <- list()
for(atac_cell_type in atac_cell_types_for_snME_analysis) {
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
  cell_type_overlap_df <- list()
  matching_snME_cell_types <- find_matching_cell_types_snME(atac_cell_type)
  for(matching_cell_type in matching_snME_cell_types) {
    print(matching_cell_type)
    snME_dms_subset <- snME_dms[snME_dms$celltype == matching_cell_type,]
    snME_dms_subset$start <- snME_dms_subset$start - 100
    snME_dms_subset$end <- snME_dms_subset$end + 100
    snME_dms_subset_granges <- makeGRangesFromDataFrame(df = snME_dms_subset, keep.extra.columns = TRUE)
    overlap <- as.data.frame(findOverlaps(atac_cell_type_peaks_granges, snME_dms_subset_granges))
    if(nrow(overlap) > 0) {
      current_cell_type_overlap_df <- combineRows(atac_cell_type_peaks, snME_dms_subset, overlap)
      current_cell_type_overlap_df$dms_cell_type <- matching_cell_type
      cell_type_overlap_df[[length(cell_type_overlap_df) + 1]] <- current_cell_type_overlap_df
    }
  }
  cell_type_overlap_df <- do.call(rbind, cell_type_overlap_df)
  colnames(cell_type_overlap_df)[11] <- "seqnames_dms"
  colnames(cell_type_overlap_df)[12] <- "start_dms"
  colnames(cell_type_overlap_df)[13] <- "end_dms"
  atac_snME_tables[[length(atac_snME_tables) + 1]] <- cell_type_overlap_df
}
atac_snME_tables <- do.call(rbind, atac_snME_tables)

# Summarize overlap for different cell types and different markers
print(table(atac_snME_tables$Cell_Type))
for(cell_type in unique(atac_snME_tables$Cell_Type)) { 
  print(cell_type)
  cell_type_subset <- atac_snME_tables[atac_snME_tables$Cell_Type == cell_type,]
  print(table(cell_type_subset$methylation))
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

