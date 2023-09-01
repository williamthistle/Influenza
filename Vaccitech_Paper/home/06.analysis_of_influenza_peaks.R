# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Subset sc_peaks to only have the lenient DAS
colnames(sc_peaks)[colnames(sc_peaks) == "value"] <- "chr"
sc_peaks_lenient_subset <- merge(sc_peaks, sc_das_lenient, by = c("chr", "idx"))[, names(sc_peaks)]
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

# Find associated motifs with peaks
# In motif file, chromosome is written like "10" versus "chr10"
current_chr <- substr(overlap_ranges$chr, 4, 6)
current_peak_start <- as.numeric(overlap_ranges$ATAC_peaks_start) + 250
current_peak_end <- as.numeric(overlap_ranges$ATAC_peaks_end) - 250
associated_motifs <- sc_motifs[sc_motifs$chr %in% current_chr & sc_motifs$point1 %in% current_peak_start & sc_motifs$point2 %in% current_peak_end,]
zero_columns <- colSums(associated_motifs) == 0
associated_motifs <- associated_motifs[, !zero_columns]

desired_sum <- 14
columns_with_desired_sum <- colSums(associated_motifs) == desired_sum
desired_columns <- names(associated_motifs)[columns_with_desired_sum]
print(desired_columns)


# OVERLAPPING GENES
mintchip_validated_genes <- sc_peaks_lenient_subset[sc_peaks_lenient_subset$nearestGene %in% mintchip_table$SYMBOL,]
# Run pos and neg genes through HumanBase - really good groupings
mintchip_validated_genes_pos <- mintchip_validated_genes[mintchip_validated_genes$sc_log2FC > 0,]
mintchip_validated_genes_neg <- mintchip_validated_genes[mintchip_validated_genes$sc_log2FC < 0,]
# What about HB within cell types?
# CD14 Mono
mintchip_validated_genes_cd14 <- mintchip_validated_genes[mintchip_validated_genes$Cell_Type == "CD14_Mono",]

# OVERLAPPING PEAKS





index <- 1
for(mint_cell_type_peak_annotation in mint_cell_type_peak_annotations) {
  current_cell_type <- cell_types[index]
  for(row_index in 1:nrow(mint_cell_type_peak_annotation)) {
    seqname <- mint_cell_type_peak_annotation[row_index,]$seqnames
    start <- mint_cell_type_peak_annotation[row_index,]$start
    end <- mint_cell_type_peak_annotation[row_index,]$end
    current_peak_table <- read.xlsx(paste0(single_cell_magical_dir, "scATAC_DAPs/", current_cell_type, "_D28_D1_diff.xlsx"), sheet = 1)
  }
  index <- index + 1
}












print(length(cell_type_peak_annotations[[1]]$SYMBOL))
print(length(intersect(cell_type_peak_annotations[[1]]$SYMBOL, mintchip_gene_table$SYMBOL)))

mint_cell_type_peak_annotations <- cell_type_peak_annotations[[1]][cell_type_peak_annotations[[1]]$SYMBOL %in% mintchip_gene_table$SYMBOL,]
mint_cell_type_genes <- unique(mint_cell_type_peak_annotations$SYMBOL)
cell_type_peak_annotations_subset <- cell_type_peak_annotations[[1]]
