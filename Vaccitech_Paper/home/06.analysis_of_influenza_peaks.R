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
overlapping_peaks <- sc_peaks_lenient_subset
overlapping_peaks <- overlapping_peaks[0,]

for(current_chr in unique(sc_peaks_lenient_subset_extended$chr)) {
  sc_peaks_lenient_subset_extended_chr_subset <- sc_peaks_lenient_subset_extended[sc_peaks_lenient_subset_extended$chr == current_chr,]
  mintchip_table_extended_subset <- mintchip_table_extended[mintchip_table_extended$seqnames == current_chr,]
  
  starts_sc_peaks_lenient_subset_extended_chr_subset <- matrix(sc_peaks_lenient_subset_extended_chr_subset$start, nrow = nrow(sc_peaks_lenient_subset_extended_chr_subset), ncol = 1)
  ends_sc_peaks_lenient_subset_extended_chr_subset <- matrix(sc_peaks_lenient_subset_extended_chr_subset$end, nrow = nrow(sc_peaks_lenient_subset_extended_chr_subset), ncol = 1)
  
  starts_mintchip_table_extended_subset <- matrix(mintchip_table_extended_subset$start, nrow = nrow(mintchip_table_extended_subset), ncol = 1)
  ends_mintchip_table_extended_subset <- matrix(mintchip_table_extended_subset$end, nrow = nrow(mintchip_table_extended_subset), ncol = 1)
  
  # Create logical matrix for overlaps using matrix operations
  overlap_matrix <- outer(starts_sc_peaks_lenient_subset_extended_chr_subset, ends_mintchip_table_extended_subset, "<=") & outer(ends_sc_peaks_lenient_subset_extended_chr_subset, starts_mintchip_table_extended_subset, ">=")
  
  # Convert the matrix to a data frame for better readability
  overlap_df <- data.frame(sc_peaks_lenient_subset_extended_chr_subset_row = rep(1:nrow(sc_peaks_lenient_subset_extended_chr_subset), each = nrow(mintchip_table_extended_subset)),
                           mintchip_table_extended_subset_row = rep(1:nrow(mintchip_table_extended_subset), times = nrow(sc_peaks_lenient_subset_extended_chr_subset)),
                           overlap = as.logical(overlap_matrix))
  print(c)
  overlap_df <- overlap_df[overlap_df$overlap == TRUE,]
  print(overlap_df)
}




# OVERLAPPING GENES
mintchip_validated_genes <- sc_peaks_lenient_subset[sc_peaks_lenient_subset$nearestGene %in% mintchip_table$SYMBOL,]
# Run pos and neg genes through HumanBase - really good groupings
mintchip_validated_genes_pos <- mintchip_validated_genes[mintchip_validated_genes$sc_log2FC > 0,]
mintchip_validated_genes_neg <- mintchip_validated_genes[mintchip_validated_genes$sc_log2FC < 0,]
# What about HB within cell types?

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
