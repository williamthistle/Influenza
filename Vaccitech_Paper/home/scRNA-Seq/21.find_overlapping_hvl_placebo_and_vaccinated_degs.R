overlapping_cell_types <- intersect(unique(scRNA_hvl_placebo_degs$Cell_Type), unique(scRNA_hvl_vaccinated_degs$Cell_Type))

for(cell_type in overlapping_cell_types) {
  print(cell_type)
  hvl_placebo_subset <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == cell_type,]
  hvl_vaccinated_subset <- scRNA_hvl_vaccinated_degs[scRNA_hvl_vaccinated_degs$Cell_Type == cell_type,]
  intersecting_genes <- intersect(hvl_placebo_subset$Gene_Name, hvl_vaccinated_subset$Gene_Name)
  print(paste0("Total DEGs overlapping: ", length(intersecting_genes)))
  print(paste0("Total Placebo DEGs: ", nrow(hvl_placebo_subset)))
  print(paste0("Percentage: ", length(intersecting_genes) / nrow(hvl_placebo_subset) * 100))
  pos_hvl_placebo_subset <- hvl_placebo_subset[hvl_placebo_subset$sc_log2FC > 0,]
  pos_hvl_vaccinated_subset <- hvl_vaccinated_subset[hvl_vaccinated_subset$sc_log2FC > 0,]
  neg_hvl_placebo_subset <- hvl_placebo_subset[hvl_placebo_subset$sc_log2FC < 0,]
  neg_hvl_vaccinated_subset <- hvl_vaccinated_subset[hvl_vaccinated_subset$sc_log2FC < 0,]
  overlapping_pos_genes <- intersect(pos_hvl_placebo_subset$Gene_Name, pos_hvl_vaccinated_subset$Gene_Name)
  overlapping_neg_genes <- intersect(neg_hvl_placebo_subset$Gene_Name, neg_hvl_vaccinated_subset$Gene_Name)
  print(overlapping_pos_genes)
  print(overlapping_neg_genes)
}
