overlapping_cell_types <- intersect(unique(scRNA_hvl_placebo_degs$Cell_Type), unique(scRNA_lvl_placebo_degs$Cell_Type))

for(cell_type in overlapping_cell_types) {
  print(cell_type)
  hvl_subset <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type == cell_type,]
  lvl_subset <- scRNA_lvl_placebo_degs[scRNA_lvl_placebo_degs$Cell_Type == cell_type,]
  pos_hvl_subset <- hvl_subset[hvl_subset$sc_log2FC > 0,]
  pos_lvl_subset <- lvl_subset[lvl_subset$sc_log2FC > 0,]
  neg_hvl_subset <- hvl_subset[hvl_subset$sc_log2FC < 0,]
  neg_lvl_subset <- lvl_subset[lvl_subset$sc_log2FC < 0,]
  overlapping_pos_genes <- intersect(pos_hvl_subset$Gene_Name, pos_lvl_subset$Gene_Name)
  overlapping_neg_genes <- intersect(neg_hvl_subset$Gene_Name, neg_lvl_subset$Gene_Name)
  print(overlapping_pos_genes)
  print(overlapping_neg_genes)
}
