# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))



for(cell_type in overlapping_cell_types) {
  print(cell_type)
  hvl_subset <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == cell_type,]
  lvl_subset <- lvl_sc_pseudobulk_deg_table[lvl_sc_pseudobulk_deg_table$Cell_Type == cell_type,]
  pos_hvl_subset <- hvl_subset[hvl_subset$sc_log2FC > 0,]
  pos_lvl_subset <- lvl_subset[lvl_subset$sc_log2FC > 0,]
  neg_hvl_subset <- hvl_subset[hvl_subset$sc_log2FC < 0,]
  neg_lvl_subset <- lvl_subset[lvl_subset$sc_log2FC < 0,]
  overlapping_pos_genes <- intersect(pos_hvl_subset$Gene_Name, pos_lvl_subset$Gene_Name)
  overlapping_neg_genes <- intersect(neg_hvl_subset$Gene_Name, neg_lvl_subset$Gene_Name)
  print(overlapping_pos_genes)
  print(overlapping_neg_genes)
}
