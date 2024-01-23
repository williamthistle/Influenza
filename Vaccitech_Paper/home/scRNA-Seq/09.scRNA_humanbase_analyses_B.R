### THIS FUNCTION:
### 1) ASSIGNS CELL TYPES TO HUMANBASE MODULES (UP AND DOWNREGULATED SC DEGS)
### 2) ASSIGNS CELL TYPES TO REACTOME PATHWAYS (MODULES FROM HB ANALYSIS)
### 3) CREATES REACTOME PATHWAY PLOTS FOR FIGURE (NOTE: DOESN'T INCLUDE CELL TYPES YET)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

B_dir <- paste0(sc_humanbase_dir, "B/")

# Grab upregulated and downregulated cell types for adaptive (trained immunity) cell types
sc_Bs <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == "B naive" | sc_pseudobulk_deg_table$Cell_Type == "B memory",]
sc_Bs_upregulated <- sc_Bs[sc_Bs$sc_log2FC > 0,]
sc_Bs_downregulated <- sc_Bs[sc_Bs$sc_log2FC < 0,]

# Assign cell types to upregulated genes from scRNA-seq (blood network)
HB_B_upregulated_genes_B <- assign_cell_types_to_humanbase_results(paste0(B_dir, "HB_B_upregulated_genes_B-lymphocyte.tsv"), sc_Bs_upregulated)
HB_B_upregulated_genes_B_modules <- HB_B_upregulated_genes_B[[1]]
HB_B_upregulated_genes_B_go_terms <- HB_B_upregulated_genes_B[[2]]

output_file <- paste0(B_dir, "HB_B_upregulated_genes_B_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_B_upregulated_genes_B_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_B_upregulated_genes_B_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Assign cell types to downregulated genes from scRNA-seq (B network)
HB_B_downregulated_genes_B <- assign_cell_types_to_humanbase_results(paste0(B_dir, "HB_B_downregulated_genes_B-lymphocyte.tsv"), sc_Bs_downregulated)
HB_B_downregulated_genes_B_modules <- HB_B_downregulated_genes_B[[1]]
HB_B_downregulated_genes_B_go_terms <- HB_B_downregulated_genes_B[[2]]

output_file <- paste0(B_dir, "HB_B_downregulated_genes_B_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_B_downregulated_genes_B_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_B_downregulated_genes_B_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
