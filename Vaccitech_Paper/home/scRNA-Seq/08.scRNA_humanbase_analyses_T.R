### THIS FUNCTION:
### 1) ASSIGNS CELL TYPES TO HUMANBASE MODULES (UP AND DOWNREGULATED SC DEGS)
### 2) ASSIGNS CELL TYPES TO REACTOME PATHWAYS (MODULES FROM HB ANALYSIS)
### 3) CREATES REACTOME PATHWAY PLOTS FOR FIGURE (NOTE: DOESN'T INCLUDE CELL TYPES YET)

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

t_naive_dir <- paste0(sc_humanbase_dir, "T_Naive/")

# Grab upregulated and downregulated cell types for adaptive (trained immunity) cell types
sc_t_naive <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == "CD4 Naive" | sc_pseudobulk_deg_table$Cell_Type == "CD8 Naive",]
sc_t_naive_upregulated <- sc_t_naive[sc_t_naive$sc_log2FC > 0,]
sc_t_naive_downregulated <- sc_t_naive[sc_t_naive$sc_log2FC < 0,]

# Assign cell types to downregulated genes from scRNA-seq (t_naive network)
HB_t_naive_downregulated_genes_t_naive <- assign_cell_types_to_humanbase_results(paste0(t_naive_dir, "HB_T_Naive_downregulated_genes_T-lymphocyte.tsv"), sc_t_naive_downregulated)
HB_t_naive_downregulated_genes_t_naive_modules <- HB_t_naive_downregulated_genes_t_naive[[1]]
HB_t_naive_downregulated_genes_t_naive_go_terms <- HB_t_naive_downregulated_genes_t_naive[[2]]

output_file <- paste0(t_naive_dir, "HB_t_naive_downregulated_genes_t_naive_with_cell_types.tsv")
cat(paste0("# GENES ASSOCIATED WITH EACH MODULE \n"), file=output_file)
utils::write.table(HB_t_naive_downregulated_genes_t_naive_modules, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("# GENES ASSOCIATED WITH EACH GO TERM WITHIN EACH MODULE \n"), file=output_file, append = TRUE)
utils::write.table(HB_t_naive_downregulated_genes_t_naive_go_terms, file = output_file, append=TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
