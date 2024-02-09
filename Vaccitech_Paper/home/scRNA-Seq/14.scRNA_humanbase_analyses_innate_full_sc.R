### THIS FUNCTION:
### 1) ASSIGNS CELL TYPES TO HUMANBASE MODULES (UP AND DOWNREGULATED SC DEGS)
### 2) ASSIGNS CELL TYPES TO REACTOME PATHWAYS (MODULES FROM HB ANALYSIS)
### 3) CREATES REACTOME PATHWAY PLOTS FOR FIGURE (NOTE: DOESN'T INCLUDE CELL TYPES YET)

run_fmd_on_scRNA <- function(gene_list) {
  if(length(gene_list) > 1 & length(gene_list) < 2000) {
    current_fmd_result <- SPEEDI::RunFMD_RNA(gene_list, "blood")
  } else { 
    current_fmd_result <- "EMPTY OR OVER 2000 GENES (TOO MANY)"
  }
  return(current_fmd_result)
}

# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

pos_innate_genes <- c()
neg_innate_genes <- c()

pos_innate_genes_by_cell_type <- list()
neg_innate_genes_by_cell_type <- list()

for(innate_cell_type in innate_cell_types) {
  innate_cell_type_for_file_name <- sub(" ", "_", innate_cell_type)
  current_deg_table <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", innate_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc.tsv"), sep = "\t", header = TRUE)
  pos_deg_table <- current_deg_table[current_deg_table$avg_log2FC > 0,]
  neg_deg_table <- current_deg_table[current_deg_table$avg_log2FC < 0,]
  pos_innate_genes <- c(pos_innate_genes, rownames(pos_deg_table))
  neg_innate_genes <- c(neg_innate_genes, rownames(neg_deg_table))
  pos_innate_genes_by_cell_type[[innate_cell_type]] <- pos_innate_genes
  neg_innate_genes_by_cell_type[[innate_cell_type]] <- neg_innate_genes
}

pos_innate_genes <- unique(pos_innate_genes)
neg_innate_genes <- unique(neg_innate_genes)

pos_innate_fmd <- run_fmd_on_scRNA(pos_innate_genes)
neg_innate_fmd <- run_fmd_on_scRNA(neg_innate_genes)


pos_adaptive_genes <- c()
neg_adaptive_genes <- c()

for(adaptive_cell_type in adaptive_cell_types) {
  adaptive_cell_type_for_file_name <- sub(" ", "_", adaptive_cell_type)
  current_deg_table <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", adaptive_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc.tsv"), sep = "\t", header = TRUE)
  pos_deg_table <- current_deg_table[current_deg_table$avg_log2FC > 0,]
  neg_deg_table <- current_deg_table[current_deg_table$avg_log2FC < 0,]
  pos_adaptive_genes <- c(pos_adaptive_genes, rownames(pos_deg_table))
  neg_adaptive_genes <- c(neg_adaptive_genes, rownames(neg_deg_table))
}

pos_adaptive_genes <- unique(pos_adaptive_genes)
neg_adaptive_genes <- unique(neg_adaptive_genes)

pos_adaptive_fmd <- run_fmd_on_scRNA(pos_adaptive_genes)
neg_adaptive_fmd <- run_fmd_on_scRNA(neg_adaptive_genes)



