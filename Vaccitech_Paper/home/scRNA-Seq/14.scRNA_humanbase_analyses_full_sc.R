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

# Innate
all_pos_innate_genes <- c()
all_neg_innate_genes <- c()

pos_innate_genes_by_cell_type <- list()
pos_innate_fmd_by_cell_type <- list()
neg_innate_genes_by_cell_type <- list()
neg_innate_fmd_by_cell_type <- list()

for(innate_cell_type in innate_cell_types) {
  innate_cell_type_for_file_name <- sub(" ", "_", innate_cell_type)
  current_deg_table <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", innate_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc.tsv"), sep = "\t", header = TRUE)
  pos_deg_table <- current_deg_table[current_deg_table$avg_log2FC > 0,]
  neg_deg_table <- current_deg_table[current_deg_table$avg_log2FC < 0,]
  all_pos_innate_genes <- c(all_pos_innate_genes, rownames(pos_deg_table))
  all_neg_innate_genes <- c(all_neg_innate_genes, rownames(neg_deg_table))
  pos_innate_genes_by_cell_type[[innate_cell_type]] <- rownames(pos_deg_table)
  pos_innate_fmd <- run_fmd_on_scRNA(rownames(pos_deg_table))
  pos_innate_fmd_by_cell_type[[innate_cell_type]] <- pos_innate_fmd
  neg_innate_genes_by_cell_type[[innate_cell_type]] <- rownames(neg_deg_table)
  neg_innate_fmd <- run_fmd_on_scRNA(rownames(neg_deg_table))
  neg_innate_fmd_by_cell_type[[innate_cell_type]] <- neg_innate_fmd
}

all_pos_innate_genes <- unique(all_pos_innate_genes)
all_neg_innate_genes <- unique(all_neg_innate_genes)

all_pos_innate_fmd <- run_fmd_on_scRNA(all_pos_innate_genes)
all_neg_innate_fmd <- run_fmd_on_scRNA(all_neg_innate_genes)

# Compare counts from pre-pseudobulk corrected to post-pseudobulk corrected for each cell type
for(innate_cell_type in innate_cell_types) {
  print(paste0("Number of uncorrected positive genes for ", innate_cell_type))
  print(length(pos_innate_genes_by_cell_type[[innate_cell_type]]))
  print(paste0("Number of corrected positive genes for ", innate_cell_type))
  pos_corrected_genes <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == innate_cell_type,]
  pos_corrected_genes <- nrow(pos_corrected_genes[pos_corrected_genes$sc_log2FC > 0,])
  print(pos_corrected_genes)
  
  print(paste0("Number of uncorrected negative genes for ", innate_cell_type))
  print(length(neg_innate_genes_by_cell_type[[innate_cell_type]]))
  print(paste0("Number of corrected negative genes for ", innate_cell_type))
  neg_corrected_genes <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == innate_cell_type,]
  neg_corrected_genes <- nrow(neg_corrected_genes[neg_corrected_genes$sc_log2FC < 0,])
  print(neg_corrected_genes)
}



# Adaptive
all_pos_adaptive_genes <- c()
all_neg_adaptive_genes <- c()

pos_adaptive_genes_by_cell_type <- list()
pos_adaptive_fmd_by_cell_type <- list()
neg_adaptive_genes_by_cell_type <- list()
neg_adaptive_fmd_by_cell_type <- list()

adaptive_cell_types_for_hb <- adaptive_cell_types[-c(8,9)]

for(adaptive_cell_type in adaptive_cell_types_for_hb) {
  adaptive_cell_type_for_file_name <- sub(" ", "_", adaptive_cell_type)
  current_deg_table <- read.table(paste0(sc_deg_dir, "D28-vs-D_minus_1-degs-", adaptive_cell_type_for_file_name, "-time_point-controlling_for_subject_id_sc.tsv"), sep = "\t", header = TRUE)
  pos_deg_table <- current_deg_table[current_deg_table$avg_log2FC > 0,]
  neg_deg_table <- current_deg_table[current_deg_table$avg_log2FC < 0,]
  all_pos_adaptive_genes <- c(all_pos_adaptive_genes, rownames(pos_deg_table))
  all_neg_adaptive_genes <- c(all_neg_adaptive_genes, rownames(neg_deg_table))
  pos_adaptive_genes_by_cell_type[[adaptive_cell_type]] <- rownames(pos_deg_table)
  pos_adaptive_fmd <- run_fmd_on_scRNA(rownames(pos_deg_table))
  pos_adaptive_fmd_by_cell_type[[adaptive_cell_type]] <- pos_adaptive_fmd
  neg_adaptive_genes_by_cell_type[[adaptive_cell_type]] <- rownames(neg_deg_table)
  neg_adaptive_fmd <- run_fmd_on_scRNA(rownames(neg_deg_table))
  neg_adaptive_fmd_by_cell_type[[adaptive_cell_type]] <- neg_adaptive_fmd
}

all_pos_adaptive_genes <- unique(all_pos_adaptive_genes)
all_neg_adaptive_genes <- unique(all_neg_adaptive_genes)

all_pos_adaptive_fmd <- run_fmd_on_scRNA(all_pos_adaptive_genes)
all_neg_adaptive_fmd <- run_fmd_on_scRNA(all_neg_adaptive_genes)




