# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

# Bulk RNA-seq will be the same regardless of cell type, so we can just grab these once
# Grab bulk RNA-seq
bulk_degs <- high_placebo_period_2_D28_vs_D_minus_1_results[[1]]
pos_bulk_degs <- bulk_degs[bulk_degs$log2FoldChange > 0,]
pos_bulk_degs <- rownames(pos_bulk_degs)
neg_bulk_degs <- bulk_degs[bulk_degs$log2FoldChange < 0,]
neg_bulk_degs <- rownames(neg_bulk_degs)
# Grab mintchip
pos_mintchip_das <- list()
neg_mintchip_das <- list()
for(mintchip_marker in mintchip_markers) {
  current_pos_mintchip_das <- read.table(paste0(mintchip_das_dir, mintchip_marker, "/upregulated/", mintchip_marker, 
                                             "_DESeq2_FC_0.1_upregulated_annotated.tsv"),
                                      sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  current_neg_mintchip_das <- read.table(paste0(mintchip_das_dir, mintchip_marker, "/downregulated/", mintchip_marker, 
                                                "_DESeq2_FC_0.1_downregulated_annotated.tsv"),
                                         sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  pos_mintchip_das[[mintchip_marker]] <- unique(current_pos_mintchip_das$SYMBOL)
  neg_mintchip_das[[mintchip_marker]] <- unique(current_neg_mintchip_das$SYMBOL)
}

# All comparisons are for closest gene in the context of ATAC / snME / Mintchip

# {ATAC, snME}
# {ATAC, Mintchip}
# {ATAC, scRNA-seq}
# {ATAC, bulk RNA-seq}
# {ATAC, snME, Mintchip}
# {ATAC, snME, scRNA-seq}
# {ATAC, snME, bulk RNA-seq}
# {ATAC, Mintchip, scRNA-seq}
# {ATAC, Mintchip, bulk RNA-seq}
# {ATAC, scRNA-seq, bulk RNA-seq}
# {ATAC, snME, Mintchip, scRNA-seq}
# {ATAC, snME, Mintchip, bulk RNA-seq}
# {ATAC, snME, scRNA-seq, bulk RNA-seq}
# {ATAC, Mintchip, scRNA-seq, bulk RNA-seq}
# {ATAC, snME, Mintchip, scRNA-seq, bulk RNA-seq}
# {snME, Mintchip}
# {snME, scRNA-seq}
# {snME, bulk RNA-seq}
# {snME, Mintchip, scRNA-seq}
# {snME, Mintchip, bulk RNA-seq}
# {snME, scRNA-seq, bulk RNA-seq}
# {snME, Mintchip, scRNA-seq, bulk RNA-seq}
# {Mintchip, scRNA-seq}
# {Mintchip, bulk RNA-seq}
# {Mintchip, scRNA-seq, bulk RNA-seq}
# {scRNA-seq, bulk RNA-seq} - already done before

ATAC_and_snME_gene_overlap <- list()
ATAC_and_mintchip_gene_overlap <- list()
ATAC_and_scRNA_gene_overlap <- list()
ATAC_and_bulkRNA_gene_overlap <- list()
ATAC_and_snME_and_mintchip_gene_overlap <- list()
ATAC_and_snME_and_scRNA_gene_overlap <- list()
ATAC_and_snME_and_bulkRNA_gene_overlap <- list()
ATAC_and_mintchip_and_scRNA_gene_overlap <- list()
ATAC_and_mintchip_and_bulkRNA_gene_overlap <- list()
ATAC_and_scRNA_and_bulkRNA_gene_overlap <- list()
ATAC_and_snME_and_mintchip_and_scRNA_gene_overlap <- list()
ATAC_and_snME_and_mintchip_and_bulkRNA_gene_overlap <- list()
ATAC_and_snME_and_scRNA_and_bulkRNA_gene_overlap <- list()
ATAC_and_mintchip_and_scRNA_and_bulkRNA_gene_overlap <- list()
ATAC_and_snME_and_mintchip_and_scRNA_and_bulkRNA_gene_overlap <- list()
# Instead of this, why not use DF?
# Row names are genes from ATAC
# I guess I need to check whether a given gene is 
# Columns are cell type, ATAC (always upregulated or downregulated), 
# mintchip (list of markers followed by upregulated, downregulated, or NA), 
# snME (cell type followed by upregulated, downregulated, or NA),
# scRNA-seq (cell type followed by upregulated, downregulated, or NA),
# bulkRNA-seq (upregulated, downregulated, or NA)
# Go gene by gene 

for(atac_cell_type in atac_cell_types) {
  atac_cell_type_for_file_name <- sub(" ", "_", atac_cell_type)
  pos_atac_peaks <- read.table(paste0(sc_das_annotated_dir, atac_cell_type_for_file_name, "_FC_0.1_upregulated_annotated.tsv"),
                               sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  pos_atac_genes <- unique(pos_atac_peaks$SYMBOL)
  neg_atac_peaks <- read.table(paste0(sc_das_annotated_dir, atac_cell_type_for_file_name, "_FC_0.1_downregulated_annotated.tsv"),
                               sep = "\t", header = TRUE, comment.char = "", quote = "\"")
  neg_atac_genes <- unique(neg_atac_peaks$SYMBOL)
  # Set associated cell types
  if(atac_cell_type == "B") {
    snME_cell_types <- c("B-Mem", "B-Naive")
    scRNA_cell_types <- c("B memory", "B naive")
  } else if(atac_cell_type == "CD14 Mono" | atac_cell_type == "CD16 Mono") {
    snME_cell_types <- "Monocyte"
    scRNA_cell_types <- c("CD14 Mono", "CD16 Mono")
  } else if(atac_cell_type == "NK") {
    snME_cell_types <- "NK-cell2"
    scRNA_cell_types <- c("NK", "NK_CD56bright")
  } else if(atac_cell_type == "T Naive") {
    snME_cell_types <- c("Tc-Naive", "Th-Naive")
    scRNA_cell_types <- c("CD4 Naive", "CD8 Naive")
  } else if(atac_cell_type == "CD4 Memory") {
    snME_cell_types <- c("Th-Mem")
    scRNA_cell_types <- c("CD4 Memory")
  } else if(atac_cell_type == "CD8 Memory") {
    snME_cell_types <- c("Tc-Mem")
    scRNA_cell_types <- c("CD8 Memory")
  }
  # Grab snME DMRs
  pos_snME_dmrs <- list()
  neg_snME_dmrs <- list() 
  for(snME_cell_type in snME_cell_types) {
    current_pos_snME_dmrs <- read.table(paste0(snME_results_dir, "Closest_Gene_Analysis/", snME_cell_type, "_D28_hypermethylated_annotated_genes.tsv"),
                                        sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    current_neg_snME_dmrs <- read.table(paste0(snME_results_dir, "Closest_Gene_Analysis/", snME_cell_type, "_D_minus_1_hypermethylated_annotated_genes.tsv"),
                                        sep = "\t", header = TRUE, comment.char = "", quote = "\"")
    pos_snME_dmrs[[snME_cell_type]] <- unique(current_pos_snME_dmrs$SYMBOL)
    neg_snME_dmrs[[snME_cell_type]] <- unique(current_neg_snME_dmrs$SYMBOL)
  }
  # Grab scRNA-seq
  pos_scRNA_degs <- list()
  neg_scRNA_degs <- list()
  for(scRNA_cell_type in scRNA_cell_types) {
    current_scRNA_degs <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == scRNA_cell_type,]
    current_pos_scRNA_degs <- current_scRNA_degs[current_scRNA_degs$sc_log2FC > 0,]
    current_neg_scRNA_degs <- current_scRNA_degs[current_scRNA_degs$sc_log2FC < 0,]
    pos_scRNA_degs[[scRNA_cell_type]] <- unique(current_pos_scRNA_degs$Gene_Name)
    neg_scRNA_degs[[scRNA_cell_type]] <- unique(current_neg_scRNA_degs$Gene_Name)
  }
  # Overlap between positive ATAC and positive snME
  pos_ATAC_pos_snME_genes <- c()
  pos_ATAC_pos_snME_scRNA_genes <- c()
  pos_ATAC_pos_snME_bulkRNA_genes <- c()
  for(snME_peak_list in pos_snME_dmrs) {
    pos_ATAC_pos_snME_genes <- c(pos_ATAC_pos_snME_genes, intersect(pos_atac_peaks$SYMBOL, snME_peak_list$SYMBOL))
  }
  pos_ATAC_pos_snME_genes <- unique(pos_ATAC_pos_snME_genes)
  # Overlap with scRNA
  for(scRNA_cell_type in scRNA_cell_types) {
    current_sc_pseudobulk_deg_table <- sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == scRNA_cell_type,]
    pos_ATAC_pos_snME_scRNA_genes <- c(pos_ATAC_pos_snME_scRNA_genes, intersect(unique(current_sc_pseudobulk_deg_table$Gene_Name))
  }
  pos_ATAC_pos_snME_scRNA_genes <- unique(pos_ATAC_pos_snME_scRNA_genes)
  # Overlap with bulk
  pos_ATAC_pos_snME_scRNA_genes <- c(pos_ATAC_pos_snME_scRNA_genes, unique(current_sc_pseudobulk_deg_table$Gene_Name))
  pos_ATAC_pos_snME_scRNA_genes <- unique(pos_ATAC_pos_snME_scRNA_genes)
  # Overlap between positive ATAC and negative snME
  pos_ATAC_neg_snME_genes <- c()
  for(snME_peak_list in neg_snME_dmrs) {
    pos_ATAC_neg_snME_genes <- c(pos_ATAC_neg_snME_genes, intersect(pos_atac_peaks$SYMBOL, snME_peak_list$SYMBOL))
  }
  pos_ATAC_neg_snME_genes <- unique(pos_ATAC_neg_snME_genes)
  # Overlap between negative ATAC and positive snME
  neg_ATAC_pos_snME_genes <- c()
  for(snME_peak_list in pos_snME_dmrs) {
    neg_ATAC_pos_snME_genes <- c(neg_ATAC_pos_snME_genes, intersect(neg_atac_peaks$SYMBOL, snME_peak_list$SYMBOL))
  }
  neg_ATAC_pos_snME_genes <- unique(neg_ATAC_pos_snME_genes)
  # Overlap between negative ATAC and negative snME
  neg_ATAC_neg_snME_genes <- c()
  for(snME_peak_list in neg_snME_dmrs) {
    neg_ATAC_neg_snME_genes <- c(neg_ATAC_neg_snME_genes, intersect(neg_atac_peaks$SYMBOL, snME_peak_list$SYMBOL))
  }
  neg_ATAC_neg_snME_genes <- unique(neg_ATAC_neg_snME_genes)
  list_of_overlap <- list(pos_ATAC_pos_snME_genes, pos_ATAC_neg_snME_genes, neg_ATAC_pos_snME_genes, neg_ATAC_neg_snME_genes)
  ATAC_and_snME_gene_overlap[[atac_cell_type]] <- list_of_overlap
}
  