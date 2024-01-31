# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))




txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

f(atac_cell_type == "B") {
  snME_cell_types <- c("B-Mem", "B-Naive")
} else if(atac_cell_type == "CD14 Mono" | atac_cell_type == "CD16 Mono") {
  snME_cell_types <- "Monocyte"
} else if(atac_cell_type == "NK") {
  snME_cell_types <- "NK-cell2"
} else if(atac_cell_type == "T Naive") {
  snME_cell_types <- c("Tc-Naive", "Th-Naive")
} else if(atac_cell_type == "CD4 Memory") {
  snME_cell_types <- c("Tc-Mem", "Th-Mem")
} else if(atac_cell_type == "CD8 Memory") {
  snME_cell_types <- c("Tc-Mem", "Th-Mem")
}

# Find overlap of genes between snME and ATAC
snME_and_ATAC_gene_overlap <- list()
for(snME_cell_type in snME_cell_types) {
  
}