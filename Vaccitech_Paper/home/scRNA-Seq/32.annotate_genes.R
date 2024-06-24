# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

innate_scRNA_hvl_placebo_degs <- scRNA_hvl_placebo_degs[scRNA_hvl_placebo_degs$Cell_Type %in% innate_cell_types,]
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c('hgnc_symbol', 'description', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'hgnc_symbol',
  values = unique(innate_scRNA_hvl_placebo_degs$Gene_Name),
  mart = ensembl
)

# Search for histone related genes
annotations[grep("histone", annotations$description),]