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

# annotate MAGICAL genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c('hgnc_symbol', 'description', 'chromosome_name', 'start_position', 'end_position'),
  filters = 'hgnc_symbol',
  values = sort(unique(overall_magical_df$Gene_symbol)) ,
  mart = ensembl
)

# ABHD17A, ASH1L, BRD1, CASP1, CSF1R, CSK, DNAJC3, DUSP1, DUSP2, DUSP6, ETS1, FOSB, FOSL2, GBP5, IFNG, IL1RAP, IRAK3, IRF2, JDP2, KDM2B, KDM3A, MAP3K11, MAP4K3, MAPK7, METTL23, MINK1, PSMB9, RIPK1, SETD2, STAT4

annotations_final <- annotations %>%
  filter(
    str_detect(hgnc_symbol, paste(gene_terms, collapse = "|")) |
      str_detect(description, paste(search_terms, collapse = "|"))
  )

important_magical_genes <- c("ABHD17A", "ASH1L", "BRD1", "CASP1", "CSF1R", "CSK", "DNAJC3", "DUSP1", "DUSP2", "DUSP6", "ETS1", "FOSB", "FOSL2", "GBP5", "IFNG", "IL1RAP", "IRAK3", "IRF2", "ITGAL", "JDP2", "KDM2B", "KDM3A", "MAP3K11", "MAP4K3", "MAPK7", "METTL23", "MINK1", "PSMB9", "RIPK1", "SETD2", "STAT4")
important_overall_magical_df <- overall_magical_df[overall_magical_df$Gene_symbol %in% important_magical_genes,]
important_magical_gene_overlap_df <- magical_gene_overlap_df[magical_gene_overlap_df$Gene_Name %in% important_magical_genes,]
write.table(important_magical_gene_overlap_df, file = "C:/Users/wat2/Desktop/important_magical_gene_overlap.tsv", sep = "\t", row.names = FALSE)