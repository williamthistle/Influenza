# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

source(paste0(base_dir, "extra_functions/MAGICAL_functions.R"))
source(paste0(base_dir, "extra_functions/my_MAGICAL_helper_functions.R"))

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

overall_magical_df <- read.table(paste0(MAGICAL_hvl_placebo_output_dir, "MAGICAL_overall_output_with_stage_1_info.tsv"), sep = "\t", header = TRUE)
magical_gene_overlap_df <- create_magical_gene_overlap_df(overall_magical_df, hvl_full_time_series_placebo_period_2_D5_vs_D_minus_1_results, hvl_full_time_series_placebo_period_2_D8_vs_D_minus_1_results)
#magical_site_overlap_df <- create_magical_site_overlap_df(overall_magical_df)



# ABHD17A, ASH1L, BRD1, CASP1, CSF1R, CSK, DNAJC3, DUSP1, DUSP2, DUSP6, ETS1, FOSB, FOSL2, GBP5, IFNG, IL1RAP, IRAK3, IRF2, JDP2, KDM2B, KDM3A, MAP3K11, MAP4K3, MAPK7, METTL23, MINK1, NFIL3, PSMB9, RIPK1, SETD2, STAT4

annotations_final <- annotations %>%
  filter(
    str_detect(hgnc_symbol, paste(gene_terms, collapse = "|")) |
      str_detect(description, paste(search_terms, collapse = "|"))
  )

important_magical_genes <- c("IL32", "CASP1", "CSF1R", "IRAK3", "IL1RAP", "NFIL3", "PTGES", "IL34", "CSF1R", "IL1B", "IL15RA", "IL21R", "IL1RN", "IL6R", "IL6R", "IL5RA", "IL4R", "IL6", "IL37", "IL17F", "IL20RB", "IRAK4", "IL10RA", "IL17B", "IL11", "IRF2", "IRF7", "IFNG", "OAS1", "PSMB9", "MNDA", "DNAJC3", "IFNGR1", "GBP5", "USP38", "IRF2", "IRF1", "IFNGR1", "IFI35", "IRF8", "IFITM3", "IRF5", "IFI6", "IRF1", "STING1", "IRF5", "IFI27", "IFNAR2", "LUARIS", "STING1", "ABHD17A", "CSK", "MAP3K11", "DUSP1", "DUSP2", "DUSP6", "MAPK7", "MAPK8", "MAP2K1", "MAP3K8", "MAP3K20", "MAP4K3", "MAPKAPK2", "PTK2B", "RELL1", "MINK1", "BRAF", "RIPK1", "MAPK13", "MAP2K5", "MAP3K11", "DUSP16", "DUSP16", "MAP3K15", "MAPK3", "MAP3K11", "MAP4K4", "MAP3K6", "MAPKBP1", "MAP3K3", "MAPK8IP1", "MINK1", "MAP4K4", "MAP4K3", "RIPK1", "MAP2K3", "MAP4K4", "MAP4K1", "MAP2K6", "MAP2K3", "MAPK15", "MAP2K2", "MAP4K2", "MAP2K1", "MAPKAP1", "MAPK10", "MAP4K1", "FOSL2", "FOSL1", "JUND", "JUN", "JUNB", "JUND", "JDP2", "FOSL2", "FOS", "JAK1", "STAT3", "STAT4", "JAK1", "STAT2", "STAT4", "STAT6", "JAK2", "JAK3", "STAT4", "JAK2", "JAK2", "JAK1", "STAT5A", "CCL3", "CX3CR1", "CCL3L1", "CXCL16 ", "CXCR2", "CCL22", "CCL26", "CXCR3", "CCL17", "CCR7", "CXCL12", "XCR1", "HIST1H1C", "HIST1H1D", "HIST1H1E", "H2AFZ", "H2BC5", "H2BC18", "H4C13", "H2BC15", "H3C2", "H2AC14", "H2AC4", "H4C8", "H2AC13", "H3-3A", "H2BC10", "H4C5", "H3-3B", "H2AC25", "H1-4", "KAT6A", "BRD1", "MSL2", "KAT2B", "KAT2A", "KAT2B", "BRD1", "MSL2", "KAT5", "HDAC5", "HDAC7", "HDAC9", "HDAC10", "HDAC4", "HDAC4", "HDAC11", "METTL23", "ASH1L", "PRDM2", "SETD2", "ASH1L", "SUV39H1", "ASH2L", "KMT2B", "SETD7", "PRMT2", "KMT2C", "KDM3A", "KDM2B", "KDM2B", "KDM5B", "KDM2A", "KDM7A", "KDM2A", "KDM8", "KDM4A", "KDM4B", "KDM2B", "KDM6A", "KDM4C", "KDM4A", "KDM4B", "KDM4E", "KDM5A")
important_overall_magical_df <- overall_magical_df[overall_magical_df$Gene_symbol %in% important_magical_genes,]
important_magical_gene_overlap_df <- magical_gene_overlap_df[magical_gene_overlap_df$Gene_Name %in% important_magical_genes,]
write.table(important_magical_gene_overlap_df, file = "C:/Users/willi/Desktop/important_magical_gene_overlap.tsv", sep = "\t", row.names = FALSE)