# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))
magical_results <- read.table(file = paste0(magical_output_dir, "MAGICAL_overall_output.tsv"), sep = "\t", header = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

atac_cell_types_for_snME_analysis <- c("B", "CD4 Memory", "CD8 Memory", "CD14 Mono", "CD16 Mono", "NK", "T Naive")
for(atac_cell_type in atac_cell_types_for_snME_analysis) {
  print(atac_cell_type)
  atac_cell_type_for_file_name <- sub(" ", "_", atac_cell_type)
  atac_cell_type_peaks <- read.table(paste0(sc_das_dir, "diff_peaks/D28-vs-D_minus_1-degs-", atac_cell_type_for_file_name, 
                                            "-time_point-controlling_for_subject_id_overlapping_peak_pct_0.01.tsv"),
                                     sep = "\t", header = TRUE)
  components <- strsplit(atac_cell_type_peaks$Peak_Name, "[:-]")
  chr <- c()
  start <- c()
  end <- c()
  for(component in components) {
    chr <- c(chr, component[1])
    start <- c(start, component[2])
    end <- c(end, component[3])
  }
  atac_cell_type_peaks$seqnames <- chr
  atac_cell_type_peaks$start <- as.integer(start)
  atac_cell_type_peaks$end <- as.integer(end)
  atac_cell_type_peaks_with_gene_info <- sc_peaks %>%
    inner_join(atac_cell_type_peaks, by = c("seqnames", "start", "end"))
  # TODO: Limit to 50kb distance from gene or TSS?
  if(atac_cell_type == "B") {
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
  for(snME_cell_type in snME_cell_types) {
    snME_dms_for_current_cell_type <- snME_dms[snME_dms$celltype == snME_cell_type,]
    upregulated_snME_dms_for_current_cell_type <- snME_dms_for_current_cell_type[snME_dms_for_current_cell_type$methylation == "upregulated",]
    downregulated_snME_dms_for_current_cell_type <- snME_dms_for_current_cell_type[snME_dms_for_current_cell_type$methylation == "downregulated",]
    
    upregulated_differential_analysis_results <- upregulated_snME_dms_for_current_cell_type[,c(2,3,4)]
    upregulated_differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(upregulated_differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
    upregulated_differential_analysis_results_df <- as.data.frame(upregulated_differential_analysis_results)
    upregulated_differential_analysis_results_df <- subset(upregulated_differential_analysis_results_df, distanceToTSS >= -2000 & distanceToTSS <= 500)
    
    downregulated_differential_analysis_results <- downregulated_snME_dms_for_current_cell_type[,c(2,3,4)]
    downregulated_differential_analysis_results <- annotatePeak(makeGRangesFromDataFrame(downregulated_differential_analysis_results), TxDb = txdb, annoDb = "org.Hs.eg.db")
    downregulated_differential_analysis_results_df <- as.data.frame(downregulated_differential_analysis_results)
    downregulated_differential_analysis_results_df <- subset(downregulated_differential_analysis_results_df, distanceToTSS >= -2000 & distanceToTSS <= 500)
    
    print(intersect(upregulated_differential_analysis_results_df$SYMBOL, atac_cell_type_peaks_with_gene_info$nearestGene))
    coregulated_up_snME_and_snATAC_genes <- intersect(upregulated_differential_analysis_results_df$SYMBOL, atac_cell_type_peaks_with_gene_info$nearestGene)
    print(intersect(downregulated_differential_analysis_results_df$SYMBOL, atac_cell_type_peaks_with_gene_info$nearestGene))
    coregulated_down_snME_and_snATAC_genes <- intersect(downregulated_differential_analysis_results_df$SYMBOL, atac_cell_type_peaks_with_gene_info$nearestGene)
    
    if(atac_cell_type == "B") {
      rna_cell_types <- c("B memory", "B naive")
    } else if(atac_cell_type == "CD14 Mono") {
      rna_cell_types <- "CD14 Mono"
    } else if(atac_cell_type == "CD16 Mono") {
      rna_cell_types <- "CD16 Mono"
    } else if(atac_cell_type == "NK") {
      rna_cell_types <- c("NK", "NK_CD56bright")
    } else if(atac_cell_type == "T Naive") {
      rna_cell_types <- c("CD4 Naive", "CD8 Naive")
    } else if(atac_cell_type == "CD4 Memory") {
      rna_cell_types <- c("CD4 Memory")
    } else if(atac_cell_type == "CD8 Memory") {
      rna_cell_types <- c("CD8 Memory")
    }
    print("Intersection between snME, snATAC, and scRNA-seq")
    for(rna_cell_type in rna_cell_types) {
      relevant_degs <- unique(sc_pseudobulk_deg_table[sc_pseudobulk_deg_table$Cell_Type == rna_cell_type,]$Gene_Name)
      print(intersect(coregulated_up_snME_and_snATAC_genes, relevant_degs))
      print(intersect(coregulated_down_snME_and_snATAC_genes, relevant_degs))
    }
    print("Intersection between snME, snATAC, and bulk")
    print(intersect(coregulated_up_snME_and_snATAC_genes, rownames(raw_high_placebo_period_2_D28_vs_D_minus_1_results)))
    print(intersect(coregulated_down_snME_and_snATAC_genes, rownames(raw_high_placebo_period_2_D28_vs_D_minus_1_results)))
  }
}