find_sex_associated_genes=function(sex_associated_dir, padj_threshold = 0.05, log2fc_threshold = 0.1) {
  sex_associated_gene_files <- list.files(sex_associated_dir, pattern = ".csv")
  sex_associated_genes <- c()
  number_of_studies <- length(sex_associated_gene_files)
  for (sex_associated_gene_file in sex_associated_gene_files) {
    print(sex_associated_gene_file)
    sex_associated_gene_file_contents <- read.csv(paste0(sex_associated_dir, sex_associated_gene_file))
    row.names <- as.character(sex_associated_gene_file_contents$X)
    sex_associated_gene_file_contents <- sex_associated_gene_file_contents[,2:ncol(sex_associated_gene_file_contents)]
    sex_associated_gene_file_contents <- as.data.frame(sex_associated_gene_file_contents)
    rownames(sex_associated_gene_file_contents) <- row.names
    sex_associated_gene_file_contents <- subset(sex_associated_gene_file_contents, padj < padj_threshold & 
                                                  abs(log2FoldChange) > log2fc_threshold)
    for (entry in rownames(sex_associated_gene_file_contents)) {
      if(entry %in% genemap$ensembl_gene_id) {
        sex_associated_genes <- c(sex_associated_genes, unique(genemap[genemap$ensembl_gene_id == entry,]$hgnc_symbol))
      } else{
        print(entry)
      }
    }
  }
  # Remove repeated genes and the "" entries (no HGNC symbol)
  #sex_associated_genes <- unique(sex_associated_genes[sex_associated_genes != ""])
  return(sex_associated_genes)
}
