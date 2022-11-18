library(DESeq2)
library(data.table)
library(tximport)

base_dir <- "C:/Users/wat2/Documents/GitHub/Influenza/"
data_dir <- "C:/Users/wat2/Documents/local_data_files/"
counts <- fread(paste0(data_dir, "rsem_genes_count.processed.txt"), header = T, sep = "\t")
row.names <- as.character(counts$gene_id)
counts <- counts[,2:ncol(counts)]
counts <- as.data.frame(counts)
rownames(counts) <- row.names

placebo_metadata_file <- paste0(base_dir, "placebo_metadata_sheet.tsv")
placebo_metadata <- read.table(placebo_metadata_file, header = TRUE, sep = "\t")
placebo_metadata <- placebo_metadata[placebo_metadata$bulkRNA_seq == TRUE,]
placebo_metadata <- placebo_metadata[placebo_metadata$subject_id %in% names(table(placebo_metadata$subject_id)[table(placebo_metadata$subject_id) == 10]),]
placebo_metadata$time_point <- paste0(placebo_metadata$period, "_", placebo_metadata$time_point)
placebo_metadata = subset(placebo_metadata, select = -c(period))
placebo_metadata$time_point[placebo_metadata$time_point == '1_D1 predose'] <- '1_D-1'
kept_aliquots <- unique(placebo_metadata$aliquot_id)
counts <- counts[kept_aliquots]
colnames(counts) <- sort(colnames(counts))
rownames(placebo_metadata) <- placebo_metadata$aliquot_id
rownames(placebo_metadata) <- sort(rownames(placebo_metadata))
placebo_metadata = subset(placebo_metadata, select = -c(aliquot_id))
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = placebo_metadata,
                              design = ~ time_point)
dds
