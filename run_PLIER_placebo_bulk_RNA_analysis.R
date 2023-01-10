library(PLIER)
library(data.table)
library(pheatmap)

##### SETUP #####
# Read in count and metadata files
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
data_dir <- "C:/Users/willi/Documents/local_data_files/"
#load(paste0(data_dir, "placebo_bulk_RNA_obj.RData"))
counts <- fread(paste0(data_dir, "rsem_genes_count.processed.txt"), header = T, sep = "\t")
all_metadata_file <- paste0(base_dir, "all_metadata_sheet.tsv")
all_metadata <- read.table(all_metadata_file, header = TRUE, sep = "\t")
viral_load_file <- paste0(base_dir, "bulk_RNA_viral_load.tsv")
viral_load <- read.table(viral_load_file, header = TRUE, sep = "\t")
viral_load_primary <- viral_load[viral_load$PARAMCD == "QPCRAUC",]
viral_load_primary <- viral_load_primary[viral_load_primary$TRT01A == "PLACEBO",]
viral_load_primary$AVAL <- as.numeric(viral_load_primary$AVAL)
# Organize by viral load (high to low) and grab top 13 subjects - they will be high
viral_load_primary <- viral_load_primary[order(viral_load_primary$AVAL, decreasing = TRUE),]
high_viral_load_subjects <- viral_load_primary$SUBJID[1:13]
# Take gene_id column from counts and use contents as the rownames of counts
row.names <- as.character(counts$gene_id)
counts <- counts[,2:ncol(counts)]
counts <- as.data.frame(counts)
rownames(counts) <- row.names
# Only keep metadata for bulk RNA-seq aliquots
all_metadata <- all_metadata[all_metadata$bulkRNA_seq == TRUE,]
# Add period into time point (by itself, time point isn't unique - we need period information
# to distinguish between D-1 in period 1 vs D-1 in period 2, for example)
all_metadata$time_point <- paste0(all_metadata$period, "_", all_metadata$time_point)
# Make time point names safe for DESeq2
all_metadata$time_point[all_metadata$time_point == '1_D1 predose'] <- '1_D_minus_1'
all_metadata$time_point[all_metadata$time_point == '2_D-2'] <- '2_D_minus_2'
all_metadata$time_point[all_metadata$time_point == '2_D-1'] <- '2_D_minus_1'
# Divide metadata into placebo
placebo_metadata <- all_metadata[all_metadata$treatment == "PLACEBO",]
viral_load_vector <- c()
for (subject_id in placebo_metadata$subject_id) {
  if(subject_id %in% high_viral_load_subjects) {
    viral_load_vector <- c(viral_load_vector, "HIGH")
  } else {
    viral_load_vector <- c(viral_load_vector, "LOW")
  }
}
placebo_metadata$viral_load <- viral_load_vector
# Find placebo-associated counts
kept_aliquots <- placebo_metadata$aliquot_id
placebo_counts <- counts[kept_aliquots]
# Sort columns in counts and rows for each so they're in same order (for DESeq2)
colnames(counts) <- sort(colnames(counts))
rownames(all_metadata) <- all_metadata$aliquot_id
rownames(all_metadata) <- sort(rownames(all_metadata))
colnames(placebo_counts) <- sort(colnames(placebo_counts))
rownames(placebo_metadata) <- placebo_metadata$aliquot_id
rownames(placebo_metadata) <- sort(rownames(placebo_metadata))
# Drop aliquot ID column (it's stored in rownames)
all_metadata <- subset(all_metadata, select = -c(aliquot_id))
placebo_metadata = subset(placebo_metadata, select = -c(aliquot_id))
# Probably OK to round expected counts from RSEM data. DESeq2 expects integers
counts <- round(counts)
placebo_counts <- round(placebo_counts)

