library(data.table)
library(pheatmap)
library(PLIER)

##### SETUP #####
# Read in count and metadata files
base_dir <- "C:/Users/wat2/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
data_dir <- "C:/Users/wat2/Documents/local_data_files/"
#load(paste0(data_dir, "placebo_bulk_RNA_obj.RData"))
counts <- fread(paste0(data_dir, "rsem_genes_count.processed.normalized.txt"), header = T, sep = "\t")
all_metadata_file <- paste0(base_dir, "all_metadata_sheet.tsv")
all_metadata <- read.table(all_metadata_file, header = TRUE, sep = "\t")
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
# Make time point names safe
all_metadata$time_point[all_metadata$time_point == '1_D1 predose'] <- '1_D_minus_1'
all_metadata$time_point[all_metadata$time_point == '2_D-2'] <- '2_D_minus_2'
all_metadata$time_point[all_metadata$time_point == '2_D-1'] <- '2_D_minus_1'
# Divide metadata into placebo
placebo_metadata <- all_metadata[all_metadata$treatment == "PLACEBO",]
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
# ALL TIME POINTS
full_time_placebo_metadata <- placebo_metadata[placebo_metadata$subject_id 
                                               %in% names(table(placebo_metadata$subject_id)
                                                          [table(placebo_metadata$subject_id) == 10]),]
full_time_placebo_counts <- placebo_counts[rownames(full_time_placebo_metadata)]
# PERIOD 1
period_1_placebo_metadata <- full_time_placebo_metadata[full_time_placebo_metadata$period == 1,]
period_1_placebo_counts <- placebo_counts[rownames(period_1_placebo_metadata)]
# Run PLIER on Period 1 and see what happens
data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)
data(svmMarkers)
data(vacData)
allPaths <- combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, canonicalPathways, vacData)
plierResult <- PLIER(as.matrix(period_1_placebo_counts), allPaths, scale = F)
plotU(plierResult, auc.cutoff = 0.70, fdr.cutoff = 0.05, top = 3)

