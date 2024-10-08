library(data.table)
library(pheatmap)
library(PLIER)

##### SETUP #####
# Read in count and metadata files
base_dir <- "C:/Users/willi/Documents/GitHub/Influenza/"
source(paste0(base_dir, "bulk_RNA_analysis_helper.R"))
data_dir <- "C:/Users/willi/Documents/local_data_files/"
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
# Load background pathway info for PLIER
data(bloodCellMarkersIRISDMAP)
data(canonicalPathways)
data(svmMarkers)
allPaths <- combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, canonicalPathways)
# ALL TIME POINTS
full_time_placebo_metadata <- placebo_metadata[placebo_metadata$subject_id 
                                               %in% names(table(placebo_metadata$subject_id)
                                                          [table(placebo_metadata$subject_id) == 10]),]
full_time_placebo_counts <- placebo_counts[rownames(full_time_placebo_metadata)]
# Remove bad subject (0 qPCR)
removed_low_value_aliquots <- rownames(full_time_placebo_metadata[full_time_placebo_metadata$subject_id == "f18c54d93cef4a4e",])
full_time_placebo_metadata <- full_time_placebo_metadata[!(full_time_placebo_metadata$subject_id %in% "f18c54d93cef4a4e"),]
full_time_placebo_counts <- full_time_placebo_counts[,!(colnames(full_time_placebo_counts) %in% removed_low_value_aliquots)]
# Run PLIER on Period 1 (including pre-vaccination)
period_1_placebo_metadata <- full_time_placebo_metadata[full_time_placebo_metadata$period == 1,]
period_1_placebo_counts <- placebo_counts[rownames(period_1_placebo_metadata)]
period_1_placebo_plier <- PLIER(as.matrix(period_1_placebo_counts), allPaths, scale = F)
plotU(period_1_placebo_plier, auc.cutoff = 0.70, fdr.cutoff = 0.05, top = 3)
# This is a clean heatmap. Seems like LV22 = NK, LV8 = B Cell, LV7 = something not cell type related,
# LV6 = T Cell, LV4 = Erythoid cell, LV3 = Type 1 interferon response (not cell type), 
# LV2 = Neutrophil, LV21 = Megakaryocyte
period_1_placebo_markers <- plierResToMarkers(period_1_placebo_plier, allPaths, index = c(2, 3, 4, 6, 7, 8, 21, 22))
colnames(period_1_placebo_markers) # Why does 3 not match the heatmap?
indexToPlot <- which(apply(period_1_placebo_plier$Uauc*(period_1_placebo_plier$Up<0.001),2,max)>0.75)
plotTopZ(period_1_placebo_plier, period_1_placebo_counts, allPaths, top = 5, index = indexToPlot)
# Run PLIER on Period 1 (after vaccination)
period_1_after_vaccination_placebo_metadata <- period_1_placebo_metadata[period_1_placebo_metadata$time_point != "1_D_minus_1",]
period_1_after_vaccination_placebo_counts <- placebo_counts[rownames(period_1_after_vaccination_placebo_metadata)]
period_1_after_vaccination_placebo_plier <- PLIER(as.matrix(period_1_after_vaccination_placebo_counts), allPaths, scale = F)
plotU(period_1_after_vaccination_placebo_plier, auc.cutoff = 0.70, fdr.cutoff = 0.05, top = 3)
# Run PLIER on Period 2 (including pre-treatment)
period_2_placebo_metadata <- full_time_placebo_metadata[full_time_placebo_metadata$period == 2,]
period_2_placebo_metadata <- period_2_placebo_metadata[period_2_placebo_metadata$time_point != "2_D_minus_2",]
period_2_placebo_counts <- placebo_counts[rownames(period_2_placebo_metadata)]
period_2_placebo_plier <- PLIER(as.matrix(period_2_placebo_counts), allPaths, scale = F)
plotU(period_2_placebo_plier, auc.cutoff = 0.70, fdr.cutoff = 0.05, top = 3)
# Run PLIER on Period 2 (after treatment)
period_2_after_treatment_placebo_metadata <- period_2_placebo_metadata[period_2_placebo_metadata$time_point != "2_D_minus_1",]
period_2_after_treatment_placebo_counts <- placebo_counts[rownames(period_2_after_treatment_placebo_metadata)]
period_2_after_treatment_placebo_plier <- PLIER(as.matrix(period_2_after_treatment_placebo_counts), allPaths, scale = F)
plotU(period_2_after_treatment_placebo_plier, auc.cutoff = 0.70, fdr.cutoff = 0.05, top = 3)
