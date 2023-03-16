library(ArchR)
library(stringr)
library(pheatmap)
library(mclust)
library(writexl)
library(ggplot2)
library(hexbin)
library(SeuratDisk)
library(dplyr)
library(openxlsx)
library(BSgenome.Hsapiens.UCSC.hg38)
library(outliers)

set.seed(1)
home_dir <- "~/SPEEDI"

# naming_token is used to select analysis and name output files
naming_token <- "high_vs_low_viral_load_D28"
date <- Sys.Date()
# data_path is where input data are stored
data_path <- paste0("/data/home/wat2/", naming_token, "/ATAC_seq_data/")
# Get list of samples that will be processed
sample_id_list <- list.dirs(data_path, recursive = FALSE)
sample_id_list <- strsplit(sample_id_list, "/")
sample_id_list <- unlist(lapply(sample_id_list, tail, n = 1L))
sample_count <- length(sample_id_list)
output_dir <- paste0("/data/home/wat2/", naming_token, "/ATAC_seq_data_output/")
viral_load_info <- read.table(paste0(home_dir, "/viral_load_info.tsv"), sep = "\t", header = TRUE)

image_dir <- paste0(output_dir, "images/")
if (!dir.exists(image_dir)) {dir.create(image_dir)}

# Identify high / low viral load samples
high_viral_load <- c()
low_viral_load <- c()
for(sample_id in sample_id_list) {
  if(sample_id %in% viral_load_info$aliquot) {
    current_info <- viral_load_info[viral_load_info$aliquot == sample_id,]
    if(current_info$viral_load == "high") {
      high_viral_load <- c(high_viral_load, sample_id)
    } else {
      low_viral_load <- c(low_viral_load, sample_id)
    }
  } else {
    stop("You have a sample that is not labeled as high or low viral load. Not currently supported.")
  }
}

high_viral_load <- sort(high_viral_load)
low_viral_load <- sort(low_viral_load)
all_viral_load <- c(high_viral_load, low_viral_load)
sample_id_list <- sample_id_list[order(match(sample_id_list, all_viral_load))]

# Label each sample as high or low viral load in input file paths and metadata
inputFiles <- paste0(data_path, sample_id_list, "/outs/atac_fragments.tsv.gz")
names(inputFiles) <- names(metadata) <- all_viral_load
metadata <- c(rep("HVL", length(high_viral_load)), rep("LVL", length(low_viral_load)))
names(metadata) <- all_viral_load

addArchRGenome("hg38")

#fragment data processing - read in data
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Don't set this too high because you can always increase later
  minFrags = 3000,
  maxFrags = 30000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# Fake cell check (find doublets)
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

# Create ArchR project based on processed input files (arrow files)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = paste0(output_dir, "ArchR/"),
  copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.
)

# Filter out doublets (fake cells)
proj <- filterDoublets(ArchRProj = proj)

load(paste0(image_dir, "atac_after_filtering_doublets.RData"))

# List information about available matrices in project
getAvailableMatrices(proj)
#------------------------------------------------------------------------------------------------
# add condition metadata (HVL or LVL) to cells in project
aa<-proj$Sample
idxSample <- which(proj$Sample %in% high_viral_load)
aa[idxSample]<-'HVL'

idxSample <- which(proj$Sample %in% low_viral_load)
aa[idxSample]<-'LVL'

proj <- addCellColData(ArchRProj = proj, data = aa, cells = proj$cellNames,name = "Conditions", force = TRUE)

# Plot out TSS Enrichment / Doublet Enrichment / Nucleosome Ratio for each sample
p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
p2 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "DoubletEnrichment", plotAs = "ridges")
p3 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges")

plotPDF(p1,p2,p3, name = "Integrated_Scores_Prefiltering.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)

# Filter out cells that don't meet TSS enrichment / doublet enrichment / nucleosome ratio criteria
idxPass <- which(proj$TSSEnrichment >= 8 & proj$NucleosomeRatio < 2 & proj$DoubletEnrichment < 5) 
cellsPass <- proj$cellNames[idxPass]
proj<-proj[cellsPass, ]

write.table(proj$cellNames, file = paste0(output_dir, "ATAC_all_cells.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
# List number of D1 and D28 cells and list number of cells remaining for each sample
table(proj$Conditions)
table(proj$Sample)

# Subset based on snRNA-seq cells and label cell types
curated_snRNA_seq_cells <- read.csv(paste0(output_dir, "high_vs_low_viral_load_D28_V3_final_cell_names_curated_2023-03-07.csv"), comment.char = "")
uncurated_snRNA_seq_cells <- read.csv(paste0(output_dir, "high_vs_low_viral_load_D28_V3_final_cell_names_uncurated_2023-03-06.csv"), comment.char = "")

# Perform dimensionality reduction on cells (addIterativeLSI), create UMAP embedding (addUMAP), 
# and add cluster information (addClusters)
addArchRThreads(threads = 8)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                        clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                        varFeatures = 15000, dims = 2:30)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 5, knnAssign = 30, maxClusters = NULL, force = TRUE)

save.image(paste0(image_dir, "atac_after_lsi_umap_clusters_1.RData"))

# Determine whether there is similarity between our clusters and clusters according to sample or condition
# (0 is random, 1 is perfect)
adjustedRandIndex(proj$Sample, proj$Clusters)
adjustedRandIndex(proj$Conditions, proj$Clusters)

# UMAP plots colored by condition, sample, cluster ID, and TSS enrichment
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)

plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_snRNA.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Load reference scRNA data
reference_dir <- "~/reference/"
if (!dir.exists(reference_dir)) {dir.create(reference_dir)}
scRNA <- LoadH5Seurat(paste0(reference_dir, "multi.h5seurat"))

# Remove certain cell types we're not interested in
idx <- which(scRNA$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
scRNA <- scRNA[,-idx]

idx <- which(scRNA$celltype.l3 == "Treg Naive")
scRNA <- scRNA[,-idx]

# Rename CD4 Proliferating and CD8 Proliferating to T Proliferating
idx <- which(scRNA$celltype.l2 %in% c("CD4 Proliferating", "CD8 Proliferating"))
scRNA$celltype.l2[idx] <- "T Proliferating"

save.image(paste0(image_dir, "atac_before_gene_integration_matrix.RData"))

# Integrate scRNA-seq reference data into scATAC-seq data
addArchRThreads(threads = 6)
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",#Harmony
  seRNA = scRNA,
  addToArrow = FALSE,
  groupRNA = "celltype.l2",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  normalization.method = "SCT",
  force = TRUE
)

rm(scRNA)
save.image(paste0(image_dir, "atac_after_gene_integration_matrix_with_filtering.RData"))

# Plot integrated data with cell type predictions 
pal <- paletteDiscrete(values = proj$predictedGroup)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE)

plotPDF(p1, name = "Integrated_annotated_gene_integration_matrix.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Add in snRNA-seq predictions from clusters
predicted_snATAC_cells <- data.frame(cell_name = proj$cellNames, voted_type = proj$predictedGroup)
for(current_row in 1:nrow(predicted_snATAC_cells)) {
  current_snATAC_cell <- predicted_snATAC_cells[current_row,]$cell_name
  if(current_snATAC_cell %in% curated_snRNA_seq_cells$cells) {
    current_voted_type <- curated_snRNA_seq_cells[curated_snRNA_seq_cells$cells == current_snATAC_cell,]$voted_type
    predicted_snATAC_cells[current_row,]$voted_type <- current_voted_type
  }
}

proj <- addCellColData(ArchRProj = proj, data = predicted_snATAC_cells$voted_type, cells = proj$cellNames, name = "predictedGroup", force = TRUE)

# Combine cell types
Cell_type_combined = proj$predictedGroup
idx <- grep("CD4 T", Cell_type_combined)
Cell_type_combined[idx] <- "CD4 Memory"

idx <- grep("CD8 T", Cell_type_combined)
Cell_type_combined[idx] <- "CD8 Memory"

idx <- grep("cDC", Cell_type_combined)
Cell_type_combined[idx] <- "cDC"

idx <- grep("Proliferating", Cell_type_combined)
Cell_type_combined[idx] <- "Proliferating"

idx <- grep("B", Cell_type_combined)
Cell_type_combined[idx] <- "B"

proj <- addCellColData(ArchRProj = proj, data = Cell_type_combined, cells = proj$cellNames, name = "predictedGroup", force = TRUE)



# Plot integrated data with cell type predictions (combined)
pal <- paletteDiscrete(values = proj$predictedGroup)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE)

plotPDF(p1, name = "Integrated_annotated_combined_gene_integration_matrix.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Combine more cell types
proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "CD4 Naive", "T Naive")
proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "CD8 Naive", "T Naive")
proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "NK_CD56bright", "NK")
proj$Cell_type_combined <- replace(proj$Cell_type_combined, proj$Cell_type_combined == "ASDC", "CD14 Mono")
proj$Cell_type_combined <- replace(proj$Cell_type_combined, proj$Cell_type_combined == "cDC", "CD14 Mono")
proj$Cell_type_combined <- replace(proj$Cell_type_combined, proj$Cell_type_combined == "Eryth", "CD14 Mono")
proj$Cell_type_combined <- replace(proj$Cell_type_combined, proj$Cell_type_combined == "HSPC", "CD14 Mono")
proj$Cell_type_combined <- replace(proj$Cell_type_combined, proj$Cell_type_combined == "pDC", "CD14 Mono")
proj$Cell_type_combined <- replace(proj$Cell_type_combined, proj$Cell_type_combined == "Plasmablast", "CD14 Mono")
proj$Cell_type_combined <- replace(proj$Cell_type_combined, proj$Cell_type_combined == "Platelet", "CD14 Mono")
proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "Treg", "T Naive")

table(proj$Conditions)
table(proj$Sample)
table(proj$predictedGroup)



# First voting scheme
cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup))
pre_cluster <- rownames(cM)
max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]

Cell_type_voting <- proj$Clusters
for (m in c(1:length(pre_cluster))){
  idxSample <- which(proj$Clusters == pre_cluster[m])
  Cell_type_voting[idxSample] <- max_celltype[m]
}
proj <- addCellColData(ArchRProj = proj, data = Cell_type_voting, cells = proj$cellNames, name = "Cell_type_voting", force = TRUE)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)


plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_Gene_Integration_Voting_1_With_snRNA_5_Combined_T.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


# Alternative voting scheme
cluster.dump <- unique(proj$Clusters)
proj$Cell_type_voting <- proj$Clusters
for (i in unique(proj$predictedGroup)) {
  print(i)
  idxPass <- which(proj$predictedGroup %in% i)
  cellsPass <- proj$cellNames[idxPass]
  sub_proj <- proj[cellsPass, ]
  freq.table <- as.data.frame(table(sub_proj$Clusters))
  freq.table <- freq.table[order(freq.table$Freq, decreasing = TRUE),]
  freq.table$diff <- abs(c(diff(freq.table$Freq), 0))
  if(nrow(freq.table) > 30) {
    freq.table <- freq.table[1:30,]
  }
  if(length(unique(freq.table$diff)) == 1) {
    next
  }
  p.values <- dixon.test(freq.table$diff)$p.value[[1]]
  print(p.values)
  max.index <- which.max(freq.table$diff)
  clusters <- as.character(freq.table$Var1[1:max.index])
  proj$Cell_type_voting[proj$Cell_type_voting %in% clusters] <- i
  cluster.dump <- cluster.dump[!cluster.dump %in% clusters]
}

if (length(cluster.dump) > 0) {
  for (i in cluster.dump) {
    idxPass <- which(proj$Cell_type_voting %in% i)
    cellsPass <- proj$cellNames[idxPass]
    sub_proj <- proj[cellsPass, ]
    freq.table <- as.data.frame(table(sub_proj$predictedGroup))
    proj$Cell_type_voting[proj$Cell_type_voting %in% i] <- as.vector(freq.table$Var1)[which.max(freq.table$Freq)]
  }
}

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)


plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_Gene_Integration_Voting_2_With_snRNA_5_Combined_T.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# See how clusters are distributed
cluster_cell_type_predictions <- vector()
cluster_cell_type_distributions <- list()
cluster_sample_distributions <- list()
cluster_conditions_distributions <- list()
idx <- 1
unique_cluster_ids <- unique(proj$Clusters)
unique_cluster_ids <- unique_cluster_ids[order(nchar(unique_cluster_ids), unique_cluster_ids)]
for (cluster in unique_cluster_ids) {
  idxPass <- which(proj$Clusters %in% cluster)
  cellsPass <- proj$cellNames[idxPass]
  filtered_cluster <-proj[cellsPass,]
  cluster_cell_type_predictions <- append(cluster_cell_type_predictions, table(filtered_cluster$Cell_type_voting))
  cluster_cell_type_distributions[[idx]] <- table(filtered_cluster$predictedGroup)
  cluster_sample_distributions[[idx]] <- table(filtered_cluster$Sample)
  cluster_conditions_distributions[[idx]] <- table(filtered_cluster$Conditions)
  idx <- idx + 1
}
names(cluster_cell_type_predictions) <- paste(unique_cluster_ids, "-", names(cluster_cell_type_predictions))

# Override vote for 3 to make it CD16 Mono
idxPass <- which(proj$Clusters %in% c("C3"))
proj$Cell_type_voting[idxPass] <- "CD16 Mono"

#Remove the messy clusters (determined through visual inspection and seeing distribution of cells in each cluster)
#idxPass <- which(proj$Clusters %in% c("C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25")) - very thorough
#idxPass <- which(proj$Clusters %in% c("C1", "C2", "C3", "C4", "C6", "C24", "C25", "C26", "C27", "C29", "C30", "C31", "C32")) # Res 4
idxPass <- which(proj$Clusters %in% c("C2", "C12", "C13", "C14", "C15", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30", "C31", "C32", "C34", "C35", "C41")) # Res 5
cellsPass <- proj$cellNames[-idxPass]
proj_minus_clusters <- proj[cellsPass, ]




save.image(paste0(output_dir, "atac_after_removing_clusters.RData"))

for (cell_type in unique(proj_minus_clusters$Cell_type_voting)) {
  print(cell_type)
  idxPass <- which(proj_minus_clusters$Cell_type_voting %in% cell_type)
  cellsPass <- proj_minus_clusters$cellNames[idxPass]
  filtered_cluster <-proj[cellsPass,]
  print(length(cellsPass))
  print(table(filtered_cluster$Conditions))
}

saveArchRProject(ArchRProj = proj, outputDirectory = "final-HVL-vs-LVL-multiome-female_unfiltered", load = FALSE)


# Print final integrated data with cell type predictions
p1 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)




plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_Gene_Integration_Voting_1_With_snRNA_minus_clusters_combined_T_5_CD16Mono_Override_FINAL_test.pdf", ArchRProj = proj_minus_clusters, addDOC = FALSE, width = 5, height = 5)
plotEmbedding(ArchRProj = proj_minus_clusters, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE) + 
  labs(title = "snATAC-seq Data Integration \n (7 Samples, 42K Cells)") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.text=element_text(size=12)) +
  theme(plot.title = element_text(face="bold", size=16)) +
  theme(axis.title=element_text(size=14))
ggsave(paste0(output_dir, "Integrated_ATAC_plot.png"), device = "png", dpi = 300)

# Set up for peak calling and other MAGICAL tasks
proj <- proj_minus_clusters

# Create cell type proportion file for MAGICAL
cell_type_proportions_df <- data.frame("Condition" = metadata, "Sample_name" = names(metadata))
total_cell_counts_df <- data.frame("Sample_name" = names(metadata))
cell_counts <- vector()
# Find total cell counts for each sample
for (sample_id in names(metadata)) {
  idxPass <- which(proj$Sample %in% sample_id)
  print(length(idxPass))
  cellsPass <- proj$cellNames[idxPass]
  sample_subset <- subsetCells(proj, cellsPass)
  cell_counts <- append(cell_counts, nCells(sample_subset))
}
total_cell_counts_df <- cbind(total_cell_counts_df, cell_counts)

for (cell_type in unique(proj$Cell_type_voting)) {
  cell_type_proportions <- vector()
  print(cell_type)
  # Grab cells associated with cell type
  idxPass <- which(proj$Cell_type_voting %in% cell_type)
  print(length(idxPass))
  cellsPass <- proj$cellNames[idxPass]
  cells_subset <- subsetCells(proj, cellsPass)
  for (sample_id in names(metadata)) {
    # Subset further based on cells associated with sample ID
    idxPass <- which(cells_subset$Sample %in% sample_id)
    print(length(idxPass))
    cellsPass <- cells_subset$cellNames[idxPass]
    sample_subset <- subsetCells(cells_subset, cellsPass)
    cell_counts <- nCells(sample_subset)
    cell_type_proportions <- append(cell_type_proportions, cell_counts / total_cell_counts_df[total_cell_counts_df$Sample_name == sample_id,]$cell_counts)
  }
  temp_df <- data.frame(cell_type_proportions)
  names(temp_df)[names(temp_df) == "cell_type_proportions"] <- cell_type
  cell_type_proportions_df <- cbind(cell_type_proportions_df, temp_df)
}
write.csv(cell_type_proportions_df, file = paste0(output_dir, "ATAC_cell_type_proportion.csv"), quote = FALSE, row.names = FALSE)

# Apparently info about ArchR Genome is not stored in the .RData file, so if we load an .RData file above,
# We need to re-add our hg38 genome
addArchRGenome("hg38")

# Peak calling - first, call addGroupCoverages to find pseudo-bulk replicates, then call peaks using MACS2
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Cell_type_voting")
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Cell_type_voting", 
  pathToMacs2 = pathToMacs2
)
# Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
# calculations
proj.2 <- addPeakMatrix(proj)
save.image(paste0(image_dir, "atac_after_peak_matrix.RData"))
# Calculate differential accessible peaks for each cell type
differential_peaks_dir <- paste0(output_dir, "diff_peaks/")
if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir)}
for (cell_type in unique(proj.2$Cell_type_voting)) {
  print(cell_type)
  # Grab cells associated with cell type
  idxPass <- which(proj.2$Cell_type_voting %in% cell_type)
  print(length(idxPass))
  cellsPass <- proj.2$cellNames[idxPass]
  cells_subset <- subsetCells(proj.2, cellsPass)
  # Find DAPs
  marker_D28_D1 <- getMarkerFeatures(ArchRProj = cells_subset, useMatrix = "PeakMatrix", groupBy = "Conditions",
                                     testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), maxCells = 15000,
                                     useGroups = "HVL", bgdGroups = "LVL")
  # Grab relevant stats
  marker_log_2_fc <- assays(marker_D28_D1)$Log2FC
  marker_mean <- assays(marker_D28_D1)$Mean
  marker_fdr <- assays(marker_D28_D1)$FDR
  marker_pval <- assays(marker_D28_D1)$Pval
  marker_mean_diff <- assays(marker_D28_D1)$MeanDiff
  marker_auc <- assays(marker_D28_D1)$AUC
  marker_mean_bgd <- assays(marker_D28_D1)$MeanBGD
  # Print stats to Excel spreadsheet (used by MAGICAL)
  cell_type <- sub(" ", "_", cell_type)
  list_of_datasets <- list("Log2FC" = marker_log_2_fc, "Mean" = marker_mean, "FDR" = marker_fdr, "Pval" = marker_pval, 
                           "MeanDiff" = marker_mean_diff, "AUC" = marker_auc, "MeanBGD" = marker_mean_bgd)
  write.xlsx(list_of_datasets, file = paste0(differential_peaks_dir, cell_type, "_", "HVL_LVL_diff.xlsx"), colNames = FALSE)
}

# Create Peaks.txt file for MAGICAL
peaks <- getPeakSet(proj.2)
seq_names <- as.data.frame(peaks@seqnames)
ranges <- as.data.frame(peaks@ranges)
peak_txt_file <- cbind(seq_names, ranges)
colnames(peak_txt_file)[1] <- "chr"
peak_txt_file <- peak_txt_file[, !(names(peak_txt_file) %in% c("width", "names"))]
peak_txt_file[,1] <- sub("chr", "", peak_txt_file[,1])
peak_txt_file[,1] <- sub("X", "23", peak_txt_file[,1])
write.table(peak_txt_file, file = paste0(output_dir, "Peaks.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

# Create peak_motif_matches.txt file for MAGICAL
proj.2 <- addMotifAnnotations(ArchRProj = proj.2, motifSet = "cisbp", name = "Motif")
peak_motif_matches <- getMatches(proj.2, name = "Motif")
peak_motif_txt_file <- as.data.frame(peak_motif_matches@assays@data$matches)
# Remove _.* from ends of column names
for (col in 1:ncol(peak_motif_txt_file)){
  colnames(peak_motif_txt_file)[col] <-  sub("_.*", "", colnames(peak_motif_txt_file)[col])
}
peak_motif_txt_file <- cbind(peak_txt_file, peak_motif_txt_file)
colnames(peak_motif_txt_file)[2] <- "point1"
colnames(peak_motif_txt_file)[3] <- "point2"
cols <- sapply(peak_motif_txt_file, is.logical)
peak_motif_txt_file[,cols] <- lapply(peak_motif_txt_file[,cols], as.numeric)
write.table(peak_motif_txt_file, file = paste0(output_dir, "peak_motif_matches.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

# Create text file for each cell type containing pseudobulk counts for peaks
Cell_types <- unique(proj$Cell_type_voting)
sample.names <- unique(proj$Sample)
peak_count <- getMatrixFromProject(ArchRProj = proj.2, useMatrix = "PeakMatrix", useSeqnames = NULL, verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),logFile = createLogFile("getMatrixFromProject"))
pseudo_bulk_dir <- paste0(output_dir, "pseudo_bulk/")
if (!dir.exists(pseudo_bulk_dir)) {dir.create(pseudo_bulk_dir)}
for (i in c(1:length(Cell_types))){
  pseudo_bulk <- matrix(nrow = length(peaks), ncol = length(sample.names), 0)
  colnames(pseudo_bulk)<-sample.names
  rownames(pseudo_bulk)<-peaks@elementMetadata$idx
  
  for (s in c(1:length(sample.names))){
    idxMatch <- which(str_detect(peak_count$Cell_type_voting,Cell_types[i]) & str_detect(as.character(peak_count$Sample),sample.names[s]))
    if (length(idxMatch)>1){
      pseudo_bulk[,s] = Matrix::rowSums(peak_count@assays@data$PeakMatrix[,idxMatch])
    }
  }
  
  write.table(pseudo_bulk, file = paste0(pseudo_bulk_dir, "pseudo_bulk_ATAC_count_", Cell_types[i], ".txt"), quote = FALSE, sep = "\t", col.names = NA)
}









# Subset based on predictions we're fully confident in
idxPass <- which(proj$cellNames %in% curated_snRNA_seq_cells$cells) 
cellsPass <- proj$cellNames[idxPass]
proj <- proj[cellsPass, ]

curated_snRNA_seq_cells <- curated_snRNA_seq_cells[curated_snRNA_seq_cells$cells %in% proj$cellNames,]
curated_snRNA_seq_cells <- curated_snRNA_seq_cells[order(match(curated_snRNA_seq_cells$cells,proj$cellNames)),]
snRNA_seq_cell_votes <- curated_snRNA_seq_cells$voted_type

proj <- addCellColData(ArchRProj = proj, data = snRNA_seq_cell_votes, cells = proj$cellNames, name = "predictedGroup", force = TRUE)


# Use both predictions we're fully confident in and those we're less confident in
predicted_snATAC_cells <- data.frame(cell_name = proj$cellNames, voted_type = proj$predictedGroup)
for(current_row in 1:nrow(predicted_snATAC_cells)) {
  current_snATAC_cell <- predicted_snATAC_cells[current_row,]$cell_name
  if(current_snATAC_cell %in% curated_snRNA_seq_cells$cells) {
    current_voted_type <- curated_snRNA_seq_cells[curated_snRNA_seq_cells$cells == current_snATAC_cell,]$voted_type
  }
  predicted_snATAC_cells[current_row,]$voted_type <- current_voted_type
}

# Only keep cells which we have a voted type for (no N/A)
predicted_snATAC_cells <- predicted_snATAC_cells[predicted_snATAC_cells$voted_type != "N/A",]
kept_cells <- predicted_snATAC_cells$cell_name
proj <- proj[kept_cells, ]

# Add predicted type
proj <- addCellColData(ArchRProj = proj, data = predicted_snATAC_cells$voted_type, cells = proj$cellNames, name = "predictedGroup", force = TRUE)

