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

set.seed(1)

# Set project dir - used to organize different projects
project_dir <- "~/PLACEBO_FLU_1/"
if (!dir.exists(project_dir)) {dir.create(project_dir)}
# Location of scATAC data - data are organized by aliquot ID
# Note that aliquot ID and sample ID are different since multiple sample types
# (scRNA-seq, scATAC-seq) can come from the same aliquot
# However, in the context of analyzing samples from a single sample type,
# you can consider aliquots and samples basically equivalent
data_dir <- paste0(project_dir, "scATAC_seq_data/")
if (!dir.exists(data_dir)) {dir.create(data_dir)}
output_dir <- paste0(project_dir, "scATAC_seq_data_output/")
if (!dir.exists(output_dir)) {dir.create(output_dir)}
sample_metadata <- read.csv(paste0(project_dir, "current_sample_metadata_minus_8d5be1a4937a7ad3.csv"))
sample_assay_types <- read.csv(paste0(project_dir, "current_set_of_scRNA_and_scATAC_seq_samples.txt"), sep = "\t")
# Look in base scATAC dir to get list of all potential aliquots
# We will only use aliquots that are paired (D-1 and D28)
aliquot_list <- list.dirs(data_dir, recursive = FALSE)
aliquot_list <- strsplit(aliquot_list, "/")
aliquot_list <- unlist(lapply(aliquot_list, tail, n = 1L))
# D1.id stores aliquot IDs for day -1 samples
# D28.id stores aliquot IDs for day 28 samples
print("Grabbing Day -1 and D28 aliquot IDs from metadata file")
D1.id <- c()
D28.id <- c()
for (aliquot in aliquot_list) {
  # Grab current sample metadata, subject associated with sample, and then check to see whether subject has two samples
  # (D-1 and D28)
  current_sample = sample_metadata[sample_metadata$X_aliquot_id == aliquot,]
  current_subject = current_sample$SUBJECT_ID
  all_samples_associated_with_current_subject = sample_metadata[sample_metadata$SUBJECT_ID == current_subject ,]
  if (nrow(all_samples_associated_with_current_subject) == 2) {
    # Now, we grab our D-1 and D28 aliquot names.
    # If D-1 is not already in D1.id AND we have scRNA-seq data from both D-1 and D28 aliquots AND we have scRNA-seq data from both D-1 and D28 aliquots,  
    # then add D-1 to D1.id and D28 to D28.id
    d_negative_1_aliquot <- all_samples_associated_with_current_subject[all_samples_associated_with_current_subject$Time_Point == "D-1",]$X_aliquot_id
    d_28_aliquot <- all_samples_associated_with_current_subject[all_samples_associated_with_current_subject$Time_Point == "D28",]$X_aliquot_id
    d_negative_1_aliquot_sample_assays <- sample_assay_types[sample_assay_types$aliquot_name == d_negative_1_aliquot,]
    presence_of_d_negative_1_RNA <- d_negative_1_aliquot_sample_assays$scRNA_seq
    d_28_aliquot_sample_assays <- sample_assay_types[sample_assay_types$aliquot_name == d_28_aliquot,]
    presence_of_d_28_RNA <- d_28_aliquot_sample_assays$scRNA_seq
    if (d_negative_1_aliquot %in% D1.id == FALSE & d_negative_1_aliquot %in% aliquot_list & d_28_aliquot %in% aliquot_list & presence_of_d_negative_1_RNA == "Yes" & presence_of_d_28_RNA == "Yes") {
      D1.id <- append(D1.id, d_negative_1_aliquot)
      D28.id <- append(D28.id, d_28_aliquot)
    }
  }
}

# Read in input files and label each sample as D1 or D28 in metadata
inputFiles <- paste0(data_dir, c(D1.id, D28.id), "/outs/fragments.tsv.gz")
names(inputFiles) <- names(metadata) <- c(paste0("Sample_", 1:length(D1.id), "_D1"), paste0("Sample_", 1:length(D28.id), "_D28"))
metadata <- c(rep("D1", length(D1.id)), rep("D28", length(D28.id)))
names(metadata) <- c(paste0("Sample_", 1:length(D1.id), "_D1"), paste0("Sample_", 1:length(D28.id), "_D28"))

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
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)

# Filter out doublets (fake cells)
proj <- filterDoublets(ArchRProj = proj)

image_dir <- paste0(output_dir, "images/")
if (!dir.exists(image_dir)) {dir.create(image_dir)}
save.image(paste0(image_dir, "atac_after_filtering_doublets.RData"))

# List available matrices in project
getAvailableMatrices(proj)
#------------------------------------------------------------------------------------------------
# add condition metadata (D1 or D28) to cells in project
aa<-proj$Sample
idxSample <- which(str_detect(proj$Sample, "D1"))
aa[idxSample]<-'D1'

idxSample <- which(str_detect(proj$Sample, "D28"))
aa[idxSample]<-'D28'

proj <- addCellColData(ArchRProj = proj, data = aa, cells = proj$cellNames,name = "Conditions", force = TRUE)

# Plot out TSS Enrichment / Doublet Enrichment / Nucleosome Ratio for each sample
p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
p2 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "DoubletEnrichment", plotAs = "ridges")
p3 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges")

plotPDF(p1,p2,p3, name = "Integrated_Scores_Prefiltering.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)

# Filter out cells that don't meet TSS enrichment / double enrichment / nucleosome ratio criteria
idxPass <- which(proj$TSSEnrichment >= 6 & proj$NucleosomeRatio < 2 & proj$DoubletEnrichment < 5) 
cellsPass <- proj$cellNames[idxPass]
proj<-proj[cellsPass, ]
# List number of D1 and D28 cells and list number of cells remaining for each sample
table(proj$Conditions)
table(proj$Sample)


# Perform dimensionality reduction on cells (addIterativeLSI), create UMAP embedding (addUMAP), 
# and add cluster information (addClusters)
addArchRThreads(threads = 8)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                        clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                        varFeatures = 15000, dims = 2:30)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 4, knnAssign = 30, maxClusters = NULL, force = TRUE)

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

plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

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

# Integrate scATAC-seq and scRNA-seq data
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
save.image(paste0(image_dir, "atac_after_gene_integration_matrix.RData"))

# Plot integrated data with cell type predictions 
pal <- paletteDiscrete(values = proj$predictedGroup)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE)

plotPDF(p1, name = "Integrated_annotated.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

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

proj <- addCellColData(ArchRProj = proj, data = Cell_type_combined, cells = proj$cellNames, name = "Cell_type_combined", force = TRUE)

# Plot integrated data with cell type predictions (combined)
pal <- paletteDiscrete(values = proj$Cell_type_combined)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Cell_type_combined", embedding = "UMAP", pal = pal, force = TRUE)

plotPDF(p1, name = "Integrated_annotated_combined.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


#####remove cells using TSS
idxPass <- which(proj$TSSEnrichment >= 15) 
cellsPass <- proj$cellNames[idxPass]
proj.filtered<-proj[cellsPass, ]

table(proj.filtered$Conditions)
table(proj.filtered$Sample)
table(proj.filtered$Cell_type_combined)

# Redo dimensionality reduction, UMAP, clustering
addArchRThreads(threads = 8)
proj.filtered <- addIterativeLSI(ArchRProj = proj.filtered, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                        clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                        varFeatures = 15000, dims = 2:30)
proj.filtered <- addUMAP(ArchRProj = proj.filtered, reducedDims = "IterativeLSI", force = TRUE)
proj.filtered <- addClusters(input = proj.filtered, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 3, knnAssign = 30, maxClusters = NULL, force = TRUE)

save.image(paste0(image_dir, "atac_after_lsi_umap_clusters_2.RData"))

# Recalculate Rand index (see above)
adjustedRandIndex(proj.filtered$Sample, proj.filtered$Clusters)
adjustedRandIndex(proj.filtered$Conditions, proj.filtered$Clusters)

# Plot integrated data with cell type predictions (combined) after TSS filtering
p1 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Cell_type_combined", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)

plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_Filtered_TSS15.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Get cell_type_voting info and add to ArchR project
cM <- as.matrix(confusionMatrix(proj.filtered$Clusters, proj.filtered$Cell_type_combined))
pre_cluster <- rownames(cM)
max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]

Cell_type_voting <- proj.filtered$Clusters
for (m in c(1:length(pre_cluster))){
  idxSample <- which(proj.filtered$Clusters == pre_cluster[m])
  Cell_type_voting[idxSample] <- max_celltype[m]
}


proj.filtered <- addCellColData(ArchRProj = proj.filtered, data = Cell_type_voting, cells = proj.filtered$cellNames, name = "Cell_type_voting", force = TRUE)
save.image(paste0(output_dir, "atac_after_cell_type_voting.RData"))

# See how clusters are distributed
cluster_predictions <- vector()
cluster_distributions <- list()
idx <- 1
for (cluster in unique(proj.filtered$Clusters)) {
  idxPass <- which(proj.filtered$Clusters %in% cluster)
  cellsPass <- proj.filtered$cellNames[idxPass]
  filtered_cluster <-proj.filtered[cellsPass,]
  cluster_predictions <- append(cluster_predictions, table(filtered_cluster$Cell_type_voting))
  cluster_distributions[[idx]] <- table(filtered_cluster$Cell_type_combined)
  idx <- idx + 1
}
names(cluster_predictions) <- paste(unique(proj.filtered$Clusters), "-", names(cluster_predictions))

#Remove the messy clusters (determined through visual inspection and seeing distribution of cells in each cluster)
idxPass <- which(proj.filtered$Clusters %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C12", "C16", "C25", "C30", "C47"))
cellsPass <- proj.filtered$cellNames[-idxPass]
proj.filtered<-proj.filtered[cellsPass, ]

# Combine more cell types
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "CD4 Naive", "T Naive")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "CD8 Naive", "T Naive")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "NK_CD56bright", "NK")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "ASDC", "CD14 Mono")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "cDC", "CD14 Mono")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "Eryth", "CD14 Mono")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "HSPC", "CD14 Mono")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "pDC", "CD14 Mono")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "Plasmablast", "CD14 Mono")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "Platelet", "CD14 Mono")
proj.filtered$Cell_type_combined <- replace(proj.filtered$Cell_type_combined, proj.filtered$Cell_type_combined == "Treg", "T Naive")

table(proj.filtered$Conditions)
table(proj.filtered$Sample)
table(proj.filtered$Cell_type_combined)

# Print final integrated data with cell type predictions
p1 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Cell_type_combined", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)


plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_Filtered_TSS15_Final.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Apparently info about ArchR Genome is not stored in the .RData file, so if we load an .RData file above,
# We need to re-add our hg38 genome
addArchRGenome("hg38")

# Peak calling - first, call addGroupCoverages to find pseudo-bulk replicates, then call peaks using MACS2
proj.filtered <- addGroupCoverages(ArchRProj = proj.filtered, groupBy = "Cell_type_combined")
pathToMacs2 <- findMacs2()
proj.filtered <- addReproduciblePeakSet(
  ArchRProj = proj.filtered, 
  groupBy = "Cell_type_combined", 
  pathToMacs2 = pathToMacs2
)
# Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
# calculations
proj.filtered.2 <- addPeakMatrix(proj.filtered)
save.image(paste0(image_dir, "atac_after_peak_matrix.RData"))
# Calculate differential accessible peaks for each cell type
differential_peaks_dir <- paste0(output_dir, "diff_peaks/")
if (!dir.exists(differential_peaks_dir)) {dir.create(differential_peaks_dir)}
for (cell_type in unique(proj.filtered.2$Cell_type_combined)) {
  print(cell_type)
  # Grab cells associated with cell type
  idxPass <- which(proj.filtered.2$Cell_type_combined %in% cell_type)
  print(length(idxPass))
  cellsPass <- proj.filtered.2$cellNames[idxPass]
  cells_subset <- subsetCells(proj.filtered.2, cellsPass)
  # Find DAPs
  marker_D28_D1 <- getMarkerFeatures(ArchRProj = cells_subset, useMatrix = "PeakMatrix", groupBy = "Conditions",
                                     testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), maxCells = 15000,
                                     useGroups = "D28", bgdGroups = "D1")
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
  write.xlsx(list_of_datasets, file = paste0(differential_peaks_dir, cell_type, "_", "D28_D1_diff.xlsx"), colNames = FALSE)
}

# Create Peaks.txt file for MAGICAL
peaks <- getPeakSet(proj.filtered.2)
seq_names <- as.data.frame(peaks@seqnames)
ranges <- as.data.frame(peaks@ranges)
peak_txt_file <- cbind(seq_names, ranges)
colnames(peak_txt_file)[1] <- "chr"
peak_txt_file <- peak_txt_file[, !(names(peak_txt_file) %in% c("width", "names"))]
peak_txt_file[,1] <- sub("chr", "", peak_txt_file[,1])
peak_txt_file[,1] <- sub("X", "23", peak_txt_file[,1])
write.table(peak_txt_file, file = paste0(output_dir, "Peaks.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

# Create peak_motif_matches.txt file for MAGICAL
proj.filtered.2 <- addMotifAnnotations(ArchRProj = proj.filtered.2, motifSet = "cisbp", name = "Motif")
peak_motif_matches <- getMatches(proj.filtered.2, name = "Motif")
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
Cell_types <- unique(proj.filtered$Cell_type_voting)
sample.names <- unique(proj.filtered$Sample)
peak_count <- getMatrixFromProject(ArchRProj = proj.filtered.2, useMatrix = "PeakMatrix", useSeqnames = NULL, verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),logFile = createLogFile("getMatrixFromProject"))
pseudo_bulk_dir <- paste0(output_dir, "pseudo_bulk/")
if (!dir.exists(pseudo_bulk_dir)) {dir.create(pseudo_bulk_dir)}
for (i in c(1:length(Cell_types))){
  pseudo_bulk <- matrix(nrow = length(peaks), ncol = length(sample.names), 0)
  colnames(pseudo_bulk)<-sample.names
  rownames(pseudo_bulk)<-peaks@elementMetadata$idx
  
  for (s in c(1:length(sample.names))){
    idxMatch <- which(str_detect(peak_count$Cell_type_voting,Cell_types[i]) & str_detect(peak_count$Sample,sample.names[s]))
    if (length(idxMatch)>1){
      pseudo_bulk[,s] = Matrix::rowSums(peak_count@assays@data$PeakMatrix[,idxMatch])
    }
  }
  
  write.table(pseudo_bulk, file = paste0(pseudo_bulk_dir, "pseudo_bulk_ATAC_count_", Cell_types[i], ".txt"), quote = FALSE, sep = "\t", col.names = NA)
}

