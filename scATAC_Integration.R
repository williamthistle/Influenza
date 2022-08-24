library(ArchR)
library(stringr)
library(pheatmap)
library(mclust)
library(writexl)
library(ggplot2)
library(hexbin)
library(SeuratDisk)

set.seed(1)

# Location of scATAC data - data are organized by aliquot ID
# Note that aliquot ID and sample ID are different since multiple sample types
# (scRNA-seq, scATAC-seq) can come from the same aliquot
# However, in the context of analyzing samples from a single sample type,
# you can consider aliquots and samples basically equivalent
base_scATAC_dir <- "~/scATAC_seq_data/"
output_dir <- "~/scATAC_seq_data_output/"
sample_metadata <- read.csv("~/current_sample_metadata_minus_8d5be1a4937a7ad3.csv")
sample_assay_types <- read.csv("~/current_set_of_scRNA_and_scATAC_seq_samples.txt", sep = "\t")
# Look in base scATAC dir to get list of all potential aliquots
# We will only use aliquots that are paired (D-1 and D28)
aliquot_list <- list.dirs(base_scATAC_dir, recursive = FALSE)
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


# TODO: Make 1:6 generic
# Read in input files and label each sample as D1 or D28 in metadata
inputFiles <- paste0(base_scATAC_dir, c(D1.id, D28.id), "/outs/fragments.tsv.gz")
names(inputFiles) <- names(metadata) <- c(paste0("Sample_", 1:6, "_D1"), paste0("Sample_", 1:6, "_D28"))
metadata <- c(rep("D1", length(D1.id)), rep("D28", length(D28.id)))
names(metadata) <- c(paste0("Sample_", 1:6, "_D1"), paste0("Sample_", 1:6, "_D28"))

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

#fake cell check (check for doublets)
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

save.image(paste0(output_dir, "atac_after_filtering_doublets.RData"))

# List available matrices in project
getAvailableMatrices(proj)
#------------------------------------------------------------------------------------------------
# add condition metadata (D1 and D28) to cells in project
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

save.image(paste0(output_dir, "atac_after_lsi_umap_clusters_1.RData"))

# Determine whether there is similarity in how cells are clustered by sample or condition and by the clusters
# determined above (0 is random, 1 is perfect)
adjustedRandIndex(proj$Sample, proj$Clusters)
adjustedRandIndex(proj$Conditions, proj$Clusters)

# UMAP plots colored by condition, sample, cluster ID, and TSS enrichment
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)

plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Load reference scRNA data
scRNA <- LoadH5Seurat("reference/multi.h5seurat")

# Remove certain cell types we're not interested in
idx <- which(scRNA$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
scRNA <- scRNA[,-idx]

idx <- which(scRNA$celltype.l3 == "Treg Naive")
scRNA <- scRNA[,-idx]

# Rename CD4 Proliferating and CD8 Proliferating to T Proliferating
idx <- which(scRNA$celltype.l2 %in% c("CD4 Proliferating", "CD8 Proliferating"))
scRNA$celltype.l2[idx] <- "T Proliferating"
# Step required in new version of ArchR to make addGeneIntegrationMatrix run successfully
#scRNA <- RenameAssays(scRNA, SCT = "RNA")
save.image(paste0(output_dir, "atac_before_gene_integration_matrix.RData"))
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
save.image(paste0(output_dir, "atac_after_gene_integration_matrix.RData"))

pal <- paletteDiscrete(values = proj$predictedGroup)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE)

plotPDF(p1, name = "Integrated_annotated.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

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

addArchRThreads(threads = 8)
proj.filtered <- addIterativeLSI(ArchRProj = proj.filtered, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                        clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                        varFeatures = 15000, dims = 2:30)
proj.filtered <- addUMAP(ArchRProj = proj.filtered, reducedDims = "IterativeLSI", force = TRUE)
proj.filtered <- addClusters(input = proj.filtered, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 3, knnAssign = 30, maxClusters = NULL, force = TRUE)

save.image(paste0(output_dir, "atac_after_lsi_umap_clusters_2.RData"))

adjustedRandIndex(proj.filtered$Sample, proj.filtered$Clusters)
adjustedRandIndex(proj.filtered$Conditions, proj.filtered$Clusters)

p1 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Cell_type_combined", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)

plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_Filtered_TSS15.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

cM <- as.matrix(confusionMatrix(proj.filtered$Clusters, proj.filtered$Cell_type_combined))
pre_cluster <- rownames(cM)
max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]

Cell_type_voting = proj.filtered$Clusters
for (m in c(1:length(pre_cluster))){
  idxSample <- which(proj.filtered$Clusters == pre_cluster[m])
  Cell_type_voting[idxSample] <- max_celltype[m]
}


proj.filtered <- addCellColData(ArchRProj = proj.filtered, data = Cell_type_voting, cells = proj.filtered$cellNames, name = "Cell_type_voting", force = TRUE)
save.image(paste0(output_dir, "atac_after_cell_type_voting.RData"))
saveArchRProject(ArchRProj = proj.filtered, paste0(output_dir, "ArchR_Filtered/"))

#Remove the messy clusters
idxPass <- which(proj.filtered$Clusters %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C12", "C16", "C24", "C25", "C30", "C36", "C37", "C45", "C47"))
cellsPass <- proj.filtered$cellNames[-idxPass]
proj.filtered<-proj.filtered[cellsPass, ]

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

p1 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "Cell_type_combined", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj.filtered, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)


plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_Filtered_TSS15_Final.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
#save.image(paste0(output_dir, "atac_final.RData"))






