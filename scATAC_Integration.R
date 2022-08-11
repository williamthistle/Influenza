library(ArchR)
library(stringr)
library(pheatmap)
library(mclust)
library(writexl)
library(ggplot2)

set.seed(1)

D1.id <- c()
D28.id <- c()

inputFiles <- paste0("path to folder/", c(D1.id, D28.id), "/outs/fragments.tsv.gz")
names(inputFiles) <- names(metadata) <- c(paste0("Sample_", 1:6, "_D1"), paste0("Sample_", 1:6, "_D28"))

metadata <- c(rep("D1", length(D1.id)), rep("D28", length(D28.id)))
names(metadata) <- c(paste0("Sample_", 1:6, "_D1"), paste0("Sample_", 1:6, "_D28"))

addArchRGenome("hg38")


#fragment data processing 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 3000,
  maxFrags = 30000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#fake cell check
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

#key step for sample selection before integration
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "output path/ArchR/",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)


proj <- filterDoublets(ArchRProj = proj)

getAvailableMatrices(proj)
#------------------------------------------------------------------------------------------------
# add conditions
aa<-proj$Sample
idxSample <- which(str_detect(proj$Sample, "D1"))
aa[idxSample]<-'D1'

idxSample <- which(str_detect(proj$Sample, "D28"))
aa[idxSample]<-'D28'

proj <- addCellColData(ArchRProj = proj, data = aa, cells = proj$cellNames,name = "Conditions", force = TRUE)


p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
p2 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "DoubletEnrichment", plotAs = "ridges")
p3 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges")

plotPDF(p1,p2,p3, name = "Integrated_Scores_Prefiltering.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)


idxPass <- which(proj$TSSEnrichment >= 6 & proj$NucleosomeRatio < 2 & proj$DoubletEnrichment < 5) 
cellsPass <- proj$cellNames[idxPass]
proj<-proj[cellsPass, ]
table(proj$Conditions)
table(proj$Sample)



addArchRThreads(threads = 8)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                        clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                        varFeatures = 15000, dims = 2:30)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 4, knnAssign = 30, maxClusters = NULL, force = TRUE)

adjustedRandIndex(proj$Sample, proj$Clusters)
adjustedRandIndex(proj$Conditions, proj$Clusters)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)

plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

scRNA <- LoadH5Seurat("path to reference/multi.h5seurat")

idx <- which(scRNA$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
scRNA <- scRNA[,-idx]

idx <- which(scRNA$celltype.l3 == "Treg Naive")
scRNA <- scRNA[,-idx]

idx <- which(scRNA$celltype.l2 %in% c("CD4 Proliferating", "CD8 Proliferating"))
scRNA$celltype.l2[idx] <- "T Proliferating"

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

#Remove the messy cluster
idxPass <- which(proj.filtered$Clusters %in% c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11")) 
cellsPass <- proj.filtered$cellNames[-idxPass]
proj.filtered<-proj.filtered[cellsPass, ]

table(proj.filtered$Conditions)
table(proj.filtered$Sample)
table(proj.filtered$Cell_type_combined)









