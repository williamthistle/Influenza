setwd("~/Desktop/ECHO/Staph/")

library(ArchR)
library(stringr)
library(pheatmap)
library(mclust)
library(Seurat)
library(SeuratDisk)
library(writexl)
library(ggplot2)

set.seed(1)

sample_meta<-read.table('MSSA_MRSA_sample_selected.txt')
inputFiles <- as.character(sample_meta[,4])
names(inputFiles) <- sample_meta[,3]

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
  outputDirectory = "Integration",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj <- filterDoublets(ArchRProj = proj)

getAvailableMatrices(proj)

save.image(file='scATAC_TSS_control_MSSA_MRSA_4_min3000_max30000.RData')
# add conditions
aa<-proj$Sample
idxSample <- which(str_detect(proj$Sample, "MSSA"))
aa[idxSample]<-'MSSA'

idxSample <- which(str_detect(proj$Sample, "MRSA"))
aa[idxSample]<-'MRSA'

idxSample <- which(str_detect(proj$Sample, "HI_control"))
aa[idxSample]<-'Control'

idxSample <- which(str_detect(proj$Sample, "ANTH_control"))
aa[idxSample]<-'Control'

proj <- addCellColData(ArchRProj = proj, data = aa, cells = proj$cellNames,name = "Conditions", force = TRUE)


#************************************************************************************************************
#************************************************************************************************************
#************************************************************************************************************


#idxSample <- which(proj$Conditions == "MRSA")
#cellsPass <- proj$cellNames[idxSample]
#proj<-proj[idxSample, ]

#filter cells with low TSS enrichment

idxPass <- which(proj$TSSEnrichment >= 12 & proj$NucleosomeRatio < 2) #12K cells > TSS4; 90K cells > TSS10
cellsPass <- proj$cellNames[idxPass]
proj<-proj[cellsPass, ]

p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
p1


#work on all
addArchRThreads(threads = 8)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                        clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                        varFeatures = 20000, dims = 1:30)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 2, knnAssign = 30, maxClusters = NULL, force = TRUE)

adjustedRandIndex(proj$Sample, proj$Clusters)
adjustedRandIndex(proj$Conditions, proj$Clusters)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)

ggAlignPlots(p1, p2, p3,type = "h")

scRNA <- LoadH5Seurat("pbmc_multimodal.h5seurat")

addArchRThreads(threads = 4)
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

pal <- paletteDiscrete(values = scRNA@meta.data$celltype.l2)
p4 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE)


#cell types from the scRNA-seq data are most abundant in each of our scATAC-seq clusters
cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup))
pre_cluster <- rownames(cM)
max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]

Cell_type_voting = proj$Clusters
for (m in c(1:length(pre_cluster))){
  idxSample <- which(proj$Clusters == pre_cluster[m])
  Cell_type_voting[idxSample] <- max_celltype[m]
}

idxSample <- which(proj$Clusters == 'C40' | proj$Clusters == 'C41')
Cell_type_voting[idxSample] <- 'CD14 Mono'

idxSample <- which(proj$Clusters == 'C12')
Cell_type_voting[idxSample] <- 'B memory'

idxSample <- which(proj$Clusters == 'C14')
Cell_type_voting[idxSample] <- 'B memory'

proj <- addCellColData(ArchRProj = proj, data = Cell_type_voting, cells = proj$cellNames, name = "Cell_type_voting", force = TRUE)

cM <- confusionMatrix(paste0(proj$Cell_type_voting), paste0(proj$Sample))
cM <- t(t(cM) / Matrix::rowSums(t(cM)))
#cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), border_color = "black")
#write.csv(as.matrix(cM ), file = 'All_sample_cell_type_proportion_20000peak_TSS12_Nu2_reso2_L2_manual.csv')
p5 <- plotEmbedding(proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE, size = 0.01, keepAxis = TRUE)
p5

idxPass <- which(proj$Clusters != 'C13') 
cellsPass <- proj$cellNames[idxPass]
proj<-proj[cellsPass, ]

p6 <- plotEmbedding(proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", rastr = FALSE, force = TRUE)
p7 <- plotEmbedding(proj, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", rastr = FALSE, force = TRUE)
ggAlignPlots(p6, p7, type = "h")


idxPass <- which(str_detect(proj$Conditions,'Control'))
cellsPass <- proj$cellNames[idxPass]
proj_control<-proj[cellsPass, ]

idxPass <- which(str_detect(proj$Conditions,'MSSA'))
cellsPass <- proj$cellNames[idxPass]
proj_MSSA<-proj[cellsPass, ]

idxPass <- which(str_detect(proj$Conditions,'MRSA'))
cellsPass <- proj$cellNames[idxPass]
proj_MRSA<-proj[cellsPass, ]

p8 <- plotEmbedding(proj_control, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)
p9 <- plotEmbedding(proj_MSSA, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)
p10 <- plotEmbedding(proj_MRSA, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)

ggAlignPlots(p8, p9, p10, type = "h")
rm(scRNA)

markersGS <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "Cell_type_voting",bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")

markerGenes  <- c(
#  'MS4A1',	'TNFRSF13B',			'AIM2',	'CD79A',	'LINC01857',	'RALGPS2',	'BANK1',	'CD79B',#	'B intermediate'
  'MS4A1',	'COCH',	'AIM2',	'BANK1',	'SSPN',	'CD79A',	'TEX9',	'RALGPS2',	'TNFRSF13C',	'LINC01781',#	'B memory'
  		'CD79A',	'IL4R',	'MS4A1',	'CXCR4',	'BTG1',	'TCL1A',	'CD79B',	'YBX3',#	'B naive'
#  'IGHA2',	'MZB1',	'TNFRSF17',	'DERL3',	'TXNDC5',	'TNFRSF13B',	'POU2AF1',	'CPNE5',	'HRASLS2',	'NT5DC2',#	'Plasmablast'
#  'GZMH',	'CD4',	'FGFBP2',	'ITGB1',	'GZMA',	'CST7',	'GNLY',	'B2M',	'IL32',	'NKG7',#	'CD4 CTL'
  'TCF7',	'CD4',	'CCR7',	'IL7R',	'FHIT',	'LEF1',	'MAL',	'NOSIP',	'LDHB',	'PIK3IP1',#	'CD4 Naive'
#  'MKI67',	'TOP2A',	'PCLAF',	'CENPF',	'TYMS',	'NUSAP1',	'ASPM',	'PTTG1',	'TPX2',	'RRM2',#	'CD4 Proliferating'
  'IL7R',	'TMSB10',	'CD4',	'ITGB1',	'LTB',		'AQP3',	'LDHB',	'IL32',	'MAL',#	'CD4 TCM'
  'IL7R',	'CCL5',	'FYB1',	'GZMK',	'IL32',	'GZMA',	'KLRB1',		'LTB',	'AQP3',#	'CD4 TEM'
  'RTKN2',	'FOXP3',		'CD4',	'IL2RA',	'TIGIT',	'CTLA4',	'FCRL3',	'LAIR2',	'IKZF2',#	'Treg'
  'CD8B',	'S100B',	'CCR7',	'RGS10',	'NOSIP',		'LEF1',	'CRTAM',	'CD8A',	'OXNAD1',#	'CD8 Naive'
#  'MKI67',	'CD8B',	'TYMS',		'PCLAF',	'CD3D',	'CLSPN',	'CD3G',	'TK1',	'RRM2',#	'CD8 Proliferating'
  'CD8B',	'ANXA1',	'CD8A',	'KRT1',		'YBX3',	'IL7R',		'NELL2',	'LDHB',#	'CD8 TCM'
  'CCL5',	'GZMH',	'CD8A',		'KLRD1',	'NKG7',	'GZMK',	'CST7',	'CD8B',	#	'CD8 TEM'
#  'PPP1R14A',	'LILRA4',	'AXL',		'SCT',	'SCN9A','LGMN',	'DNASE1L3',	'CLEC4C',	'GAS6',#	'ASDC'
  'CLEC9A',	'DNASE1L3',	'C1orf54',	'IDO1',	'CLNK',	'CADM1',	'FLT3',	'ENPP1',	'XCR1',	'NDRG2',#	'cDC1'
  'FCER1A',	'HLA-DQA1',	'CLEC10A',	'CD1C',	'ENHO',	'PLD4',	'GSN',	'SLC38A1',	'NDRG2',	'AFF3',#	'cDC2'
  'ITM2C',	'PLD4',	'SERPINF1',	'LILRA4',		'TPM2',	'MZB1',	'SPIB',	'IRF4',	'SMPD3',#	'pDC'
  'S100A9',	'CTSS',	'S100A8',	'LYZ',	'VCAN',	'S100A12',	'IL1B',	'CD14',	'G0S2',	'FCN1',#	'CD14 Mono'
  'CDKN1C',	'FCGR3A',	'PTPRC',	'LST1',	'IER5',	'MS4A7',	'RHOC',	'IFITM3',	'AIF1',	'HES4',#	'CD16 Mono'
  'GNLY',	'TYROBP',	'NKG7',	'FCER1G',	'GZMB',		'PRF1',	'FGFBP2',	'SPON2',	'KLRF1',#	'NK'
  'MKI67',	'KLRF1',	'TYMS',		'TOP2A',	'FCER1G',	'PCLAF',	'CD247',	'CLSPN',	'ASPM',#	'NK Proliferating'
  'XCL2',	'FCER1G',	'SPINK2',		'KLRC1',	'XCL1',	'SPTSSB',	'PPP1R9A',	'NCAM1',	'TNFRSF11A',#	'NK_CD56bright'
#  'HBD',	'HBM',	'AHSP',	'ALAS2',	'CA1',	'SLC4A1',	'IFIT1B',	'TRIM58',	'SELENBP1',	'TMCC2',#	'Eryth'
#  'SPINK2',	'PRSS57',	'CYTL1',	'EGFL7',	'GATA2',	'CD34',	'SMIM24',	'AVP',	'MYB',	'LAPTM4B',#	'HSPC'
#  'KIT',		'TTLL10',	'LINC01229',	'SOX4',	'KLRB1',	'TNFRSF18',	'TNFRSF4',	'IL1R1',	'HPGDS',#	'ILC'
  'PPBP',	'PF4',	'NRGN',	'GNG11',	'CAVIN2',	'TUBB1',	'CLU',	'RGS18',	'GP9',#	'Platelet'
#  'PTPN3',	'MIR4422HG',	'NUCB2',	'CAV1',	'DTHD1',	'GZMA',	'MYB',	'FXYD2',	'GZMK',	'AC004585.1',#	'dnT'
#  	'TRGC1',		'KLRC1',	'NKG7',	'TRDV2',	'CD7',	'TRGV9',	'KLRD1',	'KLRG1',#	'gdT'
  'KLRB1',	'NKG7',	'GZMK',	'IL7R',	'SLC4A10',	'GZMA',	'CXCR6',	'PRSS35',	'RBM24',	'NCR3'#	'MAIT'
)

heatmapGS <- markerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1", labelMarkers = markerGenes,transpose = TRUE)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
markergene_UMAP <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAP", quantCut = c(0.05, 0.95), imputeWeights = NULL, pal = paletteContinuous(set = "whiteBlue"))
markergene_UMAP$NKG7  

#peak calling
pathToMacs2<-"/Library/Frameworks/Python.framework/Versions/2.7/bin/macs2"

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Cell_type_voting", minCells = 50, maxCells = 500,
                          minReplicates = 3, maxReplicates = 32, force = TRUE)
proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "Cell_type_voting", pathToMacs2 = pathToMacs2, cutOff = 0.05, force = TRUE)
proj <- addPeakMatrix(proj, force = TRUE)

PromoterCount<-getMatrixFromProject(ArchRProj = proj, useMatrix = "GeneScoreMatrix", useSeqnames = NULL,
                                    verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),logFile = createLogFile("getMatrixFromProject"))
saveRDS(PromoterCount, file = "scATAC_promoter_cell_readcount_0513.rds")
rm(PromoterCount)

PeakCount<-getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix", useSeqnames = NULL,
                                     verbose = TRUE,binarize = FALSE,threads = getArchRThreads(),logFile = createLogFile("getMatrixFromProject"))
saveRDS(PeakCount, file = "scATAC_peak_cell_readcount_0513.rds")


Peaks <- getPeakSet(proj)

write_xlsx(as.data.frame(Peaks@seqnames), "1peak_chr.xlsx")
write_xlsx(as.data.frame(Peaks@ranges), "1peak_start_end.xlsx")
write_xlsx(as.data.frame(Peaks@elementMetadata), "1peak_meta.xlsx")


library(Matrix)
Cell_types <- unique(proj$Cell_type_voting)

for (i in c(1:length(Cell_types))){
  pseudo_bulk = matrix(nrow = length(Peaks), ncol = length(sample_meta[,3]), 0)
  colnames(pseudo_bulk)<-sample_meta[,3]
  rownames(pseudo_bulk)<-Peaks@elementMetadata$idx
  
  for (s in c(1:length(sample_meta[,3]))){
    idxPass <- which(str_detect(PeakCount$Cell_type_voting,Cell_types[i]) & str_detect(PeakCount$Sample,sample_meta[s,3]))
    if (length(idxPass)>1){
      pseudo_bulk[,s] = Matrix::rowSums(PeakCount@assays@data$PeakMatrix[,idxPass])
    }
  }
  write.table(pseudo_bulk, file = paste('pseudo_bulk_ATAC_count_', Cell_types[i], '.txt', sep=''), quote = FALSE, sep = "\t")
}


write_xlsx(as.data.frame(aa@seqnames), "1peak_chr.xlsx")
write_xlsx(as.data.frame(aa@ranges), "1peak_start_end.xlsx")
write_xlsx(as.data.frame(Peaks@elementMetadata), "1peak_meta.xlsx")

i=4
aa=colSums(PeakCount@assays@data$PeakMatrix)
for (i in c(1:length(Cell_types))){
  idxPass_1 <- which(str_detect(PeakCount$Cell_type_voting,Cell_types[i]) & str_detect(PeakCount$Conditions,'MSSA'))
  idxPass_2 <- which(str_detect(PeakCount$Cell_type_voting,Cell_types[i]) & str_detect(PeakCount$Conditions,'Control'))
  mean_MRSA=rowMeans(t(t(PeakCount@assays@data$PeakMatrix[, idxPass_1]*4e4)/aa[idxPass_1]))
  mean_MSSA=rowMeans(t(t(PeakCount@assays@data$PeakMatrix[, idxPass_2]*4e4)/aa[idxPass_2]))
  write.table(mean_MRSA, file = paste('MSSA_mean_', Cell_types[i], '.txt', sep=''), quote = FALSE, sep = "\t")
  write.table(mean_MSSA, file = paste('Control_mean_', Cell_types[i], '.txt', sep=''), quote = FALSE, sep = "\t")
}
mean(PeakCount@assays@data$PeakMatrix[10000, idxPass_1])
mean(PeakCount@assays@data$PeakMatrix[10000, idxPass_2])
library(chromVARmotifs)
library(dplyr)
library(tidyr)
library(data.table)

markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Cell_type_voting", bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",maxCells = 3000)
write_xlsx(as.data.frame(markersPeaks@elementMetadata@listData), paste("Differential_peaks/Cell_type_peaks.xlsx"))
write_xlsx(markersPeaks@assays@data@listData, paste("Differential_peaks/Cell_type_diff.xlsx", sep=''), col_names = FALSE)


proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE, returnMatrix = FALSE)
#pheatmap::pheatmap(mat = as.matrix(heatmapEM), color = paletteContinuous("whiteBlue"), border_color = "black")
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "left")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation = "Motif",force = TRUE)


motifPositions <- getPositions(proj)

motifs <- c("JUN_143")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs
markermotif_UMAP <- plotEmbedding(ArchRProj = proj, colorBy = "MotifMatrix", name = markerMotifs, embedding = "UMAP", imputeWeights = NULL, pal = paletteContinuous(set = "whiteBlue"))
markermotif_UMAP

motifs <- c("JUN_143")
seFoot <- getFootprints(ArchRProj = proj, positions = motifPositions[motifs], groupBy = "Cell_type_voting")
plotFootprints(seFoot = seFoot, ArchRProj = proj, normMethod = "None", plotName = "Footprints",addDOC = FALSE,smoothWindow = 5)



#**************************************************
setwd("~/Desktop/pilot/scATAC_Xi/Staph")
library(ArchR)
library(stringr)
library(writexl)
proj<-loadArchRProject(path = "Save-MRSA-0518", force = FALSE, showLogo = TRUE)

selected_cell_type <- unique(proj$Cell_type_voting)
table(proj$Cell_type_voting) 

Community_samples=c('01S0003453', '01S0003462', '01S0003464', '01S0003507', '01S0004017', '01S0003549', '01S0003989', '01S0003994', '01S0004011', '01S0004013');
Healthcare_samples=c('01S0003431', '01S0003482', '01S0004015', '01S0003466', '01S0003492', '01S0003509', '01S0003515', '01S0003527', '01S0003542', '01S0003978', '01S0003987','01S0003992');

aa<-proj$Sample
for (s in c(1:length(Community_samples))){
  idxSample <- which(str_detect(proj$Sample, Community_samples[s]))
  aa[idxSample]<-'Community'
}

for (s in c(1:length(Healthcare_samples))){
  idxSample <- which(str_detect(proj$Sample, Healthcare_samples[s]))
  aa[idxSample]<-'Healthcare'
}

idxSample <- which(str_detect(proj$Sample, "control"))
aa[idxSample]<-'Control'

table(aa)

proj <- addCellColData(ArchRProj = proj, data = aa, cells = proj$cellNames,name = "Route", force = TRUE)


for (j in c(1:length(selected_cell_type))){
  dir.create(file.path('MRSA_Differential_peaks', selected_cell_type[j]), showWarnings = FALSE)
  idxPass <- which(str_detect(proj$Cell_type_voting, selected_cell_type[j]) & str_detect(proj$Conditions, 'MSSA'))
  cellsPass <- proj$cellNames[idxPass]
  proj_temp<-proj[cellsPass, ]
  
  ## Staph vs. control
  marker_MRSA_Control <- getMarkerFeatures(ArchRProj = proj_temp, useMatrix = "PeakMatrix", groupBy = "Route",
                                           testMethod = "wilcoxon", bias = c("log10(nFrags)"), normBy = "ReadsInPeaks", maxCells = 15000,
                                           useGroups = "Healthcare", bgdGroups = "Community")
  p_MRSA_Control <- plotMarkers(seMarker = marker_MRSA_Control, name = "Healthcare", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
  plotPDF(p_MRSA_Control, name = paste(selected_cell_type[j], "MSSA_Healthcare_Community_Volcano_V2", sep='_'), width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
  write_xlsx(as.data.frame(marker_MRSA_Control@elementMetadata@listData), paste("MRSA_Differential_peaks/", selected_cell_type[j], "/MSSA_Healthcare_Community_peaks.xlsx", sep=''))
  write_xlsx(marker_MRSA_Control@assays@data@listData, paste("MRSA_Differential_peaks/", selected_cell_type[j], "/MSSA_Healthcare_Community_diff.xlsx", sep=''), col_names = FALSE)
}




for (j in c(1:length(selected_cell_type))){
  dir.create(file.path('MRSA_Differential_peaks', selected_cell_type[j]), showWarnings = FALSE)
  idxPass <- which(str_detect(proj$Cell_type_voting, selected_cell_type[j]))
  cellsPass <- proj$cellNames[idxPass]
  proj_temp<-proj[cellsPass, ]
  
  ## MRSA vs. MSSA
  marker_MRSA_MSSA <- getMarkerFeatures(ArchRProj = proj_temp, useMatrix = "PeakMatrix", groupBy = "Conditions",
                                        testMethod = "wilcoxon", bias = c("log10(nFrags)"), normBy = "ReadsInPeaks", maxCells = 15000,
                                        useGroups = "MRSA", bgdGroups = "MSSA")
  p_MRSA_MSSA <- plotMarkers(seMarker = marker_MRSA_MSSA, name = "MRSA", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
  plotPDF(p_MRSA_MSSA, name = paste(selected_cell_type[j], "MRSA_MSSA_Volcano_V2", sep='_'), width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
  write_xlsx(as.data.frame(marker_MRSA_MSSA@elementMetadata@listData), paste("MRSA_Differential_peaks/", selected_cell_type[j], "/MRSA_MSSA_peaks_V2.xlsx", sep=''))
  write_xlsx(marker_MRSA_MSSA@assays@data@listData, paste("MRSA_Differential_peaks/", selected_cell_type[j], "/MRSA_MSSA_diff_V2.xlsx", sep=''), col_names = FALSE)
  
  ## MRSA vs. control
  marker_MRSA_Control <- getMarkerFeatures(ArchRProj = proj_temp, useMatrix = "PeakMatrix", groupBy = "Conditions",
                                           testMethod = "wilcoxon", bias = c("log10(nFrags)"), normBy = "ReadsInPeaks", maxCells = 15000,
                                           useGroups = "MRSA", bgdGroups = "Control")
  p_MRSA_Control <- plotMarkers(seMarker = marker_MRSA_Control, name = "MRSA", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
  plotPDF(p_MRSA_Control, name = paste(selected_cell_type[j], "MRSA_Control_Volcano_V2", sep='_'), width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
  write_xlsx(as.data.frame(marker_MRSA_Control@elementMetadata@listData), paste("MRSA_Differential_peaks/", selected_cell_type[j], "/MRSA_Control_peaks_V2.xlsx", sep=''))
  write_xlsx(marker_MRSA_Control@assays@data@listData, paste("MRSA_Differential_peaks/", selected_cell_type[j], "/MRSA_Control_diff_V2.xlsx", sep=''), col_names = FALSE)
  
  ## MSSA vs. control
  marker_MSSA_Control <- getMarkerFeatures(ArchRProj = proj_temp, useMatrix = "PeakMatrix", groupBy = "Conditions",
                                           testMethod = "wilcoxon", bias = c("log10(nFrags)"), normBy = "ReadsInPeaks", maxCells = 15000,
                                           useGroups = "MSSA", bgdGroups = "Control")
  p_MSSA_Control <- plotMarkers(seMarker = marker_MSSA_Control, name = "MSSA", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
  plotPDF(p_MSSA_Control, name = paste(selected_cell_type[j], "MSSA_Control_Volcano_V2", sep='_'), width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
  write_xlsx(as.data.frame(marker_MSSA_Control@elementMetadata@listData), paste("MRSA_Differential_peaks/", selected_cell_type[j], "/MSSA_Control_peaks_V2.xlsx", sep=''))
  write_xlsx(marker_MSSA_Control@assays@data@listData, paste("MRSA_Differential_peaks/", selected_cell_type[j], "/MSSA_Control_diff_V2.xlsx", sep=''), col_names = FALSE)
}



proj <- addCoAccessibility(ArchRProj = proj,reducedDims = "IterativeLSI", maxDist =1e6)
Loop_plot <- plotBrowserTrack(ArchRProj = proj, groupBy = "Cell_type_voting", geneSymbol = markerGenes,
                              upstream = 50000, downstream = 50000, loops = getCoAccessibility(proj))

save.image(file='scATAC_celltype_peak_motif_coaccess.RData')

# refined peaks
proj_3 <- addIterativeLSI(ArchRProj = proj, useMatrix = "PeakMatrix", name = "IterativeLSI", iterations = 3,  force = TRUE,
                        clusterParams = list(resolution = c(4), sampleCells = 10000, n.start = 30),
                        varFeatures = 7000, dims = 1:30)
proj_3 <- addUMAP(ArchRProj = proj_3, reducedDims = "IterativeLSI", force = TRUE)
proj_3 <- addClusters(input = proj_3, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 4, knnAssign = 30, maxClusters = NULL, force = TRUE)
p1 <- plotEmbedding(ArchRProj = proj_3, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", force = TRUE)
p2 <- plotEmbedding(ArchRProj = proj_3, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
p3 <- plotEmbedding(ArchRProj = proj_3, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)

ggAlignPlots(p1, p2, p3,type = "h")

addArchRThreads(threads = 4)
proj_3 <- addGeneIntegrationMatrix(
  ArchRProj = proj_3, 
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

#cell types from the scRNA-seq data are most abundant in each of our scATAC-seq clusters
cM <- as.matrix(confusionMatrix(proj_3$Clusters, proj_3$predictedGroup))
pre_cluster <- rownames(cM)
max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]

Cell_type_voting = proj_3$Clusters
for (m in c(1:length(pre_cluster))){
  idxSample <- which(proj_3$Clusters == pre_cluster[m])
  Cell_type_voting[idxSample] <- max_celltype[m]
}

proj_3 <- addCellColData(ArchRProj = proj_3, data = Cell_type_voting, cells = proj_3$cellNames, name = "Cell_type_voting", force = TRUE)

cM <- confusionMatrix(paste0(proj_3$Cell_type_voting), paste0(proj_3$Sample))
cM <- t(t(cM) / Matrix::rowSums(t(cM)))
#cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), border_color = "black")
p6 <- plotEmbedding(proj_3, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE)
p6




proj<-loadArchRProject(path = "Staph_scATAC_integration", force = FALSE, showLogo = TRUE)

idxPass <- which(str_detect(proj$Cell_type_voting,'CD14 Mono'))
proj_CD14mono<-proj[idxPass, ]
proj_CD14mono <- addUMAP(ArchRProj = proj_CD14mono, reducedDims = "IterativeLSI", force = TRUE)

idxPass <- which(str_detect(proj$Cell_type_voting,'B memory'))
proj_B <- proj[idxPass, ]
proj_B <- addUMAP(ArchRProj = proj_B, reducedDims = "IterativeLSI", force = TRUE)

idxPass <- which(str_detect(proj$Cell_type_voting,'cDC2'))
proj_cDC2 <- proj[idxPass, ]
proj_cDC2 <- addUMAP(ArchRProj = proj_cDC2, reducedDims = "IterativeLSI", force = TRUE)

idxPass <- which(str_detect(proj$Cell_type_voting,'CD4 TCM'))
proj_T <- proj[idxPass, ]
proj_T <- addUMAP(ArchRProj = proj_T, reducedDims = "IterativeLSI", force = TRUE)


idxPass <- which(str_detect(proj$Cell_type_voting,'CD8 TEM'))
proj_CD8 <- proj[idxPass, ]
proj_CD8 <- addUMAP(ArchRProj = proj_CD8, reducedDims = "IterativeLSI", force = TRUE)

p6 <- plotEmbedding(proj_CD14mono, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", force = TRUE)
p6


markergene_UMAP_mono <- plotEmbedding(ArchRProj = proj_CD14mono, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAP", quantCut = c(0.05, 0.95), imputeWeights = NULL, pal = paletteContinuous(set = "whiteBlue"))
markergene_UMAP_mono$S100A12
markergene_UMAP_mono$S100A9
markergene_UMAP_mono$S100A8

p <- plotBrowserTrack(ArchRProj = proj_CD14mono, groupBy = "Conditions", geneSymbol = markerGenes, upstream = 10000,downstream = 500)
grid::grid.newpage()
grid::grid.draw(p$S100A12)
p_2 <- plotBrowserTrack(ArchRProj = proj_CD14mono, groupBy = "Conditions", geneSymbol = 'HLA-DRB1', upstream = 3000, downstream = 3000)
grid::grid.newpage()
grid::grid.draw(p_2$`HLA-DRB1`)


p_2 <- plotBrowserTrack(ArchRProj = proj_CD14mono, groupBy = "Conditions", geneSymbol = 'VIM', upstream = 1500, downstream = 1500)
grid::grid.newpage()
grid::grid.draw(p_2$`IER3`)


CFL1_UMAP_mono <- plotEmbedding(ArchRProj = proj_CD14mono, colorBy = "GeneScoreMatrix", name = 'CFL1', embedding = "UMAP", quantCut = c(0.65, 0.99), imputeWeights = NULL, pal = paletteContinuous(set = "whiteBlue"))
CFL1_UMAP_mono

HLADRB1_UMAP_mono <- plotEmbedding(ArchRProj = proj_CD14mono, colorBy = "GeneScoreMatrix", name = 'HLA-DRB1', embedding = "UMAP", quantCut = c(0.6, 0.95), imputeWeights = NULL, pal = paletteContinuous(set = "whiteBlue"))
HLADRB1_UMAP_mono



p_2 <- plotBrowserTrack(ArchRProj = proj_T, groupBy = "Conditions", geneSymbol = 'HSP90AB1', upstream = 10000, downstream = 10000)
grid::grid.newpage()
grid::grid.draw(p_2$HSP90AB1)

trajectory <- c("Control", "MSSA", "MRSA")
proj_CD14mono <- addTrajectory(ArchRProj = proj_CD14mono, name = "MyeloidU", groupBy = "Conditions", trajectory = trajectory, embedding = "UMAP", force = TRUE)
p_3 <- plotTrajectory(proj_CD14mono, trajectory = "MyeloidU", colorBy = "cellColData", name = "MyeloidU")
p_3








