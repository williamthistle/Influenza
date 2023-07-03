# Create Seurat object with important info for QC
load_archR_from_input_files <- function(inputFiles, analysis_dir) {
  # Create arrow files from raw input fragment files
  ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    minTSS = 4, #Don't set this too high because you can always increase later
    minFrags = 3000,
    maxFrags = 30000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
  )
  # Add doublet scores for each sample (will take a while)
  doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
    LSIMethod = 1
  )
  # Create and save ArchR project based on arrow files
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = paste0(analysis_dir, "ArchR/"),
    copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.
  )
  # Filter out doublets
  proj <- filterDoublets(ArchRProj = proj)
  # Save our project
  saveArchRProject(ArchRProj = proj, load = FALSE)
  return(proj)
}

# Add metadata for each sample to ArchR object
add_sample_metadata_atac <- function(proj, high_viral_load_samples, low_viral_load_samples,
                                d28_samples, d_minus_1_samples, male_samples, female_samples) {
  # Add viral load metadata to ArchR object
  viral_load_vector <- proj$Sample
  idxSample <- which(proj$Sample %in% high_viral_load_samples)
  viral_load_vector[idxSample]<-'HVL'
  idxSample <- which(proj$Sample %in% low_viral_load_samples)
  viral_load_vector[idxSample]<-'LVL'
  proj <- addCellColData(ArchRProj = proj, data = viral_load_vector, cells = proj$cellNames, name = "viral_load", force = TRUE)
  # Add day metadata to ArchR object
  day_vector <- proj$Sample
  idxSample <- which(proj$Sample %in% d28_samples)
  day_vector[idxSample]<-'D28'
  idxSample <- which(proj$Sample %in% d_minus_1_samples)
  day_vector[idxSample]<-'D_MINUS_1'
  proj <- addCellColData(ArchRProj = proj, data = day_vector, cells = proj$cellNames, name = "day", force = TRUE)
  # Add sex metadata to ArchR object
  sex_vector <- proj$Sample
  idxSample <- which(proj$Sample %in% male_samples)
  sex_vector[idxSample]<-'MALE'
  idxSample <- which(proj$Sample %in% female_samples)
  sex_vector[idxSample]<-'FEMALE'
  proj <- addCellColData(ArchRProj = proj, data = sex_vector, cells = proj$cellNames, name = "sex", force = TRUE)
  return(proj)
}

# Create object with sample names and associated metadata value
# This is probably not an important function for anyone else
parse_metadata_for_samples <- function(proj, group, high_viral_load_samples, low_viral_load_samples,
                                       d28_samples, d_minus_1_samples, male_samples, female_samples) {
  current_metadata <- c()
  metadata_names <- unique(proj$Sample)
  for(name in metadata_names) {
    if(group == "viral_load") {
      if(name %in% high_viral_load_samples) {
        current_metadata <- c(current_metadata, "HVL")
      } else {
        current_metadata <- c(current_metadata, "LVL")
      }
    } else if(group == "day") {
      if(name %in% d28_samples) {
        current_metadata <- c(current_metadata, "D28")
      } else {
        current_metadata <- c(current_metadata, "D_MINUS_1")
      }
    } else if(group == "sex") {
      if(name %in% male_samples) {
        current_metadata <- c(current_metadata, "MALE")
      } else {
        current_metadata <- c(current_metadata, "FEMALE")
      }
    }
  }
  names(current_metadata) <- metadata_names
  return(current_metadata)
}

# Plot QC metrics for ATAC
plot_qc_atac <- function(proj, date) {
  # Run dimensionality reduction and plot UMAP of data before any filtering is done
  addArchRThreads(threads = 8)
  proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                          clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                          varFeatures = 25000, dimsToUse = 2:30)
  proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
  plotPDF(p1, name = paste0("Dataset_Prefiltering", date), ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)
  # Plot out TSS Enrichment / Doublet Enrichment / Nucleosome Ratio for each sample to help us decide filtering thresholds
  p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
  p2 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "DoubletEnrichment", plotAs = "ridges")
  p3 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges")
  plotPDF(p1,p2,p3, name = paste0("Integrated_Scores_Prefiltering_", date), ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)
}