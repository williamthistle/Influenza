library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(clustree)
library(dplyr)
library(ggplot2)
library(hexbin)
library(mclust)
library(openxlsx)
library(outliers)
library(pheatmap)
library(scDblFinder)
library(SeuratDisk)
library(BiocParallel)
library(stringr)
library(writexl)
library(SoupX)
library(SingleR)

################## SETUP ##################
date <- Sys.Date()
home_dir <- "~/"
# Load SPEEDI (for RNA-seq analyses)
SPEEDI_dir <- paste0(home_dir, "SPEEDI")
source(paste0(SPEEDI_dir, "/prototype_API.R"))

# Load information about samples
sample_metadata <- read.table(paste0(SPEEDI_dir, "/sample_metadata.tsv"), sep = "\t", header = TRUE)

# Declare data and analysis type
data_type <- "multiome" # Can be multiome or single_cell
analysis_type <- "ATAC-seq" # Can be RNA-seq or ATAC-seq

# Directory where all data are for data_type specified above
data_path <- paste0(home_dir, data_type, "/data/")
# List of samples that can potentially be processed for the data type
all_sample_id_list <- list.dirs(data_path, recursive = FALSE)
all_sample_id_list <- strsplit(all_sample_id_list, "/")
all_sample_id_list <- unlist(lapply(all_sample_id_list, tail, n = 1L))

# data_token is used to choose subset of data that we want to analyze (pre-defined in flu_data_tokens.tsv)
data_token <- "all_multiome"
# Create directory for this particular data in the analysis directory if it doesn't exist
base_analysis_dir <- paste0(home_dir, data_type, "/analysis/", data_token, "/")
if (!dir.exists(base_analysis_dir)) {dir.create(base_analysis_dir)}
# Grab samples that we want to analyze
data_tokens <- read.table(paste0(home_dir, "flu_data_tokens.tsv"), header = TRUE)
samples <- data_tokens[data_tokens$token == data_token,]$samples
samples <-  unlist(strsplit(samples, ","))
sample_id_list <- all_sample_id_list[all_sample_id_list %in% samples]
if(length(sample_id_list) != length(samples)) {
  print(paste0("Requested samples: ", samples))
  print(paste0("Available samples (from your list): ", sample_id_list))
  stop("You have at least one invalid sample in your list")
}
sample_count <- length(sample_id_list)

# analysis_token is used to define a specific analysis
analysis_token <- "default"
analysis_dir <- paste0(home_dir, data_type, "/analysis/", data_token, "/", analysis_token, "/")
if (!dir.exists(analysis_dir)) {dir.create(analysis_dir)}
# Directory for analysis plots
plot_dir <- paste0(analysis_dir, "plots/")
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
# Directory for markers
marker_dir <- paste0(analysis_dir, "markers/")
if (!dir.exists(marker_dir)) {dir.create(marker_dir)}

# Identify high / low viral load samples, D28 / D-1 samples, and M / F samples
# Subset metadata to only contain our sample_ids
# We will add these metadata to our Seurat object later
sample_metadata <- sample_metadata[sample_metadata$aliquot %in% sample_id_list,]
high_viral_load_samples <- sort(sample_metadata[sample_metadata$viral_load == "high",]$aliquot)
low_viral_load_samples <- sort(sample_metadata[sample_metadata$viral_load == "low",]$aliquot)
all_viral_load_samples <- c(high_viral_load_samples, low_viral_load_samples)
d28_samples <- sort(sample_metadata[sample_metadata$time_point == "2_D28",]$aliquot)
d_minus_1_samples <- sort(sample_metadata[sample_metadata$time_point == "2_D_minus_1",]$aliquot)
all_day_samples <- c(d28_samples, d_minus_1_samples)
male_samples <- sort(sample_metadata[sample_metadata$sex == "M",]$aliquot)
female_samples <- sort(sample_metadata[sample_metadata$sex == "F",]$aliquot)
all_sex_samples <- c(male_samples, female_samples)

# Organize sample list by high viral load followed by low viral load
sample_id_list <- sample_id_list[order(match(sample_id_list, all_viral_load_samples))]
# Sort other types of samples according to viral load order
#all_day_samples <- all_day_samples[order(match(all_day_samples, all_viral_load_samples))]
#all_sex_samples <- all_sex_samples[order(match(all_sex_samples, all_viral_load_samples))]

# Parameters for processing both RNA-seq and ATAC data
# save_progress: If you want to save your progress
save_progress <- FALSE
# Parameters for processing RNA-seq
# record_doublets: If you want to run scDblFinder and record which cells are doublets
# run_qc: If you want to plot QC metrics (and not run the rest of the pipeline)
# run_soup: If you want to remove ambient RNA using SoupX
record_doublets <- FALSE
run_qc <- FALSE
run_soup <- FALSE
# Parameters for processing ATAC-seq data
# use_rna_labels: Uses labels from the RNA data on the ATAC-data (we assume these labels are more accurate than what the GeneIntegrationMatrix finds)
use_rna_labels <- TRUE
# subset_to_rna: Subset to only those cells we kept from RNA
subset_to_rna <- TRUE
# TODO
################## ANALYSIS ##################
if(analysis_type == "RNA_seq") {
  # Step 1 - grab matrices
  if(run_soup) {
    # TODO: Make this parallel
    all_sc_exp_matrices <- process_matrices_through_soup(data_path, sample_id_list)
  } else {
    all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
  }
  # Run QC (optional)
  if(run_qc) {
    sc_obj <- generate_qc_plots(all_sc_exp_matrices, plot_dir, date, high_viral_load_samples, low_viral_load_samples,
                    d28_samples, d_minus_1_samples, male_samples, female_samples)
    # If you run QC, you probably don't want to proceed with the rest of the pipeline, so we stop here
  } else {
    # Step 2 - filter raw data (either with adaptive thresholds or strict thresholds)
    sc_obj <- FilterRawData(all_sc_exp_matrices, human = TRUE, record_doublets = FALSE, adaptive_QC_thresholds = FALSE)
    rm(all_sc_exp_matrices)
    sc_obj <- add_sample_metadata(sc_obj, high_viral_load_samples, low_viral_load_samples,
                        d28_samples, d_minus_1_samples, male_samples, female_samples)
    # Step 3 - normalize data
    sc_obj <- InitialProcessing(sc_obj, human = TRUE)
    if(save_progress) {
      save(sc_obj, file = paste0(analysis_dir, "3_sc_obj.rds"))
    }
    # Step 4 - find batches in data
    sc_obj <- InferBatches(sc_obj)
    # Step 5 - integrate data by batch
    sc_obj <- IntegrateByBatch(sc_obj)
    if(save_progress) {
      save(sc_obj, file = paste0(analysis_dir, "5_sc_obj.rds"))
    }
    # Step 6 - process integrated assay and potentially prepare to find markers via PrepSCTFindMarkers
    sc_obj <- VisualizeIntegration(sc_obj, prep_sct_find_markers = FALSE)
    if(save_progress) {
      save(sc_obj, file = paste0(analysis_dir, "6_sc_obj.rds"))
    }
    # Load reference and remove cell types we don't like
    reference <- LoadReference("PBMC", human = TRUE)
    idx <- which(reference$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
    reference <- reference[,-idx]
    # Map cell types from reference to our data
    sc_obj <- MapCellTypes(sc_obj, reference, data_type = "snRNA")
    rm(reference)
    # Combine cell types and re-do majority vote
    sc_obj <- combine_cell_types_initial(sc_obj, resolution = 3)
    # Print UMAP by cell type (majority vote) and by cluster number - it will currently be messy
    print_UMAP(sc_obj, sample_count, "predicted_celltype_majority_vote", plot_dir, paste0("pre.clusters_by_cell_type_majority_vote_", date, ".png"))
    print_UMAP(sc_obj, sample_count, "seurat_clusters", plot_dir, paste0("pre.clusters_by_cluster_num_", date, ".png"))
    print_UMAP(sc_obj, sample_count, "predicted.id", plot_dir, paste0("pre.clusters_by_cell_type_", date, ".png"))
    print_UMAP(sc_obj, sample_count, "viral_load", plot_dir, paste0("pre.clusters_by_viral_load_", date, ".png"))
    # We always want to save our sc_obj after processing data through SPEEDI
    save(sc_obj, file = paste0(analysis_dir, "7_sc_obj_sct_markers.rds"))
    # Load sc_obj
    load(file = paste0(analysis_dir, "7_sc_obj_sct_markers.rds"))
    load(paste0(analysis_dir, "singler_labels.rds"))
    # We can use clustree to help us figure out the best resolution
    print_clustree_plot(sc_obj, plot_dir, date)
    # Re-run majority vote with best resolution
    best_res <- 3
    sc_obj <- MajorityVote(sc_obj, best_res)
    # To decide which clusters we need to remove, we will capture information about clusters
    # We will also run DE for each cluster to find cell type markers
    raw_cluster_info <- capture_cluster_info(sc_obj)
    run_differential_expression_cluster(sc_obj, marker_dir)
    # Remove messy clusters
    messy_clusters <- c(0,18,20,22,25,28,33,35,36,38,39,40,41,42,44,45,46,49,50,52,53,59)
    idxPass <- which(Idents(sc_obj) %in% messy_clusters)
    cellsPass <- names(sc_obj$orig.ident[-idxPass])
    sc_obj.minus.messy.clusters <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
    # Override labels manually where necessary
    sc_obj.minus.messy.clusters <- override_cluster_label(sc_obj.minus.messy.clusters, c(15), "CD16 Mono")
    sc_obj.minus.messy.clusters <- override_cluster_label(sc_obj.minus.messy.clusters, c(34), "pDC")
    # Print info about sample representation and breakdown of categories per cell type
    print(table(sc_obj.minus.messy.clusters$sample))
    print_celltype_counts(sc_obj.minus.messy.clusters)
    # Print plots
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted_celltype_majority_vote", plot_dir, paste0("post.clusters_by_cell_type_majority_vote_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "seurat_clusters", plot_dir, paste0("post.clusters_by_cluster_num_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "predicted.id", plot_dir, paste0("post.clusters_by_cell_type_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "viral_load", plot_dir, paste0("post.clusters_by_viral_load_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "sample", plot_dir, paste0("post.clusters_by_sample_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "day", plot_dir, paste0("post.clusters_by_day_", date, ".png"))
    print_UMAP(sc_obj.minus.messy.clusters, sample_count, "sex", plot_dir, paste0("post.clusters_by_sex_", date, ".png"))
    # write cells and associated cell type majority predictions to file (for ATAC-seq labeling)
    cells_for_ATAC <- data.frame("cells" = sc_obj.minus.messy.clusters$cell_name, voted_type = sc_obj.minus.messy.clusters$predicted_celltype_majority_vote)
    write.csv(cells_for_ATAC, file = paste0(analysis_dir, "rna_seq_labeled_cells_", date, ".csv"), quote = FALSE, row.names = FALSE)
    # Combine cell types for MAGICAL and other analyses that require snATAC-seq (granularity isn't as good for ATAC-seq)
    sc_obj.minus.messy.clusters <- combine_cell_types_magical(sc_obj.minus.messy.clusters, best_res)
    #sc_obj.minus.messy.clusters$magical_cell_types <- sc_obj.minus.messy.clusters$predicted_celltype_majority_vote # This is if MAGICAL cell types == RNA-seq cell types
    # Run differential expression for each cell type within each group of interest
    run_differential_expression_group(sc_obj.minus.messy.clusters, analysis_dir, "viral_load")
    run_differential_expression_group(sc_obj.minus.messy.clusters, analysis_dir, "day")
    run_differential_expression_group(sc_obj.minus.messy.clusters, analysis_dir, "sex")
  }
} else if(analysis_type == "ATAC_seq") {
  # Label input files
  inputFiles <- paste0(data_path, sample_id_list, "/outs/atac_fragments.tsv.gz")
  names(inputFiles) <- all_viral_load_samples
  viral_load_metadata <- c(rep("HVL", length(high_viral_load_samples)), rep("LVL", length(low_viral_load_samples)))
  names(viral_load_metadata) <- all_viral_load_samples
  # Add relevant genome for ArchR
  addArchRGenome("hg38")
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
  saveArchRProject(ArchRProj = proj, load = FALSE)
  # Load ArchR project 
  proj <- loadArchRProject(path = paste0(analysis_dir, "/ArchR/"))
  # Filter out doublets (fake cells)
  proj <- filterDoublets(ArchRProj = proj)
  # Add viral load metadata to ArchR object
  viral_load_vector <- proj$Sample
  idxSample <- which(proj$Sample %in% high_viral_load_samples)
  viral_load_vector[idxSample]<-'HVL'
  idxSample <- which(proj$Sample %in% low_viral_load_samples)
  viral_load_vector[idxSample]<-'LVL'
  proj <- addCellColData(ArchRProj = proj, data = viral_load_vector, cells = proj$cellNames,name = "viral_load", force = TRUE)
  # Add day metadata to ArchR object
  day_vector <- proj$Sample
  idxSample <- which(proj$Sample %in% d28_samples)
  day_vector[idxSample]<-'D28'
  idxSample <- which(proj$Sample %in% d_minus_1_samples)
  day_vector[idxSample]<-'D_MINUS_1'
  proj <- addCellColData(ArchRProj = proj, data = day_vector, cells = proj$cellNames,name = "day", force = TRUE)
  # Add sex metadata to ArchR object
  sex_vector <- proj$Sample
  idxSample <- which(proj$Sample %in% male_samples)
  sex_vector[idxSample]<-'MALE'
  idxSample <- which(proj$Sample %in% female_samples)
  sex_vector[idxSample]<-'FEMALE'
  proj <- addCellColData(ArchRProj = proj, data = sex_vector, cells = proj$cellNames,name = "sex", force = TRUE)
  # Plot what dataset looks like before any processing
  addArchRThreads(threads = 8)
  proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                          clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                          varFeatures = 15000, dims = 2:30)
  proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
  plotPDF(p1, name = "Dataset_Prefiltering.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)
  # Plot out TSS Enrichment / Doublet Enrichment / Nucleosome Ratio for each sample
  p1 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
  p2 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "DoubletEnrichment", plotAs = "ridges")
  p3 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges")
  plotPDF(p1,p2,p3, name = "Integrated_Scores_Prefiltering.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 5)
  # Filter out cells that don't meet TSS enrichment / doublet enrichment / nucleosome ratio criteria
  idxPass <- which(proj$TSSEnrichment >= 8 & proj$NucleosomeRatio < 2 & proj$DoubletEnrichment < 3) 
  cellsPass <- proj$cellNames[idxPass]
  proj<-proj[cellsPass, ]
  # List number of cells remaining for each condition and each sample
  table(proj$viral_load)
  table(proj$day)
  table(proj$sex)
  table(proj$Sample)
  # Perform dimensionality reduction on cells (addIterativeLSI), create UMAP embedding (addUMAP), 
  # and add cluster information (addClusters)
  addArchRThreads(threads = 8)
  proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2,  force = TRUE,
                          clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                          varFeatures = 15000, dims = 2:30)
  proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
  proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat", name = "Clusters", resolution = 4, knnAssign = 30, maxClusters = NULL, force = TRUE)
  # UMAP plots colored by condition, sample, cluster ID, and TSS enrichment
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Conditions", embedding = "UMAP", force = TRUE)
  p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE)
  p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE)
  p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE)
  plotPDF(p1,p2,p3,p4, name = "Integrated_Clustering_snRNA_doublets_removed.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
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
  addArchRThreads(threads = 8)
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
  saveArchRProject(ArchRProj = proj, load = FALSE)
  pal <- paletteDiscrete(values = proj$predictedGroup)
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE)
  plotPDF(p1, name = "Integrated_annotated_gene_integration_matrix.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
  # Add labels from RNA-seq data
  if(use_rna_labels) {
    curated_snRNA_seq_cells <- read.csv(paste0(analysis_dir, "rna_seq_labeled_cells_", date, ".csv"), comment.char = "")
    if(subset_to_rna) {
      idxPass <- which(proj$cellNames %in% curated_snRNA_seq_cells$cells) 
      cellsPass <- proj$cellNames[idxPass]
      proj <- proj[cellsPass, ]
      curated_snRNA_seq_cells <- curated_snRNA_seq_cells[curated_snRNA_seq_cells$cells %in% proj$cellNames,]
      curated_snRNA_seq_cells <- curated_snRNA_seq_cells[order(match(curated_snRNA_seq_cells$cells,proj$cellNames)),]
      snRNA_seq_cell_votes <- curated_snRNA_seq_cells$voted_type
      proj <- addCellColData(ArchRProj = proj, data = snRNA_seq_cell_votes, cells = proj$cellNames, name = "predictedGroup", force = TRUE)
    } else {
      predicted_snATAC_cells <- data.frame(cell_name = proj$cellNames, voted_type = proj$predictedGroup)
      for(current_row in 1:nrow(predicted_snATAC_cells)) {
        current_snATAC_cell <- predicted_snATAC_cells[current_row,]$cell_name
        if(current_snATAC_cell %in% curated_snRNA_seq_cells$cells) {
          current_voted_type <- curated_snRNA_seq_cells[curated_snRNA_seq_cells$cells == current_snATAC_cell,]$voted_type
          predicted_snATAC_cells[current_row,]$voted_type <- current_voted_type
        }
      }
      proj <- addCellColData(ArchRProj = proj, data = predicted_snATAC_cells$voted_type, cells = proj$cellNames, name = "predictedGroup", force = TRUE)
    }
  }
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
  pal <- paletteDiscrete(values = proj$predictedGroup)
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE)
  plotPDF(p1, name = "Integrated_annotated_combined_gene_integration_matrix.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
  # Combine other cell types
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "CD4 Naive", "T Naive")
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "CD8 Naive", "T Naive")
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "NK_CD56bright", "NK")
  proj$predictedGroup <- replace(proj$predictedGroup, proj$predictedGroup == "Treg", "T Naive")
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
  plotPDF(p1,p2,p3,p4, name = paste0("Integrated_Clustering_Gene_Integration_Voting_1_subset", date), ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
} else {
  stop("Invalid analysis type")
}

