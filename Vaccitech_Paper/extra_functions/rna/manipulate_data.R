# Initial combination of cell types
combine_cell_types_initial <- function(sc_obj, resolution = 1.5) {
  sc_obj$old.predicted.id <- sc_obj$predicted.id
  Cell_type_combined <- sc_obj$predicted.id
  idx <- grep("CD4 T", Cell_type_combined)
  Cell_type_combined[idx] <- "CD4 Memory"
  idx <- grep("CD8 T", Cell_type_combined)
  Cell_type_combined[idx] <- "CD8 Memory"
  idx <- grep("cDC", Cell_type_combined)
  Cell_type_combined[idx] <- "cDC"
  idx <- grep("Proliferating", Cell_type_combined)
  Cell_type_combined[idx] <- "Proliferating"
  sc_obj$predicted.id <- Cell_type_combined
  sc_obj <- MajorityVote(sc_obj, resolution)
  return(sc_obj)
}

# Combine more cell types (for MAGICAL)
combine_cell_types_magical <- function(sc_obj) {
  Cell_type_combined <- sc_obj$predicted_celltype_majority_vote
  levels(Cell_type_combined) <- c(levels(Cell_type_combined), "T Naive", "B", "NK_MAGICAL")
  idx <- grep("CD4 Naive", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("CD8 Naive", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("Treg", Cell_type_combined)
  Cell_type_combined[idx] <- "T Naive"
  idx <- grep("NK", Cell_type_combined)
  Cell_type_combined[idx] <- "NK_MAGICAL"
  idx <- grep("B", Cell_type_combined)
  Cell_type_combined[idx] <- "B"
  sc_obj$magical_cell_types <- Cell_type_combined
  return(sc_obj)
}

# Override cluster label from majority voting
override_cluster_label <- function(sc_obj, cluster_nums, new_cluster_label) {
  for(cluster_id in cluster_nums) {
    print(cluster_id)
    if(!(new_cluster_label %in% levels(sc_obj$predicted_celltype_majority_vote))) {
      levels(sc_obj$predicted_celltype_majority_vote) <- c(levels(sc_obj$predicted_celltype_majority_vote), new_cluster_label)
    }
    majority_vote <- sc_obj$predicted_celltype_majority_vote
    idxPass <- which(Idents(sc_obj) %in% c(cluster_id))
    majority_vote[idxPass] <- new_cluster_label
    sc_obj$predicted_celltype_majority_vote <- majority_vote
  }
  return(sc_obj)
}

# Remove specific samples from Seurat object
remove_specific_samples_from_sc_obj <- function(sc_obj, samples) {
  idxPass <- which(sc_obj$sample %in% samples)
  cellsPass <- names(sc_obj$orig.ident[-idxPass])
  sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  return(sc_obj)
}

remove_cells_based_on_umap <- function(sc_obj, first_x, second_x, first_y, second_y) {
  orig.umap.coords <- as.data.frame(sc_obj[["umap"]]@cell.embeddings)
  orig.umap.coords$cells <- rownames(orig.umap.coords)
  deleted.cells.umap.coords <- orig.umap.coords[orig.umap.coords$"UMAP_1" > first_x & orig.umap.coords$"UMAP_1" < second_x,]
  deleted.cells.umap.coords <- deleted.cells.umap.coords[deleted.cells.umap.coords$"UMAP_2" > first_y & deleted.cells.umap.coords$"UMAP_2" < second_y,]
  final.umap.coords <- orig.umap.coords[!(orig.umap.coords$cells %in% deleted.cells.umap.coords$cells),]
  cellsPass <- rownames(final.umap.coords)
  sc_obj <- subset(x = sc_obj, subset = cell_name %in% cellsPass)
  return(sc_obj)
}




#' Infer batches using LISI metric
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @return A Seurat object which contains labeled batches
#' @examples
#' \dontrun{sc_obj <- InferBatches(sc_obj)}
#' @export
InferBatches_alt <- function(sc_obj, log_flag = FALSE) {
  # Find clusters in data (prior to batch correction)
  if ('lsi' %in% Seurat::Reductions(sc_obj)) {
    sc_obj <- Seurat::FindNeighbors(object = sc_obj, reduction = "lsi", dims = 1:30)
    sc_obj <- Seurat::FindClusters(object = sc_obj, resolution = 0.2, algorithm = 2)
  } else {
    sc_obj <- Seurat::FindNeighbors(object = sc_obj, reduction = "pca", dims = 1:30)
    sc_obj <- Seurat::FindClusters(object = sc_obj, resolution = 0.1, algorithm = 2)
  }
  # Use LISI metric to guess batch labels
  X <- sc_obj@reductions$umap@cell.embeddings
  meta_data <- data.frame(sc_obj$sample)
  colnames(meta_data) <- "batch"
  meta_data$cluster <- sc_obj$seurat_clusters
  lisi.res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(lisi.res) <- c("batch", "score", "cluster", "freq")
  clusters.interest <- names(table(sc_obj$seurat_clusters))[prop.table(table(sc_obj$seurat_clusters)) > 0.01]
  for (cluster in clusters.interest) { #levels(sc_obj$seurat_clusters)) {
    cells <- names(sc_obj$seurat_clusters[sc_obj$seurat_clusters == cluster])
    X.sub <- X[which(rownames(X) %in% cells),]
    meta_data.sub <- meta_data[which(rownames(meta_data) %in% cells),]
    res <- lisi::compute_lisi(X.sub, meta_data.sub, label_colnames = "batch")
    rownames(res) <- cells
    colnames(res) <- "score"
    res$batch <- meta_data.sub$batch
    agg.res <- stats::aggregate(.~batch,data=res,mean)
    agg.res$cluster <- cluster
    agg.res$freq <- data.frame(table(res$batch))$Freq[which(data.frame(table(res$batch))$Var1 %in% agg.res$batch)]
    lisi.res <- rbind(lisi.res, agg.res)
  }
  
  p.values <- list()
  used.sample.dump <- c()
  batch.assign <- list()
  for ( i in clusters.interest) {
    lisi.res.sub <- lisi.res[lisi.res$cluster == i,]
    if (max(lisi.res.sub$score) <= 1.1) {
      samples.of.batch <- lisi.res.sub$batch[1]
      if (!(samples.of.batch %in% used.sample.dump)) {
        batch.assign <- lappend(batch.assign, samples.of.batch)
      }
      used.sample.dump <- union(used.sample.dump, samples.of.batch)
    } else {
      lisi.res.sub$scaled.score <- scale_zero_one(lisi.res.sub$score * (lisi.res.sub$freq / sum(lisi.res.sub$freq)))
      lisi.res.sub <- lisi.res.sub[order(lisi.res.sub$scaled.score, decreasing = TRUE),]
      if (dim(lisi.res.sub)[1] > 30) {
        lisi.res.sub <- lisi.res.sub[1:30,]
      }
      lisi.res.sub$diff.scaled.score <- abs(c(diff(lisi.res.sub$scaled.score), 0))
      
      if (dim(lisi.res.sub)[1] >= 3) {
        p.values[[i]] <- outliers::dixon.test(lisi.res.sub$diff.scaled.score)$p.value[[1]]
      } else {
        p.values[[i]] <- 1
      }
      
      if (p.values[[i]] < 0.05) {
        max.index <- which.max(lisi.res.sub$diff.scaled.score)
        samples.of.batch <- lisi.res.sub$batch[1:max.index]
        
        if (any(samples.of.batch %in% used.sample.dump)) {
          if (!all(samples.of.batch %in% used.sample.dump)) {
            used.index <- which(samples.of.batch %in% used.sample.dump)
            samples.of.batch <- samples.of.batch[-used.index]
            if (length(samples.of.batch) > 0) {
              batch.assign <- lappend(batch.assign, samples.of.batch)
            }
          } else if (!list(samples.of.batch) %in% batch.assign) {
            if (length(samples.of.batch) == 1) {
              batch.assign <- lappend(batch.assign, samples.of.batch)
            } else {
              used.index <- which(samples.of.batch %in% unlist(batch.assign))
              samples.of.batch <- samples.of.batch[-used.index]
              if (length(samples.of.batch) > 0) {
                batch.assign <- lappend(batch.assign, samples.of.batch)
              }
            }
          }
        } else {
          batch.assign <- lappend(batch.assign, samples.of.batch)
        }
        used.sample.dump <- union(used.sample.dump, samples.of.batch)
      }
    }
  }
  
  batch <- as.factor(sc_obj$sample)
  if (length(batch.assign) > 0) {
    levels.batch <- levels(batch)
    for (i in 1:length(batch.assign)) {
      levels.batch[which(levels(batch) %in% batch.assign[[i]])] <- i
    }
    levels.batch[!levels.batch %in% c(1:length(batch.assign))] <- length(batch.assign)+1
    levels(batch) <- levels.batch
    sc_obj$batch <- as.character(batch)
  } else {
    sc_obj$batch <- "No Batch"
  }
  gc()
  return(sc_obj)
}

MajorityVote_RNA_alt <- function(sc_obj, current_resolution = 1, log_flag = FALSE) {
  if(Seurat::DefaultAssay(sc_obj) == "integrated") {
    associated_res_attribute <- paste0("integrated_snn_res.", current_resolution)
  } else {
    associated_res_attribute <- paste0("SCT_snn_res.", current_resolution)
  }
  sc_obj <- Seurat::FindNeighbors(sc_obj, reduction = "pca", dims = 1:30)
  sc_obj <- Seurat::FindClusters(sc_obj, resolution = current_resolution)
  sc_obj$predicted.id <- as.character(sc_obj$predicted.id)
  
  integrated_snn_res_df <- sc_obj[[associated_res_attribute]]
  integrated_snn_res_cell_names <- rownames(integrated_snn_res_df)
  integrated_snn_res_values <- integrated_snn_res_df[,1]
  
  cluster.dump <- as.numeric(levels(integrated_snn_res_values))
  sc_obj$predicted_celltype_majority_vote <- sc_obj$seurat_clusters
  levels(sc_obj$predicted_celltype_majority_vote) <- as.character(levels(sc_obj$predicted_celltype_majority_vote))
  for (i in unique(sc_obj$predicted.id)) {
    print(i)
    cells <- names(sc_obj$predicted.id[sc_obj$predicted.id == i])
    freq.table <- as.data.frame(table(integrated_snn_res_df[cells,]))
    freq.table <- freq.table[order(freq.table$Freq, decreasing = TRUE),]
    freq.table$diff <- abs(c(diff(freq.table$Freq), 0))
    if(nrow(freq.table) > 30) {
      freq.table <- freq.table[1:30,]
    }
    p.values <- outliers::dixon.test(freq.table$diff)$p.value[[1]]
    max.index <- which.max(freq.table$diff)
    clusters <- as.numeric(as.character(freq.table$Var1[1:max.index]))
    levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(clusters)] <- i
    cluster.dump <- cluster.dump[!cluster.dump %in% clusters]
  }
  
  if (length(cluster.dump) > 0) {
    for (i in cluster.dump) {
      cells <- rownames(subset(integrated_snn_res_df, integrated_snn_res_df[,1] == i,))
      freq.table <- as.data.frame(table(sc_obj$predicted.id[cells]))
      levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(i)] <- as.vector(freq.table$Var1)[which.max(freq.table$Freq)]
    }
  }
  return(sc_obj)
}

#' Scale a vector to range(0,1)
#'
#' @param x Numeric vector
#' @return A scaled vector ranging from 0 to 1
scale_zero_one <- function(x) {(x - min(x))/(max(x) - min(x))}

#' Append a list to a list-of-lists
#'
#' @param lst List
#' @param ... Additional lists
#' @return A list of lists
lappend <- function (lst, ...){ c(lst, list(...))}

#' Integrate batches (ATAC)
#'
#' @param proj ArchR project containing cells for all samples
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return An ArchR object which contains integrated data
IntegrateByBatch_ATAC_alt <- function(proj, log_flag = FALSE) {
  tile_sce <- ArchR::getMatrixFromProject(proj, useMatrix='TileMatrix', binarize = TRUE)
  tile_reduc <- ArchR::getReducedDims(ArchRProj = proj, reducedDims = "IterativeLSI", returnMatrix = TRUE)
  tile_reduc <- tile_reduc[match(colnames(tile_sce), rownames(tile_reduc)),]
  for (i in colnames(SummarizedExperiment::colData(tile_sce))) {
    SummarizedExperiment::colData(tile_sce)[[i]] <- as.vector(SummarizedExperiment::colData(tile_sce)[[i]])
  }
  rownames(tile_sce) <- paste0(as.character(SummarizedExperiment::rowData(tile_sce)$seqnames),
                               '-',
                               as.character(SummarizedExperiment::rowData(tile_sce)$start))
  tile_seurat <- Seurat::CreateSeuratObject(SummarizedExperiment::assays(tile_sce)$TileMatrix[, rownames(tile_reduc)],
                                            project = "peaks",
                                            assay = "tileMatrix")
  tile_seurat <- Seurat::AddMetaData(tile_seurat, data.frame(t(SummarizedExperiment::colData(tile_sce))))
  cell.embeddings <- tile_reduc
  feature.loadings <- matrix()
  assay <- "tileMatrix"
  sdev <- 0
  reduction.key <- "LSI_"
  reduction.data <- Seurat::CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key,
    misc = list()
  )
  tile_seurat@reductions$lsi <- reduction.data
  tile_umap <- ArchR::getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
  cell.embeddings <- as.matrix(tile_umap)
  feature.loadings <- matrix()
  assay <- "tileMatrix"
  sdev <- 0
  reduction.key <- "UMAP_"
  reduction.data <- Seurat::CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key,
    misc = list()
  )
  tile_seurat@reductions$umap <- reduction.data
  tile_seurat$sample <- proj$Sample
  # Step 4: Inferring batches
  tile_seurat <- InferBatches_alt(tile_seurat, log_flag)
  proj$Batch <- tile_seurat$batch
  # If we only have one batch, we don't need to integrate by batch, so we exit the function
  if(length(unique(proj$Batch)) == 1) {
  } else {
    proj <- ArchR::addHarmony(ArchRProj = proj, reducedDims = "IterativeLSI",
                              dimsToUse = 2:30, scaleDims = TRUE,
                              corCutOff = 0.75, groupBy = "Batch", force = TRUE)
    proj <- ArchR::addUMAP(ArchRProj = proj, reducedDims = "Harmony", force = TRUE)
    proj <- ArchR::addClusters(input = proj, reducedDims = "Harmony", method = "Seurat",
                               name = "Harmony_clusters", resolution = 2, knnAssign = 30,
                               maxClusters = NULL, force = TRUE)
  }
  gc()
  return(proj)
}

InitialProcessing_ATAC_alt <- function(proj, log_flag = FALSE) {
  proj <- ArchR::addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",
                                 iterations = 2,
                                 force = TRUE,
                                 clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                                 varFeatures = 25000, dimsToUse = 2:30,
                                 saveIterations = TRUE)
  proj <- ArchR::addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
  proj <- ArchR::addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat",
                             name = "Clusters", resolution = 5, knnAssign = 30,
                             maxClusters = NULL, force = TRUE)
  p1 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p2 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  p3 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
  ArchR::plotPDF(p1,p2,p3, name = "UMAPs_After_Initial_Processing_plots", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
  gc()
  return(proj)
}
