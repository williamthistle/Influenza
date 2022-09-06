library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)
library(SeuratDisk)
library(ArchR) # Just used for confusion matrix

chunk <- function(x, n) (mapply(function(a, b) (x[a:b]), seq.int(from=1, to=length(x), by=n), pmin(seq.int(from=1, to=length(x), by=n)+(n-1), length(x)), SIMPLIFY=FALSE))

# Set project dir - used to organize different projects
project_dir <- "~/PLACEBO_FLU_1/"
if (!dir.exists(project_dir)) {dir.create(project_dir)}
# We can either process pre, post, or all data - TO DO
data_processing_choices <- c("pre", "post", "all")
data_processing_choice <- data_processing_choices[3]
# Location of scRNA data - data are organized by aliquot ID
# Note that aliquot ID and sample ID are different since multiple sample types
# (scRNA-seq, scATAC-seq) can come from the same aliquot
# However, in the context of analyzing samples from a single sample type,
# you can consider aliquots and samples basically equivalent
data_dir <- paste0(project_dir, "scRNA_seq_data/")
if (!dir.exists(data_dir)) {dir.create(data_dir)}
output_dir <- paste0(project_dir, "scRNA_seq_data_output/")
if (!dir.exists(output_dir)) {dir.create(output_dir)}
sample_metadata <- read.csv(paste0(project_dir, "current_sample_metadata_minus_8d5be1a4937a7ad3.csv"))
sample_assay_types <- read.csv(paste0(project_dir, "current_set_of_scRNA_and_scATAC_seq_samples.txt"), sep = "\t")
# Look in base scRNA dir to get list of all potential aliquots
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
    # If D-1 is not already in D1.id AND we have scRNA-seq data from both D-1 and D28 aliquots AND we have scATAC-seq data from both D-1 and D28 aliquots,  
    # then add D-1 to D1.id and D28 to D28.id
    d_negative_1_aliquot <- all_samples_associated_with_current_subject[all_samples_associated_with_current_subject$Time_Point == "D-1",]$X_aliquot_id
    d_28_aliquot <- all_samples_associated_with_current_subject[all_samples_associated_with_current_subject$Time_Point == "D28",]$X_aliquot_id
    d_negative_1_aliquot_sample_assays <- sample_assay_types[sample_assay_types$aliquot_name == d_negative_1_aliquot,]
    presence_of_d_negative_1_ATAC <- d_negative_1_aliquot_sample_assays$scATAC_seq
    d_28_aliquot_sample_assays <- sample_assay_types[sample_assay_types$aliquot_name == d_28_aliquot,]
    presence_of_d_28_ATAC <- d_28_aliquot_sample_assays$scATAC_seq
    if (d_negative_1_aliquot %in% D1.id == FALSE & d_negative_1_aliquot %in% aliquot_list & d_28_aliquot %in% aliquot_list & presence_of_d_negative_1_ATAC == "Yes" & presence_of_d_28_ATAC == "Yes") {
      D1.id <- append(D1.id, d_negative_1_aliquot)
      D28.id <- append(D28.id, d_28_aliquot)
    }
  }
}

flu.list <- list()
sample.names <- vector()
# Grab data from all all Day -1 aliquots
print("Grabbing all Day -1 h5 matrices")
for (idx in 1:length(D1.id)) {
  # Grab current sample name
  i <- D1.id[idx]
  # Read in h5 matrix associated with current sample and add sample index to cell column names
  flu.list[[i]] <- Read10X_h5(paste0(data_dir, i, "/outs/filtered_feature_bc_matrix.h5"))
  prefix <- paste0("Sample_", idx, "_D1#")
  sample.names <- append(sample.names, paste0("Sample_", idx, "_D1"))
  print(prefix)
  colnames(flu.list[[i]]) <- paste0(prefix, colnames(flu.list[[i]]))
}

# Do the exact same steps for Day 28 aliquots
print("Grabbing all Day 28 h5 matrices")
for (idx in 1:length(D28.id)) {
  i <- D28.id[idx]
  flu.list[[i]] <- Read10X_h5(paste0(data_dir, i, "/outs/filtered_feature_bc_matrix.h5"))
  prefix <- paste0("Sample_", idx, "_D28#")
  sample.names <- append(sample.names, paste0("Sample_", idx, "_D28"))
  print(prefix)
  colnames(flu.list[[i]]) <- paste0(prefix, colnames(flu.list[[i]]))
}

# Sort sample names so that D1 and D28 are grouped for each sample
sample.names <- sort(sample.names)

# From now on, we'll just use "samples" instead of "aliquots"
# Combine all cells from all samples into one object
all.flu.unbias <- do.call("cbind", flu.list)
dim(all.flu.unbias)

# Create Seurat object from scRNA-seq matrix and calculate various stats
print("Create Seurat Object from all samples")
all.flu.unbias.obj <- CreateSeuratObject(counts = all.flu.unbias,
                                        assay = "RNA",
                                        project = "unbias")
all.flu.unbias.obj <- PercentageFeatureSet(all.flu.unbias.obj,
                                          pattern = "^MT-",
                                          col.name = "percent.mt")
all.flu.unbias.obj <- PercentageFeatureSet(all.flu.unbias.obj,
                                          pattern = "^RPS",
                                          col.name = "percent.rps")
all.flu.unbias.obj <- PercentageFeatureSet(all.flu.unbias.obj,
                                          pattern = "^RPL",
                                          col.name = "percent.rpl")
all.flu.unbias.obj <- PercentageFeatureSet(all.flu.unbias.obj,
                                          pattern = "^HB",
                                          col.name = "percent.hb")
all.flu.unbias.obj <- PercentageFeatureSet(all.flu.unbias.obj,
                                          pattern = "^RP[SL]",
                                          col.name = "percent.rp")

rm(flu.list)
rm(all.flu.unbias)

# Filter cells
print("Filtering cells")
all.flu.unbias.obj <- subset(all.flu.unbias.obj,nFeature_RNA > 900 & nFeature_RNA < 4000 & percent.mt < 20 & percent.hb < 0.4 & percent.rp < 50)

# Label each cell with sample name
all.flu.unbias.obj$sample <- as.vector(sapply(strsplit(colnames(all.flu.unbias.obj), "#"), "[", 1))

# Label each cell with group name (D1 or D28)
group <- all.flu.unbias.obj$sample
group <- gsub("Sample_[1-9]_", "", group)
all.flu.unbias.obj$group <- group

# Process unbiased pooled object and visualize
all.flu.unbias.obj <- NormalizeData(all.flu.unbias.obj)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.flu.unbias.obj <- CellCycleScoring(all.flu.unbias.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
all.flu.unbias.obj$CC.Difference <- all.flu.unbias.obj$S.Score - all.flu.unbias.obj$G2M.Score
all.flu.unbias.obj <- SCTransform(all.flu.unbias.obj,
                                 vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "percent.hb", "CC.Difference"), 
                                 conserve.memory = TRUE,
                                 verbose = TRUE)
all.flu.unbias.obj <- RunPCA(all.flu.unbias.obj, npcs = 50, approx = T, verbose = T)
all.flu.unbias.obj <- RunUMAP(all.flu.unbias.obj, reduction = "pca", dims = 1:30)

# Generate plots for batch detection
print("Generating panel of plots where each plot shows cells for a given sample (batch detection)")
# Each plot will be 2x3 (except the last one, which will be whatever remains)
plot_sets <- chunk(sample.names, 6)
plot_index <- 1
for (plot_set in plot_sets) {
  print(paste0("Printing plot set ", plot_index))
  plot_subset <- subset(x = all.flu.unbias.obj, subset = sample %in% plot_set)
  DimPlot(plot_subset, reduction = "umap", split.by = "sample", group.by = "sample", ncol = 3, raster = FALSE)
  ggsave(paste0(output_dir, "all.flu.unbias.obj.", plot_index, ".PDF"), device = "pdf")
  plot_index <- plot_index + 1
}

# Label each cell as male or female
female.id <- c()
# In particular, this looks at metadata, sees which samples are female,
# and then records that info in female.id
for (idx in 1:length(D1.id)) {
  current_aliquot <- D1.id[idx]
  current_sample <- sample_metadata[sample_metadata$X_aliquot_id == current_aliquot,]
  current_sex <- current_sample$SEX
  if(current_sex == "F") {
    female_sample_name <- paste0("Sample_", idx)
    female.id <- append(female.id, female_sample_name)
  }
}

# We assume by default that sex is male
sex <- rep("Male", length(all.flu.unbias.obj$sample))

# Then we traverse our cells and if a given sample is actually female,
# we change its sex to female
for(sample in female.id){
  if(any(grep(sample, all.flu.unbias.obj$sample))){
    sex[grep(sample, all.flu.unbias.obj$sample)] <- "Female"
  }
}

# Finally, we save the sex of each cell in $sex
all.flu.unbias.obj$sex <- sex

# Use plots from above to determine batches (visually inspect and group like plots)
batch1.id <- c("Sample_1_D1", "Sample_1_D28")
batch2.id <- c("Sample_4_D1")
batch3.id <- c("Sample_5_D1", "Sample_5_D28")
batch4.id <- c("Sample_6_D1", "Sample_6_D28")
# Label samples according to batch
batch <- all.flu.unbias.obj$sample
for (i in 1:length(batch)) {
  if(any(grepl(batch[i], batch1.id))) { batch[i] <- "b1" }
  else if(any(grepl(batch[i], batch2.id))) { batch[i] <- "b2" }
  else if(any(grepl(batch[i], batch3.id))) { batch[i] <- "b3" }
  else if(any(grepl(batch[i], batch4.id))) { batch[i] <- "b4" }
  else { batch[i] <- "b5" }
}
all.flu.unbias.obj$batch <- batch
# Plot all samples grouped by batch
DimPlot(all.flu.unbias.obj, group.by = "batch", reduction = "umap",# cells.highlight = cells,
        label = TRUE, repel = TRUE, shuffle = TRUE, raster = FALSE) +
  ggtitle("Flu all samples, 5 batches")
  ggsave(paste0(output_dir, "all.flu.batches.obj.PDF"), device = "pdf")

all.flu.batch.list <- SplitObject(all.flu.unbias.obj, split.by = "batch")
rm(all.flu.unbias.obj)
# Perform normalization on each batch
for (i in 1:length(all.flu.batch.list)) {
  print(i)
  all.flu.batch.list[[i]] <- NormalizeData(all.flu.batch.list[[i]])
  all.flu.batch.list[[i]] <- CellCycleScoring(all.flu.batch.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  all.flu.batch.list[[i]]$CC.Difference <- all.flu.batch.list[[i]]$S.Score - all.flu.batch.list[[i]]$G2M.Score
  all.flu.batch.list[[i]] <- SCTransform(all.flu.batch.list[[i]],
                                         vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "percent.hb", "CC.Difference"), conserve.memory = TRUE,
                                         verbose = TRUE)
}

# Integrate batches
features <- SelectIntegrationFeatures(object.list = all.flu.batch.list, nfeatures = 2000)
all.flu.batch.list <- PrepSCTIntegration(object.list = all.flu.batch.list, anchor.features = features)
all.flu.batch.list <- lapply(X = all.flu.batch.list, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = all.flu.batch.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
flu.combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(all.flu.batch.list)
DefaultAssay(flu.combined.sct) <- "integrated"
# Run PCA / UMAP / clustering on integrated data
flu.combined.sct <- ScaleData(flu.combined.sct, verbose = T)
flu.combined.sct <- RunPCA(flu.combined.sct, npcs = 30, approx = F, verbose = T)
flu.combined.sct <- RunUMAP(flu.combined.sct, reduction = "pca", dims = 1:30)
flu.combined.sct <- FindNeighbors(flu.combined.sct, reduction = "pca", dims = 1:30, prune.SNN = 1/10, k.param = 20, n.trees = 100)
flu.combined.sct <- FindClusters(flu.combined.sct, resolution = 0.4)
# Plot integrated data
DimPlot(flu.combined.sct, reduction = "umap", label = TRUE, raster = FALSE)
ggsave(paste0(output_dir, "flu.combined.sct.PDF"), device = "pdf")
# We will integrate reference data to assign cell types
scRNA_ref <- LoadH5Seurat("reference/multi.h5seurat")
# Remove certain cell types we're not interested in
idx <- which(scRNA_ref$celltype.l2 %in% c("Doublet", "B intermediate", "CD4 CTL", "gdT", "dnT", "ILC"))
scRNA_ref <- scRNA_ref[,-idx]
idx <- which(scRNA_ref$celltype.l3 == "Treg Naive")
scRNA_ref <- scRNA_ref[,-idx]
# Rename CD4 Proliferating and CD8 Proliferating to T Proliferating
idx <- which(scRNA_ref$celltype.l2 %in% c("CD4 Proliferating", "CD8 Proliferating"))
scRNA_ref$celltype.l2[idx] <- "T Proliferating"
# Add cell type predictions
# Note that k.filter = NA skips filtering anchors stage - when this stage was in place,
# virtually all T Naive (as classified later) cells were filtered out. Not sure if this is good or bad, but
# for now we are going with bad.
flu.anchors <- FindTransferAnchors(reference = scRNA_ref, query = flu.combined.sct,
                                        dims = 1:30, reference.reduction = "pca", k.filter = NA)

predictions <- TransferData(anchorset = flu.anchors, refdata = scRNA_ref$celltype.l2,
                            dims = 1:30)
flu.combined.sct <- AddMetaData(flu.combined.sct, metadata = predictions)
rm(scRNA_ref)

image_dir <- paste0(output_dir, "images/")
if (!dir.exists(image_dir)) {dir.create(image_dir)}
save.image(paste0(image_dir, "integrated_obj_after_predictions.RData"))

# Add cell names as column
cell_names <- rownames(flu.combined.sct@meta.data)
flu.combined.sct <- AddMetaData(flu.combined.sct, metadata = cell_names, col.name = "cell_name")
# Combine cell types where we can't distinguish the difference and add as metadata column
# Note that platelets may actually be save-able, but for now, we're not doing so
Cell_type_combined = flu.combined.sct$predicted.id
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
flu.combined.sct <- AddMetaData(flu.combined.sct, metadata = Cell_type_combined, col.name = 'Cell_type_combined')
flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "CD4 Naive", "T Naive")
flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "CD8 Naive", "T Naive")
flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "NK_CD56bright", "NK")
flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "ASDC", "CD14 Mono")
#flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "cDC", "CD14 Mono")
flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "Eryth", "CD14 Mono")
#flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "HSPC", "CD14 Mono")
#flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "pDC", "CD14 Mono")
#flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "Plasmablast", "CD14 Mono")
flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "Platelet", "CD14 Mono")
flu.combined.sct$Cell_type_combined <- replace(flu.combined.sct$Cell_type_combined, flu.combined.sct$Cell_type_combined == "Treg", "T Naive")

# Add Cell_type_voting column to metadata
cM <- as.matrix(confusionMatrix(Idents(flu.combined.sct), flu.combined.sct$Cell_type_combined))
pre_cluster <- rownames(cM)
max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]
Cell_type_voting <- as.vector(Idents(flu.combined.sct))
for (m in c(1:length(pre_cluster))){
  idxSample <- which(Idents(flu.combined.sct) == pre_cluster[m])
  Cell_type_voting[idxSample] <- max_celltype[m]
}
flu.combined.sct <- AddMetaData(flu.combined.sct, metadata = Cell_type_voting, col.name = "Cell_type_voting")

# Find cluster distributions for removing messy clusters if we want
cluster_distributions <- list()
cluster_predictions <- vector()
cluster_mean_S_score <- vector()
cluster_mean_G2M_score <- vector()
cluster_mean_CC_difference <- vector()
cluster_ids <- vector()
idx <- 1
for (cluster in levels(flu.combined.sct)) {
  idxPass <- which(Idents(flu.combined.sct) %in% cluster)
  cellsPass <- names(flu.combined.sct$orig.ident[idxPass])
  filtered_cluster <- subset(x = flu.combined.sct, subset = cell_name %in% cellsPass)
  cluster_predictions <- append(cluster_predictions, table(filtered_cluster$Cell_type_voting))
  cluster_distributions[[idx]] <- table(filtered_cluster$Cell_type_combined)
  cluster_mean_S_score <- append(cluster_mean_S_score, mean(filtered_cluster$S.Score))
  cluster_mean_G2M_score <- append(cluster_mean_G2M_score, mean(filtered_cluster$G2M.Score))
  cluster_mean_CC_difference <- append(cluster_mean_CC_difference, mean(filtered_cluster$CC.Difference))
  cluster_ids <- append(cluster_ids, cluster)
  idx <- idx + 1
}
names(cluster_predictions) <- paste(levels(flu.combined.sct), "-", names(cluster_predictions))
# Maybe high S score and high G2M score indicate Proliferating cluster?
cell_cycle_df <- data.frame("Cluster" = cluster_ids, "S" = cluster_mean_S_score, "G2M" = cluster_mean_G2M_score, "CC Diff" = cluster_mean_CC_difference)


# Remove the messy clusters
idxPass <- which(Idents(flu.combined.sct) %in% c("3", "11"))
cellsPass <- names(flu.combined.sct$orig.ident[-idxPass])
flu.combined.sct.minus.clusters <- subset(x = flu.combined.sct, subset = cell_name %in% cellsPass)

# Final plot of cell types clustered
DimPlot(flu.combined.sct.minus.clusters, reduction = "umap", group.by = "Cell_type_voting", label = TRUE,
              label.size = 3, repel = TRUE, raster = FALSE) + ggtitle("Cell types")
ggsave(paste0(output_dir, "flu.combined.sct.minus.clusters.PDF"), device = "pdf")

save.image(paste0(output_dir, "integrated_obj_after_final_plot.RData"))

# Create cell type proportion file for MAGICAL
cell_type_proportions_df <- data.frame("Condition" = sub(".*_", "", sample.names), "Sample_name" = sample.names)
total_cell_counts_df <- data.frame("Sample_name" = sample.names)
cell_counts <- vector()
# Find total cell counts for each sample
for (sample_id in sample.names) {
  idxPass <- which(flu.combined.sct.minus.clusters$sample %in% sample_id)
  cellsPass <- names(flu.combined.sct.minus.clusters$orig.ident[idxPass])
  sample_subset <- subset(x = flu.combined.sct.minus.clusters, subset = cell_name %in% cellsPass)
  cell_counts <- append(cell_counts, ncol(sample_subset))
}
total_cell_counts_df <- cbind(total_cell_counts_df, cell_counts)

for (cell_type in unique(flu.combined.sct.minus.clusters$Cell_type_voting)) {
  cell_type_proportions <- vector()
  print(cell_type)
  # Grab cells associated with cell type
  idxPass <- which(flu.combined.sct.minus.clusters$Cell_type_voting %in% cell_type)
  print(length(idxPass))
  cellsPass <- names(flu.combined.sct.minus.clusters$orig.ident[idxPass])
  cells_subset <- subset(x = flu.combined.sct.minus.clusters, subset = cell_name %in% cellsPass)
  for (sample_id in sample.names) {
    # Subset further based on cells associated with sample ID
    idxPass <- which(cells_subset$sample %in% sample_id)
    cellsPass <- names(cells_subset$orig.ident[idxPass])
    if (length(cellsPass) == 0) {
      cell_type_proportions <- append(cell_type_proportions, 0)
    } else {
      sample_subset <- subset(x = cells_subset, subset = cell_name %in% cellsPass)
      cell_counts <- ncol(sample_subset)
      cell_type_proportions <- append(cell_type_proportions, cell_counts / total_cell_counts_df[total_cell_counts_df$Sample_name == sample_id,]$cell_counts)
    }
  }
  temp_df <- data.frame(cell_type_proportions)
  names(temp_df)[names(temp_df) == "cell_type_proportions"] <- cell_type
  cell_type_proportions_df <- cbind(cell_type_proportions_df, temp_df)
}
write.csv(cell_type_proportions_df, file = paste0(output_dir, "RNA_cell_type_proportion.csv"), quote = FALSE, row.names = FALSE)


print("Performing differential expression between groups (D1 and D28) for each cell type")
for (cell_type in unique(flu.combined.sct.minus.clusters$Cell_type_voting)) {
  print(cell_type)
  idxPass <- which(flu.combined.sct.minus.clusters$Cell_type_voting %in% cell_type)
  print(length(idxPass))
  cellsPass <- names(flu.combined.sct.minus.clusters$orig.ident[idxPass])
  cells_subset <- subset(x = flu.combined.sct.minus.clusters, subset = cell_name %in% cellsPass)
  DefaultAssay(cells_subset) <- "SCT"
  cells_subset <- PrepSCTFindMarkers(cells_subset)
  Idents(cells_subset) <- "group"
  #diff_markers <- FindMarkers(cells_subset, group.by = "group", ident.1 = "D1", ident.2 = "D28", logfc.threshold = 0, min.pct = 0.1)
  diff_markers <- FindMarkers(cells_subset, ident.1 = "D28", ident.2 = "D1",logfc.threshold = 0, min.pct = 0, assay = "SCT", recorrect_umi = FALSE)
  #diff_markers$p_val_adj = p.adjust(diff_markers$p_val, method='fdr')
  #diff_markers <- diff_markers[diff_markers$avg_log2FC > 0.1 | diff_markers$avg_log2FC < -0.1,]
  #diff_markers <- diff_markers[diff_markers$p_val_adj < 0.05,]
  print(nrow(diff_markers))
  cell_type <- sub(" ", "_", cell_type)
  write.csv(diff_markers, paste0(output_dir, "D28-vs-D1-degs-", cell_type, ".csv"), quote = FALSE)
}

print("Computing pseudobulk counts for each cell type")
for (cell_type in unique(flu.combined.sct.minus.clusters$Cell_type_voting)) {
  print(cell_type)
  idxPass <- which(flu.combined.sct.minus.clusters$Cell_type_voting %in% cell_type)
  print(length(idxPass))
  cellsPass <- names(flu.combined.sct.minus.clusters$orig.ident[idxPass])
  # Subset again by each sample name
  # Sum by row
  # Keep as vector in list
  cells_subset <- subset(x = flu.combined.sct.minus.clusters, subset = cell_name %in% cellsPass)
  cells_pseudobulk <- list()
  # TODO: Fix this for samples that have 0 cells in cell type
  for (sample_name in sample.names) {
    samples_subset <- subset(x = cells_subset, subset = sample %in% sample_name)
    samples_data <- samples_subset@assays$RNA@counts
    samples_data <- rowSums(as.matrix(samples_data))
    cells_pseudobulk[[sample_name]] <- samples_data
  }
  final_cells_pseudobulk_df <- bind_cols(cells_pseudobulk[1])
  for (idx in 2:length(sample.names)) {
    final_cells_pseudobulk_df <- bind_cols(final_cells_pseudobulk_df, cells_pseudobulk[idx])
  }
  final_cells_pseudobulk_df <- as.data.frame(final_cells_pseudobulk_df)
  rownames(final_cells_pseudobulk_df) <- names(cells_pseudobulk[[1]])
  final_cells_pseudobulk_df <- final_cells_pseudobulk_df[rowSums(final_cells_pseudobulk_df[])>0,]
  cell_type <- sub(" ", "_", cell_type)
  write.csv(final_cells_pseudobulk_df, paste0(output_dir, "pseudo_bulk_RNA_count_", cell_type, ".csv"))
}

# Save overall data
save.image(paste0(output_dir, "integrated_obj_final.RData"))