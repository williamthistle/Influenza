library(Seurat)
library(ggplot2)
library(dplyr)
library(clustree)

chunk <- function(x, n) (mapply(function(a, b) (x[a:b]), seq.int(from=1, to=length(x), by=n), pmin(seq.int(from=1, to=length(x), by=n)+(n-1), length(x)), SIMPLIFY=FALSE))
# Location of scRNA data - data are organized by aliquot ID
# Note that aliquot ID and sample ID are different since multiple sample types
# (scRNA-seq, scATAC-seq) can come from the same aliquot
# However, in the context of analyzing samples from a single sample type,
# you can consider aliquots and samples basically equivalent
base_scRNA_dir <- "~/scRNA_seq_data/"
output_dir <- "~/scRNA_seq_data_output/"
sample_metadata <- read.csv("~/current_sample_metadata_minus_8d5be1a4937a7ad3.csv")
sample_assay_types <- read.csv("~/current_set_of_scRNA_and_scATAC_seq_samples.txt", sep = "\t")
# Look in base scRNA dir to get list of all potential aliquots
# We will only use aliquots that are paired (D-1 and D28)
aliquot_list <- list.dirs(base_scRNA_dir, recursive = FALSE)
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
  flu.list[[i]] <- Read10X_h5(paste0(base_scRNA_dir, i, "/outs/filtered_feature_bc_matrix.h5"))
  prefix <- paste0("Sample_", idx, "_D1#")
  sample.names <- append(sample.names, paste0("Sample_", idx, "_D1"))
  print(prefix)
  colnames(flu.list[[i]]) <- paste0(prefix, colnames(flu.list[[i]]))
}

# Do the exact same steps for Day 28 aliquots
print("Grabbing all Day 28 h5 matrices")
for (idx in 1:length(D28.id)) {
  i <- D28.id[idx]
  flu.list[[i]] <- Read10X_h5(paste0(base_scRNA_dir, i, "/outs/filtered_feature_bc_matrix.h5"))
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
for (i in 1:length(all.flu.batch.list)) {
  print(i)
  all.flu.batch.list[[i]] <- NormalizeData(all.flu.batch.list[[i]])
  all.flu.batch.list[[i]] <- CellCycleScoring(all.flu.batch.list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  all.flu.batch.list[[i]]$CC.Difference <- all.flu.batch.list[[i]]$S.Score - all.flu.batch.list[[i]]$G2M.Score
  all.flu.batch.list[[i]] <- SCTransform(all.flu.batch.list[[i]],
                                         vars.to.regress = c("percent.mt", "percent.rps", "percent.rpl", "percent.hb", "CC.Difference"), conserve.memory = TRUE,
                                         verbose = TRUE)
}

features <- SelectIntegrationFeatures(object.list = all.flu.batch.list, nfeatures = 2000)
all.flu.batch.list <- PrepSCTIntegration(object.list = all.flu.batch.list, anchor.features = features)
all.flu.batch.list <- lapply(X = all.flu.batch.list, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = all.flu.batch.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
flu.combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(all.flu.batch.list)
DefaultAssay(flu.combined.sct) <- "integrated"
flu.combined.sct <- ScaleData(flu.combined.sct, verbose = T)
flu.combined.sct <- RunPCA(flu.combined.sct, npcs = 30, approx = F, verbose = T)
flu.combined.sct <- RunUMAP(flu.combined.sct, reduction = "pca", dims = 1:30)
flu.combined.sct <- FindNeighbors(flu.combined.sct, reduction = "pca", dims = 1:30, prune.SNN = 1/10, k.param = 20, n.trees = 100)
flu.combined.sct <- FindClusters(flu.combined.sct, resolution = 0.2)
# Plot integrated data
DimPlot(flu.combined.sct, reduction = "umap", label = TRUE, raster = FALSE)
ggsave(paste0(output_dir, "flu.combined.sct.PDF"), device = "pdf")
# Save overall data
save.image(paste0(output_dir, "integrated_obj_final.RData"))
