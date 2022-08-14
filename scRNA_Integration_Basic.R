library(Seurat)
library(ggplot2)

# Location of scRNA data - data are organized by aliquot ID
base_scRNA_dir <- "~/snRNA_seq_data/"
output_dir <- "~/snRNA_seq_data_output/"
sample_metadata <- read.csv("~/current_sample_metadata_minus_8d5be1a4937a7ad3.csv")
aliquot_list <- list.dirs(base_scRNA_dir, recursive = FALSE)
aliquot_list <- strsplit(aliquot_list, "/")
aliquot_list <- unlist(lapply(aliquot_list, tail, n = 1L))
# D1.id stores sample IDs for day -1 samples
# D28.id stores sample IDs for day 28 samples
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
    # If D-1 is not already in D1.id AND we have snRNA-seq data from both D-1 and D28 aliquots,  
    # then add D-1 to D1.id and D28 to D28.id
    d_negative_1_aliquot <- all_samples_associated_with_current_subject[all_samples_associated_with_current_subject$Time_Point == "D-1",]$X_aliquot_id
    d_28_aliquot <- all_samples_associated_with_current_subject[all_samples_associated_with_current_subject$Time_Point == "D28",]$X_aliquot_id
    if (d_negative_1_aliquot %in% D1.id == FALSE & d_negative_1_aliquot %in% aliquot_list & d_28_aliquot %in% aliquot_list) {
      D1.id <- append(D1.id, d_negative_1_aliquot)
      D28.id <- append(D28.id, d_28_aliquot)
    }
  }
}

flu.list <- list()
# Traverse all day -1 samples
for (idx in 1:length(D1.id)) {
  # Grab current sample name
  i <- D1.id[idx]
  # Read in h5 matrix associated with current sample and add sample index to cell column names
  flu.list[[i]] <- Read10X_h5(paste0(base_scRNA_dir, i, "/outs/filtered_feature_bc_matrix.h5"))
  prefix <- paste0("Sample_", idx, "_D1#")
  print(prefix)
  colnames(flu.list[[i]]) <- paste0(prefix, colnames(flu.list[[i]]))
}

# Do the exact same steps for day 28 samples
for (idx in 1:length(D28.id)) {
  i <- D28.id[idx]
  flu.list[[i]] <- Read10X_h5(paste0(base_scRNA_dir, i, "/outs/filtered_feature_bc_matrix.h5"))
  prefix <- paste0("Sample_", idx, "_D28#")
  print(prefix)
  colnames(flu.list[[i]]) <- paste0(prefix, colnames(flu.list[[i]]))
}

# Combine all cells from all samples into one object
all.flu.unbias <- do.call("cbind", flu.list)
dim(all.flu.unbias)

# Create Seurat object from snRNA-seq matrix and calculate various stats
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

pdf(paste0(output_dir, "all.flu.unbias.obj.PDF"))
DimPlot(all.flu.unbias.obj, reduction = "umap")
dev.off()
