library(Seurat)
library(ggplot2)

# Location of scRNA data - data are organized by aliquot ID
base_scRNA_dir <- "C:/Users/williamthistle/Desktop/snRNA_seq_data/"
output_dir <- "C:/Users/williamthistle/Desktop/snRNA_seq_data_output/"
sample_metadata <- read.csv("C:/Users/williamthistle/Documents/GitHub/Influenza/current_sample_metadata_minus_8d5be1a4937a7ad3.csv")
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

batch1.id <- c("Sample_4_D1", "Sample_4_D28")
batch2.id <- c("Sample_5_D28", "Sample_6_D1")
batch3.id <- c("Sample_1_D1", "Sample_1_D28")
batch <- all.flu.unbias.obj$sample
for (i in 1:length(batch)) {
  if(any(grepl(batch[i], batch1.id))) { batch[i] <- "b1" }
  else if(any(grepl(batch[i], batch2.id))) { batch[i] <- "b2" }
  else if(any(grepl(batch[i], batch3.id))) { batch[i] <- "b3" }
  else { batch[i] <- "b4" }
}
all.flu.unbias.obj$batch <- batch


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
all.flu.batch.list<- lapply(X = all.flu.batch.list, FUN = RunPCA, verbose = FALSE, features = features)
anchors <- FindIntegrationAnchors(object.list = all.flu.batch.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
flu.combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(flu.combined.sct) <- "integrated"
flu.combined.sct <- ScaleData(flu.combined.sct, verbose = T)
flu.combined.sct <- RunPCA(flu.combined.sct, npcs = 30, approx = F, verbose = T)
flu.combined.sct <- RunUMAP(flu.combined.sct, reduction = "pca", dims = 1:30)
flu.combined.sct <- FindNeighbors(flu.combined.sct, reduction = "pca", dims = 1:30, prune.SNN = 1/10, k.param = 20, n.trees = 100)
flu.combined.sct <- FindClusters(flu.combined.sct, resolution = 1)

rm(list=c("all.flu.unbias.obj", "all.flu.batch.list"))

save.image(paste0(output_dir, "integrated_obj.RData"))









