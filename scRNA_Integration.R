setwd("C:/Users/wat2/Desktop/snRNA_seq_data/")

library(Seurat)
library(ggplot2)


# D1.id stores sample IDs for day 1 samples
# D28.id stores sample IDs for day 28 samples
D1.id <- c("13969fce852b59a2")
D28.id <- c()


flu.list <- list()
flu.exp.list <- list()
# Traverse all day 1 samples
for (idx in 1:length(D1.id)) {
  # Grab current sample name
  i <- D1.id[idx]
  # Read in h5 matrix associated with current sample
  flu.list[[i]] <- Read10X_h5(paste0(i, "/outs/filtered_feature_bc_matrix.h5"))
  flu.exp.list[[i]] <- flu.list[[i]]
  prefix <- paste0("Sample_", idx, "_D1#")
  print(prefix)
  colnames(flu.exp.list[[i]]) <- paste0(prefix, colnames(flu.exp.list[[i]]))
}

for (idx in 1:length(D28.id)) {
  i <- D28.id[idx]
  flu.list[[i]] <- Read10X_h5(paste0("path to folder", i, "/outs/filtered_feature_bc_matrix.h5"))
  flu.exp.list[[i]] <- flu.list[[i]]
  prefix <- paste0("Sample_", idx, "_D28#")
  print(prefix)
  colnames(flu.exp.list[[i]]) <- paste0(prefix, colnames(flu.exp.list[[i]]))
}


rm(flu.list)


all.flu.unbias <- do.call("cbind", flu.exp.list)
dim(all.flu.unbias)

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

rm(flu.exp.list)
rm(all.flu.unbias)

all.flu.unbias.obj <- subset(all.flu.unbias.obj,nFeature_RNA > 900 & nFeature_RNA < 4000 & percent.mt < 20 & percent.hb < 0.4 & percent.rp < 50)


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



f.id <- c("Sample_1", "Sample_3", "Sample_4", "Sample_6")
sex <- rep("Male", length(all.flu.unbias.obj$sample))

for(sample in f.id){
    if(any(grep(sample, all.flu.unbias.obj$sample))){
        sex[grep(sample, all.flu.unbias.obj$sample)] <- "Female"
    }
}

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

save.image("output path/integrated_obj.RData")









