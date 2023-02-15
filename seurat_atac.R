library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)

counts <- Read10X_h5("C:/Users/wat2/Desktop/ATAC/filtered_feature_bc_matrix.h5")
counts <- counts$Peaks

metadata <- read.csv(
  file = "C:/Users/wat2/Desktop/ATAC/per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'C:/Users/wat2/Desktop/ATAC/atac_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations

# Compute QC
# Nucleosome signal
pbmc <- NucleosomeSignal(object = pbmc)
# TSS Enrichment
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
# Fraction of reads in peaks (since my CellRanger metadata doesn't seem to include this info)
total_fragments <- CountFragments("C:/Users/wat2/Desktop/ATAC/atac_fragments.tsv.gz")
pbmc_fragments <- total_fragments[total_fragments$CB %in% colnames(pbmc),]
pbmc_fragments <-
pbmc$fragments <- total_fragments[colnames(pbmc), "frequency_count"]
pbmc <- FRiP(
  object = pbmc,
  assay = 'peaks',
  total.fragments = 'fragments'
)
# blacklist ratio (since my CellRanger metadata doesn't seem to include this info)
pbmc$blacklist_fraction <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg38
)

