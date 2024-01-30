addArchRGenome("hg38")
atac_proj <- pseudo_bulk_replicates_and_call_peaks(atac_proj)
# Create peak matrix (matrix containing insertion counts within our merged peak set) for differential accessibility
# calculations
atac_proj <- addPeakMatrix(atac_proj)

# Convert ArchR to Signac
peak_matrix <- getPeakMatrix(atac_proj)
annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38")

seurat_atac <- ArchR2Signac(
  ArchRProject = atac_proj,
  refversion = "hg38",
  #samples = samplelist, # list of samples in the ArchRProject (default will use ArchRProject@cellColData$Sample but another list can be provided)
  fragments_dir = data_path,
  pm = peak_matrix, # peak matrix from getPeakMatrix()
  fragments_fromcellranger = "Yes", # fragments_fromcellranger This is an Yes or No selection ("NO" | "N" | "No" or "YES" | "Y" | "Yes")
  fragments_file_extension = NULL, # Default - NULL: File_Extension for fragments files (typically they should be '.tsv.gz' or '.fragments.tsv.gz')
  annotation = annotations # annotation from getAnnotation()
)


# Returns TRUE as expected - samples are in correct order
all.equal(colnames(seurat_atac), gsub('#', '_', rownames(atac_proj@cellColData)))

# Add gene score matrix from ArchR to Seurat object (could recompute it using Signac, but we'll stick with ArchR matrix for now)
gsm <- getGeneScoreMatrix(ArchRProject = atac_proj, SeuratObject = seurat_atac)
seurat_atac[['RNA']] <- CreateAssayObject(counts = gsm)

# Add UMAP
seurat_atac <- addDimRed(
  ArchRProject = atac_proj,
  SeuratObject = seurat_atac,
  addUMAPs = "UMAP",
  reducedDims = "Harmony"
)

saveRDS(seurat_atac, file = paste0(ATAC_output_dir, "ALL_seurat.RDS"))

# Method to print a well-organized UMAP plot for our snRNA-seq data
print_UMAP_ATAC_Seurat <- function(sc_obj, file_name, group_by_category = NULL, output_dir = getwd(), log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  sample_count <- length(unique(sc_obj$Sample))
  cell_count <- length(sc_obj$Sample)
  current_title <- paste0("RNA Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  if(!is.null(group_by_category)) {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, label = TRUE,
                         label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", label = TRUE,
                         label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  ggplot2::ggsave(paste0(output_dir, file_name), plot = p, device = "png", dpi = 300)
  return(TRUE)
}

# Looks good!
print_UMAP_ATAC_Seurat(seurat_atac, "seurat_atac_test.png", group_by_category = "Cell_type_voting", output_dir = ATAC_output_dir)

cell_names <- rownames(seurat_atac@meta.data)
seurat_atac <- Seurat::AddMetaData(seurat_atac, metadata = cell_names, col.name = "cell_name")

messy_clusters <- c("4", "7", "10", "13", "15", "16", "18", "20", "24", "26", "28", "29")
idxPass <- which(seurat_atac$seurat_clusters %in% messy_clusters)
cellsPass <- names(seurat_atac$orig.ident[-idxPass])
seurat_atac_minus_clusters <- subset(x = seurat_atac, subset = cell_name %in% cellsPass)

print_UMAP_ATAC_Seurat(seurat_atac_minus_clusters, "seurat_atac_minus_clusters_test.png", group_by_category = "Cell_type_voting", output_dir = ATAC_output_dir)