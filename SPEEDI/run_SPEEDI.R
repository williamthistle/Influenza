home_dir <- "~/SPEEDI"
source(paste0(home_dir, "/prototype_utils.R"))
source(paste0(home_dir, "/prototype_API.R"))

data_path <- "/data/home/wat2/multiome/snRNA_seq_data/"
sample_id_list <- c("001c96c7a4adc441", "101f214ee96a731d")
output_dir <- "/data/home/wat2/multiome/snRNA_seq_data_output/"

all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
sc_obj <- FilterRawData(all_sc_exp_matrices, human = TRUE)
sc_obj <- InitialProcessing(sc_obj, human = TRUE)
sc_obj <- InferBatches(sc_obj)
sc_obj <- IntegrateByBatch(sc_obj)
sc_obj <- VisualizeIntegration(sc_obj)
reference <- LoadReference("PBMC", human = TRUE)
sc_obj <- MapCellTypes(sc_obj, reference)

DimPlot(sc_obj, reduction = "umap", group.by = "predicted_celltype_majority_vote", label = TRUE,
        label.size = 3, repel = TRUE, raster = FALSE) + 
  labs(title = "scRNA-seq Data Integration \n (X Samples, X Cells)") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(output_dir, "sc_obj.png"), device = "png", dpi = 300)