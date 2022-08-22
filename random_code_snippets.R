output_dir <- "~/scRNA_seq_data_output/"
load(paste0(output_dir, "integrated_obj_final.RData"))

for (res in seq(0, 1, 0.1)) {
  flu.combined.sct <- FindClusters(flu.combined.sct, resolution = res)
}

clustree(flu.combined.sct, prefix = "integrated_snn_res.")
ggsave(paste0(output_dir, "flu.combined.sct.cluster.trees.PDF"), device = "pdf")

proj <- saveArchRProject()
proj <- loadArchRProject(path = paste0(output_dir, "ArchR"))