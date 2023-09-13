# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

mintchip_data_matrix_file_paths <- list.files(paste0(mintchip_dir, "data/"), full.names = TRUE)
DAS_matrices <- list()
for(current_data_file_path in mintchip_data_matrix_file_paths) {
  load(current_data_file_path)
  current_marker <- unique(dbObj.count$samples$Factor)
  # Rename sample IDs to aliquot IDs everywhere
  dbObj.count$class[DBA_ID,] <- mintchip_metadata$Aliquot
  for (i in 1:8) {
    names(dbObj.count$masks[[i]]) <- mintchip_metadata$Aliquot
  }
  colnames(dbObj.count$class) <- mintchip_metadata$Aliquot
  dbObj.count$samples$SampleID <- mintchip_metadata$Aliquot
  dbObj.count$samples$ControlID <- paste0(dbObj.count$samples$SampleID, "_C")
  dbObj.count$class[7,] <- paste0(dbObj.count$samples$SampleID, "_C")
  colnames(dbObj.count$called) <- mintchip_metadata$Aliquot
  colnames(dbObj.count$binding) <- c("CHR", "START", "END", mintchip_metadata$Aliquot)
  # Use tissue field for our subject ID as proxy (maybe should use replicate field?)
  dbObj.count$class[2,] <- mintchip_metadata$Subject
  dbObj.count$samples$Tissue <- mintchip_metadata$Subject
  # Grab subjects that have matching pre- and post-exposure and subset DiffBind obj to these subjects
  subject_subset <- table(mintchip_metadata$Subject) == 2
  subject_subset <- names(subject_subset[subject_subset])
  aliquot_subset <- mintchip_metadata[mintchip_metadata$Subject %in% subject_subset,]$Aliquot
  all_aliquots <- dbObj.count$samples$SampleID
  full_subject_aliquot_flag <- all_aliquots %in% aliquot_subset
  dbObj.count$masks$full_subject <- full_subject_aliquot_flag
  names(dbObj.count$masks$full_subject) <- all_aliquots
  dbObj.count_full_subj <- dba(dbObj.count, mask = dbObj.count$masks$full_subject)
  # Finish subsetting
  dbObj.count_full_subj$samples <- dbObj.count_full_subj$samples[dbObj.count_full_subj$samples$SampleID %in% aliquot_subset,]
  rownames(dbObj.count_full_subj$samples) <- NULL
  # Normalize data (using "safest" settings from DiffBind manual)
  dbObj.norm <- dba.normalize(dbObj.count_full_subj,normalize=DBA_NORM_NATIVE,
                              method=DBA_DESEQ2,
                              background=TRUE)
  # Find DASs
  dbObj.norm <- dba.contrast(dbObj.norm,design="~Tissue+Condition")
  dbObj.norm <- dba.analyze(dbObj.norm)
  # Find results - contrast number is currently wrong most likely
  results_fc_0 <- dba.report(dbObj.norm, contrast = 1, fold = 0)
  results_fc_0.1 <- dba.report(dbObj.norm, contrast = 1, fold = 0.1)
  results_fc_0.585 <- dba.report(dbObj.norm, contrast = 1, fold = 0.585)
  results_fc_1 <- dba.report(dbObj.norm, contrast = 1, fold = 1)
  results_fc_2 <- dba.report(dbObj.norm, contrast = 1, fold = 2)
  DAS_matrices[[current_marker]] <- list(results_fc_0, results_fc_0.1, results_fc_0.585, 
                                         results_fc_1, results_fc_2)
}
