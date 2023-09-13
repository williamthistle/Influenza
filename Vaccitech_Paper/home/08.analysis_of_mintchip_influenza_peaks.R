# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

load(paste0(mintchip_dir, "data/k4me1_diffbind_count_run1.Rdata"))

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
# Use tissue field for subject ID as proxy (maybe should use replicate?)
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




dbObj.norm <- dba.normalize(dbObj.count,normalize=DBA_NORM_NATIVE,
                         method=DBA_DESEQ2,
                         background=TRUE)

dbObj.norm <- dba.contrast(dbObj.norm,design="~Tissue+Condition")
dbObj.norm <- dba.analyze(dbObj.norm)