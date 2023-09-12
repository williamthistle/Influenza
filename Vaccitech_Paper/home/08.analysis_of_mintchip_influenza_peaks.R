# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

load(paste0(mintchip_dir, "data/k4me1_diffbind_count_run1.Rdata"))

# Rename sample IDs to aliquot IDs everywhere
dbObj.count$class[DBA_ID,] <- mintchip_metadata$Aliquot
for (i in 1:8) {
  names(dbObj.count$masks[[i]]) <- mintchip_metadata$Aliquot
}
dbObj.count$samples$SampleID <- mintchip_metadata$Aliquot
colnames(dbObj.count$called) <- mintchip_metadata$Aliquot
colnames(dbObj.count$binding) <- c("CHR", "START", "END", mintchip_metadata$Aliquot)

dbObj.norm=dba.normalize(dbObj.count,normalize=DBA_NORM_LIB,
                         library=DBA_LIBSIZE_FULL,
                         method=DBA_EDGER,
                         background=TRUE,
                         offsets=TRUE)