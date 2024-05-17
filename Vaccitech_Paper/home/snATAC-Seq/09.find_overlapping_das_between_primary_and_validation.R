# Setup environment
base_dir <- "~/GitHub/Influenza/Vaccitech_Paper/home/"
source(paste0(base_dir, "00.setup.R"))

overlapping_das_cell_types <- c("CD14_Mono", "CD16_Mono", "cDC","NK", "CD4_Memory", "CD8_Memory", "CD8_Naive", "MAIT", "B", "Proliferating", "CD4_Naive")

for(cell_type in overlapping_das_cell_types) {
  for(analysis_type in c("sc", "final")) {
    for(min_pct in c(0.01, 0.05, 0.1)) {
      print(paste0("Overlap for ", cell_type, " (", analysis_type, ", ", min_pct, ")"))
      placebo_das <- read.table(paste0(scATAC_hvl_placebo_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", min_pct, ".tsv"),
                                     sep = "\t", header = TRUE)
      vaccinated_das <- read.table(paste0(scATAC_hvl_vaccinated_das_dir, "D28-vs-D_minus_1-degs-", cell_type, "-time_point-controlling_for_subject_id_", analysis_type, "_pct_", min_pct, ".tsv"),
                                    sep = "\t", header = TRUE)
      if(analysis_type == "sc") {
        overlapping_das <- intersect(rownames(placebo_das), rownames(vaccinated_das))
        print(paste0("Total overlapping DAS: ", length(overlapping_das)))
        print(paste0("Total placebo DAS: ", nrow(placebo_das)))
        print(paste0("Total vaccinated DAS: ", nrow(vaccinated_das)))
        print(paste0("Percent of placebo DAS that overlap: ", (length(overlapping_das) / nrow(placebo_das) * 100)))
      } else {
        overlapping_das_with_robust <- intersect(placebo_das$Peak_Name, vaccinated_das$Peak_Name)
        print(paste0("Total overlapping DAS (using robust pvalue): ", length(overlapping_das_with_robust)))
        print(paste0("Total placebo DAS (using robust pvalue): ", nrow(placebo_das)))
        print(paste0("Percent of placebo DAS that overlap (using robust pvalue): ", (length(overlapping_das_with_robust) / nrow(placebo_das) * 100)))
        placebo_das <- placebo_das[placebo_das$pseudo_bulk_pval < 0.05,]
        vaccinated_das <- vaccinated_das[vaccinated_das$pseudo_bulk_pval < 0.05,]
        overlapping_das_without_robust <- intersect(placebo_das$Peak_Name, vaccinated_das$Peak_Name)
        print(paste0("Total overlapping DAS (not using robust pvalue): ", length(overlapping_das_without_robust)))
        print(paste0("Total placebo DAS (not using robust pvalue): ", nrow(placebo_das)))
        print(paste0("Percent of placebo DAS that overlap (not using robust pvalue): ", (length(overlapping_das_without_robust) / nrow(placebo_das) * 100)))
      }
    }
  }
}
