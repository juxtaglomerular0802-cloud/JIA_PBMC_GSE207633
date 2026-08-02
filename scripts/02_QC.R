library(Seurat)
library(dplyr)
library(ggplot2)

setwd("D:/scRNA/GSE207633_JIA/GSE207633_RAW")


# Check cells before QC
qc_summary_before <- combined@meta.data %>%
  group_by(GSM, Group, Patient) %>%
  summarise(
    n_cells = n(),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    .groups = "drop"
  )

qc_summary_before


# Mitochondrial gene percentage calculation
combined[["percent.mt"]] <- PercentageFeatureSet(
  combined,
  pattern = "^MT-"
)

# Ribisomal gene percentage calculation
combined[["percent.ribo"]] <- PercentageFeatureSet(
  combined,
  pattern = "^RP[SL]"
)

# Hemoglobin gene percentage calculation
combined[["percent.hb"]] <- PercentageFeatureSet(
  combined,
  pattern = "^HB[^(P)]"
)


# Quality Check before QC_Group
VlnPlot(
  combined,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "Group",
  ncol = 3,
  pt.size = 0
)

# Quality Check before QC_Patients individual
VlnPlot(
  combined,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "Patient",
  ncol = 1,
  pt.size = 0
)

# Save as pdf
pdf(
  "QC_violin_by_patient_before_filtering.pdf",
  width = 16,
  height = 12
)

print(
  VlnPlot(
    combined,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    group.by = "Patient",
    ncol = 1,
    pt.size = 0
  )
)

dev.off()

# Check relationship between QC variables
FeatureScatter(
  combined,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

FeatureScatter(
  combined,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

# Check QC stats and save as csv
qc_summary <- combined@meta.data %>%
  group_by(GSM, Group, Patient) %>%
  summarise(
    cells = n(),
    median_feature = median(nFeature_RNA),
    q01_feature = quantile(nFeature_RNA, 0.01),
    q99_feature = quantile(nFeature_RNA, 0.99),
    median_count = median(nCount_RNA),
    q99_count = quantile(nCount_RNA, 0.99),
    median_mt = median(percent.mt),
    q95_mt = quantile(percent.mt, 0.95),
    .groups = "drop"
  )

qc_summary

write.csv(
  qc_summary,
  "QC_summary_before_filtering.csv",
  row.names = FALSE
)


# 1st QC
combined_qc <- subset(
  combined,
  subset =
    nFeature_RNA >= 100 &
    nFeature_RNA <= 6000 &
    percent.mt <= 25
)

# Comparision between before and after QC cell count
ncol(combined)
ncol(combined_qc)

before_cells <- combined@meta.data %>%
  count(GSM, Group, Patient, name = "before")

after_cells <- combined_qc@meta.data %>%
  count(GSM, Group, Patient, name = "after")

qc_cell_counts <- left_join(
  before_cells,
  after_cells,
  by = c("GSM", "Group", "Patient")
) %>%
  mutate(
    after = ifelse(is.na(after), 0, after),
    removed = before - after,
    retained_percent = after / before * 100
  )

qc_cell_counts






# save result as csv
write.csv(
  qc_cell_counts,
  "QC_cell_counts_before_after.csv",
  row.names = FALSE
)

# save Violin Plot
VlnPlot(
  combined_qc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "Group",
  ncol = 3,
  pt.size = 0
)

pdf(
  "QC_violin_by_patient_after_filtering.pdf",
  width = 16,
  height = 12
)

print(
  VlnPlot(
    combined_qc,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    group.by = "Patient",
    ncol = 1,
    pt.size = 0
  )
)

dev.off()


# SaveRDS file
saveRDS(
  combined_qc,
  file = "combined_QC.rds"
)




