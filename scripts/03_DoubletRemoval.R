if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "scDblFinder",
  "SingleCellExperiment"
))

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(dplyr)
library(ggplot2)

setwd("D:/scRNA/GSE207633_JIA/GSE207633_RAW")

combined_qc <- readRDS("combined_QC.rds")


# Transform Seurat obj as SingleCellExperiment
combined_qc <- JoinLayers(combined_qc)

sce <- as.SingleCellExperiment(
    combined_qc,
    assay = "RNA"
)

# Check whether there is "counts" matrix
assayNames(sce)

# Execute scDblFinder
set.seed(1234)

sce <- scDblFinder(
  sce,
  samples = "GSM"
)

table(sce$scDblFinder.class)


# Put result in Seruat obj
identical(
  colnames(combined_qc),
  colnames(sce)
)

combined_qc$scDblFinder.score <-
  colData(sce)$scDblFinder.score

combined_qc$scDblFinder.class <-
  colData(sce)$scDblFinder.class

head(
  combined_qc@meta.data[
    ,
    c(
      "GSM",
      "Group",
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt",
      "scDblFinder.score",
      "scDblFinder.class"
    )
  ]
)

# Calculate Doublet proportion by individual samples
doublet_summary <- combined_qc@meta.data %>%
  dplyr::count(
    GSM,
    Group,
    scDblFinder.class,
    name = "n_cells"
  ) %>%
  group_by(GSM, Group) %>%
  mutate(
    total_cells = sum(n_cells),
    percent = n_cells / total_cells * 100
  ) %>%
  ungroup()

doublet_summary

# Calculate Doublet proportion by individual samples_Simple version
doublet_by_sample <- combined_qc@meta.data %>%
  group_by(GSM, Group, Patient) %>%
  summarise(
    total_cells = n(),
    doublets = sum(scDblFinder.class == "doublet"),
    singlets = sum(scDblFinder.class == "singlet"),
    doublet_percent = doublets / total_cells * 100,
    .groups = "drop"
  )

doublet_by_sample


# save as csv
write.csv(
  doublet_by_sample,
  "scDblFinder_doublet_summary_by_sample.csv",
  row.names = FALSE
)

# Visualization
ggplot(
  doublet_by_sample,
  aes(
    x = GSM,
    y = doublet_percent,
    fill = Group
  )
) +
  geom_col() +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  ) +
  labs(
    x = NULL,
    y = "Doublet percentage"
  )

# Save Visualized graph
pdf(
  "scDblFinder_doublet_rate_by_sample.pdf",
  width = 12,
  height = 6
)

print(
  ggplot(
    doublet_by_sample,
    aes(
      x = GSM,
      y = doublet_percent,
      fill = Group
    )
  ) +
    geom_col() +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1
      )
    ) +
    labs(
      x = NULL,
      y = "Doublet percentage"
    )
)

dev.off()


# Remove Doublet
combined_singlet <- subset(
  combined_qc,
  subset = scDblFinder.class == "singlet"
)

ncol(combined_qc)
ncol(combined_singlet)

before_doublet <- combined_qc@meta.data %>%
  count(GSM, Group, Patient, name = "before_doublet_removal")

after_doublet <- combined_singlet@meta.data %>%
  count(GSM, Group, Patient, name = "after_doublet_removal")

doublet_cell_counts <- left_join(
  before_doublet,
  after_doublet,
  by = c("GSM", "Group", "Patient")
) %>%
  mutate(
    after_doublet_removal = ifelse(
      is.na(after_doublet_removal),
      0,
      after_doublet_removal
    ),
    removed_doublets =
      before_doublet_removal - after_doublet_removal,
    retained_percent =
      after_doublet_removal /
      before_doublet_removal *
      100
  )

doublet_cell_counts

# Save as csv file
write.csv(
  doublet_cell_counts,
  "cell_counts_after_doublet_removal.csv",
  row.names = FALSE
)


# remove unnecessary obj
rm(sce)
gc()

# Save RDS
saveRDS(
  combined_singlet,
  file = "GSE207633_QC_doublet_removed_before_normalization.rds"
)




