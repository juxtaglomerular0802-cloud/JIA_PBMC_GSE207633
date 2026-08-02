library(Seurat)
library(stringr)
library(dplyr)

# Set Directory
setwd("D:/scRNA/GSE207633_JIA/GSE207633_RAW")

# Check Files
files <- list.files(
  pattern = "\\.h5$",
  full.names = TRUE
)


# Read Files (may take long time)
sample.list <- lapply(files, function(f){

  counts <- Read10X_h5(f)

  obj <- CreateSeuratObject(
    counts = counts,
    min.cells = 3,
    min.features = 200
  )

  obj
})


# Extract sample names
sample.names <- basename(files)

sample.names

# Add Metadata (ex. sample GSM000000, Group MAS, Patient 000)
for(i in seq_along(sample.list)){

  fname <- basename(files[i])

  parts <- str_split(fname, "_", simplify = TRUE)

  sample.list[[i]]$GSM <- parts[1]
  sample.list[[i]]$Group <- parts[2]
  sample.list[[i]]$Patient <- parts[3]

}

# Modify cell name for prevent duplication
for(i in seq_along(sample.list)){

  sample.list[[i]] <- RenameCells(
    sample.list[[i]],
    add.cell.id = sample.list[[i]]$Patient[1]
  )

}

# Merge Samples
combined <- merge(
  sample.list[[1]],
  y = sample.list[-1]
)


# Save as Seurat object using in R
saveRDS(
  combined,
  file = "combined_raw.rds"
)


# IF want to load file again
combined <- readRDS("combined_raw.rds")
