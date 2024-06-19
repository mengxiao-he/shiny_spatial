# global.R

# Check if 'Seurat' package is loaded and load it if necessary
if (!require(Seurat, quietly = TRUE)) {
  library(Seurat)
}

# Check if 'mb_hd' is already loaded to avoid reloading
if (!exists("mb_hd")) {
  mb_hd <- readRDS("path/to/your/initial_file.rds")
}