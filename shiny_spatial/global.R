# global.R

set_data_path <- "/Users/mengxiaohe/10x_Visium_HD_data/mouse_brain/"

avail_datasets <- NULL

# Check if 'Seurat' package is loaded and load it if necessary
if (!require(Seurat, quietly = TRUE)) {
  library(Seurat)
}

# Check if 'mb_hd' is already loaded to avoid reloading
if (!exists("mb_hd")) {
  mb_hd <- readRDS(file = paste0(set_data_path, "seurat_hd_mb_banksy_rctd.rds"))
  avail_datasets <- c(avail_datasets, "mb_hd")
} else {
  avail_datasets <- c(avail_datasets, "mb_hd")
}