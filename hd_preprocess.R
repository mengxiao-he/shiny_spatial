# Load packages
library(hdf5r)
library(arrow)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(Banksy)

# load data
localdir <- "/Users/mengxiaohe/10x_Visium_HD_data/mouse_brain/"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))

# setting default assay changes between 8um and 16um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"

# normalize data
object <- NormalizeData(object)

# unsupervised clustering with sketch
object <- FindVariableFeatures(object)
object <- ScaleData(object)

# Select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 3)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

# project back the cluster labels to the full data
object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# # plot umaps with cluster results
# DefaultAssay(object) <- "sketch"
# Idents(object) <- "seurat_cluster.sketched"
# p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)")
#
# # switch to full dataset
# DefaultAssay(object) <- "Spatial.008um"
# Idents(object) <- "seurat_cluster.projected"
# p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)")
#
# plot cluster results spatially
# SpatialDimPlot(object, label = T, repel = T, label.size = 4)
#
# Idents(object) <- "seurat_cluster.projected"
# cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
# SpatialDimPlot(object,
#                cells.highlight = cells[setdiff(names(cells), "NA")],
#                cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
# ) + NoLegend()

# saveRDS(object, file = paste0(localdir, "seurat_hd_mb.rds"))

# Identifying spatially-defined tissue domains
# uses banksy
# run on a server to satisfy memory requirments

# object <- RunBanksy(object,
#                     lambda = 0.8, verbose = TRUE,
#                     assay = "Spatial.008um", slot = "data", features = "variable",
#                     k_geom = 50
# )
# 
# DefaultAssay(object) <- "BANKSY"
# object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
# object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
# object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)

# replace object with banksy results object
object <- readRDS(paste0(localdir, "seurat_hd_mb_banksy_label.rds"))
Idents(object) <- "banksy_cluster"

# plots of banksy results
# SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)

# banksy_cells <- CellsByIdentities(object)
# SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], 
#                cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()

# Integration with single-cell
# load allen brain atlas ref
ref <- readRDS(file = paste0(localdir, "allen_sc_ref/allen_cortex.Rds"))

# remove celltypes with low counts
ref <- subset(ref, idents = "CR", invert = TRUE)

# create the RCTD reference object
Idents(ref) <- "subclass"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)
reference <- Reference(counts, cluster, nUMI)

# create the RCTD query object
counts_hd <- object[["sketch"]]$counts
cells_hd <- colnames(object[["sketch"]])
coords <- GetTissueCoordinates(object)[cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 12)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
object <- AddMetaData(object, metadata = RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
object$first_type <- as.character(object$first_type)
object$first_type[is.na(object$first_type)] <- "Unknown"
object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)

# plot results
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "full_first_type"

# SpatialDimPlot(object, group.by = "full_first_type", label = T, repel = T, label.size = 4)

cells <- CellsByIdentities(object)
# sort and put unknown at end
sorted_names <- sort(names(cells))
sorted_names_without_unknown <- sorted_names[sorted_names != "Unknown"]
sorted_cells <- cells[sorted_names_without_unknown]
# Append the "Unknown" entry to the end if it exists
if ("Unknown" %in% cell_names) {
  sorted_cells[["Unknown"]] <- cells[["Unknown"]]
}
cells <- sorted_cells

# SpatialDimPlot(object, cells.highlight = cells[1:16], cols.highlight = c("#FFFF00", "grey50"), 
#                facet.highlight = T, combine = T, ncol = 4)


saveRDS(object, file = paste0(localdir, "seurat_hd_mb_banksy_rctd.rds"))
