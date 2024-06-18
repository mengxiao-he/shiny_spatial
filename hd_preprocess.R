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

saveRDS(object, file = "/Users/mengxiaohe/10x_Visium_HD_data/mouse_brain/seurat_hd_mb.rds")

# Identifying spatially-defined tissue domains
# uses banksy
# run on a server to have memory

object <- RunBanksy(object,
                    lambda = 0.8, verbose = TRUE,
                    assay = "Spatial.008um", slot = "data", features = "variable",
                    k_geom = 50
)

DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)

# plots of banksy results
Idents(object) <- "banksy_cluster"
SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)

banksy_cells <- CellsByIdentities(object)
SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()

