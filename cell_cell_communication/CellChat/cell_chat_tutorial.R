# Note: here we're going to use 3k PBMCs from 10x Genomics in order to create a CellChat object and run the basic analysis.
# The data consists in 3k PBMCs (Peripheral Blood Mononuclear Cells)from a Healthy Donor and is freely available from 10x Genomics.
# PBMCs include lymphocytes (T cells, B cells, and NK cells), monocytes, and dendritic cells.
# Tutorial: https://www.youtube.com/watch?v=mOSG31aaj14&list=PLOLdjuxsfI4M15_EeIJY7S63_-Wce-3QY 
pak::pak("sqjin/CellChat")
library(CellChat)
library(Seurat)
library(reticulate)

########--Part I: Create a CellChat object 
# 1. Create a CellChat object from a Seurat RDS file
seurat_object <- readRDS("seurat_tutorial_R/output/pbmc3k_final.rds")

data.input <- GetAssayData(seurat_object, assay = "RNA", layer = "data") 

labels <- Idents(seurat_object)

meta <- data.frame(group = labels, row.names = names(labels)) 

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# Note: you can create. CellChat object from other Scanpy h5ad file as well as from .rda file,
# as long as you have the normalized data matrix and the meta data with cell labels.

saveRDS(cellchat, file = "cell_cell_communication/CellChat/output/cellchat_object_pbmc3k.rds")