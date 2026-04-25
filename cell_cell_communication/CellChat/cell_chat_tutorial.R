# CellChat - to use after you already have processed SCRNASEq data + labeled cells for Cell-Cell interactions


# Note: here we're going to use 3k PBMCs from 10x Genomics in order to create a CellChat object and run the basic analysis.
# The data consists in 3k PBMCs (Peripheral Blood Mononuclear Cells)from a Healthy Donor and is freely available from 10x Genomics.
# PBMCs include lymphocytes (T cells, B cells, and NK cells), monocytes, and dendritic cells.
# Tutorial: https://www.youtube.com/watch?v=VGjddh2yNsQ
pak::pak("sqjin/CellChat")
library(CellChat)
library(Seurat)
library(reticulate)
library(dplyr)

########--Part I: Create a CellChat object 
# 1. Create a CellChat object from a Seurat RDS file
seurat_object <- readRDS("seurat_tutorial_R/output/pbmc3k_final.rds") # load Seurat object (if it's not labeled you can let CellChat label it for you- unrecommended)

data.input <- GetAssayData(seurat_object, assay = "RNA", layer = "data") 

labels <- Idents(seurat_object)

meta <- data.frame(group = labels, row.names = names(labels)) 

cellChat <- createCellChat(object = data.input, meta = meta, group.by = "group") # convert to CellChat

#cellChat <- createCellChat(object = data.input, group.by = "cell.labels") # convert to CellChat

# Note: you can create. CellChat object from other Scanpy h5ad file as well as from .rda file,
# as long as you have the normalized data matrix and the meta data with cell labels.

saveRDS(cellChat, file = "cell_cell_communication/CellChat/output/cellchat_object_pbmc3k.rds")

# 2. Initialize CellChat DB
CellChatDB_obj <- CellChatDB.human

showDatabaseCategory(CellChatDB_obj)
dplyr::glimpse(CellChatDB_obj$interaction) # look at the DB structure

# 3. Add CellChatDB to CellChat object 
CellChatDB_obj.use <- subsetDB(CellChatDB_obj, search = "Secreted Signaling", key = "annotation") # Secreted Signaling = cell-cell comm like hormones, neurotransmitters, or growth factors.
cellChat@DB <- CellChatDB_obj
cellChat <- subsetData(cellChat)

# 4. Identify highly variable ligand/receptor pairs (differentially expressed across prelabeled celltypes)
future::plan("multisession", workers = 4) # run in parallel

# takes a looong time
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)

# The num of highly variable ligand-receptor pairs used for signaling inference is <>

# 5. Calculate communication probablity and infer cellular communication network
cellChat <- computeCommunProb(cellChat, type = "triMean") # trimean, is a measure of a probability distribution's central tendency defined as a weighted average of the distribution's quartiles
saveRDS(cellchat, file = "cell_cell_communication/CellChat/output/cellchat_trimean_pbmc3k.rds")


########--Part II: Analysis of the cellular communication network
# 1. Load RDS
seurat_object <- readRDS("cell_cell_communication/CellChat/output/cellchat_trimean_pbmc3k.rds") 

# 2. Filter-out cell-cell communications where the prob of 2 celltypes interacting is low
cellChat <- filterCommunication(cellChat, min.cells = 10) # the min num of cells in each cell group is 10

# 3. Recalculate cell-cell communication after filtering 
cellChat <- computeCommunProbPathway(cellChat)

# 4. Count the number of edges (cell-cell communications) on a per-network basis
cellChat <- aggregateNet(cellChat)

########--Part III: Vizualizations and pathway selection
# Note: the thicker the edge is, the more likely that the two cells are communicating with each other (high communication probability)

groupSize <- as.numeric(table(cellChat@idents)) # we have set the group to represent all cell types

# plot circle plot
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, lable.edge = F, title.name = "Interactions weights.strength")

# We can look at a single interaction profile for a given pathway - but how can we know which pathway to look at??
# Compute the network centrality scores
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")

# Outgoing heatmap
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing") # "sending" cells

#Incoming heatmap
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming") # "recieving" cells

# show both heatmaps next to each other (the darker the color is - the stringer the signal is)

ht1 + ht2



# Analyze specific pathway