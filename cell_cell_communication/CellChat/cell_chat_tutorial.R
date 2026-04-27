# CellChat - to use after you already have processed SCRNASEq data + labeled cells for Cell-Cell interactions
# scRNA-seq data of human skin from patients with atopic dermatitis
# from https://figshare.com/projects/Example_data_for_cell-cell_communication_analysis_using_CellChat/157272

pak::pak("sqjin/CellChat")
library(CellChat)
library(Seurat)
library(reticulate)
library(dplyr)
library(ComplexHeatmap)

# #### (A) Generate data input starting from a count data matrix:
# # Upload the count data matrix in a .rda or other format:
# load("data/data_humanSkin_CellChat.rda")
# # (i) Obtain the normalized data matrix:
# data.input = data_humanSkin$data
# # (ii) Generate a data frame with row names containing cell metadata:
# meta = data_humanSkin$meta
# # (iii)Subset the data from one condition for further analysis:
# cell.use = rownames(meta)[meta$condition == "LS"]
# data.input = data.input[, cell.use]
# meta = meta[cell.use,]


# #(A) Create a CellChat object from the digital gene expression matrix and cell label
# # (you can also create a CellChat object from a Seurat object, SingleCellExperiment object, or Anndata object)
# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")


# # select the L–R interaction database (CellChatDB.human, in this case)
# CellChatDB <- CellChatDB.human
# showDatabaseCategory(CellChatDB)
# dplyr::glimpse(CellChatDB$interaction)


# #(A) Select a specific subset of the CellChatDB database
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # cell-cell where cells release molecules
# # #(B) Select all CellChatDB database except for nonprotein signaling
# # CellChatDB.use <- subsetDB(CellChatDB)
# # #(C) Select all CellChatDB database
# # CellChatDB.use <- CellChatDB


# # Set the selected DB in the CellChat object and then subset the expression data matrix
# # using genes relevant to the selected L–R pairs
# cellchat@DB <- CellChatDB.use
# cellchat <- subsetData(cellchat)

# # Identify over-expressed ligands or receptors in each cell group to infer the cell state-specific communications.
# future::plan("multisession", workers = 4)
# cellchat <- identifyOverExpressedGenes(cellchat)

# # For each overexpressed ligand and receptor obtained identify over-expressed L–R interactions if either its associated ligand or receptor is over expressed:
# cellchat <- identifyOverExpressedInteractions(cellchat)

# # Infer cell–cell communication at a L–R pair level
# cellchat <- computeCommunProb(cellchat, type = "triMean", trim = NULL, raw.use = TRUE)


# # Filter the cell–cell communication, based on the number of cells in each group. By default,
# # the minimum number of cells required in each cell group for cell–cell communication is 10.
# cellchat <- filterCommunication(cellchat, min.cells = 10)


# # CellChat computes the communication probability at the signaling pathway level by summarizing the
# # communication probabilities of all L–R pairs associated with each signaling pathway
# cellchat <- computeCommunProbPathway(cellchat)

# # Calculate the aggregated cell–cell communication network. 
# #(A) Perform calculation across all the cell groups
# cellchat <- aggregateNet(cellchat)
# # #(B) Perform calculation across a subset of cell groups
# # sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB")
# # targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC")
# # cellchat <- aggregateNet(cellchat, sources.use = sources.use, targets.use = targets.use)


# # Export the CellChat object together with the inferred cell–cell communication networks
# # and save them as a .rds file.
# saveRDS(cellchat, file = "data/cellchat_humanSkin_LS.rds")


############################################################
# CellChat Tutorial – Visualizations Only
############################################################
# Load the saved CellChat object
cellchat <- readRDS("data/cellchat_humanSkin_LS.rds")

############################
# 0. Available pathways
############################
pathways.show.all <- cellchat@netP$pathways
pathways.show <- c("MIF")   # choose one valid pathway

############################################################
# 1. Visualize communication network of a signaling pathway
############################################################

# Circle plot: shows the overall strength and direction of signaling cell-cell for a selected pathway.
netVisual_aggregate(
  cellchat,
  signaling = pathways.show,
  layout = "circle"
)

# Hierarchy plot: shows directional communication from sender cell groups to specified receiver groups in a layered source-to-target structure.
vertex.receiver <- c(1,2,3)   # target groups on left panel
netVisual_aggregate(
  cellchat,
  signaling = pathways.show,
  layout = "hierarchy",
  vertex.receiver = vertex.receiver
)

# Chord diagram
netVisual_aggregate(
  cellchat,
  signaling = pathways.show,
  layout = "chord"
)

# Heatmap (use ONE pathway at a time)
netVisual_heatmap(
  cellchat,
  signaling = pathways.show,
  color.heatmap = "Reds"
)

############################################################
# 2. Visualize individual ligand-receptor pairs
############################################################

# Contribution of LR pairs within pathway
netAnalysis_contribution(
  cellchat,
  signaling = pathways.show
)

# Extract enriched LR pairs
pairLR <- extractEnrichedLR(
  cellchat,
  signaling = pathways.show,
  geneLR.return = FALSE
)

# Visualize first LR pair
LR.show <- pairLR[1,]

netVisual_individual(
  cellchat,
  signaling = pathways.show,
  pairLR.use = LR.show,
  layout = "circle"
)

############################################################
# 3. Bubble plots – multiple pathways / interactions
############################################################

# All interactions from source cluster 4 to targets 5:8
netVisual_bubble(
  cellchat,
  sources.use = 4,
  targets.use = 5:8,
  remove.isolate = FALSE
)

# Restrict to selected pathways
netVisual_bubble(
  cellchat,
  sources.use = 4,
  targets.use = 5:8,
  signaling = c("MIF","GALECTIN"),
  remove.isolate = FALSE
)

# Restrict to LR pairs
pairLR.use <- extractEnrichedLR(
  cellchat,
  signaling = c("MIF","GALECTIN")
)

netVisual_bubble(
  cellchat,
  sources.use = 4,
  targets.use = 5:8,
  pairLR.use = pairLR.use,
  remove.isolate = TRUE
)

############################################################
# 4. Chord diagram of LR interactions
############################################################

netVisual_chord_gene(
  cellchat,
  sources.use = 4,
  targets.use = 5:8
)

# Pathway-level chord diagram
netVisual_chord_gene(
  cellchat,
  sources.use = 4,
  targets.use = 5:8,
  slot.name = "netP"
)

############################################################
# 5. Gene expression visualization
############################################################

# Built-in CellChat plot
plotGeneExpression(
  cellchat,
  signaling = "MIF",
  enriched.only = TRUE,
  type = "violin"
)

# Using Seurat
genes.use <- extractEnrichedLR(
  cellchat,
  signaling = "MIF",
  geneLR.return = TRUE
)$geneLR

Seurat::VlnPlot(
  seu_obj,
  features = genes.use
)

############################################################
# 6. Signaling role / centrality analysis
############################################################

cellchat <- netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)

# Heatmap of sender/receiver roles
netAnalysis_signalingRole_network(
  cellchat,
  signaling = pathways.show
)

# Scatter plot of dominant senders/receivers
netAnalysis_signalingRole_scatter(cellchat)

# Outgoing signaling heatmap
ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing"
)

# Incoming signaling heatmap
ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming"
)

ht1 + ht2

############################################################
# 7. Global communication patterns
############################################################

library(NMF)
library(ggalluvial)

# Determine K
selectK(cellchat, pattern = "outgoing")

# Learn patterns
cellchat <- identifyCommunicationPatterns(
  cellchat,
  pattern = "outgoing",
  k = 4
)

# River plot
netAnalysis_river(
  cellchat,
  pattern = "outgoing"
)

# Dot plot
netAnalysis_dot(
  cellchat,
  pattern = "outgoing"
)

############################################################
# 8. Network similarity / pathway clustering
############################################################

cellchat <- computeNetSimilarity(
  cellchat,
  type = "functional"
)

cellchat <- netEmbedding(
  cellchat,
  type = "functional"
)

cellchat <- netClustering(
  cellchat,
  type = "functional"
)

netVisual_embedding(
  cellchat,
  type = "functional"
)

############################################################
# 9. Save object
############################################################

saveRDS(cellchat, file = "data/cellchat_results.rds")








# ########--Part I: Create a CellChat object 
# # 1. Create a CellChat object from a Seurat RDS file
# seurat_object <- readRDS("seurat_tutorial_R/output/pbmc3k_final.rds") # load Seurat object (if it's not labeled you can let CellChat label it for you- unrecommended)

# data.input <- GetAssayData(seurat_object, assay = "RNA", layer = "data") 

# labels <- Idents(seurat_object)

# meta <- data.frame(group = labels, row.names = names(labels)) 

# cellChat <- createCellChat(object = data.input, meta = meta, group.by = "group") # convert to CellChat

# #cellChat <- createCellChat(object = data.input, group.by = "cell.labels") # convert to CellChat

# # Note: you can create. CellChat object from other Scanpy h5ad file as well as from .rda file,
# # as long as you have the normalized data matrix and the meta data with cell labels.

# saveRDS(cellChat, file = "cell_cell_communication/CellChat/output/cellchat_object_pbmc3k.rds")

# # 2. Initialize CellChat DB
# CellChatDB_obj <- CellChatDB.human

# showDatabaseCategory(CellChatDB_obj)
# dplyr::glimpse(CellChatDB_obj$interaction) # look at the DB structure

# # 3. Add CellChatDB to CellChat object 
# CellChatDB_obj.use <- subsetDB(CellChatDB_obj, search = "Secreted Signaling", key = "annotation") # Secreted Signaling = cell-cell comm like hormones, neurotransmitters, or growth factors.
# cellChat@DB <- CellChatDB_obj
# cellChat <- subsetData(cellChat)

# # 4. Identify highly variable ligand/receptor pairs (differentially expressed across prelabeled celltypes)
# future::plan("multisession", workers = 4) # run in parallel

# # takes a long time
# cellChat <- identifyOverExpressedGenes(cellChat)
# cellChat <- identifyOverExpressedInteractions(cellChat)

# # The num of highly variable ligand-receptor pairs used for signaling inference is <>

# # 5. Calculate communication probablity and infer cellular communication network (takes a loooong time)
# cellChat <- computeCommunProb(cellChat, type = "triMean") # trimean, is a measure of a probability distribution's central tendency defined as a weighted average of the distribution's quartiles


# ########--Part II: Analysis of the cellular communication network
# # 1. Load
# #### load saved comm network

# # 2. Filter-out cell-cell communications where the prob of 2 celltypes interacting is low
# cellChat <- filterCommunication(cellChat, min.cells = 10) # the min num of cells in each cell group is 10

# # 3. Recalculate cell-cell communication after filtering 
# cellChat <- computeCommunProbPathway(cellChat)

# # 4. Count the number of edges (cell-cell communications) on a per-network basis
# cellChat <- aggregateNet(cellChat)

# ########--Part III: Vizualizations and pathway selection
# # Note: the thicker the edge is, the more likely that the two cells are communicating with each other (high communication probability)

# groupSize <- as.numeric(table(cellChat@idents)) # we have set the group to represent all cell types

# # plot circle plot
# pdf(
#   file = "cell_cell_communication/CellChat/output/circle_plot.pdf",
#   width = 10,
#   height = 10
# )

# netVisual_circle(
#   cellChat@net$weight,
#   vertex.weight = groupSize,
#   weight.scale = TRUE,
#   label.edge = FALSE,
#   title.name = "Interactions weights.strength"
# )

# dev.off()

# # We can look at a single interaction profile for a given pathway - but how can we know which pathway to look at??
# # Compute the network centrality scores
# cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")

# # Outgoing heatmap
# ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing") # "sending" cells

# #Incoming heatmap
# ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming") # "recieving" cells

# # show both heatmaps next to each other (the darker the color is - the stringer the signal is)

# pdf("cell_cell_communication/CellChat/output/signaling_heatmaps.pdf", 14, 8)

# ComplexHeatmap::draw(ht1 + ht2)

# dev.off()



# # plot each pathway separately AND remove unused cell identity levels

# cellChat@idents <- droplevels(cellChat@idents)

# # rebuild network after dropping empty groups
# cellChat <- aggregateNet(cellChat)

# # Select signaling pathways to visualize
# # "MIF", "GALECTIN", "APP", "ANNEXIN"
# pathways.show <- c("MIF")

# # Set plotting layout to a single panel
# par(mfrow = c(1, 1))

# # Plot aggregated communication network for selected pathways
# # Circle layout shows overall interaction strength between cell groups
# netVisual_aggregate(
#   cellChat,
#   signaling = pathways.show,
#   layout = "circle"
# )

# # Reset plotting layout to a single panel
# par(mfrow = c(1, 1))

# # Plot heatmap of communication probabilities/intensities
# # Rows/columns typically represent sender and receiver cell groups
# # Darker red = stronger signaling interactions
# netVisual_heatmap(
#   cellChat,
#   signaling = pathways.show,
#   color.heatmap = "Reds"
# )

# # Show contribution of individual ligand-receptor pairs
# # within each selected signaling pathway
# # Helps identify which interactions drive pathway activity
# netAnalysis_contribution(
#   cellChat,
#   signaling = pathways.show
# )

# # Display cell identity labels / cluster names
# # Used to determine which numeric indices correspond to which cell groups
# levels(cellChat@idents)

# # Bubble plot of interactions from source cluster 4
# # to target clusters 5 and 7 across all pathways
# # Bubble size/color usually reflect communication probability/significance
# # remove.isolate = FALSE keeps groups with no interactions visible
# netVisual_bubble(
#   cellChat,
#   sources.use = 4,
#   targets.use = c(5,7),
#   remove.isolate = FALSE
# )

# # Bubble plot restricted to CCL and CXCL signaling only
# # Useful for chemokine-specific communication patterns
# netVisual_bubble(
#   cellChat,
#   sources.use = 4,
#   targets.use = c(5,7),
#   signaling = c("MIF"),
#   remove.isolate = FALSE
# )