# Single Cell RNA-seq Analysis with Seurat - PBMC 3K Tutorial


# Clean up objects and memory
rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors=FALSE, dplyr.summarise.inform = FALSE)


# import libraries
library(dplyr)
library(Seurat)
library(patchwork)


# Load the PBMC dataset (from Single Cell RNA sequencing)
pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized) data.
# Chooses features (genes) that are expressed in at least 200 cells, and cells that have at least 3 features expressed.
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Quantity control (QC) - remove cells with low/high counts or high mitochondrial content
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# n_count = the total number of molecules detected within a cell
# n_feature = the total number of genes detected within a cell

# Visualize QC metrics as a violin plot (visual representation to identify outliers)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# filter out cells that have unique feature counts over 2,500 or less than 200, or >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
# Scale the data (all genes) - remove unwanted sources of variation, e.g. batch effects, cell cycle, etc.
# ==standardizing data to a mean of zero and a standard deviation of one (z-score normalization)
pbmc <- ScaleData(pbmc, features = all.genes)

# Dimentionality reduction - PCA
# For each principal component, every gene has a loading (a weight within the linear transformation) that indicates how much that gene contributes to the PC
# PC1 explains more of the variance in the data than PC2, which explains more than PC3, etc.
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# visualize the first few principal components. Each point on the plot is a cell.
# here we're representing PC1 against PC2 - we can see that cells group in certain clusters according to these PCs
DimPlot(pbmc, reduction = "pca") + NoLegend()

# The heatmap below shows the loadings of the genes associated with PC1 across a random subset of 500 cells.
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset = How many PCs do we need to explain the variation in the data? 
# we tried to find the "elbow"=  the point in which the variation explained by each additional PC isn't as big anymore
ElbowPlot(pbmc)

# find the clusters by computing a k-nearest neighbor graph based on the euclidean distance in PCA 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# lower resolution parameter leads to fewer clusters, higher resolution leads to more clusters
# now in the metadata we'll find a number indicating to which cluster each of the cells belongs
pbmc <- FindClusters(pbmc, resolution = c(0.1, 0.5, 1))
# Show pca plot with clusters (+colors)
DimPlot(pbmc, group.by = "RNA_snn_res.0.5", label = TRUE)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# dimensional reduction technique (non linear. another option is tSNE) to visualize clusters in 2D space
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

saveRDS(pbmc, file = "seurat_tutorial_R/output/pbmc_tutorial.rds")


#----- DONE WITH CLUSTERING, NOW FIND GENETIC MARKERS -----#
# Giving identity to clusters according to the expression of canonical markers (genes that are known to be expressed in certain cell types)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "seurat_tutorial_R/output/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(pbmc, file = "seurat_tutorial_R/output/pbmc3k_final.rds")