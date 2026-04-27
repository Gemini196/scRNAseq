# Lignature is a computational tool and curated database that predicts active ligands by matching transcriptomic changes in receiver cells against a library of experimentally derived signatures for 362 human ligands.
# Unlike network-based methods (such as NicheNet) that rely on incomplete or "noisy" signaling models, Lignature utilizes actual gene expression "fingerprints" to capture the direction and magnitude of cellular responses
# This data-driven approach allows it to identify biologically relevant communication events with higher accuracy and fewer false positives in both in vitro and single-cell datasets
# Would you like me to create a comparison report that highlights the specific performance rankings of Lignature and NicheNet across different datasets?

install.packages(c("tidyverse", "Seurat", "RColorBrewer", "pheatmap", "lsa", "data.table", "circlize", "randomcoloR", "hrbrthemes", "ggrepel", "caret", "glmnet"))
bio_pkgs = c("fgsea", "graphite")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(bio_pkgs)

pak::pak("yingxinac/Lignature")

############################
# Load packages
############################

library(Seurat)
library(RColorBrewer)
library(pheatmap)
library(fgsea)
library(graphite)
library(lsa)
library(devtools)
library(data.table)
library(circlize)
library(tidyverse)
library(randomcoloR)
library(hrbrthemes)
library(ggrepel)
library(caret)
library(glmnet)
library(Lignature)

`%!in%` <- Negate(`%in%`)
`%notin%` <- Negate(`%in%`)

############################
# Load Lignature data
############################

load("siglist.RData")
load("sigmeta.RData")
load("lr_network_hg.RData")

# lr_network columns:
# L      = ligand name
# Lgene  = ligand gene symbol(s)
# R      = receptor name
# Rgene  = receptor gene symbol(s)
# Multi-subunit genes are joined with "_"

lr_network = unique(as.matrix(lr_network)[, c("L", "Lgene", "R", "Rgene")])

load("mykegg_10500.RData")

# mykegg = list of gene nodes in KEGG pathways

############################
# Group signatures by treatment ligands
############################

siggroupsoutput = groupSigs(siglist = siglist, sigmeta = sigmeta)
sigLs = siggroupsoutput$sigLs
siggroups = siggroupsoutput$siggroups

############################
# Load input data
############################

load("seuratobj.RData")

# Input Seurat object should be LogNormalized
# Metadata must contain:
# cond = sample condition
# cl   = cell type / cluster

############################
# Get expression summary and DE information
############################

conds = c("C", "NC")

condpairs = matrix(c("C", "NC"), ncol = 2)
colnames(condpairs) = c("cond1", "cond2")

# Each row of condpairs defines cond1 vs cond2 comparison

allcls = unique(seuratobj@meta.data$cl)

exprinfo = getExprInfo(
  seuratobj = seuratobj,
  cls = allcls,
  conds = conds,
  condpairs = condpairs
)

############################
# Get signatures of input data
############################

mysiglist = getMySigs(
  DElist = exprinfo$DElist,
  cls = allcls,
  condpairs = condpairs,
  keggnodes = mykegg,
  pctcut = 0.1,
  belowpctcut = "lowpct.remove"
)

# Genes detected in <10% cells are removed

############################
# Get LR interaction scores
############################

LRIscorelist = getLRIScore(
  lr_network = lr_network,
  exprinfo = exprinfo,
  LRI.method = "scsigr"
)

# LRI.method options:
# "cpdb", "natmi", "scsigr"

############################
# Get ligand scores for receiver cells
############################

cl_to = "Myeloid"
condpair = "C_vs_NC"

sigscorelist = getSigScore(
  siglist = siglist,
  testsig = mysiglist[[cl_to]][[condpair]],
  whichsig = "lfc",
  method = "cor.pearson",
  nump = 1000,
  corpby.permu = FALSE
)

Lscoredf = getLScore(
  siggroups = siggroups,
  lr_network = lr_network,
  sigscorelist = sigscorelist,
  pvalcut = 0.1
)

# Only signatures with p < 0.1 are used

############################
# Optional confidence scores
############################

whichLscore = "max"

Lscoredf$confscore = getconfscore(
  myscores = Lscoredf[, whichLscore],
  siglist = siglist,
  siggroups = siggroups,
  whichsig = "lfc",
  method = "cor.pearson"
)

############################
# Optional combinatorial ligand analysis
############################

Lcomb = getLcomb(
  siglist = siglist,
  testsig = mysiglist[[cl_to]][[condpair]],
  siggroups = siggroups,
  sigLs = sigLs,
  whichsig = "lfc",
  combmethod = "lasso"
)

# combmethod:
# "lasso", "ridge", "elasticnet"

Lscoredf$Lweight = Lcomb$Lscores[rownames(Lscoredf)]

############################
# Organize LR information
############################

LRinfolist = getLRinfosummary(
  exprinfo = exprinfo,
  LRIscorelist = LRIscorelist,
  Lscoredf = Lscoredf,
  lr_network = lr_network,
  cl_to = cl_to,
  condpair = condpair,
  otherLscore = c("confscore", "Lweight")
)

############################
# Scatter plot of ligand vs LR scores
############################

getScatterLRpairs(
  LRinfolist = LRinfolist,
  cls_from = allcls,
  cl_to = cl_to,
  Lscore.by = "Lscore_signed_maxabs",
  SLRI.by = "SLRI_condmax",
  point.size = 3,
  Lscore.cut = 0.1,
  SLRI.cut = 0.3
)

############################
# Filter LR interactions
############################

cut.by = c(
  "Lscore_max_abs",
  "SLRI_condmax",
  "SLRIspec_condmax",
  "Lpct_condmax",
  "Rpct_condmax"
)

cut.bounds = list(
  c(0.1, Inf),
  c(0.3, Inf),
  c(0.75, Inf),
  c(0.25, Inf),
  c(0.25, Inf)
)

names(cut.bounds) = cut.by

LRinfo_fil = filterLRinfo(
  LRinfolist = LRinfolist,
  cut.bounds = cut.bounds
)

############################
# Write tables
############################

mydir = "outputtables/"

writeLRtable(
  tablelist = LRinfo_fil$LRinfolist_fil,
  mydir = mydir
)

############################
# Plot DE ligand/receptor heatmaps
############################

whichLs = "in.Lignature.only"

col = rev(colorRampPalette(brewer.pal(3, "RdBu"))(1000))

plotDELRheatmap(
  exprinfo = exprinfo,
  siggroups = siggroups,
  lr_network = lr_network,
  LorR = "L",
  whichLs = whichLs,
  cls = allcls,
  condpair = condpair,
  pct.cut = 0.25,
  padj.cut = 0.01,
  lfc.cut = 0.2,
  col = col
)

plotDELRheatmap(
  exprinfo = exprinfo,
  siggroups = siggroups,
  lr_network = lr_network,
  LorR = "R",
  whichLs = whichLs,
  cls = allcls,
  condpair = condpair,
  pct.cut = 0.25,
  padj.cut = 0.01,
  lfc.cut = 0.2,
  col = col
)

############################
# KEGG NES barplot
############################

getKEGGNESbarplot(
  mysig = mysiglist[[cl_to]][[condpair]],
  padj.cut = 0.1
)

############################
# Scatter plot of ligand scores vs expression
############################

getScatter(
  LRinfolist = LRinfolist,
  LRs = LRinfo_fil$lrs_fil_union,
  cls_from = allcls,
  cl_to = cl_to,
  condpair = condpair,
  Lscore.by = "Lscore_signed_maxabs",
  point.size = 3
)

############################
# Plot LR heatmaps
############################

col = rev(colorRampPalette(brewer.pal(3, "RdBu"))(1000))

plotLRs(
  LRinfolist = LRinfolist,
  LRs = LRinfo_fil$lrs_fil_union,
  cls_from = allcls,
  cl_to = cl_to,
  condpair = condpair,
  plotwhich = "Lavgexpr",
  Lscore.by = "Lscore_signed_maxabs",
  order.by = "Lscore_signed_maxabs",
  order.decreasing = TRUE,
  col = col,
  scale = "row",
  show_rownames = TRUE
)

col = rev(colorRampPalette(brewer.pal(3, "PuOr"))(1000))

plotLRs(
  LRinfolist = LRinfolist,
  LRs = LRinfo_fil$lrs_fil_union,
  cls_from = allcls,
  cl_to = cl_to,
  condpair = condpair,
  plotwhich = "Ravgexpr",
  Lscore.by = "Lscore_signed_maxabs",
  order.by = "Lscore_signed_maxabs",
  order.decreasing = TRUE,
  col = col,
  scale = "none",
  show_rownames = TRUE
)

col = colorRampPalette(brewer.pal(3, "PiYG"))(1000)

plotLRs(
  LRinfolist = LRinfolist,
  LRs = LRinfo_fil$lrs_fil_union,
  cls_from = allcls,
  cl_to = cl_to,
  condpair = condpair,
  plotwhich = "LScore",
  Lscore.by = "Lscore_signed_maxabs",
  order.by = "Lscore_signed_maxabs",
  order.decreasing = TRUE,
  col = col,
  scale = "none",
  show_rownames = TRUE
)

############################
# Circos plot
############################

cut.by = c(
  "Lscore_max",
  "SLRI_condmax",
  "SLRIspec_condmax",
  "Lpct_condmax",
  "Rpct_condmax"
)

cut.bounds = list(
  c(0.1, Inf),
  c(0.3, Inf),
  c(0.9, Inf),
  c(0.25, Inf),
  c(0.25, Inf)
)

names(cut.bounds) = cut.by

LRinfo_fil2 = filterLRinfo(
  LRinfolist = LRinfolist,
  cut.bounds = cut.bounds
)

mycols = c(
  "paleturquoise3", "#D098EE", "lightgoldenrod1",
  "#2CA030", "thistle2", "salmon",
  "#3CB7CC", "#FFB900", "#AC613C",
  "#84BD00", "#FF5A5A", "lightsteelblue",
  "#F9D23C", "lightpink"
)

getCircosPlot(
  plotLRinfolist = LRinfo_fil2$LRinfolist_fil,
  plotcond = "condmax",
  cls_from = allcls,
  cl_to = cl_to,
  mycols = mycols,
  width_same_ligand_cluster = 0.2,
  width_different_ligand_cluster = 0.5,
  width_same_receptor_cluster = 0.2,
  width_different_receptor_cluster = 0.5,
  width_ligand_receptor = 1,
  cplotthresh = 0.1,
  cex = 1,
  transupperbase = 1.5,
  transidx = 0.5
)