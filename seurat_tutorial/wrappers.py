"""Lightweight wrappers around a few Seurat functions using rpy2.

This module lives outside the normal package structure because it is a
standalone tutorial demonstrating how Python code can call into R/Seurat.
"""

from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

# enable automatic conversion between pandas DataFrames and R data.frames
# Note: pandas2ri.activate() is deprecated in newer rpy2 versions
# Use pandas2ri.converter for context-based conversion

# Load Seurat library in R and access functions via robjects
robjects.r("library(Seurat)")


def read_10x(data_dir):
    """Wrap Seurat::Read10X.

    Parameters
    ----------
    data_dir : str
        Path to the directory containing the 10X Genomics output files
        (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz).

    Returns
    -------
    R object
        A sparse matrix (dgCMatrix) with genes as rows and cells as columns.
    """
    return robjects.r.Read10X(data_dir=data_dir)


def create_seurat_object(counts_df, project="PythonSeurat", min_cells=0, min_features=0):
    """Wrap Seurat::CreateSeuratObject.

    Parameters
    ----------
    counts_df : pandas.DataFrame or R object
        Genes x cells counts matrix. Can be a pandas DataFrame or an R matrix/dgCMatrix.
    project : str
        Project name passed to the Seurat object.
    min_cells : int
        Include features detected in at least this many cells. Defaults to 0.
    min_features : int
        Include cells where at least this many features are detected. Defaults to 0.

    Returns
    -------
    R object
        The Seurat object returned by the R function. Python code can
        continue to manipulate it via rpy2.
    """
    # convert pandas to R if needed; otherwise assume it's already an R object
    if hasattr(counts_df, '__class__') and 'DataFrame' in str(counts_df.__class__):
        with pandas2ri.converter.context():
            r_counts = pandas2ri.py2rpy(counts_df)
    else:
        r_counts = counts_df
    
    return robjects.r.CreateSeuratObject(**{
        'counts': r_counts,
        'project': project,
        'min.cells': min_cells,
        'min.features': min_features,
    })


def normalize_data(seurat_obj, normalization_method="LogNormalize", scale_factor=10000):
    """Wrap Seurat::NormalizeData.

    Parameters are the same as the R function.

    Returns
    -------
    R object
        The modified Seurat object (returned invisibly by the R function).
    """
    # rpy2 does not allow Python identifiers with dots, so build a dict
    return robjects.r.NormalizeData(seurat_obj, **{
        'normalization.method': normalization_method,
        'scale.factor': scale_factor,
    })


def find_variable_features(seurat_obj, selection_method="vst", nfeatures=2000):
    """Wrap Seurat::FindVariableFeatures.

    Parameters mirror the R function documentation.

    Returns
    -------
    R object
        The Seurat object with variable features identified (returned invisibly).
    """
    return robjects.r.FindVariableFeatures(seurat_obj, **{
        'selection.method': selection_method,
        'nfeatures': nfeatures,
    })


def scale_data(seurat_obj, features=None, vars_to_regress=None, split_by=None, **kwargs):
    """Wrap Seurat::ScaleData.

    Any additional keyword arguments are forwarded to the R function.
    """
    args = {}
    if features is not None:
        args["features"] = features
    if vars_to_regress is not None:
        args["vars.to.regress"] = vars_to_regress
    if split_by is not None:
        args["split.by"] = split_by
    args.update(kwargs)
    return robjects.r.ScaleData(seurat_obj, **args)


def run_pca(seurat_obj, features=None, npcs=50, verbose=False, **kwargs):
    """Wrap Seurat::RunPCA.

    Common parameters are provided; additional args may be passed.
    """
    args = {"npcs": npcs, "verbose": verbose}
    if features is not None:
        args["features"] = features
    args.update(kwargs)
    return robjects.r.RunPCA(seurat_obj, **args)


def run_tsne(seurat_obj, dims=1, perplexity=30, **kwargs):
    """Wrap Seurat::RunTSNE.

    """
    return robjects.r.RunTSNE(seurat_obj, dims=dims, perplexity=perplexity, **kwargs)


def run_umap(seurat_obj, dims=1, **kwargs):
    """Wrap Seurat::RunUMAP.

    """
    return robjects.r.RunUMAP(seurat_obj, dims=dims, **kwargs)


def find_neighbors(seurat_obj, dims=1, k_param=20, **kwargs):
    """Wrap Seurat::FindNeighbors.

    Parameters
    ----------
    k_param : int
        The `k.param` argument in Seurat (captured here with an underscore).
    """
    # `k.param` uses a dot, so pass via dict expansion
    return robjects.r.FindNeighbors(seurat_obj, **{
        'dims': dims,
        'k.param': k_param,
        **kwargs,
    })


def find_clusters(seurat_obj, resolution=0.8, **kwargs):
    """Wrap Seurat::FindClusters.

    """
    return robjects.r.FindClusters(seurat_obj, resolution=resolution, **kwargs)



