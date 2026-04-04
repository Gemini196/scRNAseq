import os
import tempfile
import subprocess


import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
import torch
from scvi.external import CellAssign
import numpy as np

# Tutorial: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/cellassign_tutorial.html#create-and-fit-cellassign-model


def main():
    scvi.settings.seed = 0
    print("Last run with scvi-tools version:", scvi.__version__)
    
    sc.set_figure_params(figsize=(6, 6), frameon=False)
    sns.set_theme()
    torch.set_float32_matmul_precision("high")
    save_dir = tempfile.TemporaryDirectory()

    sce_follicular_path = os.path.join(save_dir.name, "sce_follicular.h5ad")
    sce_hgsc_path = os.path.join(save_dir.name, "sce_hgsc.h5ad")
    fl_celltype_path = os.path.join(save_dir.name, "fl_celltype.csv")
    hgsc_celltype_path = os.path.join(save_dir.name, "hgsc_celltype.csv")

    # Download the data
    downloads = [
        ("https://ndownloader.figshare.com/files/27458798", sce_follicular_path),
        ("https://ndownloader.figshare.com/files/27458822", sce_hgsc_path),
        ("https://ndownloader.figshare.com/files/27458828", hgsc_celltype_path),
        ("https://ndownloader.figshare.com/files/27458831", fl_celltype_path),
    ]

    for url, path in downloads:
        subprocess.run(["wget", "-q", url, "-O", path], check=True)

    # Follicular Lymphoma Data
    follicular_adata = sc.read(sce_follicular_path)
    fl_celltype_markers = pd.read_csv(fl_celltype_path, index_col=0)

    follicular_adata.obs.index = follicular_adata.obs.index.astype("str")
    follicular_adata.var.index = follicular_adata.var.index.astype("str")
    follicular_adata.var_names_make_unique()
    follicular_adata.obs_names_make_unique()

    follicular_adata

    follicular_bdata = follicular_adata[:, fl_celltype_markers.index].copy()
    lib_size = follicular_bdata.X.sum(1)
    follicular_bdata.obs["size_factor"] = lib_size / np.mean(lib_size)
    scvi.external.CellAssign.setup_anndata(follicular_bdata, size_factor_key="size_factor")
    follicular_model = CellAssign(follicular_bdata, fl_celltype_markers)
    follicular_model.train()

    follicular_model.history["elbo_validation"].plot()

    # Predict and plot assigned cell types
    predictions = follicular_model.predict()
    predictions.head()

    sns.clustermap(predictions, cmap="viridis")

    follicular_bdata.obs["cellassign_predictions"] = predictions.idxmax(axis=1).values

    # celltype is the original CellAssign prediction
    sc.pl.umap(
        follicular_bdata,
        color=["celltype", "cellassign_predictions"],
        frameon=False,
        ncols=1,
    )
    
    # Model reproducibility
    df = follicular_bdata.obs
    confusion_matrix = pd.crosstab(
        df["cellassign_predictions"],
        df["celltype"],
        rownames=["cellassign_predictions"],
        colnames=["Original predictions"],
    )
    confusion_matrix /= confusion_matrix.sum(1).ravel().reshape(-1, 1)
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.heatmap(
        confusion_matrix,
        cmap=sns.diverging_palette(245, 320, s=60, as_cmap=True),
        ax=ax,
        square=True,
        cbar_kws={"shrink": 0.4, "aspect": 12},
    )



if __name__ == "__main__":
    main()
