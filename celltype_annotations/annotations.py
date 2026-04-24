from __future__ import annotations
import anndata as ad
import pooch
import scanpy as sc
import scvi
import torch
import pandas as pd
import numpy as np
from scvi.external import CellAssign
import celltypist as ct
import decoupler as dc
import celltypeai as cta

sc.set_figure_params(dpi=50, facecolor="white")

# Tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html#

def main():
    # =========================
    # Download example dataset
    # =========================
    # Location:
    #    Linux	~/.cache/scverse_tutorials/
    #    macOS	~/Library/Caches/scverse_tutorials/
    #    Windows	C:\Users\<user>\AppData\Local\scverse_tutorials\Cache\
    EXAMPLE_DATA = pooch.create(
        path=pooch.os_cache("scverse_tutorials"),
        base_url="doi:10.6084/m9.figshare.22716739.v1/",
    )
    EXAMPLE_DATA.load_registry_from_doi()


    # =========================
    # Load 10x datasets
    # =========================
    samples = {
        "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
        "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
    }
    adatas = {}

    for sample_id, filename in samples.items():
        path = EXAMPLE_DATA.fetch(filename)
        sample_adata = sc.read_10x_h5(path)
        sample_adata.var_names_make_unique()
        adatas[sample_id] = sample_adata


    # =========================
    # Merge samples into one AnnData - combines multiple single-cell datasets into one dataset and ensures that every cell has a unique identifier.
    # =========================
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()


    # =========================
    # Annotate gene categories for QC (quality control)
    # =========================
    adata.var["mt"] = adata.var_names.str.startswith("MT-")          # mitochondrial genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL")) # ribosomal genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")       # hemoglobin genes


    # =========================
    # Calculate quality control metrics
    # =========================
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        log1p=True,
    )

    # Inspect violin plots of some of the computed QC metrics (visualize the distribution of QC metrics across cells))
    sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,)

    # Inspect scatter plot of total counts vs. percentage of mitochondrial counts
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")


    # Filter out low-quality cells and genes (remove cells with fewer than 100 genes detected and genes that are expressed in fewer than 3 cells)
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)

    # Run a doublet detection algorithm
    sc.pp.scrublet(adata, batch_key="sample")
    # One can now either filter directly on predicted_doublet or use the doublet_score later during clustering to filter clusters with high doublet scores.

    # =========================
    # Store raw counts
    # =========================
    adata.layers["counts"] = adata.X.copy()


    # =========================
    # Normalization and log transform
    # =========================
    sc.pp.normalize_total(adata)  # normalize counts per cell
    sc.pp.log1p(adata)            # log-transform. log(1+x), so that zero counts remain zero after transformation


    # =========================
    # Identify highly variable genes
    # =========================
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")


    # =========================
    # Dimensionality reduction (PCA)
    # =========================
    sc.tl.pca(adata)

    sc.pl.pca_variance_ratio(
        adata,
        n_pcs=50,
        log=True,
    )

    sc.pl.pca(
        adata,
        color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
    )



    # =========================
    # UMAP visualization
    # =========================
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
    )

    # =========================
    # Clustering (Leiden)
    # =========================
    # Using the igraph implementation and a fixed number of iterations can be significantly faster,
    # especially for larger datasets
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)


    for res in [0.02, 0.5, 2.0]:
        sc.tl.leiden(
            adata,
            key_added=f"leiden_res_{res:4.2f}",
            resolution=res,
            flavor="igraph",
        )

    sc.pl.umap(adata, color=["leiden"])

    # =========================
    # Re-assess quality control and cell filtering
    # =========================
    sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,)

    sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,)

    # =========================
    # Manual cell-type annotation
    # =========================
    for res in [0.02, 0.5, 2.0]:
        sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph")

    sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",)

    sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",)



    # =========================
    # Known marker genes for annotation
    # =========================
    marker_genes = {
        "CD14+ Mono": ["FCN1", "CD14"],
        "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
        "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
        "Erythroblast": ["MKI67", "HBA1", "HBB"],
        "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
        "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
        "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
        "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
        "B cells": [
            "MS4A1",
            "ITGB1",
            "COL4A4",
            "PRDM1",
            "IRF4",
            "PAX5",
            "BCL11A",
            "BLK",
            "IGHD",
            "IGHM",
        ],
        "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
        "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
        "CD4+ T": ["CD4", "IL7R", "TRBC2"],
        "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
        "T naive": ["LEF1", "CCR7", "TCF7"],
        "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
    }

    # sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")

    adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Lymphocytes",
        "1": "Monocytes",
        "2": "Erythroid",
        "3": "B Cells",
    })

    #sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")

    # =========================
    # Differentially-expressed Genes as Markers
    # =========================

    # Obtain cluster-specific differentially expressed genes
    sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

    # sc.pl.rank_genes_groups_dotplot(adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5)

    sc.get.rank_genes_groups_df(adata, group="7").head(5)

    dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
    sc.pl.umap(
        adata,
        color=[*dc_cluster_genes, "leiden_res_0.50"],
        legend_loc="on data",
        frameon=False,
        ncols=3,
    )

    # =========================
    # Automatic label prediction with CellTypeAI
    # =========================
    cta.cell_annotator(
                "Human",
                "Immune system",
                adata,
                model="phi4:14b",
                num_genes=200,
                n_iterations=1, #recommended:3
            )

    sc.pl.umap(
        adata,
        color="cell_type_ai",
        ncols=1,
        size=3,
        )

    # =========================
    # Automatic label prediction with CellTypist
    # =========================
    ct.models.download_models(model=["Immune_All_Low.pkl"], force_update=True)

    model = ct.models.Model.load(model="Immune_All_Low.pkl")
    predictions = ct.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True, over_clustering="leiden_res_0.50")
    # convert back to anndata||
    adata = predictions.to_adata()

    sc.pl.umap(adata, color="majority_voting", ncols=1)

    # =========================
    # CellAssign Annotation
    # =========================
    scvi.settings.seed = 0
    torch.set_float32_matmul_precision("high")

    # -------------------------
    # Prepare marker gene matrix
    # -------------------------
    all_genes = sorted(set(g for genes in marker_genes.values() for g in genes))
    marker_df = pd.DataFrame(0, index=all_genes, columns=marker_genes.keys())

    for celltype, genes in marker_genes.items():
        marker_df.loc[genes, celltype] = 1

    # Keep only genes
    marker_df = marker_df.loc[marker_df.index.intersection(adata.var_names)]

    # -------------------------
    # Prepare AnnData for CellAssign
    # -------------------------
    adata_ca = adata[:, marker_df.index].copy()

    # Use RAW counts
    adata_ca.X = adata_ca.layers["counts"]

    # Compute size factors
    lib_size = np.array(adata_ca.X.sum(axis=1)).flatten()
    adata_ca.obs["size_factor"] = lib_size / np.mean(lib_size)

    # Setup scvi
    CellAssign.setup_anndata(adata_ca, size_factor_key="size_factor")

    # -------------------------
    # Train model
    # -------------------------
    model = CellAssign(adata_ca, marker_df)
    model.train()

    # -------------------------
    # Predict cell types
    # -------------------------
    predictions = model.predict()

    adata.obs["cellassign_labels"] = predictions.idxmax(axis=1).values

    
    # -------------------------
    # Annotation with enrichment analysis
    # -------------------------
    # Query Omnipath and get PanglaoDB
    markers = dc.op.resource(name="PanglaoDB", organism="human")
    # Keep canonical cell type markers alone
    markers = markers[markers["canonical_marker"]]

    # Remove duplicated entries
    markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]
    markers.head()

    dc.mt.mlm(adata, net=markers.rename(columns=dict(cell_type="source", genesymbol="target")), verbose=True)

    acts = dc.pp.get_obsm(adata, "score_mlm")
    sc.pl.umap(
        acts,
        color=[
            "majority_voting",
            "B cells",
            "T cells",
            "Monocytes",
            "Erythroid-like and erythroid precursor cells",
            "NK cells",
        ],
        wspace=0.5,
        ncols=3,
    )

    enr = dc.tl.rankby_group(acts, groupby="leiden_res_0.50")
    annotation_dict = enr[enr["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()
    adata.obs["dc_anno"] = [annotation_dict[clust] for clust in adata.obs["leiden_res_0.50"]]
    

    # -------------------------
    # Visualization (FINAL)
    # -------------------------
    sc.pl.umap(
        adata,
        color=["majority_voting", "cellassign_labels", "dc_anno", "cell_type_ai"],
        wspace=0.4,   # space between plots
        ncols=2,      # ensures side-by-side layout
        size=3,
    )


if __name__ == "__main__":
    main()
