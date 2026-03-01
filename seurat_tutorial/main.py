"""Seurat tutorial: PBMC 3K dataset analysis.

This module demonstrates using Python to interact with Seurat via rpy2 wrappers.
"""

import sys
from pathlib import Path

# Add parent directory to path to allow imports from sibling modules
sys.path.insert(0, str(Path(__file__).parent))

from wrappers import read_10x, create_seurat_object
from wrappers import percentage_feature_set, add_metadata_column, vln_plot, save_plot


def main():
    """Entry point for the scRNAseq application.

    This loads the PBMC 3K dataset and creates a Seurat object with basic filtering.
    """
    # Construct the path to the data directory relative to this script
    script_dir = Path(__file__).parent
    data_dir = script_dir.parent / "filtered_gene_bc_matrices" / "hg19"
    data_dir_str = str(data_dir)
    
    print("Loading PBMC 10X data...")
    print(f"Data directory: {data_dir_str}")
    pbmc_data = read_10x(data_dir_str)
    
    # Initialize the Seurat object with the raw (non-normalized) data
    print("Creating Seurat object...")
    pbmc = create_seurat_object(
        pbmc_data,
        project="pbmc3k",
        min_cells=3,
        min_features=200
    )
    
    # Display the Seurat object summary
    print("\nSeurat object summary:")
    print(pbmc)

    # The [[ operator can add columns to object metadata.
    # This is a great place to stash QC stats
    print("\nCalculating QC metrics (percent.mt)...")
    percent_mt = percentage_feature_set(pbmc, pattern="^MT-")
    pbmc = add_metadata_column(pbmc, "percent.mt", percent_mt)

    print("Seurat object after QC metrics:")
    print(pbmc)

    # Visualize QC metrics as a violin plot and save to file
    print("\nVisualizing QC metrics (violin plots)...")
    plot_obj = vln_plot(pbmc, features=["nFeature_RNA", "nCount_RNA", "percent.mt"], ncol=3)
    out_dir = script_dir.parent / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "qc_vln.png"
    save_plot(plot_obj, str(out_file), width=10, height=4, dpi=300)
    print(f"Saved violin plot to: {out_file}")


if __name__ == "__main__":
    main()

