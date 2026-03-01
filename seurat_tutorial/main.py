"""Seurat tutorial: PBMC 3K dataset analysis.

This module demonstrates using Python to interact with Seurat via rpy2 wrappers.
"""

import sys
from pathlib import Path

# Add parent directory to path to allow imports from sibling modules
sys.path.insert(0, str(Path(__file__).parent))

from wrappers import read_10x, create_seurat_object


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


if __name__ == "__main__":
    main()

