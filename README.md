# scRNAseq
scRNASeq tool usage

## Setup Instructions 📦

This repository primarily contains Python code, but some workflows leverage R. To get the environment running correctly:

1. **Python dependencies** are listed in `requirements.txt` and can be installed with:
   ```sh
   pip install -r requirements.txt
   ```

2. **R dependencies** are managed separately. At a minimum, install:
   - **Seurat** – the main package used for single-cell RNA-seq analysis
   - **languageserver** – for editor support (R diagnostics/autocomplete)

   Run these commands in an R session (or RStudio):

   ```r
   install.packages("Seurat")
   install.packages("languageserver")
   library(Seurat)
   library(languageserver)
   ```

   Or as a one-liner from the shell:

   ```sh
   R -e 'install.packages(c("Seurat", "languageserver"))'
   ```

   > You can also add R package requirements to a `DESCRIPTION` file or use `renv`/`packrat` for
   > reproducible R environments.

3. Document any additional setup in this README as your project grows.
