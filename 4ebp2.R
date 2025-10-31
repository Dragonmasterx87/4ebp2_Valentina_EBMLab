# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 03/23/2022
# R version 4.2.1 (2019-12-12) 'Funny-Looking Kid'

#!/usr/bin/env Rscript
# 01_install_packages.R - Install ALL required packages
# Run ONCE after fresh R installation with Rtools45

cat("Installing packages...\n\n")

# CRAN packages
cran_pkgs <- c(
  "Seurat", "SeuratObject", "harmony", "Signac", "qs", "arrow",
  "ggplot2", "ggrepel", "cowplot", "patchwork", "viridis", "viridisLite",
  "RColorBrewer", "circlize", "dplyr", "tidyr", "readr", "tidyverse",
  "data.table", "stringr", "Matrix", "ggridges", "plotly", "clustree",
  "future", "hdf5r", "R.utils", "devtools", "ragg", "reticulate",
  "SoupX", "leiden", "gprofiler2", "Rcpp"
)

install.packages(cran_pkgs, repos = "https://cran.rstudio.com", Ncpus = 4)

# Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)

bioc_pkgs <- c(
  "scDblFinder", "EnsDb.Hsapiens.v86", "GenomeInfoDb", "glmGamPoi",
  "GOSemSim", "org.Hs.eg.db", "AnnotationHub", "MeSHDbi",
  "clusterProfiler", "DOSE", "dittoSeq", "escape", "EnrichmentBrowser",
  "ComplexHeatmap", "Nebulosa", "DropletUtils", "EnhancedVolcano",
  "JASPAR2020", "TFBSTools", "BSgenome.Hsapiens.UCSC.hg38",
  "motifmatchr", "chromVAR", "BiocParallel", "DESeq2", "limma",
  "edgeR", "SingleCellExperiment", "SpatialExperiment"
)

BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

# GitHub packages
if (!require("remotes")) install.packages("remotes")

remotes::install_github('cole-trapnell-lab/monocle3', upgrade = "never")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade = "never")
remotes::install_github('yanlinlin82/ggvenn', upgrade = "never")
remotes::install_github('gaospecial/ggVennDiagram', upgrade = "never")
remotes::install_github('quadbiolab/Pando', upgrade = "never")
remotes::install_github('cole-trapnell-lab/cicero-release', ref = "monocle3", upgrade = "never")
remotes::install_github('satijalab/seurat-wrappers', upgrade = "never")

cat("\n✓ Installation complete\n")

# 02_load_libraries.R
# Source this at the start of every analysis: source("02_load_libraries.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(leiden)
  library(stringr)
  library(hdf5r)
  library(SoupX)
  library(Rcpp)
  library(cowplot)
  library(Matrix)
  library(ggridges)
  library(dplyr)
  library(tidyverse)
  library(data.table)
  library(reticulate)
  library(Seurat)
  library(monocle3)
  library(harmony)
  library(Signac)
  library(scDblFinder)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  library(plotly)
  library(clustree)
  library(patchwork)
  library(future)
  library(DoubletFinder)
  library(EnhancedVolcano)
  library(glmGamPoi)
  library(GOSemSim)
  library(org.Hs.eg.db)
  library(AnnotationHub)
  library(MeSHDbi)
  library(clusterProfiler)
  library(DOSE)
  library(dittoSeq)
  library(escape)
  library(EnrichmentBrowser)
  library(viridisLite)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(Nebulosa)
  library(DropletUtils)
  library(ggvenn)
  library(ggVennDiagram)
  library(devtools)
  library(R.utils)
  library(qs)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(motifmatchr)
  library(chromVAR)
  library(cicero)
  library(BiocParallel)
  library(DESeq2)
  library(Pando)
  library(ragg)
  library(gprofiler2)
  library(arrow)
})

options(future.globals.maxSize = Inf)
set.seed(1234)

cat("✓ All libraries loaded\n")


