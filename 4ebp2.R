# CODING COMPENDIUM ####
# The following set of code is a description of the analysis performed in the 
# paper entitled "enter name of paper here"
# Author Fahd Qadir FMJ Lab Tulane University, Schoool of Medicine
# Date code was written: 31/10/2025
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

# Explore scRNA-seq directory structure

# ============================================
# COMPLETE WORKFLOW: DecontX + Full Analysis
# Ad Astra Per Aspera - Never Give Up!
# ============================================
library(Seurat)
library(celda)
library(SingleCellExperiment)
library(harmony)
library(ggplot2)
library(dplyr)
library(qs)
library(gprofiler2)
library(RColorBrewer)

# Install celda if needed
if (!requireNamespace("celda", quietly = TRUE)) {
  BiocManager::install("celda")
  library(celda)
}

setwd("C:/Users/mqadir/Box/scRNAseq NOD pancreas immune cells/Analysis")

cat("✓ All libraries loaded\n")
cat("✓ Working directory set to:\n  ", getwd(), "\n")

# ============================================
# PART 1: Define Sample Information
# ============================================
dir1 <- "C:/Users/mqadir/Box/scRNAseq NOD pancreas immune cells/eMizrachi_10X3primev4OCM_07302025"
dir2 <- "C:/Users/mqadir/Box/scRNAseq NOD pancreas immune cells/eMizrachi_10X3primev4OCM_08212025"

samples <- data.frame(
  sample = c("Sample1", "Sample2", "Sample3", "Sample4", 
             "Sample5", "Sample6", "Sample7", "Sample8"),
  genotype = c("BP2_KO", "Control", "BP2_KO", "Control",
               "BP2_KO", "BP2_KO", "Control", "Control"),
  batch = c("Batch1", "Batch1", "Batch1", "Batch1",
            "Batch2", "Batch2", "Batch2", "Batch2"),
  directory = c(rep(dir1, 4), rep(dir2, 4))
)

print(samples)

# ============================================
# PART 2: DecontX Ambient RNA Removal
# ============================================
cat("\n=== PART 2: DecontX Ambient RNA Removal ===\n\n")

decontx_corrected_list <- list()

for(i in 1:nrow(samples)) {
  cat("\n--- Processing", samples$sample[i], "with DecontX ---\n")
  
  sample_path <- file.path(samples$directory[i], samples$sample[i])
  
  # Read filtered matrix
  filt_matrix <- Read10X_h5(file.path(sample_path, "sample_filtered_feature_bc_matrix.h5"))
  
  cat("Filtered matrix dims:", dim(filt_matrix), "\n")
  
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(assays = list(counts = filt_matrix))
  
  # Run decontX
  cat("Running decontX (this may take a few minutes)...\n")
  sce <- decontX(sce, verbose = TRUE)
  
  # Get corrected counts
  corrected_matrix <- decontXcounts(sce)
  
  # Get contamination estimates
  contamination <- colData(sce)$decontX_contamination
  cat("Mean contamination:", round(mean(contamination), 4), "\n")
  cat("Median contamination:", round(median(contamination), 4), "\n")
  cat("Max contamination:", round(max(contamination), 4), "\n")
  
  cat("Corrected matrix dims:", dim(corrected_matrix), "\n")
  
  # Store results
  decontx_corrected_list[[samples$sample[i]]] <- list(
    matrix = corrected_matrix,
    contamination = contamination
  )
}

cat("\n✓ DecontX correction completed for all samples!\n")

# ============================================
# PART 3: Create Seurat Objects
# ============================================
cat("\n=== PART 3: Creating Seurat Objects ===\n\n")

seurat_list <- list()

for(i in 1:nrow(samples)) {
  cat("\nCreating Seurat object for", samples$sample[i], "...\n")
  
  corrected_counts <- decontx_corrected_list[[samples$sample[i]]]$matrix
  contamination <- decontx_corrected_list[[samples$sample[i]]]$contamination
  
  seurat_obj <- CreateSeuratObject(counts = corrected_counts, 
                                   project = samples$sample[i],
                                   min.cells = 3, 
                                   min.features = 200)
  
  seurat_obj$sample <- samples$sample[i]
  seurat_obj$genotype <- samples$genotype[i]
  seurat_obj$batch <- samples$batch[i]
  seurat_obj$decontX_contamination <- contamination
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  seurat_list[[samples$sample[i]]] <- seurat_obj
  
  cat("Cells:", ncol(seurat_obj), "| Features:", nrow(seurat_obj), "\n")
}

# ============================================
# PART 4: Merge and QC
# ============================================
cat("\n=== PART 4: Merging & QC ===\n\n")

seurat <- merge(seurat_list[[1]], y = seurat_list[2:8], 
                add.cell.ids = samples$sample)

cat("Before QC:\n")
print(seurat)
print(table(seurat$sample))
print(table(seurat$genotype))

seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                   percent.mt < 20)

cat("\nAfter QC:\n")
print(seurat)
cat("Mean contamination across all cells:", round(mean(seurat$decontX_contamination), 4), "\n")

# ============================================
# PART 5: Normalization & Harmony
# ============================================
cat("\n=== PART 5: Normalization & Harmony ===\n\n")

seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- JoinLayers(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunHarmony(object = seurat, group.by.vars = "batch")
seurat <- RunUMAP(seurat, reduction = "harmony", dims = 1:30)
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.5)

qsave(seurat, "seurat_decontx_harmony_integrated.qs")

# Basic plots
p1 <- DimPlot(seurat, reduction = "umap", group.by = "sample")
p2 <- DimPlot(seurat, reduction = "umap", group.by = "genotype")
p3 <- DimPlot(seurat, reduction = "umap", group.by = "batch")
p4 <- DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p5 <- FeaturePlot(seurat, features = "decontX_contamination", cols = c("lightgrey", "red"))

ggsave("umap_by_sample_decontx.pdf", p1, width = 10, height = 8)
ggsave("umap_by_genotype_decontx.pdf", p2, width = 10, height = 8)
ggsave("umap_by_batch_decontx.pdf", p3, width = 10, height = 8)
ggsave("umap_clusters_decontx.pdf", p4, width = 10, height = 8)
ggsave("umap_contamination_decontx.pdf", p5, width = 10, height = 8)

print(p1); print(p2); print(p3); print(p4); print(p5)

cat("\n✓ Integration complete!\n")

# ============================================
# PART 6: Find Cluster Markers
# ============================================
cat("\n=== PART 6: Finding Cluster Markers ===\n\n")

cat("Total clusters:", length(unique(seurat$seurat_clusters)), "\n\n")

markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25)

top10 <- markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(markers, "all_markers_decontx.csv", row.names = FALSE)
write.csv(top10, "top10_markers_per_cluster_decontx.csv", row.names = FALSE)

for(i in sort(unique(top10$cluster))) {
  cat("\n========== Cluster", i, "==========\n")
  cluster_genes <- top10 %>% filter(cluster == i)
  print(cluster_genes[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])
}

cat("\n\nSaved: top10_markers_per_cluster_decontx.csv\n")

# ============================================
# PART 7: Cell Type Annotation (WITH T CELL RE-BINNING)
# ============================================
cat("\n=== PART 7: Cell Type Annotation ===\n\n")

# Initial cell type annotations
cell_types <- c(
  "0" = "Activated T cells",
  "1" = "Naive T cells",
  "2" = "Effector T cells",
  "3" = "CD8+ T cells",
  "4" = "Monocyte-derived DCs",
  "5" = "Plasma cells",
  "6" = "T cells",              # Will re-annotate this!
  "7" = "Macrophages",
  "8" = "cDC1",
  "9" = "cDC2",
  "10" = "Proliferating",
  "11" = "Stromal",
  "12" = "Gamma delta T cells",
  "13" = "Tregs",
  "14" = "pDC",
  "15" = "Eosinophils",
  "16" = "NK cells",
  "17" = "Basophils",
  "18" = "ISG+ cells",
  "19" = "B cells",
  "20" = "Acinar cells",
  "21" = "Neutrophils",
  "22" = "M2 Macrophages",
  "23" = "M2 Macrophages",       # Merged with 22
  "24" = "cDC2",
  "25" = "Neutrophils",
  "26" = "Endothelial"
)

seurat@meta.data$cell_type_base <- cell_types[as.character(seurat@meta.data$seurat_clusters)]

# ============================================
# RE-ANNOTATE GENERIC "T CELLS" (Cluster 6)
# ============================================
cat("\n=== Re-annotating generic 'T cells' cluster ===\n\n")

tcell_generic <- seurat@meta.data$cell_type_base == "T cells"
cat("Generic 'T cells' to re-annotate:", sum(tcell_generic), "\n")

if(sum(tcell_generic) > 0) {
  
  seurat_tcells <- subset(seurat, cells = colnames(seurat)[tcell_generic])
  data <- GetAssayData(seurat_tcells, layer = "data")
  
  # Calculate comprehensive T cell subtype scores
  tcell_scores <- data.frame(
    # CD8+ vs CD4+
    CD8_score = colMeans(data[c("Cd8a", "Cd8b1"), , drop = FALSE]),
    CD4_score = colMeans(data[c("Cd4"), , drop = FALSE]),
    
    # Naive markers
    Naive_score = colMeans(data[c("Ccr7", "Sell", "Lef1", "Tcf7"), , drop = FALSE]),
    
    # Activated markers
    Activated_score = colMeans(data[c("Cd69", "Cd44", "Tnfrsf9"), , drop = FALSE]),
    
    # Effector markers
    Effector_score = colMeans(data[c("Gzma", "Gzmb", "Prf1", "Nkg7", "Ccl5"), , drop = FALSE]),
    
    # Treg markers
    Treg_score = colMeans(data[c("Foxp3", "Il2ra", "Ctla4", "Ikzf2"), , drop = FALSE]),
    
    # Gamma-delta markers
    GD_score = colMeans(data[c("Trdc", "Trgc1", "Trgc2"), , drop = FALSE])
  )
  
  # Assign T cell subtypes based on marker expression
  new_annotations <- character(nrow(tcell_scores))
  
  for(i in 1:nrow(tcell_scores)) {
    scores <- tcell_scores[i, ]
    
    # Check for Tregs first (highest specificity)
    if(scores$Treg_score > 0.5) {
      new_annotations[i] <- "Tregs"
    }
    # Check for gamma-delta
    else if(scores$GD_score > 0.3) {
      new_annotations[i] <- "Gamma delta T cells"
    }
    # Check CD8+ vs CD4+
    else if(scores$CD8_score > scores$CD4_score) {
      # CD8+ T cells
      if(scores$Effector_score > 1.0) {
        new_annotations[i] <- "CD8+ Effector T cells"
      } else if(scores$Naive_score > 0.5) {
        new_annotations[i] <- "CD8+ Naive T cells"
      } else if(scores$Activated_score > 0.5) {
        new_annotations[i] <- "CD8+ Activated T cells"
      } else {
        new_annotations[i] <- "CD8+ T cells"
      }
    }
    else if(scores$CD4_score > scores$CD8_score) {
      # CD4+ T cells
      if(scores$Effector_score > 1.0) {
        new_annotations[i] <- "CD4+ Effector T cells"
      } else if(scores$Naive_score > 0.5) {
        new_annotations[i] <- "CD4+ Naive T cells"
      } else if(scores$Activated_score > 0.5) {
        new_annotations[i] <- "CD4+ Activated T cells"
      } else {
        new_annotations[i] <- "CD4+ T cells"
      }
    }
    # Double negative or unclear
    else {
      if(scores$Effector_score > 1.0) {
        new_annotations[i] <- "Effector T cells"
      } else if(scores$Activated_score > 0.5) {
        new_annotations[i] <- "Activated T cells"
      } else if(scores$Naive_score > 0.5) {
        new_annotations[i] <- "Naive T cells"
      } else {
        new_annotations[i] <- "T cells"
      }
    }
  }
  
  # Update the annotations
  seurat@meta.data$cell_type_base[tcell_generic] <- new_annotations
  
  cat("\nRe-annotated T cell breakdown:\n")
  print(table(new_annotations))
}

# ============================================
# ASSIGN PARENT TYPES FOR PROLIFERATING & ISG+
# ============================================

assign_parent_type <- function(seurat_subset) {
  data <- GetAssayData(seurat_subset, layer = "data")
  
  scores <- data.frame(
    T_cell = colMeans(data[c("Cd3e", "Cd3d"), , drop = FALSE]),
    CD8_T = colMeans(data[c("Cd8a", "Cd8b1"), , drop = FALSE]),
    CD4_T = colMeans(data[c("Cd4"), , drop = FALSE]),
    Myeloid = colMeans(data[c("Lyz2", "Cd68", "Csf1r"), , drop = FALSE]),
    NK = colMeans(data[c("Ncr1", "Klrb1c", "Nkg7"), , drop = FALSE]),
    B_cell = colMeans(data[c("Cd19", "Ms4a1", "Cd79a"), , drop = FALSE]),
    DC = colMeans(data[c("Itgax", "H2-Ab1"), , drop = FALSE]),
    Plasma = colMeans(data[c("Jchain", "Mzb1"), , drop = FALSE])
  )
  
  parent_types <- apply(scores, 1, function(x) names(which.max(x)))
  
  parent_map <- c(
    "T_cell" = "T cells",
    "CD8_T" = "CD8+ T cells",
    "CD4_T" = "CD4+ T cells",
    "Myeloid" = "Myeloid",
    "NK" = "NK cells",
    "B_cell" = "B cells",
    "DC" = "DCs",
    "Plasma" = "Plasma cells"
  )
  
  return(parent_map[parent_types])
}

# Proliferating cells
cat("\nAssigning parent types to proliferating cells...\n")
prolif_cells <- seurat@meta.data$cell_type_base == "Proliferating"
if(sum(prolif_cells) > 0) {
  seurat_prolif <- subset(seurat, cells = colnames(seurat)[prolif_cells])
  parent_types_prolif <- assign_parent_type(seurat_prolif)
  seurat@meta.data$cell_type[prolif_cells] <- paste0("Proliferating ", parent_types_prolif)
  
  cat("Proliferating cell breakdown:\n")
  print(table(parent_types_prolif))
}

# ISG+ cells
cat("\nAssigning parent types to ISG+ cells...\n")
isg_cells <- seurat@meta.data$cell_type_base == "ISG+ cells"
if(sum(isg_cells) > 0) {
  seurat_isg <- subset(seurat, cells = colnames(seurat)[isg_cells])
  parent_types_isg <- assign_parent_type(seurat_isg)
  seurat@meta.data$cell_type[isg_cells] <- paste0("ISG+ ", parent_types_isg)
  
  cat("ISG+ cell breakdown:\n")
  print(table(parent_types_isg))
}

# All other cells
other_cells <- !(prolif_cells | isg_cells)
seurat@meta.data$cell_type[other_cells] <- seurat@meta.data$cell_type_base[other_cells]
seurat@meta.data$cell_type_base <- NULL

cat("\n=== FINAL cell type counts ===\n")
print(sort(table(seurat@meta.data$cell_type), decreasing = TRUE))

qsave(seurat, "seurat_decontx_annotated.qs")
cat("\nSaved: seurat_decontx_annotated.qs\n")

# ============================================
# PART 8: Visualization with Custom Colors
# ============================================
cat("\n=== PART 8: Creating Visualizations ===\n\n")

custom_colors <- c(
  # T cells - all subtypes
  "Activated T cells" = "#4E79A7",
  "Naive T cells" = "#A0CBE8",
  "Effector T cells" = "#59A14F",
  "CD8+ T cells" = "#8CD17D",
  "CD8+ Naive T cells" = "#B6E3C4",
  "CD8+ Activated T cells" = "#6FB587",
  "CD8+ Effector T cells" = "#4A9C6E",
  "CD4+ T cells" = "#B6992D",
  "CD4+ Naive T cells" = "#D9C78F",
  "CD4+ Activated T cells" = "#C4AB5C",
  "CD4+ Effector T cells" = "#A88F2A",
  "T cells" = "#86BCB6",
  "Gamma delta T cells" = "#499894",
  "Tregs" = "#F28E2B",
  
  # Proliferating
  "Proliferating T cells" = "#5A8AC6",
  "Proliferating CD8+ T cells" = "#76B7B2",
  "Proliferating CD4+ T cells" = "#B07AA1",
  "Proliferating Myeloid" = "#D4A373",
  "Proliferating DCs" = "#C49A00",
  "Proliferating NK cells" = "#79706E",
  "Proliferating B cells" = "#86BCB6",
  "Proliferating Plasma cells" = "#8DBDD6",
  
  # ISG+
  "ISG+ T cells" = "#9C755F",
  "ISG+ CD8+ T cells" = "#BAB0AC",
  "ISG+ CD4+ T cells" = "#D37295",
  "ISG+ Myeloid" = "#CA6680",
  "ISG+ DCs" = "#F1CE63",
  "ISG+ NK cells" = "#A0CBE8",
  "ISG+ B cells" = "#9D7660",
  "ISG+ Plasma cells" = "#B5A7D5",
  
  # Myeloid
  "Macrophages" = "#E15759",
  "M2 Macrophages" = "#FF9DA7",
  "Monocyte-derived DCs" = "#F28E2B",
  "Myeloid" = "#FFBE7D",
  
  # DCs
  "cDC1" = "#EDC948",
  "cDC2" = "#F1CE63",
  "pDC" = "#FABFD2",
  
  # Granulocytes
  "Neutrophils" = "#D37295",
  "Eosinophils" = "#FABFD2",
  "Basophils" = "#B07AA1",
  
  # NK cells
  "NK cells" = "#59A14F",
  
  # B cells
  "B cells" = "#A0CBE8",
  "Plasma cells" = "#76B7B2",
  
  # Stromal/Other
  "Stromal" = "#BAB0AC",
  "Acinar cells" = "#D4A373",
  "Endothelial" = "#FABFD2"
)

# Handle any missing cell types
present_types <- sort(unique(seurat@meta.data$cell_type))
missing_types <- setdiff(present_types, names(custom_colors))
if (length(missing_types) > 0) {
  cat("Adding colors for:", paste(missing_types, collapse = ", "), "\n")
  additional_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(missing_types))
  names(additional_colors) <- missing_types
  custom_colors <- c(custom_colors, additional_colors)
}

ordered_types <- names(custom_colors)[names(custom_colors) %in% present_types]
seurat@meta.data$cell_type <- factor(seurat@meta.data$cell_type, levels = ordered_types)

# Create UMAP plots
p1 <- DimPlot(seurat, reduction = "umap", group.by = "cell_type", label = TRUE, 
              repel = TRUE, cols = custom_colors, label.size = 2.5) + 
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 6)) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2)))

p2 <- DimPlot(seurat, reduction = "umap", group.by = "genotype", shuffle = TRUE,
              cols = c("BP2_KO" = "#E15759", "Control" = "#4E79A7")) +
  theme_classic()

p3 <- DimPlot(seurat, reduction = "umap", group.by = "sample", shuffle = TRUE) +
  theme_classic()

ggsave("umap_cell_types_decontx.pdf", p1, width = 16, height = 10)
ggsave("umap_genotype_decontx.pdf", p2, width = 10, height = 8)
ggsave("umap_sample_decontx.pdf", p3, width = 12, height = 8)

print(p1); print(p2); print(p3)

# Cell type counts
cat("\n=== Cell type counts ===\n")
print(table(seurat@meta.data$cell_type))

cat("\n=== By genotype ===\n")
ct_by_geno <- table(seurat@meta.data$cell_type, seurat@meta.data$genotype)
print(ct_by_geno)

summary_df <- as.data.frame(ct_by_geno)
colnames(summary_df) <- c("Cell_Type", "Genotype", "Count")
write.csv(summary_df, "cell_type_by_genotype_decontx.csv", row.names = FALSE)

# ============================================
# PART 9: Cell Type Distribution Plots
# ============================================
cat("\n=== PART 9: Cell Type Distribution Plots ===\n\n")

# By sample
counts_sample <- as.data.frame(table(seurat@meta.data$sample, seurat@meta.data$cell_type))
colnames(counts_sample) <- c("Sample", "Cell_Type", "Count")

counts_sample <- counts_sample %>%
  group_by(Sample) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

sample_info <- unique(seurat@meta.data[, c("sample", "genotype")])
counts_sample <- merge(counts_sample, sample_info, by.x = "Sample", by.y = "sample", all.x = TRUE)
counts_sample$Cell_Type <- factor(counts_sample$Cell_Type, levels = ordered_types)

sample_order <- counts_sample %>%
  distinct(Sample, genotype) %>%
  arrange(genotype, Sample) %>%
  pull(Sample)
counts_sample$Sample <- factor(counts_sample$Sample, levels = unique(sample_order))

p_sample <- ggplot(counts_sample, aes(x = Sample, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors[levels(counts_sample$Cell_Type)], drop = FALSE) +
  labs(title = "Cell Type Distribution by Sample", y = "Percentage (%)", x = "", fill = "Cell Type") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100))

ggsave("celltype_distribution_by_sample_decontx.pdf", p_sample, width = 16, height = 8)
ggsave("celltype_distribution_by_sample_decontx.png", p_sample, width = 16, height = 8, dpi = 300)
print(p_sample)

# By genotype
counts_genotype <- as.data.frame(table(seurat@meta.data$genotype, seurat@meta.data$cell_type))
colnames(counts_genotype) <- c("Genotype", "Cell_Type", "Count")

counts_genotype <- counts_genotype %>%
  group_by(Genotype) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

counts_genotype$Cell_Type <- factor(counts_genotype$Cell_Type, levels = ordered_types)

p_genotype <- ggplot(counts_genotype, aes(x = Genotype, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors[levels(counts_genotype$Cell_Type)], drop = FALSE) +
  labs(title = "Cell Type Distribution by Genotype", y = "Percentage (%)", x = "", fill = "Cell Type") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100))

ggsave("celltype_distribution_by_genotype_decontx.pdf", p_genotype, width = 14, height = 8)
ggsave("celltype_distribution_by_genotype_decontx.png", p_genotype, width = 14, height = 8, dpi = 300)
print(p_genotype)

write.csv(counts_sample, "celltype_percentages_by_sample_decontx.csv", row.names = FALSE)
write.csv(counts_genotype, "celltype_percentages_by_genotype_decontx.csv", row.names = FALSE)

cat("\n✓ Distribution plots saved\n")

# ============================================
# PART 10: Final UMAP Visualization
# ============================================
cat("\n=== PART 10: Final UMAP Visualization ===\n\n")

p_umap_final <- DimPlot(
  seurat,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE, repel = TRUE, label.size = 2.5,
  cols = custom_colors[levels(seurat@meta.data$cell_type)]
) +
  ggtitle("UMAP: Cell Types (DecontX-corrected)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 6)
  ) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

ggsave("UMAP_cell_type_decontx_FINAL.pdf", p_umap_final, width = 16, height = 10)
ggsave("UMAP_cell_type_decontx_FINAL.png", p_umap_final, width = 16, height = 10, dpi = 300)
print(p_umap_final)

# Faceted by genotype
p_umap_geno <- DimPlot(
  seurat,
  reduction = "umap",
  group.by = "cell_type",
  split.by = "genotype",
  label = FALSE, 
  cols = custom_colors[levels(seurat@meta.data$cell_type)],
  ncol = 2
) +
  ggtitle("UMAP: Cell Types by Genotype") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 6)) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2)))

ggsave("UMAP_cell_type_by_genotype_decontx_FINAL.pdf", p_umap_geno, width = 18, height = 7)
ggsave("UMAP_cell_type_by_genotype_decontx_FINAL.png", p_umap_geno, width = 18, height = 7, dpi = 300)
print(p_umap_geno)

# ============================================
# PART 11: DGE and ORA Analysis
# ============================================
cat("\n=== PART 11: DGE & ORA Analysis ===\n\n")

dge_dir <- "DGE_analysis_decontx"
ora_dir <- "ORA_analysis_decontx"

dir.create(dge_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ora_dir, "UP"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ora_dir, "DOWN"), showWarnings = FALSE, recursive = TRUE)

cell_types_unique <- unique(seurat@meta.data$cell_type)
cell_types_unique <- cell_types_unique[!is.na(cell_types_unique)]

cat("Analyzing", length(cell_types_unique), "cell types\n")

DefaultAssay(seurat) <- "RNA"

for(ct in cell_types_unique) {
  
  cat("\n========================================\n")
  cat("Processing:", ct, "\n")
  
  seurat_subset <- subset(seurat, subset = cell_type == ct)
  
  n_ko <- sum(seurat_subset@meta.data$genotype == "BP2_KO")
  n_ctrl <- sum(seurat_subset@meta.data$genotype == "Control")
  
  cat("BP2_KO cells:", n_ko, "| Control cells:", n_ctrl, "\n")
  
  if(n_ko < 3 || n_ctrl < 3) {
    cat("Skipping - not enough cells in both groups\n")
    next
  }
  
  cat("Running DGE...\n")
  Idents(seurat_subset) <- seurat_subset@meta.data$genotype
  
  tryCatch({
    markers <- FindMarkers(seurat_subset, 
                           ident.1 = "BP2_KO", 
                           ident.2 = "Control",
                           logfc.threshold = 0,
                           min.pct = 0.1,
                           test.use = "wilcox")
    
    markers$gene <- rownames(markers)
    
    clean_name <- gsub(" ", "_", ct)
    clean_name <- gsub("/", "_", clean_name)
    clean_name <- gsub("\\+", "plus", clean_name)
    
    dge_file <- file.path(dge_dir, paste0(clean_name, "_BP2KOvsControl.csv"))
    write.csv(markers, dge_file, row.names = FALSE)
    cat("Saved DGE:", dge_file, "\n")
    
    up_genes <- filter(markers, p_val_adj < 0.1 & avg_log2FC > 0)$gene
    down_genes <- filter(markers, p_val_adj < 0.1 & avg_log2FC < 0)$gene
    
    cat("Upregulated:", length(up_genes), "| Downregulated:", length(down_genes), "\n")
    
    # ORA for upregulated genes
    if(length(up_genes) > 0) {
      cat("Running ORA for UP genes...\n")
      go_up <- gost(up_genes, organism = "mmusculus", significant = TRUE, 
                    user_threshold = 0.2, correction_method = "fdr", 
                    domain_scope = "annotated", sources = "GO", evcodes = TRUE)
      
      if(!is.null(go_up$result)) {
        go_up_res <- as.data.frame(go_up$result)
        
        if("p_value" %in% names(go_up_res)) {
          names(go_up_res)[names(go_up_res) == "p_value"] <- "hypergeometric_FDR"
        }
        
        drop_cols <- c("parents", "source_order", "effective_domain_size", 
                       "query", "precision", "recall", "evidence_codes")
        go_up_res <- select(go_up_res, -any_of(drop_cols))
        
        ora_up_file <- file.path(ora_dir, "UP", paste0(clean_name, "_BP2KOvsControl.csv"))
        write.csv(go_up_res, ora_up_file, row.names = FALSE)
        cat("Saved ORA UP:", ora_up_file, "\n")
      } else {
        cat("No significant GO terms for UP genes\n")
      }
    }
    
    # ORA for downregulated genes
    if(length(down_genes) > 0) {
      cat("Running ORA for DOWN genes...\n")
      go_down <- gost(down_genes, organism = "mmusculus", significant = TRUE,
                      user_threshold = 0.2, correction_method = "fdr",
                      domain_scope = "annotated", sources = "GO", evcodes = TRUE)
      
      if(!is.null(go_down$result)) {
        go_down_res <- as.data.frame(go_down$result)
        
        if("p_value" %in% names(go_down_res)) {
          names(go_down_res)[names(go_down_res) == "p_value"] <- "hypergeometric_FDR"
        }
        
        drop_cols <- c("parents", "source_order", "effective_domain_size",
                       "query", "precision", "recall", "evidence_codes")
        go_down_res <- select(go_down_res, -any_of(drop_cols))
        
        ora_down_file <- file.path(ora_dir, "DOWN", paste0(clean_name, "_BP2KOvsControl.csv"))
        write.csv(go_down_res, ora_down_file, row.names = FALSE)
        cat("Saved ORA DOWN:", ora_down_file, "\n")
      } else {
        cat("No significant GO terms for DOWN genes\n")
      }
    }
    
  }, error = function(e) {
    cat("Error processing", ct, ":", e$message, "\n")
  })
}

# ============================================
# FINAL SUMMARY
# ============================================
cat("\n========================================\n")
cat("🚀 AD ASTRA PER ASPERA! 🚀\n")
cat("✓✓✓ COMPLETE ANALYSIS FINISHED! ✓✓✓\n")
cat("========================================\n")
cat("DecontX-corrected & annotated object: seurat_decontx_annotated.qs\n")
cat("Total cells:", ncol(seurat), "\n")
cat("Mean contamination removed:", round(mean(seurat$decontX_contamination), 4), "\n")
cat("DGE results in:", dge_dir, "\n")
cat("ORA results in:", ora_dir, "\n")
cat("\nBy genotype:\n")
print(table(seurat$genotype))
cat("\nBy cell type (Top 20):\n")
print(head(sort(table(seurat$cell_type), decreasing = TRUE), 20))
cat("\n========================================\n")
cat("🎖️ NEVER GIVE UP, NEVER SURRENDER! 🎖️\n")
cat("========================================\n")


