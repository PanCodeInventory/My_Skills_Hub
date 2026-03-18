---
name: r-visualization
description: R-based single-cell and spatial data visualization using scplotter. Use for creating publication-quality plots for scRNA-seq, scTCR-seq, scBCR-seq, and spatial transcriptomics data. Supports Seurat, Giotto, and anndata formats with LLM-assisted plot generation.
license: GPL-3.0 license
metadata:
    skill-author: K-Dense Inc.
    source: https://github.com/pwwang/scplotter
---

# R Visualization with scplotter

## Overview

scplotter is an R package built upon [`plotthis`](https://github.com/pwwang/plotthis) that provides a comprehensive set of functions to visualize single-cell sequencing and spatial data. Use this skill when working with:

- scRNA-seq data (dimensionality reduction, clustering, marker genes)
- scTCR-seq/scBCR-seq data (clonal analysis, repertoire diversity)
- Spatial transcriptomics data (Visium, Xenium, SlideSeq, CODEX, etc.)
- Cell-cell communication analysis
- LLM-assisted plot generation

## When to Use This Skill

This skill should be used when:
- Creating publication-quality visualizations for single-cell RNA-seq data
- Analyzing TCR/BCR repertoire sequencing data
- Visualizing spatial transcriptomics from 10x Visium, Xenium, Nanostring CosMx, etc.
- Generating cell-cell communication plots from CellPhoneDB/LIANA results
- Working with Seurat or Giotto objects
- Loading and visualizing anndata (.h5ad) files in R
- Using LLM assistance to generate visualization code automatically

## Installation

```r
# From GitHub
remotes::install_github("pwwang/scplotter")
# or
devtools::install_github("pwwang/scplotter")

# Using conda
# conda install pwwang::r-scplotter

# Load the package
library(scplotter)
```

## Quick Start Templates

### Template 1: scRNA-seq Visualization

```r
library(scplotter)
library(Seurat)

# Load data (example with Seurat object)
data(ifnb_sub)  # Built-in dataset

# Dimensionality reduction plot (UMAP/t-SNE)
CellDimPlot(ifnb_sub, reduction = "umap", group.by = "cell_type")

# Cell statistics plot
CellStatPlot(ifnb_sub, group.by = "cell_type", 
             metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# Feature expression on embedding
FeatureStatPlot(ifnb_sub, features = c("CD3D", "CD14", "MS4A1"),
                plot_type = "violin")

# Marker genes visualization (volcano plot)
MarkersPlot(marker_results, top_n = 20)

# GSEA results
GSEASummaryPlot(gsea_results)
GSEAPlot(gsea_results, term = "INTERFERON_ALPHA_RESPONSE")
```

### Template 2: TCR/BCR Repertoire Analysis

```r
library(scplotter)
library(scRepertoire)

# Load TCR/BCR data
data(tcr_data)  # Your tcr object

# Clonal volume distribution
ClonalVolumePlot(tcr_data, group.by = "sample")

# Clonal abundance
ClonalAbundancePlot(tcr_data, top_n = 20)

# Clonal residency across tissues
ClonalResidencyPlot(tcr_data, group.by = "tissue")

# Clonal composition
ClonalCompositionPlot(tcr_data, group.by = "cell_type")

# Clonal overlap between samples
ClonalOverlapPlot(tcr_data, samples = c("sample1", "sample2"))

# V/J gene usage
ClonalGeneUsagePlot(tcr_data, gene_type = "v")

# Clonal diversity analysis
ClonalDiversityPlot(tcr_data, group.by = "condition")

# CDR3 length distribution
ClonalLengthPlot(tcr_data, chain = "TRA")

# Rarefaction curve
ClonalRarefactionPlot(tcr_data)

# K-mer analysis
ClonalKmerPlot(tcr_data, k = 3)

# Positional analysis
ClonalPositionalPlot(tcr_data)
```

### Template 3: Spatial Transcriptomics

```r
library(scplotter)
library(Giotto)  # or Seurat

# Load spatial data (Giotto object)
data(spatial_obj)

# Spatial dimension reduction plot
SpatDimPlot(spatial_obj, reduction = "umap", group.by = "cell_type")

# Spatial feature expression
SpatFeaturePlot(spatial_obj, features = c("Gene1", "Gene2"),
                plot_type = "heatmap")

# For Seurat spatial objects
CellDimPlot(seurat_spatial, reduction = "spatial", group.by = "cluster")
```

### Template 4: Cell-Cell Communication

```r
library(scplotter)
library(LIANA)

# Load CellPhoneDB/LIANA results
data(cellphonedb_res)

# Dot plot of interactions
CCCPlot(cellphonedb_res, plot_type = "dot",
        top_n = 20, group.by = "cell_type")

# Heatmap of interactions
CCCPlot(cellphonedb_res, plot_type = "heatmap",
        ligand_target = "CD40", receptor_target = "CD40LG")

# Circle plot
CCCPlot(cellphonedb_res, plot_type = "circle")
```

## Core Visualization Functions

### scRNA-seq Functions

| Function | Description | Use Case |
|----------|-------------|----------|
| `CellDimPlot()` | Dimensionality reduction plot | UMAP/t-SNE/PCA visualization |
| `CellStatPlot()` | Cell statistics plot | QC metrics distribution |
| `CellVelocityPlot()` | RNA velocity plot | Trajectory/directionality |
| `FeatureStatPlot()` | Feature statistic plot | Gene expression violin/box plots |
| `ClustreePlot()` | Clustering tree | Compare clustering resolutions |
| `EnrichmentPlot()` | Enrichment analysis plot | GO/KEGG enrichment |
| `MarkersPlot()` | Marker gene visualization | Volcano plots for DEGs |
| `GSEASummaryPlot()` | GSEA summary | Enrichment map/dotplot |
| `GSEAPlot()` | GSEA enrichment plot | Running enrichment score |
| `CCCPlot()` | Cell-cell communication | Ligand-receptor interactions |

### scTCR-seq/scBCR-seq Functions

| Function | Description | Use Case |
|----------|-------------|----------|
| `ClonalVolumePlot()` | Clonal volume | Clone size distribution |
| `ClonalAbundancePlot()` | Clonal abundance | Top clones by frequency |
| `ClonalLengthPlot()` | CDR3 length | CDR3 length distribution |
| `ClonalResidencyPlot()` | Clonal residency | Tissue residency patterns |
| `ClonalCompositionPlot()` | Clonal composition | Clone makeup by group |
| `ClonalOverlapPlot()` | Clonal overlap | Shared clones between samples |
| `ClonalStatPlot()` | Clonal statistics | Summary statistics |
| `ClonalDiversityPlot()` | Clonal diversity | Diversity indices |
| `ClonalGeneUsagePlot()` | V/J gene usage | V(D)J gene frequency |
| `ClonalPositionalPlot()` | Positional analysis | CDR3 positional patterns |
| `ClonalKmerPlot()` | K-mer analysis | Motif enrichment |
| `ClonalRarefactionPlot()` | Rarefaction curve | Sampling depth assessment |

### Spatial Data Functions

| Function | Description | Use Case |
|----------|-------------|----------|
| `SpatDimPlot()` | Spatial dimension plot | Clusters on tissue |
| `SpatFeaturePlot()` | Spatial feature plot | Gene expression on tissue |

## Working with Different Data Formats

### Seurat Objects

```r
library(scplotter)
library(Seurat)

# Load Seurat object
seurat_obj <- readRDS("seurat_object.rds")

# Most scplotter functions work directly with Seurat objects
CellDimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")

# Feature plots
FeatureStatPlot(seurat_obj, features = c("CD3D", "CD8A", "CD4"),
                plot_type = "box")
```

### Giotto Objects

```r
library(scplotter)
library(Giotto)

# Load Giotto object
giotto_obj <- loadGiottoObject("giotto_data")

# Spatial plots
SpatDimPlot(giotto_obj, reduction = "umap", group.by = "cluster")
SpatFeaturePlot(giotto_obj, features = c("Gene1", "Gene2"))
```

### AnnData (.h5ad) Files

```r
library(scplotter)
library(reticulate)

# Load anndata file
adata <- read_h5ad("data.h5ad")

# Convert to Seurat or work directly
seurat_obj <- ConvertToSeurat(adata)
CellDimPlot(seurat_obj, reduction = "umap")
```

See `articles/Working_with_anndata_h5ad_files.html` for detailed workflow.

## Clone Selectors (for TCR/BCR analysis)

Helper functions to select clones based on various criteria:

```r
# Select top N clones
top_clones <- top(tcr_data, n = 10)

# Select clones with specific criteria
selected <- sel(tcr_data, count > 5)

# Unique clones
unique_clones <- uniq(tcr_data)

# Shared clones between samples
shared_clones <- shared(tcr_data, samples = c("A", "B"))

# Comparison operators
gt(tcr_data, count, 10)    # greater than
ge(tcr_data, count, 10)    # greater or equal
lt(tcr_data, count, 10)    # less than
le(tcr_data, count, 10)    # less or equal
eq(tcr_data, count, 10)    # equal
ne(tcr_data, count, 10)    # not equal

# Logical operators
and(gt(tcr_data, count, 5), lt(tcr_data, count, 100))
or(eq(tcr_data, tissue, "blood"), eq(tcr_data, tissue, "tumor"))
```

## LLM-Assisted Visualization

scplotter includes built-in LLM support for automatic plot generation:

```r
library(scplotter)
library(tidyprompt)

# Setup LLM provider
provider <- tidyprompt::llm_provider_openai(
    api_key = Sys.getenv("OPENAI_API_KEY")
)

# Create chat instance
chat <- SCPlotterChat$new(provider = provider)

# Ask for visualization
chat$ask("Generate a UMAP plot colored by cell type for the ifnb_sub data")
# Automatically identifies: CellDimPlot(ifnb_sub, group.by = "cell_type")

# Iterate on visualization
chat$ask("Make it a dot plot instead")
chat$ask("Show only T cells and B cells")
chat$ask("Add gene expression overlay for CD3D")
```

See `articles/Visualizing_data_with_LLMs.html` for complete workflow.

## Publication-Quality Customization

### General Plot Styling

```r
library(scplotter)
library(ggplot2)

# Set global theme
theme_set(theme_minimal(base_size = 12))

# Custom color palette
custom_colors <- c("B cell" = "#1f77b4", "T cell" = "#ff7f0e", 
                   "Monocyte" = "#2ca02c", "NK" = "#d62728")

CellDimPlot(seurat_obj, group.by = "cell_type",
            palette = custom_colors,
            point_size = 2,
            alpha = 0.7,
            legend.title = "Cell Type",
            title = "Single Cell Atlas")

# Save publication quality
ggsave("figure1.pdf", width = 8, height = 6, dpi = 300)
```

### Multi-Panel Figures

```r
library(patchwork)

# Create multiple plots
p1 <- CellDimPlot(obj, group.by = "cell_type")
p2 <- CellStatPlot(obj, group.by = "cell_type", metrics = "nFeature_RNA")
p3 <- FeatureStatPlot(obj, features = top_markers, plot_type = "heatmap")

# Combine
combined <- p1 + p2 + p3 + 
            plot_layout(ncol = 2, heights = c(2, 1))

ggsave("multi_panel_figure.pdf", combined, width = 12, height = 10)
```

## Built-in Datasets

```r
# Load example datasets
data(ifnb_sub)        # Subsetted IFNB-stimulated PBMCs
data(pancreas_sub)    # Mouse pancreas dataset
data(cellphonedb_res) # CellPhoneDB results example

# Explore dataset
?ifnb_sub
str(ifnb_sub)
```

## Common Workflows

### Workflow 1: Complete scRNA-seq Analysis

```r
library(scplotter)
library(Seurat)

# 1. Load and QC
seurat_obj <- readRDS("processed.rds")

# 2. Dimensionality reduction
CellDimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")

# 3. QC metrics visualization
CellStatPlot(seurat_obj, group.by = "seurat_clusters",
             metrics = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

# 4. Marker gene visualization
markers <- FindAllMarkers(seurat_obj)
MarkersPlot(markers, top_n = 20, group = "0")

# 5. Feature expression
FeatureStatPlot(seurat_obj, features = c("CD3D", "CD14", "MS4A1"),
                plot_type = "violin")

# 6. Save
ggsave("umap_clusters.pdf", width = 8, height = 6)
```

### Workflow 2: Spatial Analysis

```r
library(scplotter)
library(Giotto)

# 1. Load spatial data
spatial_obj <- loadGiottoObject("visium_data")

# 2. Spatial clustering
SpatDimPlot(spatial_obj, group.by = "cluster", point_size = 1.5)

# 3. Gene expression on tissue
SpatFeaturePlot(spatial_obj, features = c("GFAP", "MBP", "NEUN"),
                plot_type = "dot", spot_size = 2)

# 4. Combined view
p1 <- SpatDimPlot(spatial_obj)
p2 <- SpatFeaturePlot(spatial_obj, features = "Gene1")
p1 + p2
```

### Workflow 3: TCR Repertoire Analysis

```r
library(scplotter)
library(scRepertoire)

# 1. Load TCR data
tcr_obj <- readRDS("tcr_data.rds")

# 2. Clonal statistics
ClonalStatPlot(tcr_obj, group.by = "sample")

# 3. Top clones
ClonalAbundancePlot(tcr_obj, top_n = 30)

# 4. Diversity analysis
ClonalDiversityPlot(tcr_obj, group.by = "condition",
                    metric = "shannon")

# 5. Clonal overlap
ClonalOverlapPlot(tcr_obj, samples = c("tumor", "blood"))

# 6. V gene usage
ClonalGeneUsagePlot(tcr_obj, gene_type = "v", top_n = 15)
```

## Key Parameters Reference

### CellDimPlot
- `reduction`: "umap", "tsne", "pca"
- `group.by`: Metadata column for coloring
- `palette`: Custom color vector
- `point_size`: Point size (default: 1)
- `alpha`: Transparency (0-1)

### CellStatPlot
- `group.by`: Grouping variable
- `metrics`: QC metrics to plot
- `plot_type`: "box", "violin", "dot"

### FeatureStatPlot
- `features`: Gene names to plot
- `plot_type`: "violin", "box", "dot", "heatmap"
- `group.by`: Grouping variable

### CCCPlot
- `plot_type`: "dot", "heatmap", "circle"
- `top_n`: Top N interactions to show
- `ligand_target`, `receptor_target`: Specific interactions

## Tips for Effective Visualization

1. **Always start with CellDimPlot**: Get overview of cell types/clusters
2. **Use CellStatPlot for QC**: Check data quality before interpretation
3. **Choose appropriate plot types**:
   - Violin plots for distribution comparison
   - Dot plots for many genes/groups
   - Heatmaps for patterns across groups
4. **Customize colors**: Use consistent color schemes across figures
5. **Use patchwork**: Combine related plots into multi-panel figures
6. **LLM assistance**: Use SCPlotterChat for rapid iteration
7. **Export high resolution**: Always use 300+ DPI for publication

## Troubleshooting

### "Object not compatible with scplotter"
- Ensure object is Seurat, Giotto, or compatible class
- Check that required slots are populated (reductions, metadata)

### "No reduction found"
- Run dimensionality reduction first: `RunUMAP()`, `RunTSNE()`
- Specify correct reduction name: `reduction = "umap"`

### "Color palette mismatch"
- Provide custom palette matching number of groups
- Use `unique()` to count group levels

### Spatial data not displaying
- Verify spatial coordinates are present in object
- Check `SpatDimPlot` vs `CellDimPlot` for spatial data

## Additional Resources

- **Official documentation**: https://pwwang.github.io/scplotter/
- **GitHub repository**: https://github.com/pwwang/scplotter
- **plotthis (core)**: https://github.com/pwwang/plotthis
- **Seurat**: https://satijalab.org/seurat/
- **Giotto**: https://drieslab.github.io/Giotto_website/
- **scRepertoire**: https://github.com/BorchLab/scRepertoire
- **LIANA**: https://github.com/saezlab/liana-py

## Articles and Vignettes

- **Visualizing data with LLMs**: `articles/Visualizing_data_with_LLMs.html`
- **Working with anndata files**: `articles/Working_with_anndata_h5ad_files.html`
- **Spatial data guide**: `articles/Knowing_your_spatial_data_and_visualization.html`

### Spatial Data Examples

**Seurat-based:**
- 10x Visium: `articles/Seurat_10x_Visium.html`
- VisiumHD: `articles/Seurat_10x_VisiumHD.html`
- Xenium: `articles/Seurat_Xenium.html`
- Nanostring CosMx: `articles/Seurat_Nanostring_CosMx.html`
- SlideSeq: `articles/Seurat_SlideSeq.html`

**Giotto-based:**
- Visium: `articles/Giotto_Visium.html`
- Xenium: `articles/Giotto_Xenium.html`
- CODEX: `articles/Giotto_CODEX.html`
- Vizgen: `articles/Giotto_vizgen.html`
- seqFISH: `articles/Giotto_seqFISH.html`
