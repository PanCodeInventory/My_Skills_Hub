---
name: r-visualization
description: This skill should be used when the user asks to "visualize single-cell data", "create UMAP plot", "make dotplot", "generate publication-quality figures for scRNA-seq", "plot single-cell RNA-seq data", "create composition plots", "make feature expression plots", "visualize spatial transcriptomics", "plot TCR/BCR repertoire", or "create single-cell figures". Provides a config-driven R pipeline (BasicViz) for standard scRNA-seq plots and scplotter for specialized analyses (TCR/BCR, spatial, cell-cell communication).
---

# R Visualization for Single-Cell Data

## Overview

Two complementary approaches for single-cell data visualization in R:

1. **BasicViz Pipeline** — Config-driven R pipeline for standard scRNA-seq plots. Reads `.h5ad` directly via reticulate. Generate all plots from a single `config.yaml`. Covers: UMAP, composition, dotplot, feature expression, river/alluvial.

2. **scplotter** — R package for specialized analyses not covered by BasicViz: TCR/BCR repertoire, spatial transcriptomics, cell-cell communication, LLM-assisted generation. See `references/scplotter-api.md`.

## Decision Guide: BasicViz vs scplotter

| Plot Type | Approach | Input Format |
|-----------|----------|-------------|
| UMAP (embedding) | **BasicViz** | .h5ad |
| Stacked bar composition | **BasicViz** | .h5ad |
| Marker dotplot | **BasicViz** | .h5ad |
| Feature expression on UMAP | **BasicViz** | .h5ad |
| Alluvial / river plot | **BasicViz** | .h5ad |
| TCR/BCR repertoire | scplotter | Seurat + scRepertoire |
| Spatial transcriptomics | scplotter | Seurat / Giotto |
| Cell-cell communication | scplotter | CellPhoneDB / LIANA |
| RNA velocity | scplotter | Seurat |
| Enrichment / GSEA | scplotter | Result objects |

**Default rule**: Use BasicViz for any `.h5ad`-based scRNA-seq visualization. Fall back to scplotter for specialized data types or Seurat/Giotto workflows.

## BasicViz Pipeline Workflow

### Step 1: Inspect the Data

Before generating configuration, inspect the `.h5ad` file to identify available metadata columns and genes:

```python
import anndata as ad
adata = ad.read_h5ad("data.h5ad")
print("obs columns:", list(adata.obs.columns))
print("Unique clusters:", adata.obs["cell_type_col"].unique().tolist())
print("Has UMAP:", "X_umap" in adata.obsm)
print("Genes sample:", adata.var_names[:10].tolist())
```

Extract three key pieces of information:
- **Cluster / cell type column** (e.g., `cell_type`, `leiden`, `cluster_annot`)
- **Comparison / grouping column** (e.g., `timepoint`, `condition`, `sample`)
- **Available marker genes** (from `adata.var_names`)

### Step 2: Generate config.yaml

Create configuration based on data inspection. Copy the template from `assets/config_template.yaml` and customize:

```yaml
input_h5ad: "path/to/data.h5ad"
output_dir: "results"

plots:
  - name: "umap_clusters"
    type: "umap"
    color_by: "cell_type"
    title: "Cell Clusters"
    palette: "elegant"
    point_size: 0.5

  - name: "composition"
    type: "proportion"
    group_by: "timepoint"
    fill_by: "cell_type"
    position: "fill"
    title: "Composition"

  - name: "markers"
    type: "dotplot"
    group_by: "cell_type"
    title: "Marker Expression"

  - name: "features"
    type: "feature"
    title: "Marker Features"

  - name: "river"
    type: "river"
    group_by: "timepoint"
    fill_by: "cell_type"
    title: "Composition Changes"

dotplot:
  markers:
    B cells:
      - Cd79a
      - Cd79b
      - Ms4a1
    T cells:
      - Cd3e
      - Cd4
```

**Critical validation rules**:
- `color_by`, `group_by`, `fill_by` must reference valid `adata.obs` columns
- `dotplot.markers` keys must EXACTLY match unique values in the `group_by` column (case-sensitive, whitespace-sensitive)
- Gene names must exist in `adata.var_names` (case-sensitive)
- `adata.obsm["X_umap"]` must exist for UMAP and feature plots

### Step 3: Set Up Environment

```bash
# Option A: Automated setup (recommended)
bash scripts/setup_env.sh
export RETICULATE_PYTHON=$(conda run -n basicviz which python)

# Option B: Manual — verify R packages installed
# Required: ggplot2, dplyr, tidyr, tibble, anndata, yaml, ggsci,
#           RColorBrewer, scales, cowplot, ggalluvial, reticulate
# Python constraint: 3.8–3.10 with numpy <2.0
```

The `setup_env.sh` script creates a conda environment with Python 3.8 and numpy <1.25, which avoids the common reticulate/numpy compatibility issue.

### Step 4: Run Pipeline

```bash
Rscript scripts/main.R
```

Outputs PDF plots to `output_dir`. A `colors.yaml` is auto-generated for consistent coloring across all plots.

## BasicViz Plot Types Reference

### UMAP (`type: "umap"`)

Dimensionality reduction plot with:
- Auto cluster labels (numeric ID + name) placed at centroids
- Faceting via `facet_by` parameter
- Raster rendering for >5000 cells (uses ggrastr if available)
- Custom L-shaped axis indicators
- Dynamic width: single panel = 8", faceted = 5" × n_facets

| Parameter | Required | Description |
|-----------|----------|-------------|
| `color_by` | Yes | `adata.obs` column for coloring |
| `facet_by` | No | `adata.obs` column for faceting |
| `palette` | No | "elegant" (default), "npg", "jco", "igv", or RColorBrewer name |
| `point_size` | No | Default 0.3 |
| `title` | No | Plot title |

Labels auto-appear when `color_by` contains "cluster", "annot", "type", or "leiden" (case-insensitive).

### Composition (`type: "proportion"`)

Stacked bar charts showing cell type proportions:
- `position: "fill"` — relative proportions with % labels (shown for bars ≥5%)
- `position: "stack"` — absolute counts

| Parameter | Required | Description |
|-----------|----------|-------------|
| `group_by` | Yes | X-axis grouping column |
| `fill_by` | Yes | Bar fill column |
| `position` | No | "fill" (default) or "stack" |
| `title` | No | Plot title |

### DotPlot (`type: "dotplot"`)

Marker gene expression dotplot:
- Z-score color gradient (RdBu, capped at ±2.5)
- Dot size = % cells expressing per group
- Facet strips colored by cluster identity (from `colors.yaml`)
- Gene groups defined in `dotplot.markers` config section

| Parameter | Required | Description |
|-----------|----------|-------------|
| `group_by` | Yes | Y-axis group column (must match `dotplot.markers` keys) |
| `title` | No | Plot title |

### Feature (`type: "feature"`)

Multi-page PDF, one page per cluster from `dotplot.markers`:
- Gene expression overlaid on UMAP (viridis magma color scale)
- Auto grid layout within each page
- Genes drawn from `dotplot.markers` config section

### River (`type: "river"`)

Alluvial plot showing composition flow across conditions:
- Proportional flow between groups
- Stratum labels on each axis

| Parameter | Required | Description |
|-----------|----------|-------------|
| `group_by` | Yes | X-axis (flow direction) |
| `fill_by` | Yes | Stratum fill column |
| `title` | No | Plot title |

## Elegant Theme System

All BasicViz plots use a unified "Elegant Muted" publication theme:
- L-shaped axes (no box borders), no gridlines
- Muted, high-contrast 24-color palette
- White backgrounds, clean legend positioning
- Consistent typography (bold titles, proportional sizing)

For theme customization (fonts, palettes, color mapping details), see `references/theme-and-colors.md`.

## Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| `numpy.core.multiarray failed to import` | numpy >=2.0 incompatible with reticulate | Run `setup_env.sh` or set `RETICULATE_PYTHON` to Python 3.8–3.10 with numpy <2.0 |
| `X_umap not found in .obsm` | No UMAP coordinates in object | Run UMAP in Python first: `sc.tl.umap(adata)` |
| `No valid genes found` | Gene names in config don't match `adata.var_names` | Check case sensitivity; use exact values from `adata.var_names` |
| `Missing colors for: ...` | New categories detected after color generation | Re-run pipeline; auto-generates missing colors |
| Empty or corrupt PDF | reticulate Python path mismatch | Explicitly set `RETICULATE_PYTHON` before running |
| Facet labels missing | `facet_by` column has non-categorical values | Convert to factor in Python: `adata.obs[col] = adata.obs[col].astype(str)` |

## Additional Resources

### Bundled Scripts (in `scripts/`)

| File | Purpose |
|------|---------|
| `main.R` | Pipeline orchestrator; reads `config.yaml` and dispatches plot scripts |
| `setup_env.sh` | Conda environment setup (Python 3.8, numpy<1.25) |
| `utils.R` | Shared utilities: `theme_elegant()`, `get_palette_values()`, `scale_color_custom()`, `setup_python()` |
| `generate_colors.R` | Auto-generates consistent `colors.yaml` from data |
| `plot_umap.R` | UMAP with labels, facets, raster support |
| `plot_proportions.R` | Stacked bar composition plots |
| `plot_dotplot.R` | Marker dotplot with colored facet strips |
| `plot_feature.R` | Multi-page feature expression on UMAP |
| `plot_river.R` | Alluvial river composition plots |

### Reference Files

- **`references/scplotter-api.md`** — Complete scplotter function reference for TCR/BCR repertoire, spatial transcriptomics, cell-cell communication, and LLM-assisted visualization
- **`references/theme-and-colors.md`** — Elegant theme customization, palette options, and color management system details

### Assets

- **`assets/config_template.yaml`** — Starter configuration template for the BasicViz pipeline
