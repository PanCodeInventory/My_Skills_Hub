---
name: CellChat Analysis
description: This skill should be used when the user asks to "analyze cell-cell communication", "run CellChat", "cell communication analysis", "ligand-receptor interaction", "cellchat visualization", "intercellular communication", or mentions CellChatDB. Provides comprehensive guidance for cell-cell communication analysis from single-cell and spatial transcriptomics data.
version: 1.0.0
---

# CellChat Analysis Skill

Guide for inference, visualization, and analysis of cell-cell communication from single-cell and spatially resolved transcriptomics using CellChat.

## Overview

CellChat is an R package that enables systematic analysis of cell-cell communication by inferring intercellular signaling networks from gene expression data. This skill provides workflows for setting up the environment, running analyses, and interpreting results.

## When to Use

- Analyzing cell-cell communication from scRNA-seq data
- Identifying ligand-receptor interactions between cell types
- Visualizing signaling networks and pathways
- Comparing communication patterns across conditions
- Analyzing spatial transcriptomics data
- Investigating cell microenvironment interactions

## Prerequisites

### Environment Setup

Activate the cellchat environment before starting:

```bash
mamba activate cellchat
R
```

### Required R Packages

Load necessary packages:

```r
library(CellChat)
library(Seurat)  # If using Seurat objects
library(dplyr)
library(ggplot2)
```

## Core Analysis Workflow

### Step 1: Data Preparation

Create CellChat object from Seurat object:

```r
data.input <- GetAssayData(seurat_obj, slot = "data")
meta <- seurat_obj@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
```

Or from expression matrix:

```r
cellchat <- createCellChat(object = expression_matrix, meta = cell_labels, group.by = "labels")
```

### Step 2: Database Configuration

Set up the ligand-receptor database:

```r
# For human data
cellchat@DB <- CellChatDB.human

# For mouse data
# cellchat@DB <- CellChatDB.mouse

# View database categories
showDatabaseCategory(CellChatDB)
```

### Step 3: Data Preprocessing

Identify overexpressed genes and interactions:

```r
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
```

### Step 4: Communication Network Inference

Calculate communication probabilities:

```r
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
```

### Step 5: Network Visualization

Generate network visualizations:

```r
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

# Circle plot
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = TRUE, label.edge = FALSE)

# Bubble plot
netVisual_bubble(cellchat, remove.isolate = TRUE)
```

### Step 6: Signaling Pathway Analysis

Analyze specific pathways:

```r
# Get available pathways
pathways <- cellchat@netP$pathways

# Visualize specific pathway
netVisual_aggregate(cellchat, signaling = "WNT", layout = "circle")

# Compute network centrality
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = "WNT")
```

## Quick Reference Commands

| Task | Command |
|------|---------|
| Create object | `createCellChat(object, meta, group.by)` |
| Set database | `cellchat@DB <- CellChatDB.human` |
| Preprocess | `identifyOverExpressedGenes(cellchat)` |
| Compute probabilities | `computeCommunProb(cellchat, raw.use = TRUE)` |
| Network circle plot | `netVisual_circle(cellchat@net$count)` |
| Bubble plot | `netVisual_bubble(cellchat)` |
| Pathway analysis | `netVisual_aggregate(cellchat, signaling = pathway)` |
| Extract network | `subsetCommunication(cellchat)` |

## Common Parameters

### computeCommunProb

- `raw.use = TRUE`: Use raw count data
- `population.size = TRUE`: Account for cell population size
- `type = "triMean"`: Use trimean for expression aggregation

### Visualization Functions

- `layout = "circle"`: Circular layout
- `layout = "hierarchy"`: Hierarchical layout
- `layout = "chord"`: Chord diagram
- `signaling`: Specific pathway name
- `sources.use`: Source cell types
- `targets.use`: Target cell types

## Output Interpretation

### Communication Network

- **Number of interactions**: Count of significant ligand-receptor pairs
- **Interaction strength**: Weighted probability of communication
- **Pathways**: Biological signaling pathways involved

### Key Metrics

- **Prob**: Communication probability (0-1)
- **Pval**: Statistical significance
- **Interaction_name**: Ligand-receptor pair
- **Pathway_name**: Signaling pathway

## Best Practices

1. **Data Quality**: Ensure proper normalization before CellChat analysis
2. **Cell Types**: Use biologically meaningful cell type annotations
3. **Database**: Choose appropriate species database (human/mouse)
4. **Filtering**: Adjust `min.cells` parameter based on dataset size
5. **Validation**: Cross-reference findings with known literature

## Troubleshooting

### Installation Issues

If NMF installation fails:
```r
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
devtools::install_github("renozao/NMF")
```

If circlize/ComplexHeatmap fails:
```r
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
```

### Memory Issues

For large datasets:
```r
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
```

### Visualization Issues

If graphics device fails:
```r
pdf("output.pdf")
netVisual_circle(cellchat@net$count)
dev.off()
```

## Additional Resources

### Reference Files

- **`references/workflow-detailed.md`**: Complete step-by-step workflow with explanations
- **`references/visualization-guide.md`**: Comprehensive visualization options and customization
- **`references/database-reference.md`**: CellChatDB structure and customization
- **`references/troubleshooting.md`**: Common issues and solutions

### Example Files

- **`examples/basic-analysis.R`**: Minimal working example
- **`examples/seurat-integration.R`**: Integration with Seurat workflow
- **`examples/spatial-analysis.R`**: Spatial transcriptomics analysis
- **`examples/comparative-analysis.R`**: Comparing multiple conditions

### Utility Scripts

- **`scripts/setup-environment.R`**: Environment setup and package installation
- **`scripts/run-analysis.R`**: Complete analysis pipeline script
- **`scripts/visualization-batch.R`**: Batch generation of visualizations

## Related Tools

- **Seurat**: For scRNA-seq preprocessing and clustering
- **SingleR**: For automated cell type annotation
- **scCATCH**: Alternative cell type annotation
- **CellPhoneDB**: Alternative cell communication tool
- **NicheNet**: Alternative ligand-target inference

## Citation

If using CellChat in research, cite:

- Jin et al., Inference and analysis of cell-cell communication using CellChat, Nature Communications 2021
- Jin et al., CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics, Nature Protocols 2024

## Version Information

- CellChat Version: 2.2.0+
- R Version: 4.3+
- Bioconductor: 3.18+
