# CellChat Analysis Skill

A comprehensive skill for cell-cell communication analysis using CellChat.

## Quick Start

### Activation

This skill automatically activates when users ask about:
- "analyze cell-cell communication"
- "run CellChat"
- "cell communication analysis"
- "ligand-receptor interaction"
- "cellchat visualization"
- "intercellular communication"
- "CellChatDB"

### Basic Usage

```r
# Load CellChat
library(CellChat)

# Create object from Seurat
cellchat <- createCellChat(
  object = GetAssayData(seurat_obj, slot = "data"),
  meta = seurat_obj@meta.data,
  group.by = "cell_type"
)

# Run analysis pipeline
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- aggregateNet(cellchat)

# Visualize
netVisual_circle(cellchat@net$count)
netVisual_bubble(cellchat)
```

## Skill Structure

```
cellchat-analysis/
├── SKILL.md                          # Main skill documentation
├── references/
│   ├── workflow-detailed.md          # Complete step-by-step guide
│   ├── visualization-guide.md        # All visualization options
│   ├── database-reference.md         # CellChatDB customization
│   └── troubleshooting.md            # Common issues & solutions
├── examples/
│   ├── basic-analysis.R              # Minimal working example
│   ├── seurat-integration.R          # Seurat workflow
│   ├── spatial-analysis.R            # Spatial transcriptomics
│   └── comparative-analysis.R        # Multi-condition comparison
└── scripts/
    └── setup-environment.R           # Environment setup
```

## Key Features

### 📊 Analysis Types
- **Single Dataset**: Basic cell-cell communication
- **Seurat Integration**: Direct Seurat object input
- **Spatial Analysis**: Spatial transcriptomics support
- **Comparative**: Multi-condition comparisons

### 🎨 Visualization Types
- Circle networks
- Hierarchical layouts
- Chord diagrams
- Bubble plots
- Heatmaps
- Scatter plots

### 🧬 Supported Databases
- CellChatDB.human (default)
- CellChatDB.mouse
- Custom database support

### 🔧 Customization
- Database subsetting
- Custom interactions
- Parameter tuning
- Multi-species support

## Environment Requirements

### System
- R >= 4.0
- Bioconductor >= 3.18
- 16GB+ RAM recommended

### Conda Environment
```bash
mamba create -n cellchat -c conda-forge r-base=4.3
mamba activate cellchat
```

### R Packages
Core packages:
- CellChat (GitHub: jinworks/CellChat)
- Seurat (optional)
- dplyr, ggplot2, igraph
- BiocManager, BiocNeighbors
- NMF, circlize, ComplexHeatmap

## Usage Examples

### Example 1: Basic Analysis

See `examples/basic-analysis.R` for a complete minimal example.

```r
# Load skill guidance
source("examples/basic-analysis.R")
```

### Example 2: Seurat Integration

See `examples/seurat-integration.R` for Seurat workflow.

```r
# Full Seurat to CellChat pipeline
cellchat <- createCellChat(
  object = GetAssayData(seurat_obj, slot = "data"),
  meta = seurat_obj@meta.data,
  group.by = "cell_type"
)
```

### Example 3: Comparative Analysis

See `examples/comparative-analysis.R` for comparing conditions.

```r
# Merge multiple datasets
cellchat.list <- list(Control = chat1, Treatment = chat2)
cellchat.merged <- mergeCellChat(cellchat.list)

# Compare
netVisual_diffInteraction(cellchat.merged)
```

## Common Workflows

### Workflow 1: Quick Analysis

```r
# 1. Load data
library(CellChat)
cellchat <- createCellChat(object = expr, meta = meta, group.by = "type")

# 2. Analyze
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat) %>%
  identifyOverExpressedGenes() %>%
  identifyOverExpressedInteractions() %>%
  computeCommunProb(raw.use = TRUE) %>%
  filterCommunication(min.cells = 10) %>%
  aggregateNet()

# 3. Visualize
netVisual_circle(cellchat@net$count)
netVisual_bubble(cellchat)
```

### Workflow 2: Pathway-Focused

```r
# Analyze specific pathways
pathways <- c("WNT", "BMP", "VEGF")
for (p in pathways) {
  netVisual_aggregate(cellchat, signaling = p, layout = "circle")
}
```

### Workflow 3: Export Results

```r
# Save results
saveRDS(cellchat, "cellchat_results.rds")
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "interactions.csv")
```

## Visualization Gallery

### Network Visualizations
- **Circle plots**: Overview of all interactions
- **Hierarchical**: Shows directionality
- **Chord diagrams**: Alternative network view

### Statistical Visualizations
- **Bubble plots**: Interaction details with p-values
- **Heatmaps**: Pathway activity across cell types
- **Scatter plots**: Signaling roles

### Comparative Visualizations
- **Differential networks**: Compare conditions
- **Stacked bar plots**: Rank pathways
- **Side-by-side**: Direct comparisons

## Troubleshooting

### Common Issues

See `references/troubleshooting.md` for detailed solutions.

**Quick fixes:**

```r
# NMF installation fails
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
devtools::install_github("renozao/NMF")

# Memory issues
cellchat <- computeCommunProb(cellchat, population.size = TRUE)

# No interactions found
cellchat <- filterCommunication(cellchat, min.cells = 5)
```

## Database Customization

See `references/database-reference.md` for detailed database information.

### Subsetting

```r
# Use only secreted signaling
cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")

# Include specific pathways
cellchat@DB$interaction <- cellchat@DB$interaction[
  cellchat@DB$interaction$pathway_name %in% c("WNT", "BMP"),
]
```

### Custom Interactions

```r
# Add custom LR pairs
custom <- data.frame(
  interaction_name = "Ligand_Receptor",
  pathway_name = "CUSTOM",
  ligand = "Ligand",
  receptor = "Receptor",
  annotation = "Secreted Signaling"
)
cellchat@DB$interaction <- rbind(cellchat@DB$interaction, custom)
```

## Performance Tips

### Large Datasets (>50,000 cells)

1. Use `population.size = TRUE`
2. Enable parallel processing
3. Process in chunks
4. Filter cell types first

```r
# Parallel processing
library(future)
plan("multisession", workers = 4)

# Memory-efficient
cellchat <- computeCommunProb(cellchat, 
                               raw.use = TRUE,
                               population.size = TRUE)
```

### Speed Optimization

1. Use fast algorithms
2. Reduce permutation tests
3. Pre-filter cell types

```r
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = TRUE)
cellchat <- computeCommunProb(cellchat, nrun = 100)
```

## Citation

When using CellChat in research:

```
Jin et al., Inference and analysis of cell-cell communication using CellChat, 
Nature Communications 2021

Jin et al., CellChat for systematic analysis of cell–cell communication 
from single-cell transcriptomics, Nature Protocols 2024
```

## Resources

- **GitHub**: https://github.com/jinworks/CellChat
- **Tutorials**: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
- **Protocol**: https://www.nature.com/articles/s41596-024-01045-4

## Support

For issues:
1. Check `references/troubleshooting.md`
2. Visit GitHub Issues: https://github.com/jinworks/CellChat/issues
3. Review examples in `examples/`

## License

This skill is for educational and research purposes. CellChat is licensed under GPL-3.0.

---

**Version**: 1.0.0  
**Last Updated**: 2025  
**Compatible with**: CellChat >= 2.2.0, R >= 4.0
