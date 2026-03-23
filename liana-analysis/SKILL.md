---
name: liana-analysis
description: This skill should be used when the user asks to "analyze cell-cell communication", "run LIANA", "ligand-receptor analysis", "CellPhoneDB", "CellChat", "NATMI", "single-cell communication", "CCC inference", "intercellular communication", "LR interaction", or needs guidance on cell-cell communication (CCC) inference methods. Provides comprehensive guidance for LIANA+ all-in-one framework for ligand-receptor interaction analysis from single-cell RNA-seq data.
version: 1.0.0
---

# LIANA+: Cell-Cell Communication Analysis

## Overview

LIANA+ (LIgand-receptor ANalysis frAmework) is a unified framework for inferring and consensus-scoring cell-cell communication (CCC) from single-cell transcriptomics data. Rather than running a single method in isolation, LIANA+ aggregates results across 8+ CCC inference methods (CellPhoneDB, CellChat, NATMI, Consensus, etc.) and produces consensus ligand-receptor interaction scores combining both communication **magnitude** and **specificity**.

Key capabilities:
- Run multiple CCC methods through a single interface
- Produce consensus scores (`magnitude_rank` + `specificity_rank`) for robust interaction ranking
- Generate publication-ready visualizations (bubble plots, dot plots, chord diagrams)
- Support for spatial transcriptomics via `liana.mt.bivar`
- Extensible resource system with curated ligand-receptor databases

## When to Use This Skill

- Inferring ligand-receptor interactions from scRNA-seq data
- Comparing CCC across multiple methods (CellPhoneDB, CellChat, NATMI, etc.)
- Generating consensus-ranked communication networks
- Analyzing spatially-resolved cell-cell communication
- Validating cell-type annotations through communication patterns
- Investigating signaling pathways between specific cell populations
- Cross-condition differential communication analysis

## Prerequisites

- Python >= 3.10
- `scanpy`, `anndata`, `pandas`, `numpy`, `matplotlib`, `seaborn`
- LIANA+ installed via pip or conda/mamba

```bash
# pip install
pip install liana

# or via conda-forge (recommended)
mamba install -c conda-forge liana
```

For spatial analysis, also install `musclex` or `spatialdata` as needed.

## Core Analysis Workflow

### Step 1: Load and Prepare Data

Load an annotated AnnData object with cell-type labels. The object should contain raw (un-normalized) counts in `adata.X` or `adata.raw.X`.

```python
import scanpy as sc
import liana as li

# Load pre-processed data with cell type annotations
adata = sc.read_h5ad("annotated_data.h5ad")

# Ensure cell type labels exist in adata.obs
# Required: adata.obs must contain a grouping column (e.g., 'cell_type')
assert 'cell_type' in adata.obs.columns

# Use raw counts for LIANA (recommended)
# liana uses adata.raw by default if available
```

### Step 2: Run LIANA+

Execute the consensus pipeline across all methods. The primary entry point is `liana.mt.rank_aggregate`.

```python
# Run the consensus rank aggregate across all methods
liana_res = li.mt.rank_aggregate(
    adata,
    groupby='cell_type',         # cell type annotation column
    resource_name='consensus',   # LR resource (default)
    expr_prop=0.1,              # minimum expression proportion
    use_raw=True,                # use adata.raw.X
    n_perms=1000,               # permutations for statistical testing (CellPhoneDB)
    verbose=True,
)

# Result is a DataFrame with columns:
# - ligand, receptor, source, target
# - magnitude_rank, specificity_rank
# - per-method scores (e.g. cellphonedb_mean, cellchat, natmi)
```

### Step 3: Visualize Results

Generate publication-ready plots to interpret communication patterns.

```python
# Top interactions bubble plot (default)
li.pl.dotplot(
    liana_res,
    uns_keys=['magnitude_rank', 'specificity_rank'],
    figure_size=(8, 10),
    source_labels=['B cells', 'CD4 T cells'],
    target_labels=['CD8 T cells', 'NK cells'],
)

# Heatmap of interactions
li.pl.heatmap(
    liana_res,
    figure_size=(8, 10),
)

# Chord diagram for inter-cellular communication network
li.pl.chord_graph(
    liana_res,
    source_labels=['B cells', 'CD4 T cells'],
    target_labels=['CD8 T cells', 'NK cells'],
)

# Tile plot for method agreement
li.pl.tileplot(
    liana_res,
    source_labels=['B cells'],
    target_labels=['CD8 T cells'],
)
```

### Step 4: Interpret Results

LIANA+ produces two core consensus metrics:

- **`magnitude_rank`**: Ranks interactions by overall communication strength (combines expression level and interaction probability across methods). Lower values = stronger interactions.
- **`specificity_rank`**: Ranks interactions by how specifically they occur between source-target pairs. Lower values = more specific interactions.

Combine both metrics: the strongest, most specific interactions appear at the top-left of dot plots (low values in both ranks).

```python
# Filter for top interactions
top_interactions = liana_res[
    (liana_res['magnitude_rank'] < 0.05) &
    (liana_res['specificity_rank'] < 0.05)
]

# Inspect individual method scores
top_interactions[['ligand', 'receptor', 'source', 'target',
                   'magnitude_rank', 'specificity_rank',
                   'cellphonedb', 'cellchat', 'natmi']]
```

### Step 5: Export Results

```python
# Save results as CSV
liana_res.to_csv("liana_results.csv", index=False)

# Save filtered top interactions
top_interactions.to_csv("liana_top_interactions.csv", index=False)

# Save plots
import matplotlib.pyplot as plt
plt.savefig("liana_dotplot.pdf", dpi=300, bbox_inches='tight')
```

## Available Methods

LIANA+ integrates multiple CCC methods through a unified interface:

| Method | Key Score | Description |
|--------|-----------|-------------|
| **CellPhoneDB** | `cellphonedb` | Permutation-based statistical testing for significant LR pairs |
| **CellChat** | `cellchat` | Probabilistic model with signaling pathway analysis |
| **NATMI** | `natmi` | Weighted expression product with specificity weighting |
| **Consensus** | (meta) | Aggregate of multiple curated resources |
| **Connectome** | `connectome` | Expression-weighted mean with edge specificity |
| **Geometric Mean** | `geometric_mean` | Geometric mean of per-method scores |
| **SingleCellSignalR** | `sca` | Network-based scoring with signaling flow |
| **log2FC** | `log2fc` | Log2 fold-change based scoring |

Consensus scores (`magnitude_rank`, `specificity_rank`) aggregate across all methods using a rank-based approach.

Refer to `references/methods_comparison.md` for a detailed comparison of methods, their scores, assumptions, and appropriate use cases.

## Key Parameters

### `liana.mt.rank_aggregate`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `groupby` | required | Column in `adata.obs` for cell type labels |
| `resource_name` | `'consensus'` | LR resource database to use |
| `expr_prop` | `0.1` | Minimum proportion of cells expressing gene in a group |
| `use_raw` | `True` | Use `adata.raw.X` for expression values |
| `n_perms` | `1000` | Number of permutations (affects CellPhoneDB) |
| `min_cells` | `5` | Minimum cells per group |
| `return_all_lrs` | `False` | Return all LR pairs regardless of expression |
| `seed` | `0` | Random seed for reproducibility |

### Resource Options

Common resources available via `resource_name`:
- `'consensus'` — aggregated from multiple curated databases (default)
- `'cellphonedb'` — CellPhoneDB v5 database
- `'cellchatdb'` — CellChat v2 database
- `'connectomedb'` — ConnectomeDB 2020
- `'icellnet'` — iCellNet database
- `'guide2pharma'` — Guide2Pharma curated interactions

## Best Practices

1. **Use raw counts**: Set `use_raw=True` and ensure `adata.raw` contains un-normalized counts. Normalized or scaled data can distort interaction scores.
2. **Require sufficient cells per group**: Filter cell types with fewer than 10-20 cells. Small groups produce unreliable statistics. Adjust `min_cells` accordingly.
3. **Choose expression proportion carefully**: The `expr_prop` parameter (default 0.1) controls the minimum fraction of cells that must express a gene. Increase to 0.2 for noisy datasets; decrease to 0.05 for rare cell types.
4. **Inspect method agreement**: Check whether top interactions are supported by multiple methods using tile plots. Interactions ranked highly by only one method may be less reliable.
5. **Validate with known biology**: Cross-reference top interactions with known signaling pathways for the tissue/cell types under study.
6. **Compare conditions separately**: For differential communication, run LIANA+ on each condition independently, then compare results.
7. **Use spatial validation when available**: Combine with `liana.mt.bivar` for spatial transcriptomics to confirm physical proximity of interacting cells.
8. **Report resource version**: Always document which LR resource and LIANA+ version were used for reproducibility.

## Common Pitfalls

- **Gene symbol mismatches**: LIANA+ expects gene symbols matching the LR resource. Verify gene names in `adata.var_names` match the resource (typically human: HGNC symbols; mouse: MGI symbols). Use `sc.pp.convert_genes_to_symbols()` if needed.
- **Expression proportion too strict**: Setting `expr_prop` too high (>0.3) filters out biologically relevant rare interactions. Start at 0.1 and adjust based on data quality.
- **Insufficient permutations**: Low `n_perms` (<100) reduces statistical power for CellPhoneDB. Use 1000+ for publication-quality results.
- **Missing cell type annotations**: LIANA+ requires a grouping variable. Ensure cell type labels exist in `adata.obs` before running.
- **Over-interpreting individual method scores**: Consensus scores are more robust than any single method. Prioritize `magnitude_rank` and `specificity_rank` over individual method outputs.
- **Forgetting to set `use_raw=True`**: Using normalized/scaled data instead of raw counts can lead to incorrect interaction scores.

## Bundled Resources

### scripts/setup_environment.sh
Bash script to create a mamba environment with LIANA+ and all dependencies:

```bash
bash scripts/setup_environment.sh liana_env
```

### scripts/lr_analysis.py
CLI tool for automated LR analysis with configurable parameters:

```bash
python scripts/lr_analysis.py input.h5ad --groupby cell_type --output results/
```

### assets/analysis_template.py
Full workflow template with configurable parameters. Copy and customize:

```bash
cp assets/analysis_template.py my_communication_analysis.py
python my_communication_analysis.py
```

### examples/basic_analysis.py
Minimal working example using pbmc68k data.

### examples/spatial_analysis.py
Spatial bivariate analysis example for spatial transcriptomics.

## Additional Resources

### Reference Files

- **`references/workflow_detailed.md`** — Complete step-by-step workflow with code, covering data preparation, running individual methods, spatial analysis, differential communication, and advanced filtering.
- **`references/api_reference.md`** — Quick API reference for `liana.mt`, `liana.pl`, and `liana.rs` modules with function signatures and parameters.
- **`references/methods_comparison.md`** — Detailed comparison of 8+ CCC methods, their scoring mechanisms, assumptions, strengths, and recommended use cases.

### Example Files

- **`examples/basic_analysis.py`** — Minimal working example with pbmc68k data
- **`examples/spatial_analysis.py`** — Spatial bivariate analysis with squidpy

## Citation

Dimitrov, D., Türei, D., Hege, N., Dugourd, A. & Saez-Rodriguez, J. Comprehensive characterization of cell-cell communication. *Nature Cell Biology* **26**, 591–603 (2024). https://doi.org/10.1038/s41556-024-0744-4

## Related Skills

- **scanpy** — For preprocessing and QC before running LIANA+ (prerequisite)
- **scientific-visualization** — For publication-quality figure customization
- **anndata** — For AnnData format questions and data structure
