---
name: scrna-de-analysis
description: This skill should be used when the user asks about "single-cell differential expression", "scRNA-seq DE analysis", "pseudobulk analysis", "MAST", "marker genes", "compare conditions in single-cell", "Wilcoxon rank-sum test for scRNA", "DESeq2 for single-cell", or needs guidance on choosing appropriate DE methods for single-cell data. Provides comprehensive code patterns and best practices for differential expression analysis in scRNA-seq.
version: 1.0.0
---

# Single-Cell RNA-seq Differential Expression Analysis

## Overview

Differential expression (DE) analysis in single-cell RNA-seq requires careful method selection to avoid inflated false discovery rates. This skill provides production-ready code patterns for four distinct DE methods, each optimized for specific experimental scenarios.

## Method Selection Guide

| Scenario | Recommended Method | Key Advantage |
|----------|-------------------|---------------|
| Quick marker gene identification | Wilcoxon (Scanpy) | Fast, built-in |
| Multi-biological replicate studies | Pseudobulk + DESeq2/edgeR | Gold standard FDR control |
| Sparse data with high dropout | MAST | Zero-inflation modeling |
| Complex batch effects | SCVI-tools | Integrated batch correction |

## Quick Reference

### Wilcoxon Rank-Sum (Scanpy Native)

```python
import scanpy as sc

# Basic marker gene identification
sc.tl.rank_genes_groups(adata, groupby="cell_type", method="wilcoxon")

# Filter by fold change and expression
sc.tl.filter_rank_genes_groups(
    adata,
    min_fold_change=1.5,
    min_in_group_fraction=0.3
)

# Extract results
df = sc.get.rank_genes_groups_df(adata, group="CD4_T")
```

### Pseudobulk + DESeq2 (Multi-replicate Studies)

```python
from decoupler import pseudobulk

# Aggregate to pseudobulk
pb_df, pb_obs = pseudobulk(
    adata,
    sample_col='donor_id',
    groupby='cell_type',
    layer='counts',
    min_cells=10
)
```

Then use DESeq2 in R or PyDESeq2 for statistical testing.

### SCVI-tools (Deep Learning with Batch Correction)

```python
import scvi

scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
model = scvi.model.SCVI(adata)
model.train()

de_df = model.differential_expression(
    groupby="cell_type",
    group1="CD4_T",
    group2="CD8_T",
    batch_correction=True
)
sig_genes = de_df[de_df["proba_de"] > 0.95]
```

### MAST (Zero-inflated Data)

```r
library(MAST)
zlmFit <- zlm(~ condition + donor, sca = sca_object)
summaryFit <- summary(zlmFit, doLRT = "condition")
```

## Critical Thresholds

| Parameter | Recommended Value | Purpose |
|-----------|------------------|---------|
| `min_fold_change` | 1.5-2.0 | Filter weak effects |
| `min_in_group_fraction` | 0.25-0.30 | Ensure expression in target |
| `max_out_group_fraction` | 0.5 | Ensure specificity |
| `pvals_adj` | < 0.05 | FDR-corrected significance |
| `proba_de` (SCVI) | > 0.95 | Bayesian probability |
| `min_cells` (pseudobulk) | >= 10 | Reliable aggregation |

## Common Pitfalls

1. **Using Wilcoxon for multi-replicate studies**: Causes pseudoreplication bias and inflated FDR
2. **Applying batch correction before DE**: May remove true biological signal
3. **Ignoring donor effects**: Persists even after batch correction
4. **Simple t-tests**: Poor variance estimation with zero-inflated data

## Decision Flow

```
                    Have biological replicates?
                           /            \
                         YES             NO
                          |               |
               Use Pseudobulk       Complex batches?
           (DESeq2/edgeR)              /        \
                                     YES          NO
                                      |            |
                              Use SCVI-tools   Use Wilcoxon
                              (batch_corr=True)  (quick markers)
```

## When to Use This Skill

Apply this skill when:
- Comparing gene expression between conditions in scRNA-seq
- Identifying marker genes for cell clusters
- Analyzing differential expression with multiple donors/patients
- Handling batch effects in DE analysis
- Choosing between DE methods for single-cell data

## Reference Files

For detailed implementations and code examples:

- **`references/wilcoxon-method.md`** - Complete Wilcoxon workflow with parameter tuning
- **`references/pseudobulk-method.md`** - Pseudobulk aggregation and DESeq2/edgeR analysis
- **`references/mast-method.md`** - MAST zero-inflated model implementation
- **`references/scvi-de-method.md`** - SCVI-tools deep learning DE with batch correction
- **`references/best-practices.md`** - 2024-2025 benchmarks and anti-patterns

## Example Files

Working code examples in `examples/`:
- **`examples/wilcoxon_markers.py`** - Marker gene identification workflow
- **`examples/pseudobulk_pipeline.py`** - End-to-end pseudobulk DE analysis
- **`examples/scvi_de_analysis.py`** - SCVI differential expression

## Key References

1. Squair et al. (2021) Nat Commun - Pseudoreplication bias in scRNA-seq DE
2. Nguyen et al. (2023) Nat Commun - Benchmarking scRNA-seq DE methods
3. Wu et al. (2025) Genome Biology - GLIMES paradigm and DE pitfalls
4. Single-cell Best Practices: https://www.sc-best-practices.org/conditions/differential_gene_expression.html
