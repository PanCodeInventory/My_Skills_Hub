---
name: ATAC-seq Chromosome Visualization
description: This skill should be used when the user asks to "可视化 ATAC-seq 染色体", "visualize ATAC-seq chromosome", "create ideogram", "plot ATAC-seq peaks on chromosome", or mentions RIdeogram for ATAC-seq visualization. Provides guidance for chromosome-scale ATAC-seq visualization using R with pre-configured r_env environment.
version: 0.1.0
---

# ATAC-seq Chromosome Visualization Skill

## Purpose

This skill provides workflows for ATAC-seq data analysis using R, with focus on genome-wide visualization using RIdeogram. All analysis runs in the pre-configured `r_env` mamba environment containing required packages.

## Environment Setup

### Activate R Environment

Before any R analysis, activate the conda/mamba environment:

```bash
mamba activate r_env
```

Verify R and packages are available:

```bash
R --version
R -e "library(RIdeogram); library(GenomicRanges); library(dplyr)"
```

## Core Workflow: ATAC-seq Visualization with RIdeogram

### Step 1: Prepare Input Data

RIdeogram requires three data structures:

| Data Type | Format | Purpose |
|-----------|--------|---------|
| Karyotype | 3 or 5 columns | Chromosome boundaries and centromeres |
| Overlaid (heatmap) | 4 columns | Continuous values (peak density) |
| Label (markers) | 6 columns | Discrete features (specific peaks) |

### Step 2: Convert ATAC-seq Peaks to RIdeogram Format

Transform BED/narrowPeak files from MACS2 or similar peak callers:

```r
# Read ATAC-seq peaks
peaks <- read.table("peaks.narrowPeak", header = FALSE)
colnames(peaks) <- c("chrom", "start", "end", "name", "score", 
                     "strand", "signal", "pval", "qval", "peak")

# Filter significant peaks (q-value < 0.1)
peaks_filtered <- peaks[peaks$qval > 10, ]

# Create label data
label_data <- data.frame(
  Type = "marker",
  Shape = "box",
  Chr = peaks_filtered$chrom,
  Start = peaks_filtered$start,
  End = peaks_filtered$end,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b73027"))(50)[
    cut(peaks_filtered$signal, breaks = 50)]
)
```

### Step 3: Calculate Peak Density (Heatmap Mode)

For genome-wide accessibility patterns:

```r
# Calculate density in windows (e.g., 10Mb)
window_size <- 10000000
peak_density <- data.frame()

for (chr in unique(karyotype$Chr)) {
  chr_end <- karyotype[karyotype$Chr == chr, "End"]
  windows <- seq(0, chr_end, by = window_size)
  
  for (i in 1:(length(windows)-1)) {
    w_start <- windows[i] + 1
    w_end <- windows[i+1]
    n_peaks <- sum(peaks_filtered$chrom == chr & 
                   peaks_filtered$start >= w_start & 
                   peaks_filtered$end <= w_end)
    
    peak_density <- rbind(peak_density, data.frame(
      Chr = chr, Start = w_start, End = w_end, Value = n_peaks
    ))
  }
}
```

### Step 4: Generate Ideogram

```r
library(RIdeogram)

# Load human karyotype (built-in) or custom
data(human_karyotype)

# Create ideogram with markers
ideogram(
  karyotype = human_karyotype,
  label = label_data,
  label_type = "marker",
  output = "atac_peaks.svg"
)

# Convert to PNG (300 DPI)
convertSVG("atac_peaks.svg", device = "png", dpi = 300)
```

## Customization Options

### Color Schemes

```r
# Blue-White-Red for heatmaps
colorset1 = c("#4575b4", "#ffffbf", "#d73027")

# Orange-White-Purple for labels
colorset2 = c("#b35806", "#f7f7f7", "#542788")
```

### Plot Dimensions

```r
ideogram(
  karyotype = karyotype,
  width = 170,    # Plot width
  Lx = 160,       # Legend X position
  Ly = 35,        # Legend Y position
  output = "custom.svg"
)
```

### Shape Types

Available shapes for labels: `box`, `triangle`, `circle`

## Best Practices

| Task | Recommendation |
|------|----------------|
| Peak filtering | q-value < 0.05 or score > 100 |
| Window size | 1-10Mb for genome-wide view |
| Label limit | Top 50-100 peaks to avoid clutter |
| Export format | SVG for publication, PNG for preview |
| Resolution | 300 DPI minimum for publication |

## Additional Resources

### Reference Files

For detailed documentation and examples:

- **`references/rideogram-guide.md`** - Complete RIdeogram API reference and ATAC-seq adaptation guide
- **`references/data-formats.md`** - Input data format specifications and conversion utilities
- **`references/color-schemes.md`** - Color palette recommendations for genomic data

### Example Files

Working code examples in `examples/`:

- **`basic-ideogram.R`** - Minimal working example with built-in human data
- **`atac-peaks-full.R`** - Complete workflow from BED file to publication figure
- **`density-heatmap.R`** - Peak density calculation and heatmap overlay

### Utility Scripts

Automation scripts in `scripts/`:

- **`convert-bed-to-rideogram.sh`** - Convert BED/narrowPeak to RIdeogram label format
- **`calculate-density.R`** - Calculate peak density in configurable windows
- **`validate-input.R`** - Validate input data formats before plotting

## Citation

When publishing work using RIdeogram:

> Hao Z, Lv D, Ge Y, Shi J, Weijers D, Yu G, Chen J. 2020. RIdeogram: drawing SVG graphics to visualize and map genome-wide data on the idiograms. *PeerJ Computer Science* 6:e251. https://doi.org/10.7717/peerj-cs.251

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| Package not found | Run `mamba activate r_env` first |
| Chromosome mismatch | Ensure Chr names match between karyotype and data |
| SVG not rendering | Check file path is writable |
| Colors not applying | Verify color hex format (6 characters, no #) |

### Get Help

For complex issues, consult:
1. `references/rideogram-guide.md` for detailed API
2. `examples/` for working code patterns
3. RIdeogram vignette: https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html
