# Color Schemes for Genomic Data Visualization

## Principles for ATAC-seq Visualization

### Color Selection Guidelines

1. **Use perceptually uniform colormaps** - Equal data changes should appear as equal color changes
2. **Consider colorblind accessibility** - ~8% of males have red-green colorblindness
3. **Use diverging schemes for centered data** - When data has meaningful midpoint (e.g., log2FC)
4. **Use sequential schemes for magnitude** - When showing intensity or density

## Recommended Palettes

### Sequential (Single Direction)

For peak density, signal intensity, or accessibility scores:

```r
# Blues (light to dark)
colorset1 = c("#f7fbff", "#6baed6", "#084594")

# Greens (light to dark)
colorset1 = c("#f7fcf5", "#74c476", "#006d2c")

# Oranges (light to dark)
colorset1 = c("#fff5eb", "#fd8d3c", "#e6550d")

# Purples (light to dark)
colorset1 = c("#fcfbfd", "#9e9ac8", "#3f007d")
```

### Diverging (Two Directions from Center)

For differential accessibility (log2 fold change):

```r
# Blue-White-Red (standard)
colorset1 = c("#4575b4", "#ffffbf", "#d73027")

# Green-White-Purple
colorset1 = c("#1a9850", "#fee08b", "#d73027")

# Brown-Blue-Green (colorblind friendly)
colorset1 = c("#a6611a", "#dfc27d", "#80cdc1", "#018571")

# RdYlBu (common in genomics)
colorset1 = c("#a50021", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6")
```

### Qualitative (Distinct Categories)

For different peak types or sample groups:

```r
# Set1 (9 colors, colorblind safe)
colors = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", 
           "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999")

# Dark2 (8 colors)
colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
           "#66a61e", "#e6ab02", "#a6761d", "#666666")

# Paired (12 colors)
colors = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
           "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
           "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
```

## Implementation in RIdeogram

### Heatmap Colors

```r
# Apply to overlaid data
ideogram(
  karyotype = karyotype,
  overlaid = peak_density,
  colorset1 = c("#4575b4", "#ffffbf", "#d73027"),
  output = "density.svg"
)
```

### Label Colors

```r
# Apply to marker data
ideogram(
  karyotype = karyotype,
  label = label_data,
  label_type = "marker",
  colorset2 = c("#b35806", "#f7f7f7", "#542788"),
  output = "markers.svg"
)
```

## Custom Gradient Generation

### Using colorRampPalette

```r
# Create 100-color gradient
gradient <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b73027"))(100)

# Apply to data based on quantiles
peaks$color_group <- cut(peaks$signal, breaks = 100, labels = FALSE)
peaks$color <- gradient[peaks$color_group]

# Format for RIdeogram (remove #)
label_data$color <- gsub("#", "", peaks$color)
```

### Using RColorBrewer

```r
library(RColorBrewer)

# Sequential
brewer.pal(n = 9, name = "Blues")
brewer.pal(n = 9, name = "Reds")
brewer.pal(n = 9, name = "Greens")

# Diverging
brewer.pal(n = 11, name = "RdBu")
brewer.pal(n = 11, name = "PiYG")
brewer.pal(n = 11, name = "BrBG")

# Qualitative
brewer.pal(n = 8, name = "Dark2")
brewer.pal(n = 9, name = "Set1")
```

### Using Viridis (Perceptually Uniform)

```r
library(viridis)

# Viridis (default)
viridis(100)

# Magma
magma(100)

# Plasma
plasma(100)

# Inferno
inferno(100)

# Cividis (optimized for colorblindness)
cividis(100)
```

## Color by Data Type

### Peak Signal Strength

```r
# Sequential blue (low to high signal)
colorset1 = c("#f7fbff", "#6baed6", "#084594")
```

### Differential Accessibility

```r
# Diverging (negative to positive log2FC)
colorset1 = c("#2c7bb6", "#ffffff", "#d7191c")
#         blue (down)    white    red (up)
```

### Peak Significance

```r
# Binary (significant vs non-significant)
label_data$color <- ifelse(peaks$qval > 10, "d73027", "fee08b")
```

### Chromosome-Specific Colors

```r
# Unique color per chromosome
chr_colors <- setNames(RColorBrewer::brewer.pal(12, "Set3"), 
                       unique(karyotype$Chr))
label_data$color <- chr_colors[label_data$Chr]
```

## Hex Color Reference

| Color | Hex Code | RGB |
|-------|----------|-----|
| Red | #d73027 | (215, 48, 39) |
| Blue | #4575b4 | (69, 117, 180) |
| Green | #1a9850 | (26, 152, 80) |
| Orange | #fd8d3c | (253, 141, 60) |
| Purple | #9e9ac8 | (158, 154, 200) |
| Yellow | #ffffbf | (255, 255, 191) |

## Colorblind Accessibility

### Safe Palettes

```r
# Okabe-Ito (8 colors, universal design)
okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
              "#0072B2", "#D55E00", "#CC79A7", "#000000")
```

### Avoid

- Red-Green combinations (deuteranopia)
- Blue-Purple with low contrast
- Fine details in yellow

### Test Your Colors

```r
# Using colorblindr package
library(colorblindr)

# Simulate colorblind vision
cvd_grid(your_plot)
```

## RIdeogram Color Format

**Important:** RIdeogram expects hex colors WITHOUT the # symbol:

```r
# Correct
color = "ff7f00"

# Wrong
color = "#ff7f00"

# Convert if needed
color <- gsub("#", "", "#ff7f00")
```
