# RIdeogram Complete Guide for ATAC-seq

## Package Overview

RIdeogram is an R package that draws SVG (Scalable Vector Graphics) to visualize and map genome-wide data on idiograms.

### Three Core Functions

| Function | Purpose | Parameters |
|----------|---------|------------|
| `ideogram()` | Create ideogram with overlays | karyotype, overlaid, label, label_type, colorset1, colorset2, width, Lx, Ly, output |
| `convertSVG()` | Convert SVG to raster/PDF | svg, device, width, height, dpi |
| `GFFex()` | Extract features from GFF | input, karyotype, feature, window |

## Installation

```r
# From CRAN
install.packages("RIdeogram")

# Or from R-universe (latest version)
install.packages('RIdeogram', repos = c('https://tickingclock1992.r-universe.dev', 'https://cloud.r-project.org'))

library(RIdeogram)
```

## Input Data Formats

### Karyotype File

Contains chromosome information. Two formats supported:

**3-column format (no centromere):**
```
Chr  Start   End
1    0       248956422
2    0       242193529
3    0       198295559
```

**5-column format (with centromere):**
```
Chr  Start   End     CE_start    CE_end
1    0       248956422 122026459 124932724
2    0       242193529 92188145  94090557
3    0       198295559 90772458  93655574
```

### Overlaid Data (Heatmap)

For showing continuous values like peak density or signal intensity:

```
Chr   Start     End     Value
1     1         1000000  65
1     1000001   2000000  76
1     2000001   3000000  35
```

### Label Data (Markers)

For showing discrete features like specific peaks or genes:

```
Type     Shape    Chr    Start      End     color
marker   box      1      86357632   86357700 ff7f00
marker   triangle 6     69204486   69204568 6a3d9a
gene     circle   3     68882967   68883091 33a02c
```

**Shape options:** `box`, `triangle`, `circle`  
**Type options:** `marker`, `heatmap`, `line`, `polygon`

## ATAC-seq Data Conversion

### From BED Format

Standard BED output from peak callers:
```
chr1    10000   10500   Peak1    50    +
chr1    15000   15800   Peak2    100   +
chr1    25000   26000   Peak3    75    +
```

Convert to RIdeogram label format:
```r
# Read BED file
peaks <- read.table("peaks.bed", header = FALSE)
colnames(peaks) <- c("chrom", "start", "end", "name", "score", "strand")

# Create label data
label_data <- data.frame(
  Type = "marker",
  Shape = "box",
  Chr = peaks$chrom,
  Start = peaks$start,
  End = peaks$end,
  color = ifelse(peaks$score > 100, "ff0000", "ff7f00")
)
```

### From narrowPeak Format

MACS2 narrowPeak output:
```
chr1    10000   10500   Peak1    50    .    100    500    10    5000
```

Columns: chrom, start, end, name, score, strand, signal, pval, qval, peak

```r
peaks <- read.table("peaks.narrowPeak", header = FALSE)
colnames(peaks) <- c("chrom", "start", "end", "name", "score", 
                     "strand", "signal", "pval", "qval", "peak")

# Filter by significance (q-value < 0.1 means qval > 10 in -log10)
peaks_sig <- peaks[peaks$qval > 10, ]

# Color by signal strength
label_data <- data.frame(
  Type = "marker",
  Shape = "box",
  Chr = peaks_sig$chrom,
  Start = peaks_sig$start,
  End = peaks_sig$end,
  color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b73027"))(50)[
    cut(peaks_sig$signal, breaks = 50)]
)
```

## Peak Density Calculation

For genome-wide accessibility patterns, calculate density in windows:

```r
calculate_peak_density <- function(peaks_df, karyotype, window_size = 1000000) {
  result <- data.frame()
  
  for (i in 1:nrow(karyotype)) {
    chr <- karyotype$Chr[i]
    chr_end <- karyotype$End[i]
    
    # Create windows
    starts <- seq(0, chr_end, by = window_size)
    ends <- pmin(starts + window_size, chr_end)
    
    # Count peaks in each window
    chr_peaks <- peaks_df[peaks_df$chrom == chr, ]
    
    for (j in 1:(length(starts) - 1)) {
      peak_count <- sum(chr_peaks$start >= starts[j] & 
                        chr_peaks$end <= ends[j])
      
      result <- rbind(result, data.frame(
        Chr = chr,
        Start = starts[j] + 1,
        End = ends[j],
        Value = peak_count
      ))
    }
  }
  return(result)
}

# Usage
peak_density <- calculate_peak_density(peaks_sig, human_karyotype, window_size = 10000000)
```

## Visualization Examples

### Basic Ideogram

```r
library(RIdeogram)
data(human_karyotype)

# Simple karyotype only
ideogram(karyotype = human_karyotype, output = "basic.svg")
convertSVG("basic.svg", device = "png", dpi = 300)
```

### With Peak Markers

```r
ideogram(
  karyotype = human_karyotype,
  label = label_data,
  label_type = "marker",
  output = "atac_markers.svg"
)
```

### With Density Heatmap

```r
ideogram(
  karyotype = human_karyotype,
  overlaid = peak_density,
  colorset1 = c("#4575b4", "#ffffbf", "#d73027"),
  output = "atac_density.svg"
)
```

### Combined View

```r
ideogram(
  karyotype = human_karyotype,
  overlaid = peak_density,
  label = label_data[1:100, ],  # Top 100 peaks only
  label_type = "marker",
  colorset1 = c("#4575b4", "#ffffbf", "#d73027"),
  output = "atac_combined.svg"
)
```

## Color Schemes

### Predefined Palettes

```r
# Blue-White-Red (diverging)
colorset1 = c("#4575b4", "#ffffbf", "#d73027")

# Orange-White-Purple (diverging)
colorset2 = c("#b35806", "#f7f7f7", "#542788")

# Green-Yellow-Red (sequential)
colorset1 = c("#1a9850", "#fee08b", "#d73027")

# Blue-Purple (sequential)
colorset1 = c("#2166ac", "#67a9cf", "#f7f7f7")
```

### Custom Color Generation

```r
# Generate gradient palette
colors <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b73027"))(100)

# Apply to data based on value ranges
label_data$color <- colors[cut(peaks$signal, breaks = 100)]
```

## Plot Customization

### Dimensions and Legend Position

```r
ideogram(
  karyotype = karyotype,
  width = 170,    # Width of plot region (mm)
  Lx = 160,       # Legend X position
  Ly = 35,        # Legend Y position
  output = "custom.svg"
)
```

### Output Format Conversion

```r
# PNG (raster, good for presentations)
convertSVG("output.svg", device = "png", dpi = 300)

# PDF (vector, good for publications)
convertSVG("output.svg", device = "pdf")

# TIFF (high-quality raster)
convertSVG("output.svg", device = "tiff", dpi = 600)

# JPG (compressed)
convertSVG("output.svg", device = "jpg", quality = 95)
```

## GFFex Function

Extract feature density from GFF/GTF annotation files:

```r
gene_density <- GFFex(
  input = "annotation.gff3",
  karyotype = karyotype,
  feature = "gene",      # Feature type to extract
  window = 1000000       # Window size in bp
)
```

Useful for comparing ATAC-seq accessibility with gene density or other annotations.

## Species Support

### Built-in Data

```r
data(human_karyotype)
data(gene_density)
data(Random_RNAs_500)
```

### Custom Species

Create karyotype from chrom.sizes file:

```r
# Download from UCSC
# wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes

chrom_sizes <- read.table("mm39.chrom.sizes", header = FALSE)
colnames(chrom_sizes) <- c("Chr", "End")
chrom_sizes$Start <- 0
karyotype <- chrom_sizes[, c("Chr", "Start", "End")]
```

For centromere positions, consult species-specific databases (e.g., UCSC Genome Browser).
