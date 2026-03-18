# Data Format Specifications

## Karyotype Format

### Required Columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| Chr | string | Chromosome ID | 1, 2, X, Y, chr1 |
| Start | integer | Chromosome start position | 0 |
| End | integer | Chromosome end position | 248956422 |

### Optional Columns (for centromere display)

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| CE_start | integer | Centromere start position | 122026459 |
| CE_end | integer | Centromere end position | 124932724 |

### Example Files

**3-column format (karyotype_minimal.txt):**
```
Chr	Start	End
1	0	248956422
2	0	242193529
3	0	198295559
X	0	156040895
Y	0	57227415
```

**5-column format (karyotype_full.txt):**
```
Chr	Start	End	CE_start	CE_end
1	0	248956422	122026459	124932724
2	0	242193529	92188145	94090557
3	0	198295559	90772458	93655574
```

### Loading in R

```r
# Tab-separated
karyotype <- read.table("karyotype.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Or from built-in data
data(human_karyotype)
```

## Overlaid Data Format (Heatmap)

### Required Columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| Chr | string | Chromosome ID | 1 |
| Start | integer | Window start position | 1 |
| End | integer | Window end position | 1000000 |
| Value | numeric | Feature value (density, signal) | 65 |

### Example File (peak_density.txt)

```
Chr	Start	End	Value
1	1	1000000	65
1	1000001	2000000	76
1	2000001	3000000	35
2	1	1000000	42
```

### Generating from ATAC-seq Peaks

```r
# Window-based density
generate_density <- function(peaks, karyotype, window_size = 1e6) {
  result <- data.frame()
  
  for (chr in unique(karyotype$Chr)) {
    chr_end <- karyotype[karyotype$Chr == chr, "End"]
    windows <- seq(1, chr_end, by = window_size)
    
    for (i in 1:(length(windows) - 1)) {
      w_start <- windows[i]
      w_end <- min(windows[i + 1] - 1, chr_end)
      
      count <- sum(peaks$chrom == chr & 
                   peaks$start >= w_start & 
                   peaks$end <= w_end)
      
      result <- rbind(result, data.frame(
        Chr = chr, Start = w_start, End = w_end, Value = count
      ))
    }
  }
  return(result)
}
```

## Label Data Format (Markers)

### Required Columns

| Column | Type | Description | Valid Values |
|--------|------|-------------|--------------|
| Type | string | Label type | marker, heatmap, line, polygon |
| Shape | string | Marker shape | box, triangle, circle |
| Chr | string | Chromosome ID | Must match karyotype |
| Start | integer | Feature start | 1-based |
| End | integer | Feature end | > Start |
| color | string | Hex color (no #) | ff0000, 00ff00, 0000ff |

### Example File (markers.txt)

```
Type	Shape	Chr	Start	End	color
marker	box	1	86357632	86357700	ff7f00
marker	triangle	6	69204486	69204568	6a3d9a
marker	circle	3	68882967	68883091	33a02c
```

### Color Encoding Strategies

**Binary (significant vs non-significant):**
```r
color <- ifelse(peaks$qval > 10, "ff0000", "ff7f00")
```

**Quantile-based:**
```r
colors <- c("2166ac", "fee08b", "d73027")
color_index <- cut(peaks$signal, breaks = 3, labels = FALSE)
color <- colors[color_index]
```

**Continuous gradient:**
```r
color_palette <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b73027"))(100)
color <- color_palette[cut(peaks$signal, breaks = 100)]
# Remove # for RIdeogram
color <- gsub("#", "", color)
```

## Input File Conversion Scripts

### BED to Label Format

```bash
#!/bin/bash
# convert-bed-to-label.sh

INPUT_BED=$1
OUTPUT_TXT=$2

awk -v OFS='\t' '{
  print "marker", "box", $1, $2, $3, "ff7f00"
}' "$INPUT_BED" > "$OUTPUT_TXT"

# Add header
sed -i '1i Type\tShape\tChr\tStart\tEnd\tcolor' "$OUTPUT_TXT"
```

### narrowPeak to Label Format (R)

```r
convert_narrowpeak <- function(input_file, output_file, color_by = "signal") {
  peaks <- read.table(input_file, header = FALSE)
  colnames(peaks) <- c("chrom", "start", "end", "name", "score", 
                       "strand", "signal", "pval", "qval", "peak")
  
  # Filter
  peaks <- peaks[peaks$qval > 10, ]
  
  # Color
  if (color_by == "signal") {
    colors <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b73027"))(100)
    color <- colors[cut(peaks$signal, breaks = 100)]
  } else {
    color <- ifelse(peaks$score > 100, "ff0000", "ff7f00")
  }
  
  label_data <- data.frame(
    Type = "marker",
    Shape = "box",
    Chr = peaks$chrom,
    Start = peaks$start,
    End = peaks$end,
    color = gsub("#", "", color)
  )
  
  write.table(label_data, output_file, sep = "\t", 
              row.names = FALSE, quote = FALSE)
}
```

## Validation

### Check Karyotype

```r
validate_karyotype <- function(karyotype) {
  required_cols <- c("Chr", "Start", "End")
  if (!all(required_cols %in% colnames(karyotype))) {
    stop("Missing required columns: ", 
         setdiff(required_cols, colnames(karyotype)))
  }
  
  if (any(karyotype$Start >= karyotype$End)) {
    stop("Invalid: Start >= End")
  }
  
  if (any(karyotype$Start < 0)) {
    stop("Invalid: Start < 0")
  }
  
  message("Karyotype valid: ", nrow(karyotype), " chromosomes")
}
```

### Check Label Data

```r
validate_label <- function(label_data, karyotype) {
  required_cols <- c("Type", "Shape", "Chr", "Start", "End", "color")
  if (!all(required_cols %in% colnames(label_data))) {
    stop("Missing required columns")
  }
  
  # Check shapes
  valid_shapes <- c("box", "triangle", "circle")
  if (!all(label_data$Shape %in% valid_shapes)) {
    stop("Invalid shape. Must be: ", paste(valid_shapes, collapse = ", "))
  }
  
  # Check chromosome match
  if (!all(label_data$Chr %in% karyotype$Chr)) {
    unmatched <- setdiff(label_data$Chr, karyotype$Chr)
    warning("Chromosomes not in karyotype: ", paste(unmatched, collapse = ", "))
  }
  
  message("Label data valid: ", nrow(label_data), " markers")
}
```

## Common Chromosome Naming Conventions

| Convention | Example | Usage |
|------------|---------|-------|
| UCSC | chr1, chr2, chrX | GRCh38, hg38 |
| Ensembl | 1, 2, X | GRCh38 |
| NCBI | NC_000001.11 | RefSeq |

**Important:** Ensure consistency between karyotype, peaks, and annotation files.
