# ATAC-seq Peak Visualization - Complete Workflow
# From BED/narrowPeak file to publication-ready ideogram

library(RIdeogram)
library(dplyr)

# ============================================================
# STEP 1: Load R Environment
# ============================================================
# Run externally: mamba activate r_env

# ============================================================
# STEP 2: Load ATAC-seq Peak Data
# ============================================================

# Option A: From narrowPeak file (MACS2 output)
peaks <- read.table("peaks.narrowPeak", header = FALSE)
colnames(peaks) <- c("chrom", "start", "end", "name", "score", 
                     "strand", "signal", "pval", "qval", "peak")

# Option B: From BED file
# peaks <- read.table("peaks.bed", header = FALSE)
# colnames(peaks) <- c("chrom", "start", "end", "name", "score", "strand")

# ============================================================
# STEP 3: Filter Significant Peaks
# ============================================================

# Filter by q-value (qval > 10 means q < 0.1 in -log10 scale)
peaks_filtered <- peaks %>%
  filter(qval > 10) %>%
  filter(signal > 50)  # Optional: filter by signal strength

cat("Total peaks:", nrow(peaks), "\n")
cat("Filtered peaks:", nrow(peaks_filtered), "\n")

# ============================================================
# STEP 4: Prepare Karyotype
# ============================================================

# Option A: Use built-in human data
data(human_karyotype)
karyotype <- human_karyotype

# Option B: Load custom karyotype
# karyotype <- read.table("karyotype.txt", sep = "\t", header = TRUE)

# Ensure chromosome naming consistency
# Convert peaks chrom to match karyotype (remove 'chr' prefix if needed)
peaks_filtered$chrom <- gsub("chr", "", peaks_filtered$chrom)

# ============================================================
# STEP 5: Create Label Data (Marker Mode)
# ============================================================

# Color peaks by signal strength using gradient
color_palette <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b73027"))(100)
peak_color_index <- cut(peaks_filtered$signal, breaks = 100, labels = FALSE)
peak_colors <- color_palette[peak_color_index]

# Remove # for RIdeogram format
peak_colors <- gsub("#", "", peak_colors)

# Create label dataframe
label_data <- data.frame(
  Type = "marker",
  Shape = "box",  # Options: box, triangle, circle
  Chr = peaks_filtered$chrom,
  Start = peaks_filtered$start,
  End = peaks_filtered$end,
  color = peak_colors
)

# Optional: Limit to top peaks to avoid clutter
top_n <- 200
label_data_top <- label_data %>%
  arrange(desc(peaks_filtered$signal[order(row_number())])) %>%
  head(top_n)

# ============================================================
# STEP 6: Create Peak Density (Heatmap Mode)
# ============================================================

window_size <- 10000000  # 10Mb windows
peak_density <- data.frame()

for (chr in unique(karyotype$Chr)) {
  chr_info <- karyotype[karyotype$Chr == chr, ]
  chr_end <- chr_info$End
  
  # Create windows
  windows <- seq(0, chr_end, by = window_size)
  
  for (i in 1:(length(windows) - 1)) {
    w_start <- windows[i] + 1
    w_end <- windows[i + 1]
    
    # Count peaks in window
    chr_peaks <- peaks_filtered[peaks_filtered$chrom == chr, ]
    n_peaks <- sum(chr_peaks$start >= w_start & chr_peaks$end <= w_end)
    
    peak_density <- rbind(peak_density, data.frame(
      Chr = chr,
      Start = w_start,
      End = w_end,
      Value = n_peaks
    ))
  }
}

# ============================================================
# STEP 7: Generate Ideograms
# ============================================================

# 7A: Marker-only view
ideogram(
  karyotype = karyotype,
  label = label_data_top,
  label_type = "marker",
  output = "atac_peaks_markers.svg"
)

# 7B: Density heatmap view
ideogram(
  karyotype = karyotype,
  overlaid = peak_density,
  colorset1 = c("#4575b4", "#ffffbf", "#d73027"),  # Blue-White-Red
  output = "atac_density_heatmap.svg"
)

# 7C: Combined view (heatmap + markers)
ideogram(
  karyotype = karyotype,
  overlaid = peak_density,
  label = label_data_top,
  label_type = "marker",
  colorset1 = c("#4575b4", "#ffffbf", "#d73027"),
  output = "atac_combined.svg"
)

# ============================================================
# STEP 8: Convert to Publication Formats
# ============================================================

convertSVG("atac_peaks_markers.svg", device = "png", dpi = 300)
convertSVG("atac_density_heatmap.svg", device = "png", dpi = 300)
convertSVG("atac_combined.svg", device = "png", dpi = 300)

# Also create PDF for publications
convertSVG("atac_combined.svg", device = "pdf")

cat("\nGenerated files:\n")
cat("- atac_peaks_markers.svg/png\n")
cat("- atac_density_heatmap.svg/png\n")
cat("- atac_combined.svg/png/pdf\n")
