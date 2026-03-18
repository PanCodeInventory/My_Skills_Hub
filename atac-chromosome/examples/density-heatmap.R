# Peak Density Heatmap Example
# Calculate and visualize ATAC-seq peak density across genome

library(RIdeogram)
library(dplyr)

# ============================================================
# Configuration
# ============================================================
window_size <- 5000000  # 5Mb windows
input_peaks <- "peaks.narrowPeak"
output_file <- "density_heatmap.svg"

# ============================================================
# Load Data
# ============================================================
peaks <- read.table(input_peaks, header = FALSE)
colnames(peaks) <- c("chrom", "start", "end", "name", "score", 
                     "strand", "signal", "pval", "qval", "peak")

# Filter significant peaks
peaks <- peaks[peaks$qval > 10, ]

# Clean chromosome names
peaks$chrom <- gsub("chr", "", peaks$chrom)

# Load karyotype
data(human_karyotype)

# ============================================================
# Density Calculation Function
# ============================================================
calculate_density <- function(peaks_df, karyotype, window_size) {
  result <- data.frame()
  
  for (i in 1:nrow(karyotype)) {
    chr <- karyotype$Chr[i]
    chr_end <- karyotype$End[i]
    
    # Create windows
    starts <- seq(0, chr_end, by = window_size)
    
    for (j in 1:(length(starts) - 1)) {
      w_start <- starts[j] + 1
      w_end <- starts[j + 1]
      
      # Count peaks in window
      chr_peaks <- peaks_df[peaks_df$chrom == chr, ]
      n_peaks <- sum(chr_peaks$start >= w_start & chr_peaks$end <= w_end)
      
      result <- rbind(result, data.frame(
        Chr = chr,
        Start = w_start,
        End = w_end,
        Value = n_peaks
      ))
    }
  }
  
  return(result)
}

# ============================================================
# Calculate Density
# ============================================================
cat("Calculating peak density with", window_size / 1e6, "Mb windows...\n")
peak_density <- calculate_density(peaks, human_karyotype, window_size)

# Summary statistics
cat("\nDensity summary:\n")
print(summary(peak_density$Value))

# ============================================================
# Visualize
# ============================================================

# Option 1: Blue-White-Red (standard)
ideogram(
  karyotype = human_karyotype,
  overlaid = peak_density,
  colorset1 = c("#4575b4", "#ffffbf", "#d73027"),
  width = 170,
  Lx = 160,
  Ly = 35,
  output = output_file
)

# Option 2: Green-Yellow-Red
# ideogram(
#   karyotype = human_karyotype,
#   overlaid = peak_density,
#   colorset1 = c("#1a9850", "#fee08b", "#d73027"),
#   output = "density_green.svg"
# )

# Option 3: Viridis-like
# ideogram(
#   karyotype = human_karyotype,
#   overlaid = peak_density,
#   colorset1 = c("#440154", "#31688e", "#35b779", "#fde725"),
#   output = "density_viridis.svg"
# )

# Convert to PNG
convertSVG(output_file, device = "png", dpi = 300)

cat("\nGenerated:", output_file, "and PNG version\n")
