#!/usr/bin/env Rscript
# calculate-density.R
# Calculate ATAC-seq peak density in configurable windows

args <- commandArgs(trailingOnly = TRUE)

# Default parameters
input_file <- "peaks.narrowPeak"
output_file <- "peak_density.txt"
window_size <- 1000000  # 1Mb default
karyotype_file <- NULL

# Parse arguments
if (length(args) >= 1) input_file <- args[1]
if (length(args) >= 2) output_file <- args[2]
if (length(args) >= 3) window_size <- as.numeric(args[3])
if (length(args) >= 4) karyotype_file <- args[4]

cat("=== Peak Density Calculator ===\n")
cat("Input:", input_file, "\n")
cat("Output:", output_file, "\n")
cat("Window size:", window_size / 1e6, "Mb\n\n")

# Load required library
suppressPackageStartupMessages(library(dplyr))

# Load peaks
cat("Loading peaks...\n")
peaks <- read.table(input_file, header = FALSE)

if (ncol(peaks) >= 10) {
  # narrowPeak format
  colnames(peaks) <- c("chrom", "start", "end", "name", "score", 
                       "strand", "signal", "pval", "qval", "peak")
  cat("Format: narrowPeak\n")
} else if (ncol(peaks) >= 4) {
  # BED format
  colnames(peaks)[1:4] <- c("chrom", "start", "end", "name")
  peaks$signal <- ifelse(ncol(peaks) >= 5, peaks[,5], 1)
  cat("Format: BED\n")
} else {
  stop("Unrecognized peak format")
}

# Filter significant peaks
if ("qval" %in% colnames(peaks)) {
  n_before <- nrow(peaks)
  peaks <- peaks[peaks$qval > 10, ]
  cat("Filtered:", n_before, "->", nrow(peaks), "peaks (qval > 10)\n")
}

# Clean chromosome names
peaks$chrom <- gsub("chr", "", peaks$chrom)

# Load karyotype
if (is.null(karyotype_file)) {
  cat("Using built-in human karyotype...\n")
  data(human_karyotype, package = "RIdeogram")
  karyotype <- human_karyotype
} else {
  cat("Loading karyotype from", karyotype_file, "\n")
  karyotype <- read.table(karyotype_file, header = TRUE, stringsAsFactors = FALSE)
}

# Calculate density
cat("\nCalculating density...\n")
peak_density <- data.frame()

for (chr in unique(karyotype$Chr)) {
  chr_info <- karyotype[karyotype$Chr == chr, ]
  chr_end <- chr_info$End
  
  windows <- seq(0, chr_end, by = window_size)
  
  for (i in 1:(length(windows) - 1)) {
    w_start <- windows[i] + 1
    w_end <- windows[i + 1]
    
    chr_peaks <- peaks[peaks$chrom == chr, ]
    n_peaks <- sum(chr_peaks$start >= w_start & chr_peaks$end <= w_end)
    
    peak_density <- rbind(peak_density, data.frame(
      Chr = chr,
      Start = w_start,
      End = w_end,
      Value = n_peaks
    ))
  }
  
  # Progress indicator
  cat("  ", chr, ": ", length(windows) - 1, " windows\n", sep = "")
}

# Summary
cat("\nDensity summary:\n")
print(summary(peak_density$Value))
cat("\nWindows with peaks:", sum(peak_density$Value > 0), "/", nrow(peak_density), "\n")

# Save output
write.table(peak_density, output_file, sep = "\t", 
            row.names = FALSE, quote = FALSE)

cat("\nSaved to:", output_file, "\n")
cat("Ready for use with RIdeogram overlaid parameter\n")
