#!/usr/bin/env Rscript
# validate-input.R
# Validate input data formats before RIdeogram plotting

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: validate-input.R <file1.txt> [file2.txt ...]\n")
  cat("\nValidates:\n")
  cat("  - Karyotype files (3 or 5 columns)\n")
  cat("  - Overlaid/heatmap files (4 columns)\n")
  cat("  - Label/marker files (6 columns)\n")
  quit(status = 0)
}

validate_karyotype <- function(df, filename) {
  errors <- character()
  warnings <- character()
  
  # Check columns
  required <- c("Chr", "Start", "End")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing columns:", paste(missing, collapse = ", ")))
  }
  
  # Check Start < End
  if ("Start" %in% colnames(df) && "End" %in% colnames(df)) {
    invalid <- which(df$Start >= df$End)
    if (length(invalid) > 0) {
      errors <- c(errors, paste("Rows with Start >= End:", length(invalid)))
    }
  }
  
  # Check Start >= 0
  if ("Start" %in% colnames(df)) {
    negative <- which(df$Start < 0)
    if (length(negative) > 0) {
      errors <- c(errors, paste("Rows with Start < 0:", length(negative)))
    }
  }
  
  # Check centromere if present
  if (all(c("CE_start", "CE_end") %in% colnames(df))) {
    invalid_ce <- which(df$CE_start >= df$CE_end)
    if (length(invalid_ce) > 0) {
      warnings <- c(warnings, paste("Centromere issues (CE_start >= CE_end):", length(invalid_ce)))
    }
  }
  
  # Report
  cat("\n=== Karyotype:", filename, "===\n")
  cat("Rows:", nrow(df), "\n")
  cat("Columns:", paste(colnames(df), collapse = ", "), "\n")
  
  if (length(errors) > 0) {
    cat("ERRORS:\n")
    for (e in errors) cat("  -", e, "\n")
    return(FALSE)
  }
  
  if (length(warnings) > 0) {
    cat("WARNINGS:\n")
    for (w in warnings) cat("  -", w, "\n")
  }
  
  cat("STATUS: VALID\n")
  return(TRUE)
}

validate_overlaid <- function(df, filename) {
  errors <- character()
  
  required <- c("Chr", "Start", "End", "Value")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing columns:", paste(missing, collapse = ", ")))
  }
  
  if ("Value" %in% colnames(df)) {
    if (!is.numeric(df$Value)) {
      errors <- c(errors, "Value column must be numeric")
    }
    if (any(df$Value < 0, na.rm = TRUE)) {
      warnings <- c(warnings, "Negative values found")
    }
  }
  
  cat("\n=== Overlaid/Heatmap:", filename, "===\n")
  cat("Rows:", nrow(df), "\n")
  cat("Value range:", min(df$Value, na.rm = TRUE), "-", max(df$Value, na.rm = TRUE), "\n")
  
  if (length(errors) > 0) {
    cat("ERRORS:\n")
    for (e in errors) cat("  -", e, "\n")
    return(FALSE)
  }
  
  cat("STATUS: VALID\n")
  return(TRUE)
}

validate_label <- function(df, filename) {
  errors <- character()
  warnings <- character()
  
  required <- c("Type", "Shape", "Chr", "Start", "End", "color")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    errors <- c(errors, paste("Missing columns:", paste(missing, collapse = ", ")))
  }
  
  # Check shapes
  if ("Shape" %in% colnames(df)) {
    valid_shapes <- c("box", "triangle", "circle")
    invalid <- setdiff(unique(df$Shape), valid_shapes)
    if (length(invalid) > 0) {
      errors <- c(errors, paste("Invalid shapes:", paste(invalid, collapse = ", "), 
                                "- valid: box, triangle, circle"))
    }
  }
  
  # Check color format (should be hex without #)
  if ("color" %in% colnames(df)) {
    bad_colors <- grep("^[^0-9a-fA-F]", df$color, value = TRUE)
    if (length(bad_colors) > 0) {
      warnings <- c(warnings, paste("Colors may have invalid format (should be hex without #):",
                                    paste(head(bad_colors, 3), collapse = ", ")))
    }
  }
  
  cat("\n=== Label/Markers:", filename, "===\n")
  cat("Rows:", nrow(df), "\n")
  cat("Shapes:", paste(unique(df$Shape), collapse = ", "), "\n")
  cat("Types:", paste(unique(df$Type), collapse = ", "), "\n")
  
  if (length(errors) > 0) {
    cat("ERRORS:\n")
    for (e in errors) cat("  -", e, "\n")
    return(FALSE)
  }
  
  if (length(warnings) > 0) {
    cat("WARNINGS:\n")
    for (w in warnings) cat("  -", w, "\n")
  }
  
  cat("STATUS: VALID\n")
  return(TRUE)
}

# Main
cat("=== RIdeogram Input Validator ===\n")

all_valid <- TRUE

for (file in args) {
  if (!file.exists(file)) {
    cat("\nERROR: File not found:", file, "\n")
    all_valid <- FALSE
    next
  }
  
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  ncols <- ncol(df)
  
  if (ncols == 3 || ncols == 5) {
    valid <- validate_karyotype(df, file)
  } else if (ncols == 4) {
    valid <- validate_overlaid(df, file)
  } else if (ncols == 6) {
    valid <- validate_label(df, file)
  } else {
    cat("\n=== Unknown format:", file, "===\n")
    cat("Columns:", ncols, "\n")
    cat("Expected: 3/5 (karyotype), 4 (overlaid), or 6 (label)\n")
    valid <- FALSE
  }
  
  all_valid <- all_valid && valid
}

cat("\n=== Summary ===\n")
if (all_valid) {
  cat("All files: VALID\n")
  quit(status = 0)
} else {
  cat("Some files: INVALID\n")
  quit(status = 1)
}
