# Basic RIdeogram Example
# Minimal working example with built-in human karyotype data

# Load library
library(RIdeogram)

# Load built-in human karyotype data
data(human_karyotype)

# View data structure
head(human_karyotype)
#   Chr Start       End  CE_start    CE_end
# 1   1     0 248956422 122026459 124932724
# 2   2     0 242193529  92188145  94090557

# Create basic ideogram (karyotype only)
ideogram(
  karyotype = human_karyotype,
  output = "basic_ideogram.svg"
)

# Convert to PNG (300 DPI for publication)
convertSVG("basic_ideogram.svg", device = "png", dpi = 300)

# Result: basic_ideogram.svg and basic_ideogram.png in working directory
