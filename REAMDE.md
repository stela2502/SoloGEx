# SoloGEx

**Solo Gene Expression Explorer** (SoloGEx) is a minimal R package for exploratory analysis 
of single-sample (“solo”) gene expression data against reference datasets.

## Features

- Wrap expression matrices and sample metadata into a `ExpressionDataSet` object
- Compute per-group statistics (mean, SD, CV) for reference groups
- Compare singleton samples to reference and calculate Z-scores
- Identify closest baseline group per gene
- Calculate overall deviation per sample to rank unusual profiles
- Export a combined table with expression, group stats, and singleton analysis
- Pure base R, no Bioconductor dependencies

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# From local folder
devtools::install("/path/to/SoloGEx")