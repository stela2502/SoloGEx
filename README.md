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

# or directly from github
devtools::install_github("stela2502/SoloGEx")
```

## Quck Start

```r
library(SoloGEx)


# Reference
ref_expr <- matrix(rnorm(10*4), nrow = 10)
rownames(ref_expr) = paste(sep=".", "Gene", 1:10);
ref_colData <- data.frame(SampleID = paste0("R", 1:4),
      BiologicalGroup = rep(c("A","B"), each = 2))
ref_edat <- SoloGEx(ref_expr, ref_colData)
  
# Singleton
singleton_expr <- matrix(rnorm(10*2), nrow = 10)
rownames(singleton_expr) = paste(sep=".", "Gene", 1:10);
singleton_colData <- data.frame(SampleID = paste0("S", 1:2))
singleton_edat <- SoloGEx(singleton_expr, singleton_colData)
  
singleton_edat <- analyze_singletons(singleton_edat, ref_edat, group_col="BiologicalGroup")
  
# Export the results to a file:

write_combined_file ( singleton_edat, "your_file.csv")

# Create singleton dataset
singleton_expr
```
