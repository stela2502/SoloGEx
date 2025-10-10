# SoloGEx

**Solo Gene Expression Explorer** (SoloGEx) is a minimal R package for exploratory analysis 
of single-sample (“solo”) gene expression data against reference datasets.

## Features

- Wrap expression matrices and sample metadata into a `SoloGEx` object (expression and colData)
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

## For a real dataset you should nortmalize the data!
ref_edat <- SoloGEx(ref_expr, ref_colData)
  
# Singleton
singleton_expr <- matrix(rnorm(10*2), nrow = 10)
rownames(singleton_expr) = paste(sep=".", "Gene", 1:10);
singleton_colData <- data.frame(SampleID = paste0("S", 1:2))

## For a real dataset you should nortmalize the data!
singleton_edat <- SoloGEx(singleton_expr, singleton_colData)

# One stop create all the values in one go  
singleton_edat <- analyze_singletons(singleton_edat, ref_edat, group_col="BiologicalGroup")
  
# Export the results + logFC to a file:
write_combined_file ( singleton_edat, "your_file.csv")

singleton_edat = compute_corr_to_reference( singleton_edat, ref_edat )

# write the correlation heatmap to SVG
svg("correlation_heatmap.svg", width = 6, height = 6)  # open SVG device
plot_corr(singleton_edat)                      # generate the plot
dev.off()                                              # close device and save file

# Create singleton dataset
singleton_expr
```

## 🧠 What Happens Inside `analyze_singletons()`

The `analyze_singletons()` function is the analytical core of SoloGEx.  
It takes one or more **singleton samples** (for example, individual patients or experimental outliers) and compares their expression profiles to a set of **reference groups** (e.g., healthy controls, treated vs untreated conditions, or known tissue types).

The goal is to understand:
- How similar each singleton is to the existing groups,  
- Which group each gene behaves most like,  
- And how overall distinct the sample appears.

---

### 🔍 Step 1: Define the Reference Landscape

Before any comparisons are made, SoloGEx establishes the **expression characteristics** of the reference groups.

For every gene, it summarizes:
- The **average expression level** (group mean)  
- The **variation within the group** (standard deviation, SD)  
- Optionally, the **coefficient of variation** (CV = SD / mean)

This creates a statistical fingerprint of what “normal” expression looks like for each group.  

These fingerprints are used to judge how far a singleton’s expression differs from expected group behavior.

---

### ⚖️ Step 2: Measure Gene-Level Deviations

Each gene in a singleton sample is then compared to every reference group.  
SoloGEx asks:  
> “If this sample belonged to Group A, how unusual would this gene’s expression be?”

This is captured by a **Z-score**, which expresses how many standard deviations a gene’s expression lies above or below the group average.  

- A **Z-score near 0** means the gene behaves like that group’s typical pattern.  
- A **high absolute Z-score (|Z|)** means the gene is strongly over- or under-expressed relative to that group.

The process repeats for each reference group, producing a map of possible similarities for every gene.

---

### 🧭 Step 3: Identify the Most Similar Group per Gene

For each gene, SoloGEx determines which reference group shows the **smallest deviation** — that is, where the expression looks most at home.

The result is:
- A **closest group** label (e.g., “Healthy”, “Tumor_A”, “Control”)  
- And a **minimal Z-score**, describing how well that gene fits in there.

This is biologically informative:  
- Some genes may clearly align with a specific group,  
- Others may deviate from all known patterns — suggesting deregulation or an unknown state.

### ⚠️ Important Note About “Closest Group”

If the Z-scores across all reference groups are very similar for a gene, the assigned closest group is **essentially arbitrary**.

* In this case, do not over-interpret the closest group.
* Focus on genes with clear deviations and patterns across multiple genes.

---

### 🌐 Step 4: Summarize Sample-Level Similarity

Once each gene has found its “closest home,” SoloGEx looks across all genes in the singleton sample and asks:

> “On average, how far does this sample sit from the nearest reference patterns?”

It calculates the mean of all minimal deviations — a single metric called **overall deviation**.  

- A **low overall deviation** means the sample globally resembles one of the reference groups.  
- A **high deviation** suggests a transcriptomic shift — such as treatment response, disease progression, or a unique patient-specific phenotype.

This gives a compact summary of how typical or exceptional a sample appears in its biological context.

---

### 📦 Step 5: Store and Organize the Results

The function stores all findings within the SoloGEx object, under `@singleton_analysis`:

~~~
singleton_analysis:
  ├─ Zscores        → per-gene deviations to all groups
  ├─ closest_group  → which group each gene resembles
  ├─ min_Zscore     → how far it still deviates
  └─ overall_dev    → sample-wide average deviation
~~~

These layers allow you to interpret results from different perspectives:
- **Gene-level:** which pathways or markers behave unexpectedly  
- **Group-level:** which biological condition each sample most resembles  
- **Sample-level:** how extreme or normal the sample is overall

---

### 🧬 Step 6: Export for Interpretation

Once analysis is complete, you can summarize and share results using:

~~~
write_combined_file(single_edat, "SoloGEx_results.csv")
~~~

This export merges:
- The original expression values,  
- Reference group statistics, and  
- The computed Z-scores, closest groups, and deviation summaries.

The resulting table can be directly used for:
- Heatmaps or cluster visualization,  
- Pathway enrichment of deregulated genes,  
- Or sample classification analyses.

---

### 🧩 Intuitive Summary

| Step | Concept | What It Reveals |
|------|----------|-----------------|
| **1. Reference statistics** | Define “normal” gene behavior in each group | Establishes baseline variation |
| **2. Gene-level comparison** | Compare singleton gene expression vs each group | Detects specific deregulations |
| **3. Closest group mapping** | Identify best-fitting group per gene | Highlights partial similarity patterns |
| **4. Sample-level deviation** | Summarize global transcriptomic distance | Measures overall normality vs divergence |
| **5. Export** | Combine all into interpretable table | Ready for downstream biology |

---

### 💡 Example Biological Interpretation

Imagine comparing a **new tumor biopsy** against a panel of known tumor types:

- Most genes align with the *melanoma* profile → **closest_group = Melanoma**  
- Some immune-related genes match the *lymphoma* pattern → mixed signature  
- Overall deviation is **high**, suggesting an atypical or hybrid state  

This result immediately guides hypotheses:
- Tumor might show lineage plasticity  
- Immune infiltration could be driving mixed signatures  
- Or the sample could represent a new subtype worth further validation


