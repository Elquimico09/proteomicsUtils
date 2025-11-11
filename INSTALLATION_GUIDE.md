# Installation Guide for proteomicsUtils

## What's Been Done

Your `utilities.R` file has been successfully converted into an R package called `proteomicsUtils`. The package is now structured properly with:

- **Package directory**: `proteomicsUtils/`
- **All 29 functions** from your utilities.R file
- **DESCRIPTION file** with metadata and dependencies
- **NAMESPACE file** that exports all functions
- **Built package**: `proteomicsUtils_0.1.0.tar.gz`

## Installation Options

### Option 1: Full Installation (Recommended)

This installs all dependencies and then the package:

```r
# 1. Install CRAN dependencies
install.packages(c("pROC", "e1071", "truncnorm", "ggplot2", "dplyr",
                   "tidyr", "readr", "readxl", "stringr", "ggprism",
                   "ggrepel", "ggpubr", "matrixStats", "preprocessCore"))

# 2. Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "ComplexHeatmap"))

# 3. Install the package
install.packages("proteomicsUtils_0.1.0.tar.gz", repos = NULL, type = "source")
```

### Option 2: Quick Source Installation (Alternative)

If you want to use the functions immediately without waiting for all dependencies, you can still source the file as before, but now it's organized in the package:

```r
source("proteomicsUtils/R/utilities.R")
```

### Option 3: Manual R CMD INSTALL (From Terminal)

```bash
# From the ProteomicsDP directory
R CMD INSTALL proteomicsUtils_0.1.0.tar.gz
```

Note: This requires all dependencies to be installed first (same as Option 1, steps 1-2).

## Using the Package

Once installed, simply load it in any R session:

```r
library(proteomicsUtils)

# Now all your functions are available!
data <- importData("mydata.csv")
filtered_data <- filterMissingValues(data, 0.5)
volcano <- makeVolcano(results)
```

## Troubleshooting

### Missing Dependencies

If you get an error about missing dependencies:

1. Check which package is missing from the error message
2. Install it using:
   - For CRAN packages: `install.packages("package_name")`
   - For Bioconductor packages: `BiocManager::install("package_name")`

### clusterProfiler Installation Issues

The `clusterProfiler` package is from Bioconductor and can take several minutes to install. If it fails:

```r
# Try installing with dependencies
BiocManager::install("clusterProfiler", dependencies = TRUE)
```

## What This Means for Your Workflow

**Before**: You had to copy/paste or source `utilities.R` in every new project

**Now**: You can simply add `library(proteomicsUtils)` to any R script and all 29 functions are immediately available!

## Next Steps

1. Wait for the dependency installation to complete (currently running in background)
2. Install the package using one of the options above
3. Test it in a new R session with `library(proteomicsUtils)`
4. Delete or archive your old utilities.R sourcing code from your projects

## Package Contents

All 29 functions from your utilities.R are included:

- makeVolcano_deprecated, makeVolcano
- scale_matrix_rows
- extractGeneName
- geneOntology, performGO
- performTTest, performMWTest, performANOVA4, performANOVA6
- importData, categorizeData, compareCohorts
- performPCA, plotPCA
- makeClustermap
- filterHighAbundance, filterKeratin, filterMissingValues
- normalizeData
- makeBarplot, makeBarplot3, makeBarplot_simp
- imputeValues, impute_minprob
- calculateROC, calculateSVM
- create_grouped_boxplot
- perform_rowwise_anova_tukey_fc
