# proteomicsUtils

A collection of utility functions for proteomics data analysis.

## Description

This package provides a comprehensive set of functions for proteomics data analysis, including:

- Data import and preprocessing
- Normalization and filtering
- Statistical testing (t-tests, Mann-Whitney, ANOVA)
- Data visualization (volcano plots, heatmaps, PCA plots, barplots)
- Gene ontology analysis
- Machine learning utilities (ROC curves, SVM)
- Missing value imputation

## Installation

### From GitHub (Recommended)

```r
# Install remotes if you don't have it
install.packages("remotes")
remotes::install_github("https://github.com/Elquimico09/proteomicsUtils")
```

Or using devtools:

```r
# Install devtools if you don't have it
install.packages("devtools")
devtools::install_github("https://github.com/Elquimico09/proteomicsUtils")
```

The package will automatically install most dependencies. However, you may need to manually install Bioconductor packages first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "ComplexHeatmap"))
```

### From Local Source

```r
install.packages("proteomicsUtils_0.1.0.tar.gz", repos = NULL, type = "source")
```

Or from terminal:

```bash
R CMD INSTALL proteomicsUtils_0.1.0.tar.gz
```

## Usage

### Example Workflow

```r
library(proteomicsUtils)

# Import data
data <- importData("proteomics_data.csv")

# Filter and normalize
filtered <- filterMissingValues(data, threshold = 0.5)
normalized <- normalizeData(filtered) # log2quantile normalization
relative_normalized <- normalizeData(filtered, method = "relative")

# visualize distribution
visualizeDist(filtered)

# Visualization
volcano_plot <- makeVolcano(results, fc_cutoff = 1)
pca_results <- performPCA(imputed)
pca_plot <- plotPCA(pca_results)
```

### Example Workflow Using Pipe

```r
library(proteomicsUtils)

# Import filter and normalize all at once
normalized <- importData("proteomics_data.csv") %>%
    filterMissingValues(threshold = 0.5) %>%
    filterKeratin() %>%
    # If serum && depleted
    # filterHighAbundance() %>%
    normalizeData()

# pca
pca_plot <- normalized %>%
    performPCA() %>%
    plotPCA()
```


### Example Functions

- `importData()` - Import proteomics data from CSV or Excel files
- `normalizeData()` - Normalize your data
- `filterMissingValues()` - Filter rows with too many missing values
- `imputeValues()` - Impute missing values
- `performTTest()` - Perform t-tests
- `performANOVA4()` - Perform ANOVA with 4 cohorts
- `makeVolcano()` - Create volcano plots
- `performPCA()` - Perform PCA analysis
- `plotPCA()` - Create PCA plots
- `makeClustermap()` - Create heatmaps
- `geneOntology()` - Perform gene ontology analysis
- `calculateROC()` - Calculate and plot ROC curves
- `create_grouped_boxplot()` - Create grouped boxplots

## Available Functions

The package includes 31 functions:

1. calculateROC
2. calculateSVM
3. convertFormat
4. create_grouped_boxplot
5. extractGeneName
6. extractGeneNames
7. filterHighAbundance
8. filterKeratin
9. filterMissingValues
10. geneOntology
11. importData
12. imputeMinProb
13. makeBarplot
14. makeBarplot_simp
15. makeBarplot3
16. makeClustermap
17. makeVolcano
18. normalizeData
19. perform_rowwise_anova_tukey_fc
20. perform2WayANOVA
21. performANOVA
22. performANOVA6
23. performGO
24. performKW
25. performMWTest
26. performPCA
27. performTTest
28. plotPCA
29. scale_matrix_rows
30. testNormality
31. visualizeDist

## License

MIT License

## Author

Vishal Sandilya  
<vishal.sandilya@ttu.edu>
