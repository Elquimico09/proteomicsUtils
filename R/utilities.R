### Scale Rows Function 

#' Scale Matrix Rows
#' @title Scale Matrix Rows
#' @description
#' Scales rows of a matrix by centering and/or scaling. Similar to base::scale()
#' but operates on rows instead of columns.
#' @param x A numeric matrix
#' @param center Logical indicating whether to center the rows (default TRUE)
#' @param scale Logical indicating whether to scale the rows (default TRUE)
#' @param add_attr Logical indicating whether to add attributes for center and scale (default TRUE)
#' @param rows Optional vector of row indices to subset
#' @param cols Optional vector of column indices to subset
#' @return A matrix with scaled rows
#' @export
#' @examples
#' # mat <- matrix(rnorm(100), nrow=10)
#' # scaled_mat <- scale_matrix_rows(mat)
scale_matrix_rows <- function(x,
                             center = TRUE,
                             scale = TRUE,
                             add_attr = TRUE,
                             rows = NULL,
                             cols = NULL) {
  
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  
  # Get the row means
  cm <- rowMeans(x, na.rm = TRUE)
  
  # Get the row sd
  if (scale) {
    csd <- matrixStats::rowSds(x, center = cm, na.rm = TRUE)
  } else {
    # just divide by 1 if not
    csd <- rep(1, length = length(cm))
  }
  
  if (!center) {
    # just subtract 0
    cm <- rep(0, length = length(cm))
  }
  
  x <- (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  
  return(x)
}

# Function: performGO
# Takes a dataframe with a gene column and performs gene ontology on the
# specified ontology
#
# Args:
#   - df: dataframe containing a column labelled gene with the hgnc gene symbol
#   - ontology = "BP", other options are "CC" for cellular component and
#                "MF" for molecular function
# Returns:
#   - dataframe containing GO results
#' Perform Gene Ontology Analysis
#' @title Gene Ontology Analysis
#' @description
#' Performs gene ontology enrichment analysis on a set of genes using clusterProfiler.
#' @param df A dataframe containing a 'Gene' column with gene symbols
#' @param ontology The ontology to use: "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
#' @param database The organism database to use (default: org.Rn.eg.db for rat)
#' @param keytype The type of gene identifier (default: "SYMBOL")
#' @return A dataframe containing GO enrichment results
#' @export
#' @examples
#' # df <- data.frame(Gene = c("ACTB", "GAPDH", "TP53"))
#' # go_results <- geneOntology(df, ontology = "BP")
geneOntology <- function(df, ontology = "BP", database = org.Rn.eg.db, keytype = "SYMBOL") {
  genes <- df$Gene
  go_analysis <- as.data.frame(clusterProfiler::enrichGO(gene = genes, OrgDb = database,
                                        keyType = keytype, ont = ontology))
  return(go_analysis)
}

# Function: performGO
# Takes a list of dataframes and performs gene ontology on each dataframe and returns
# a list of gene ontology results
#
# Args:
#   - data: list of dataframes containing each dataframe to perform GO Analysis on
#
# Returns:
#   - goresult: list containing dataframes with the go_results
#' Perform Gene Ontology on Multiple Dataframes
#' @title Perform GO Analysis on List
#' @description
#' Takes a list of dataframes and performs gene ontology analysis on each,
#' across all three ontologies (BP, MF, CC).
#' @param data A list of dataframes, each containing a 'Gene' column
#' @return A list of lists containing GO results for each dataframe and ontology
#' @export
#' @examples
#' # df1 <- data.frame(Gene = c("ACTB", "GAPDH"))
#' # df2 <- data.frame(Gene = c("TP53", "MYC"))
#' # results <- performGO(list(df1, df2))
performGO <- function(data) {
  ontologies <- c("BP", "MF", "CC")
  goresult <- list() # Empty list to store the final results
  for (i in seq_along(data)) {
    gene_ontology_data <- data[[i]]
    gene_ont <- list() # Empty list to store result for each go before merging with goresult
    n <- 1
    for (onto in ontologies) {
      gene_ont[[n]] <- geneOntology(gene_ontology_data, ontology = onto, database = org.Rn.eg.db)
      n <- n + 1    # move on to the next ontology
    }
    goresult[[i]] <- gene_ont
    }
  return(goresult)
}


#' Perform T-Test on Two Cohorts
#'
#' Takes two dataframes and performs a t-test on each row of the dataframes.
#' Computes p-values, fold changes, and log2 fold changes with optional p-value adjustment.
#'
#' @param df1 Dataframe containing the first cohort. Must have the same number of rows as df2.
#' @param df2 Dataframe containing the second cohort. Must have the same number of rows as df1.
#' @param adjust Logical, whether to adjust p-values for multiple testing (default: FALSE).
#' @param p_adjust Method for p-value adjustment if adjust = TRUE. Options: "BH" (default),
#'   "bonferroni", "holm", "hochberg", "hommel", "BY", "none".
#' @param log_scale Logical, whether data is already log-transformed (default: TRUE).
#'   If TRUE: log2fc = mean(df2) - mean(df1) (difference of log values).
#'   If FALSE: foldchange = mean(df2) / mean(df1), log2fc = log2(foldchange).
#'
#' @return A dataframe containing the combined original data plus computed columns:
#'   - pvalue: raw p-value from t-test
#'   - adjusted_pvalue: adjusted p-value (if adjust = TRUE)
#'   - log2fc: log2 fold change
#'   - foldchange: fold change (if log_scale = FALSE)
#'
#' @export
#'
#' @examples
#' # Basic t-test without adjustment
#' result <- performTTest(control_data, treatment_data)
#'
#' # T-test with BH adjustment
#' result <- performTTest(control_data, treatment_data, adjust = TRUE)
#'
#' # T-test on non-log data
#' result <- performTTest(control_data, treatment_data, log_scale = FALSE)
performTTest <- function(df1, df2, adjust = FALSE, p_adjust = "BH", log_scale = TRUE) {
  # Ensure dataframes have the same number of rows
  if (nrow(df1) != nrow(df2)) {
    stop("df1 and df2 must have the same number of rows")
  }

  # Warn if adjust = TRUE but p_adjust not explicitly specified
  if (adjust && missing(p_adjust)) {
    warning("Adjustment method not specified, defaulting to BH")
  }

  # Print log scale status
  if (log_scale) {
    message("log_scale = TRUE: Data assumed to be log-transformed. Fold change will NOT be calculated.")
  } else {
    message("log_scale = FALSE: Data assumed to be linear scale. Fold change WILL be calculated.")
  }

  # Initialize result vectors
  pvalues <- numeric(nrow(df1))
  foldchange <- numeric(nrow(df1))
  log2fc <- numeric(nrow(df1))

  for (i in 1:nrow(df1)) {
    # Extract values for row i (exclude NAs)
    values_df1 <- unlist(df1[i, ], use.names = FALSE)
    values_df1 <- values_df1[!is.na(values_df1)]

    values_df2 <- unlist(df2[i, ], use.names = FALSE)
    values_df2 <- values_df2[!is.na(values_df2)]

    # Check valid sample size for t-test
    if (length(values_df1) >= 3 && length(values_df2) >= 3) {
      # Compute t-test
      test_result <- tryCatch(
        t.test(values_df1, values_df2),
        error = function(e) list(p.value = 1)
      )

      pvalues[i] <- test_result$p.value

      mean1 <- mean(values_df1, na.rm = TRUE)
      mean2 <- mean(values_df2, na.rm = TRUE)

      if (log_scale) {
        # Data is log-transformed: difference of means = log2 fold change
        log2fc[i] <- mean2 - mean1
      } else {
        # Data is linear scale: calculate ratio then take log2
        if (mean1 != 0) {
          foldchange[i] <- mean2 / mean1
          log2fc[i] <- log2(foldchange[i])
        } else {
          foldchange[i] <- NA
          log2fc[i] <- NA
        }
      }
    } else {
      pvalues[i] <- 1
      foldchange[i] <- ifelse(log_scale, NA, 0)
      log2fc[i] <- 0
    }
  }

  # Adjust p-values if requested
  if (adjust) {
    pvalues <- p.adjust(pvalues, method = p_adjust)
  }

  nlog10p <- -log10(pvalues)

  # Combine results with original data
  c_df <- cbind(df1, df2)
  c_df$pvalue <- pvalues
  if (!log_scale) {
    c_df$foldchange <- foldchange
  }
  c_df$nlog10p <- nlog10p
  c_df$log2fc <- log2fc

  return(c_df)
}


#' Perform Mann-Whitney U Test on Two Cohorts
#'
#' Takes two dataframes and performs a non-parametric Mann-Whitney U test on each row.
#' This is the non-parametric alternative to the t-test for comparing two groups.
#' Computes p-values, fold changes, and log2 fold changes with optional p-value adjustment.
#'
#' @param df1 Dataframe containing the first cohort. Must have the same number of rows as df2.
#' @param df2 Dataframe containing the second cohort. Must have the same number of rows as df1.
#' @param adjust Logical, whether to adjust p-values for multiple testing (default: FALSE).
#' @param p_adjust Method for p-value adjustment if adjust = TRUE. Options: "BH" (default),
#'   "bonferroni", "holm", "hochberg", "hommel", "BY", "none".
#' @param log_scale Logical, whether data is already log-transformed (default: TRUE).
#'   If TRUE: log2fc = mean(df2) - mean(df1) (difference of log values).
#'   If FALSE: foldchange = mean(df2) / mean(df1), log2fc = log2(foldchange).
#'
#' @return A dataframe containing the combined original data plus computed columns:
#'   - pvalue: raw p-value from Mann-Whitney U test
#'   - adjusted_pvalue: adjusted p-value (if adjust = TRUE)
#'   - log2fc: log2 fold change
#'   - foldchange: fold change (if log_scale = FALSE)
#'
#' @export
#'
#' @examples
#' # Basic Mann-Whitney test without adjustment
#' result <- performMWTest(control_data, treatment_data)
#'
#' # Mann-Whitney test with BH adjustment
#' result <- performMWTest(control_data, treatment_data, adjust = TRUE)
#'
#' # Test on non-log data
#' result <- performMWTest(control_data, treatment_data, log_scale = FALSE)
performMWTest <- function(df1, df2, adjust = FALSE, p_adjust = "BH", log_scale = TRUE) {
  # Ensure dataframes have the same number of rows
  if (nrow(df1) != nrow(df2)) {
    stop("df1 and df2 must have the same number of rows")
  }

  # Warn if adjust = TRUE but p_adjust not explicitly specified
  if (adjust && missing(p_adjust)) {
    warning("Adjustment method not specified, defaulting to BH")
  }

  # Print log scale status
  if (log_scale) {
    message("log_scale = TRUE: Data assumed to be log-transformed. Fold change will NOT be calculated.")
  } else {
    message("log_scale = FALSE: Data assumed to be linear scale. Fold change WILL be calculated.")
  }

  # Initialize result vectors
  pvalues <- numeric(nrow(df1))
  foldchange <- numeric(nrow(df1))
  log2fc <- numeric(nrow(df1))

  for (i in 1:nrow(df1)) {
    values_df1 <- as.numeric(df1[i, ])
    values_df1 <- values_df1[!is.na(values_df1)]

    values_df2 <- as.numeric(df2[i, ])
    values_df2 <- values_df2[!is.na(values_df2)]

    # Check if both rows have at least 3 non-NA values
    if (length(values_df1) >= 3 && length(values_df2) >= 3) {
      test_result <- tryCatch(
        wilcox.test(values_df1, values_df2),
        error = function(e) list(p.value = 1)
      )

      pvalues[i] <- test_result$p.value

      mean1 <- mean(values_df1, na.rm = TRUE)
      mean2 <- mean(values_df2, na.rm = TRUE)

      if (log_scale) {
        # Data is log-transformed: difference of means = log2 fold change
        log2fc[i] <- mean2 - mean1
      } else {
        # Data is linear scale: calculate ratio then take log2
        if (mean1 != 0) {
          foldchange[i] <- mean2 / mean1
          log2fc[i] <- log2(foldchange[i])
        } else {
          foldchange[i] <- NA
          log2fc[i] <- NA
        }
      }
    } else {
      pvalues[i] <- 1
      foldchange[i] <- ifelse(log_scale, NA, 0)
      log2fc[i] <- 0
    }
  }

  # Adjust p-values if requested
  if (adjust) {
    pvalues <- p.adjust(pvalues, method = p_adjust)
  }

  nlog10p <- -log10(pvalues)

  # Combine results with original data
  c_df <- cbind(df1, df2)
  c_df$pvalue <- pvalues
  if (!log_scale) {
    c_df$foldchange <- foldchange
  }
  c_df$nlog10p <- nlog10p
  c_df$log2fc <- log2fc

  return(c_df)
}


#' Perform Principal Component Analysis (PCA)
#'
#' Given a numeric dataframe, performs PCA by transposing the dataframe (so that samples
#' become rows and features become columns). Automatically filters out rows with zero
#' variance and handles missing values by replacing them with zeros. Returns PCA results
#' including transformed coordinates, variance explained, and axis labels.
#'
#' @param df Dataframe containing numeric values. Rows are features (e.g., genes, proteins)
#'   and columns are samples.
#' @param sample A character vector containing the sample names/labels for each column in df.
#'   Must have the same length as the number of columns in df.
#' @param scale Logical, whether to scale the data to unit variance (default: TRUE).
#' @param center Logical, whether to center the data to zero mean (default: TRUE).
#'
#' @return A list containing three elements:
#'   - [[1]]: dataframe with PCA coordinates (PC1, PC2, ...) and Sample column
#'   - [[2]]: character string for x-axis label (PC1 with variance explained)
#'   - [[3]]: character string for y-axis label (PC2 with variance explained)
#'
#' @export
#'
#' @examples
#' # Perform PCA on proteomics data
#' sample_names <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
#' pca_results <- performPCA(protein_data, sample = sample_names)
#'
#' # PCA without scaling
#' pca_results <- performPCA(protein_data, sample = sample_names, scale = FALSE)
performPCA <- function(df, sample, scale = TRUE, center = TRUE) {
  pca_results <- list()
  df <- df %>% select(where(is.numeric))  %>%
    replace(is.na(.), 0) %>%
    filter(apply(., 1, var, na.rm = TRUE) > 0)
  pca_df <- df %>% t()
  pca_matrix <- as.matrix(pca_df)
  pca <- prcomp(pca_matrix, scale. = scale, center = center)
  pca_result_df <- as.data.frame(pca$x)
  pca_result_df$Sample <- sample
  variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
  xlabel <- paste0("Principal Component 1", " (", round(variance_explained[1] * 100, 2), "%)")
  ylabel <- paste0("Principal Component 2", " (", round(variance_explained[2] * 100, 2), "%)")
  pca_results[[1]] <- pca_result_df
  pca_results[[2]] <- xlabel
  pca_results[[3]] <- ylabel
  return(pca_results)
}

#' Plot Principal Component Analysis (PCA) Results
#'
#' Creates a 2D scatter plot of PCA results showing PC1 vs PC2. Takes the output from
#' performPCA() and generates a publication-ready plot with variance explained in axis labels.
#' Points are colored by sample groups and styled with the ggprism theme.
#'
#' @param pca_list A list containing PCA results from performPCA(). Expected structure:
#'   - [[1]]: dataframe with PC1, PC2, and Sample columns
#'   - [[2]]: x-axis label (PC1 with variance explained)
#'   - [[3]]: y-axis label (PC2 with variance explained)
#' @param title Character, the title of the plot (default: "Principal Component Analysis").
#'
#' @return A ggplot2 object containing the PCA scatter plot.
#'
#' @export
#'
#' @examples
#' # Perform PCA and plot
#' sample_names <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
#' pca_results <- performPCA(protein_data, sample = sample_names)
#' plotPCA(pca_results)
#'
#' # Plot with custom title
#' plotPCA(pca_results, title = "PCA of Protein Expression")
plotPCA <- function(pca_list,
                    title = "Principal Component Analysis") {
  pca_df <- pca_list[[1]]
  xlabel <- pca_list[[2]]
  ylabel <- pca_list[[3]]
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Sample, color = Sample)) +
    geom_point(shape = 21, color = "black", size = 4) +
    ggprism::theme_prism(base_size = 18) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold"),
          axis.text = element_text(size = 22, face = "bold"),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.minor.ticks.length = rel(0.5),
          plot.title = element_text(size = 24, face = "bold")) +
    labs(title = title, x = xlabel, y = ylabel) + stat_ellipse(geom = "polygon", alpha = 0.2, linetype = 0)
  return(pca_plot)
}

#' Create a Clustered Heatmap
#'
#' Given a numeric dataframe, creates a clustered heatmap with optional filtering
#' by p-value and log2 fold change. Automatically filters the data and generates
#' a ComplexHeatmap object with hierarchical clustering.
#'
#' @param df Dataframe containing numeric values. Should include pvalue and log2fc columns for filtering.
#'   If a "Gene" column is present, it will be used for row names.
#' @param scale Logical, whether to scale the data (default: TRUE).
#' @param show_rownames Logical, whether to show the row names (default: FALSE).
#' @param rownames Vector, custom row names for the dataframe (default: NULL).
#' @param p_cutoff Numeric, p-value cutoff for filtering rows (default: 0.05).
#' @param log2fc_cutoff Numeric, absolute log2 fold change cutoff for filtering rows (default: 1).
#' @param col Color mapping for the heatmap. Default is green-black-red scale from -2 to 2.
#' @param width Unit, the width of the heatmap (default: NULL for automatic).
#' @param height Unit, the height of the heatmap (default: NULL for automatic).
#'
#' @return A ComplexHeatmap object containing the clustered heatmap.
#'
#' @export
#'
#' @examples
#' # Basic heatmap
#' makeClustermap(results_df)
#'
#' # Heatmap with custom cutoffs
#' makeClustermap(results_df, p_cutoff = 0.01, log2fc_cutoff = 2)
makeClustermap <- function(df, scale = TRUE, show_rownames = FALSE, rownames = NULL,
                           p_cutoff = 0.05, log2fc_cutoff = 1,
                           col = circlize::colorRamp2(c(-2, 0, 2), c("green", "black", "red")),
                           width = NULL, height = NULL) {
  # if the df has a gene column, set it as the rownames
  df <- df %>%
    filter(pvalue < p_cutoff & abs(log2fc) >= log2fc_cutoff)
  if ("Gene" %in% colnames(df)) {
    rownames <- df$Gene
  }
  
  df <- df %>% dplyr::select(where(is.numeric))  %>%
    select(-c(log2fc, pvalue, padj, foldchange, nlog10p))
  require(ComplexHeatmap)
  colnames <- colnames(df)
  if (is.null(rownames)) {
    rownames(df) <- rownames
  }
  if (scale) {
    df <- as.data.frame(scale_matrix_rows(as.matrix(df), center = TRUE, scale = TRUE))
    colnames(df) <- colnames
    df <- replace(df, is.na(df), 0)
  }
  if (show_rownames) {
    hm <- ComplexHeatmap::Heatmap(df, cluster_columns = FALSE, col = col,
                                  show_row_names = TRUE, row_labels = rownames,
                                  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                  width = width, height = height,
                                  rect_gp = gpar(col = "black"), border_gp = gpar(col = "black"),
                                  heatmap_legend_param = list(
                                    title = "Z-Score", at = c(-2, 0, 2),
                                    labels = c("-2", "0", "2"),
                                    direction = "horizontal",
                                    legend_width = width,
                                    border = TRUE,
                                    font_size = 8
                                  ))
  } else {
    hm <- ComplexHeatmap::Heatmap(df, cluster_columns = FALSE, col = col,
                                  show_row_names = FALSE,
                                  column_names_gp = gpar(fontsize = 12),
                                  width = width, height = height, border_gp = gpar(col = "black"),
                                  heatmap_legend_param = list(
                                    title = "Z-Score", at = c(-2, 0, 2),
                                    labels = c("-2", "0", "2"),
                                    direction = "horizontal",
                                    legend_width = width,
                                    border = TRUE,
                                    labels_gp = gpar(fontsize = 12, fontface = "bold")
                                  ))
  }
  return(hm)
}

#' Create Volcano Plot
#' @title Make Volcano Plot
#' @description
#' Creates a volcano plot from differential expression results.
#' @param df A dataframe containing 'log2fc' and 'nlog10p' columns
#' @param fc_cutoff Fold change cutoff (default 0)
#' @param p_cutoff P-value cutoff (default 0.05)
#' @return A ggplot2 object containing the volcano plot
#' @export
#' @examples
#' # df <- data.frame(log2fc=rnorm(100), nlog10p=runif(100, 0, 5), pvalue=runif(100))
#' # plot <- makeVolcano(df, fc_cutoff=1, p_cutoff=0.05)
makeVolcano <- function(df, fc_cutoff = 0, p_cutoff = 0.05) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  if (!("log2fc" %in% names(df))) {
    stop("Must have a columns named log2fc")
  }
  if (!("nlog10p" %in% names(df))) {
    stop("Must have a column named nlog10p")
  }


  df <- df %>% 
    mutate(
      pvalue = ifelse(is.na(pvalue), 1, pvalue),
      log2fc = ifelse(is.na(log2fc), 0, log2fc)
    )
  ymax <- ceiling(max(df$nlog10p, na.rm = TRUE) * 1.3)
  xmax <- ceiling(max(abs(df$log2fc), na.rm = TRUE) * 1.3)
  xmin <- -xmax
  x_lims <- max(xmax, abs(xmin))
  df <- df %>% 
    mutate(status = case_when(
      (nlog10p > -log10(p_cutoff)) & (log2fc >= fc_cutoff) ~ "Upregulated",
      (nlog10p > -log10(p_cutoff)) & (log2fc <= -fc_cutoff) ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))
  volcanoPlot <- ggplot() +
    geom_point(data = df %>% filter(status != "Significant"),
               aes(x = log2fc, y = nlog10p, color = factor(status)),
               alpha = 0.6, shape = 19, size = 2) +
    geom_point(data = df %>% filter(status == "Significant"),
               aes(x = log2fc, y = nlog10p, color = factor(status)),
               alpha = 0.8, shape = 19, size = 3) +
    ggprism::theme_prism(base_size = 18) +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 14, face = "bold"),
          axis.title.y = element_text(size = 24, face = "bold"),
          axis.text = element_text(size = 18, face = "bold"),
          axis.title.x = element_text(size = 24, face = "bold"),
          axis.minor.ticks.length = rel(0.5),
          plot.title = element_text(size = 24, face = "bold")) +
    scale_color_manual(values = c("Upregulated" = "red2",
                                 "Downregulated" = "forestgreen",
                                 "Not Significant" = "darkgrey")) +
    labs(x = expression(log[2]*foldchange), y = expression(-log[10]*p)) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black",
               size = 0.75) 
  num_upregulated <- sum(df$status == "Upregulated")
  num_downregulated <- sum(df$status == "Downregulated")
  volcanoPlot <- volcanoPlot +
    geom_label(aes(x = (x_lims * 0.3), y = (ymax/2)), label = paste(num_upregulated),
              size = 5, fontface = "bold") +
    geom_label(aes(x = (-x_lims * 0.3), y = (ymax/2)), label = paste(num_downregulated),
              size = 5, fontface = "bold")
  if (fc_cutoff != 0) {
    volcanoPlot <- volcanoPlot + geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                                            linetype = "dashed", color = "black",
                                            size = 0.75)
  } else {
    volcanoPlot <- volcanoPlot
  }
  return(volcanoPlot)
}


#' Perform ANOVA on Multiple Cohorts
#'
#' Takes data frames as arguments, each belonging to a cohort, and performs an ANOVA on each row.
#' Works dynamically with 3 or more groups. The function returns a data frame with the p-value
#' from the ANOVA, Tukey's HSD p-values for pairwise comparisons, and fold changes for each
#' pairwise comparison. Column names are automatically derived from the variable names passed in.
#'
#' @param ... Data frames, each containing data for one cohort.
#'   All data frames must have the same number of rows.
#' @param comparisons Which pairwise comparisons to include in output. Options:
#'   - "all" (default): include all pairwise comparisons
#'   - A list of character vectors specifying pairs, e.g.,
#'     list(c("control", "treated"), c("control", "placebo")).
#'     Group names must match the variable names passed in.
#'
#' @return A data frame with the combined original data plus computed columns:
#'   - anova_pvalue: overall ANOVA p-value for each row
#'   - tukey_pvalue_X_Y: Tukey's HSD adjusted p-values for each pairwise comparison
#'   - foldchange_X_Y: fold changes (mean difference) for each pairwise comparison
#'
#' @export
#'
#' @examples
#' # For 3 groups (all pairwise comparisons)
#' result <- performANOVA(control, treated, placebo)
#'
#' # Only specific comparisons
#' result <- performANOVA(control, treated, placebo,
#'                        comparisons = list(c("control", "treated"),
#'                                           c("control", "placebo")))
performANOVA <- function(..., comparisons = "all") {
  require(dplyr)
  require(stats)

  # Capture the call to extract variable names
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  # Remove 'comparisons' from args if present
  args$comparisons <- NULL
  arg_names <- as.character(args)

  # Get the actual data frames
  df_list <- list(...)
  # Remove comparisons from df_list if it got included
  if (!is.null(names(df_list))) {
    df_list <- df_list[names(df_list) != "comparisons"]
  }

  n_groups <- length(df_list)
  if (n_groups < 3) stop("Need at least 3 groups for ANOVA")

  # Use variable names as group names
  group_names <- arg_names

  # Number of rows (assuming all data frames have the same number of rows)
  n_rows <- nrow(df_list[[1]])

  # Generate all pairwise combinations (needed for Tukey regardless of output filter)
  all_pairs <- combn(n_groups, 2)
  all_n_pairs <- ncol(all_pairs)
  all_pair_names <- apply(all_pairs, 2, function(p) paste0(group_names[p[1]], "_", group_names[p[2]]))

  # Determine which pairs to include in output
  if (identical(comparisons, "all")) {
    output_pair_indices <- seq_len(all_n_pairs)
  } else {
    # comparisons is a list of pairs like list(c("control", "treated"), ...)
    output_pair_indices <- integer(0)
    for (comp in comparisons) {
      if (length(comp) != 2) stop("Each comparison must be a vector of 2 group names")
      # Find this pair in all_pairs (order doesn't matter)
      found <- FALSE
      for (j in seq_len(all_n_pairs)) {
        g1_name <- group_names[all_pairs[1, j]]
        g2_name <- group_names[all_pairs[2, j]]
        if ((comp[1] == g1_name && comp[2] == g2_name) ||
            (comp[1] == g2_name && comp[2] == g1_name)) {
          output_pair_indices <- c(output_pair_indices, j)
          found <- TRUE
          break
        }
      }
      if (!found) {
        stop(paste0("Comparison not found: ", comp[1], " vs ", comp[2],
                    ". Available groups: ", paste(group_names, collapse = ", ")))
      }
    }
  }

  # Initialize result storage for ALL pairs (Tukey computes all anyway)
  pvalues <- numeric(n_rows)
  tukey_pvalues <- matrix(1, nrow = n_rows, ncol = all_n_pairs)
  foldchanges <- matrix(0, nrow = n_rows, ncol = all_n_pairs)

  # Loop over each row
  for (i in seq_len(n_rows)) {
    tryCatch({
      # Extract values for each group
      values_list <- lapply(df_list, function(df) as.numeric(df[i, ]))

      # Check if all groups have >= 3 non-NA values
      enough_data <- all(sapply(values_list, function(v) sum(!is.na(v)) >= 3))

      if (enough_data) {
        # Combine values and create group factor
        values <- unlist(values_list)
        groups <- factor(rep(group_names, sapply(values_list, length)),
                         levels = group_names)

        # Perform ANOVA
        anova_result <- aov(values ~ groups)
        pvalues[i] <- summary(anova_result)[[1]][["Pr(>F)"]][1]

        if (!is.na(pvalues[i])) {
          tukey_result <- TukeyHSD(anova_result)
          comp <- tukey_result$groups

          # Extract Tukey p-values and fold changes for each pair
          for (j in seq_len(all_n_pairs)) {
            g1 <- all_pairs[1, j]
            g2 <- all_pairs[2, j]
            comp_name <- paste0(group_names[g2], "-", group_names[g1])

            tukey_pvalues[i, j] <- comp[comp_name, "p adj"]
            foldchanges[i, j] <- mean(values_list[[g2]], na.rm = TRUE) -
                                 mean(values_list[[g1]], na.rm = TRUE)
          }
        }
      } else {
        # If not enough data for ANOVA, set p-value = 1
        pvalues[i] <- 1
        # tukey_pvalues and foldchanges already initialized to 1 and 0
      }
    }, error = function(e) {
      # If an error occurs, use safe defaults (already initialized)
      pvalues[i] <<- 1
    })
  }

  # Build result data frame
  c_df <- dplyr::bind_cols(df_list)
  c_df$anova_pvalue <- pvalues

  # Add Tukey p-values and fold changes ONLY for requested comparisons
  for (j in output_pair_indices) {
    c_df[[paste0("tukey_pvalue_", all_pair_names[j])]] <- tukey_pvalues[, j]
    c_df[[paste0("foldchange_", all_pair_names[j])]] <- foldchanges[, j]
  }

  return(c_df)
}

#' Perform 2-Way ANOVA with Automatic Grouping
#'
#' Performs 2-way ANOVA on each row of a dataframe with automatic grouping based on two factors.
#' Takes a single dataframe where columns represent different combinations of two factors.
#' For example, if you have male control, female control, male disease, and female disease samples,
#' you can specify factor1 = c("Male", "Female") and factor2 = c("control", "disease"), and the
#' function will automatically group the columns accordingly.
#'
#' @param df A dataframe where each row is a variable and columns represent different samples.
#' @param factor1 A character vector of levels for the first factor (e.g., c("Male", "Female")).
#' @param factor2 A character vector of levels for the second factor (e.g., c("control", "disease")).
#' @param comparisons Either "all" for all pairwise comparisons, or a list of specific comparisons
#'   for Tukey's post-hoc test. Each comparison should be a vector of 2 group names formed by
#'   combining factor levels with a space (e.g., "Male control").
#'   Example: list(c("Male control", "Male disease"), c("Female control", "Female disease"))
#'
#' @return A dataframe with ANOVA results, including main effects, interaction, and Tukey post-hoc tests:
#'   - anova2way_pvalue_overall: minimum p-value across all effects
#'   - anova2way_pvalue_factor1: main effect p-value for first factor
#'   - anova2way_pvalue_factor2: main effect p-value for second factor
#'   - anova2way_pvalue_interaction: interaction effect p-value
#'   - tukey_pvalue_X_Y: Tukey's HSD p-values for specified comparisons
#'   - foldchange_X_Y: fold changes (mean differences) for specified comparisons
#'
#' @export
#'
#' @examples
#' # Example with male/female and control/disease
#' result <- perform2WayANOVA(data,
#'                            factor1 = c("Male", "Female"),
#'                            factor2 = c("control", "disease"),
#'                            comparisons = list(c("Male control", "Male disease"),
#'                                             c("Female control", "Female disease")))
perform2WayANOVA <- function(df, factor1, factor2, comparisons = "all") {
  require(dplyr)
  require(stats)

  # Create all combinations of the two factors
  combinations <- expand.grid(factor1 = factor1, factor2 = factor2, stringsAsFactors = FALSE)
  group_names <- paste(combinations$factor1, combinations$factor2)
  n_groups <- length(group_names)

  if (n_groups < 2) stop("Need at least 2 groups for 2-way ANOVA")

  # Verify that we have the right number of columns in the dataframe
  # Assuming each group has multiple replicates, we need to determine how many replicates per group
  n_cols <- ncol(df)
  if (n_cols %% n_groups != 0) {
    stop(paste0("Number of columns (", n_cols, ") is not evenly divisible by number of groups (",
                n_groups, "). Each group should have the same number of replicates."))
  }

  n_replicates <- n_cols / n_groups
  n_rows <- nrow(df)

  # Create group assignment vector (assumes columns are ordered by group)
  group_vector <- rep(group_names, each = n_replicates)
  factor1_vector <- rep(combinations$factor1, each = n_replicates)
  factor2_vector <- rep(combinations$factor2, each = n_replicates)

  # Generate all pairwise combinations for Tukey
  all_pairs <- combn(n_groups, 2)
  all_n_pairs <- ncol(all_pairs)
  all_pair_names <- apply(all_pairs, 2, function(p) paste0(group_names[p[1]], "_", group_names[p[2]]))

  # Determine which pairs to include in output
  if (identical(comparisons, "all")) {
    output_pair_indices <- seq_len(all_n_pairs)
  } else {
    output_pair_indices <- integer(0)
    for (comp in comparisons) {
      if (length(comp) != 2) stop("Each comparison must be a vector of 2 group names")
      found <- FALSE
      for (j in seq_len(all_n_pairs)) {
        g1_name <- group_names[all_pairs[1, j]]
        g2_name <- group_names[all_pairs[2, j]]
        if ((comp[1] == g1_name && comp[2] == g2_name) ||
            (comp[1] == g2_name && comp[2] == g1_name)) {
          output_pair_indices <- c(output_pair_indices, j)
          found <- TRUE
          break
        }
      }
      if (!found) {
        stop(paste0("Comparison not found: ", comp[1], " vs ", comp[2],
                    ". Available groups: ", paste(group_names, collapse = ", ")))
      }
    }
  }

  # Initialize result storage
  pvalues_factor1 <- numeric(n_rows)
  pvalues_factor2 <- numeric(n_rows)
  pvalues_interaction <- numeric(n_rows)
  pvalues_overall <- numeric(n_rows)
  tukey_pvalues <- matrix(1, nrow = n_rows, ncol = all_n_pairs)
  foldchanges <- matrix(0, nrow = n_rows, ncol = all_n_pairs)

  # Loop over each row
  for (i in seq_len(n_rows)) {
    tryCatch({
      # Extract values for this row
      values <- as.numeric(df[i, ])

      # Check if we have enough non-NA values
      n_valid <- sum(!is.na(values))

      if (n_valid >= n_groups * 2) {  # At least 2 replicates per group minimum
        # Create dataframe for ANOVA
        anova_data <- data.frame(
          value = values,
          factor1 = factor(factor1_vector, levels = factor1),
          factor2 = factor(factor2_vector, levels = factor2),
          group = factor(group_vector, levels = group_names)
        )

        # Remove NA values
        anova_data <- anova_data[!is.na(anova_data$value), ]

        # Perform 2-way ANOVA
        anova_result <- aov(value ~ factor1 * factor2, data = anova_data)
        anova_summary <- summary(anova_result)[[1]]

        # Extract p-values for main effects and interaction
        pvalues_factor1[i] <- anova_summary[["Pr(>F)"]][1]
        pvalues_factor2[i] <- anova_summary[["Pr(>F)"]][2]
        pvalues_interaction[i] <- anova_summary[["Pr(>F)"]][3]

        # Overall model p-value (minimum of the three)
        pvalues_overall[i] <- min(pvalues_factor1[i], pvalues_factor2[i],
                                  pvalues_interaction[i], na.rm = TRUE)

        # Perform Tukey's HSD on the group factor for post-hoc comparisons
        anova_group <- aov(value ~ group, data = anova_data)

        if (!is.na(pvalues_overall[i])) {
          tukey_result <- TukeyHSD(anova_group)
          comp <- tukey_result$group

          # Calculate group means for fold changes
          group_means <- sapply(seq_len(n_groups), function(g) {
            group_data <- values[group_vector == group_names[g]]
            mean(group_data, na.rm = TRUE)
          })

          # Extract Tukey p-values and fold changes for each pair
          for (j in seq_len(all_n_pairs)) {
            g1 <- all_pairs[1, j]
            g2 <- all_pairs[2, j]
            comp_name <- paste0(group_names[g2], "-", group_names[g1])

            tukey_pvalues[i, j] <- comp[comp_name, "p adj"]
            foldchanges[i, j] <- group_means[g2] - group_means[g1]
          }
        }
      } else {
        # Not enough data
        pvalues_factor1[i] <- 1
        pvalues_factor2[i] <- 1
        pvalues_interaction[i] <- 1
        pvalues_overall[i] <- 1
      }
    }, error = function(e) {
      # If an error occurs, use safe defaults
      pvalues_factor1[i] <<- 1
      pvalues_factor2[i] <<- 1
      pvalues_interaction[i] <<- 1
      pvalues_overall[i] <<- 1
    })
  }

  # Build result data frame
  result_df <- df
  result_df$anova2way_pvalue_overall <- pvalues_overall
  result_df[[paste0("anova2way_pvalue_", deparse(substitute(factor1)))]] <- pvalues_factor1
  result_df[[paste0("anova2way_pvalue_", deparse(substitute(factor2)))]] <- pvalues_factor2
  result_df$anova2way_pvalue_interaction <- pvalues_interaction

  # Add Tukey p-values and fold changes ONLY for requested comparisons
  for (j in output_pair_indices) {
    result_df[[paste0("tukey_pvalue_", all_pair_names[j])]] <- tukey_pvalues[, j]
    result_df[[paste0("foldchange_", all_pair_names[j])]] <- foldchanges[, j]
  }

  return(result_df)
}

#' Perform Kruskal-Wallis Test on Multiple Cohorts
#'
#' Takes data frames as arguments, each belonging to a cohort, and performs a non-parametric
#' Kruskal-Wallis test on each row. This is the non-parametric alternative to ANOVA for comparing
#' three or more groups. Works dynamically with 3 or more groups. The function returns a data frame
#' with the p-value from the Kruskal-Wallis test, Dunn's test p-values for pairwise comparisons,
#' and fold changes for each pairwise comparison. Column names are automatically derived from
#' the variable names passed in.
#'
#' @param ... Data frames, each containing data for one cohort.
#'   All data frames must have the same number of rows.
#' @param comparisons Which pairwise comparisons to include in output. Options:
#'   - "all" (default): include all pairwise comparisons
#'   - A list of character vectors specifying pairs, e.g.,
#'     list(c("control", "treated"), c("control", "placebo")).
#'     Group names must match the variable names passed in.
#' @param p_adjust Method for p-value adjustment in Dunn's test. Default is "BH" (Benjamini-Hochberg).
#'   Other options: "bonferroni", "holm", "hochberg", "hommel", "BY", "none".
#'
#' @return A data frame with the combined original data plus computed columns:
#'   - kw_pvalue: overall Kruskal-Wallis p-value for each row
#'   - dunn_pvalue_X_Y: Dunn's test adjusted p-values for each pairwise comparison
#'   - foldchange_X_Y: fold changes (mean difference) for each pairwise comparison
#'
#' @export
#'
#' @examples
#' # For 3 groups (all pairwise comparisons)
#' result <- performKW(control, treated, placebo)
#'
#' # Only specific comparisons
#' result <- performKW(control, treated, placebo,
#'                     comparisons = list(c("control", "treated"),
#'                                        c("control", "placebo")))
performKW <- function(..., comparisons = "all", p_adjust = "BH") {
  require(dplyr)
  require(stats)

  # Check if dunn.test is installed
  if (!requireNamespace("dunn.test", quietly = TRUE)) {
    stop("Package 'dunn.test' is required for Dunn's post-hoc test.\n",
         "Please install it with: install.packages('dunn.test')")
  }

  # Capture the call to extract variable names
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  # Remove named arguments from args
  args$comparisons <- NULL
  args$p_adjust <- NULL
  arg_names <- as.character(args)

  # Get the actual data frames
  df_list <- list(...)
  # Remove named arguments from df_list if they got included
  if (!is.null(names(df_list))) {
    df_list <- df_list[!names(df_list) %in% c("comparisons", "p_adjust")]
  }

  n_groups <- length(df_list)
  if (n_groups < 3) stop("Need at least 3 groups for Kruskal-Wallis")

  # Use variable names as group names
  group_names <- arg_names

  # Number of rows (assuming all data frames have the same number of rows)
  n_rows <- nrow(df_list[[1]])

  # Generate all pairwise combinations
  all_pairs <- combn(n_groups, 2)
  all_n_pairs <- ncol(all_pairs)
  all_pair_names <- apply(all_pairs, 2, function(p) paste0(group_names[p[1]], "_", group_names[p[2]]))

  # Determine which pairs to include in output
  if (identical(comparisons, "all")) {
    output_pair_indices <- seq_len(all_n_pairs)
  } else {
    output_pair_indices <- integer(0)
    for (comp in comparisons) {
      if (length(comp) != 2) stop("Each comparison must be a vector of 2 group names")
      found <- FALSE
      for (j in seq_len(all_n_pairs)) {
        g1_name <- group_names[all_pairs[1, j]]
        g2_name <- group_names[all_pairs[2, j]]
        if ((comp[1] == g1_name && comp[2] == g2_name) ||
            (comp[1] == g2_name && comp[2] == g1_name)) {
          output_pair_indices <- c(output_pair_indices, j)
          found <- TRUE
          break
        }
      }
      if (!found) {
        stop(paste0("Comparison not found: ", comp[1], " vs ", comp[2],
                    ". Available groups: ", paste(group_names, collapse = ", ")))
      }
    }
  }

  # Initialize result storage
  pvalues <- numeric(n_rows)
  dunn_pvalues <- matrix(1, nrow = n_rows, ncol = all_n_pairs)
  foldchanges <- matrix(0, nrow = n_rows, ncol = all_n_pairs)

  # Loop over each row
  for (i in seq_len(n_rows)) {
    tryCatch({
      # Extract values for each group
      values_list <- lapply(df_list, function(df) as.numeric(df[i, ]))

      # Check if all groups have >= 3 non-NA values
      enough_data <- all(sapply(values_list, function(v) sum(!is.na(v)) >= 3))

      if (enough_data) {
        # Combine values and create group factor
        values <- unlist(values_list)
        groups <- factor(rep(group_names, sapply(values_list, length)),
                         levels = group_names)

        # Remove NAs
        valid_idx <- !is.na(values)
        values <- values[valid_idx]
        groups <- groups[valid_idx]

        # Perform Kruskal-Wallis test
        kw_result <- kruskal.test(values ~ groups)
        pvalues[i] <- kw_result$p.value

        if (!is.na(pvalues[i])) {
          # Perform Dunn's test (suppress all printed output)
          sink(tempfile())
          dunn_result <- tryCatch(
            dunn.test::dunn.test(values, groups, method = p_adjust, kw = FALSE, table = FALSE, list = FALSE),
            finally = sink()
          )

          # Extract Dunn's p-values for each pair
          # dunn.test returns comparisons in format "group1 - group2"
          for (j in seq_len(all_n_pairs)) {
            g1 <- all_pairs[1, j]
            g2 <- all_pairs[2, j]

            # Try both orderings to find the comparison
            comp_name1 <- paste0(group_names[g1], " - ", group_names[g2])
            comp_name2 <- paste0(group_names[g2], " - ", group_names[g1])

            idx <- which(dunn_result$comparisons == comp_name1 | dunn_result$comparisons == comp_name2)

            if (length(idx) > 0) {
              dunn_pvalues[i, j] <- dunn_result$P.adjusted[idx[1]]
            }

            # Calculate fold change
            foldchanges[i, j] <- mean(values_list[[g2]], na.rm = TRUE) -
                                 mean(values_list[[g1]], na.rm = TRUE)
          }
        }
      } else {
        pvalues[i] <- 1
      }
    }, error = function(e) {
      pvalues[i] <<- 1
    })
  }

  # Build result data frame
  c_df <- dplyr::bind_cols(df_list)
  c_df$kw_pvalue <- pvalues

  # Add Dunn's p-values and fold changes ONLY for requested comparisons
  for (j in output_pair_indices) {
    c_df[[paste0("dunn_pvalue_", all_pair_names[j])]] <- dunn_pvalues[, j]
    c_df[[paste0("foldchange_", all_pair_names[j])]] <- foldchanges[, j]
  }

  return(c_df)
}


#' convertFormat
#' @description
#' Converts ANOVA or Kruskal-Wallis results to the format expected by makeVolcano. Extracts
#' the p-value and fold change for a specific comparison and renames them to pvalue, nlog10p,
#' and log2fc. Automatically detects whether input is from performANOVA (tukey_pvalue) or
#' performKW (dunn_pvalue).
#'
#' @param results Data frame output from performANOVA or performKW
#' @param comparison Character vector of length 2 specifying the comparison,
#'   e.g., c("control", "treated"). Order doesn't matter.
#'
#' @return A data frame with added columns:
#'   - pvalue: the Tukey/Dunn p-value for the specified comparison
#'   - nlog10p: -log10 transformed p-value
#'   - log2fc: the fold change for the specified comparison
#' @export
#'
#' @examples
#' # From ANOVA:
#' anova_res <- performANOVA(control, treated, placebo)
#' volcano_data <- convertFormat(anova_res, c("control", "treated"))
#' makeVolcano(volcano_data)
#'
#' # From Kruskal-Wallis:
#' kw_res <- performKW(control, treated, placebo)
#' volcano_data <- convertFormat(kw_res, c("control", "treated"))
#' makeVolcano(volcano_data)
convertFormat <- function(results, comparison) {
  if (length(comparison) != 2) {
    stop("comparison must be a character vector of length 2")
  }

  # Try both orderings of the comparison to find the matching column
  pair_name1 <- paste0(comparison[1], "_", comparison[2])
  pair_name2 <- paste0(comparison[2], "_", comparison[1])

  # Check for tukey_pvalue (from performANOVA) or dunn_pvalue (from performKW)
  tukey_col1 <- paste0("tukey_pvalue_", pair_name1)
  tukey_col2 <- paste0("tukey_pvalue_", pair_name2)
  dunn_col1 <- paste0("dunn_pvalue_", pair_name1)
  dunn_col2 <- paste0("dunn_pvalue_", pair_name2)
  fc_col1 <- paste0("foldchange_", pair_name1)
  fc_col2 <- paste0("foldchange_", pair_name2)

  # Find which column exists
  pvalue_col <- NULL
  fc_col <- NULL
  flip_sign <- FALSE

  if (tukey_col1 %in% names(results)) {
    pvalue_col <- tukey_col1
    fc_col <- fc_col1
  } else if (tukey_col2 %in% names(results)) {
    pvalue_col <- tukey_col2
    fc_col <- fc_col2
    flip_sign <- TRUE
  } else if (dunn_col1 %in% names(results)) {
    pvalue_col <- dunn_col1
    fc_col <- fc_col1
  } else if (dunn_col2 %in% names(results)) {
    pvalue_col <- dunn_col2
    fc_col <- fc_col2
    flip_sign <- TRUE
  } else {
    # List available comparisons from either tukey or dunn columns
    available <- grep("^(tukey_pvalue_|dunn_pvalue_)", names(results), value = TRUE)
    available <- gsub("^(tukey_pvalue_|dunn_pvalue_)", "", available)
    stop(paste0("Comparison '", comparison[1], "' vs '", comparison[2], "' not found.\n",
                "Available comparisons: ", paste(available, collapse = ", ")))
  }

  # Create output with renamed columns
  result <- results
  result$pvalue <- result[[pvalue_col]]
  result$log2fc <- if (flip_sign) -result[[fc_col]] else result[[fc_col]]
  result$nlog10p <- -log10(result$pvalue)

  return(result)
}


# Function: makeBarplot
# Takes 2 dataframes, one long dataframe to be used to make the dataframe and a
# wide dataframe containing the pvalues
#
# Args:
# @param df: dataframe containing the data
# @param pvalue_df: dataframe containing the pvalues
# @param samples: vector containing sample labels (must be in order you want the barplot)
#
# Returns:
# @return: ggplot2 object containing the barplot
#' Create Barplot with Statistical Annotations
#' @title Make Barplot
#' @description
#' Creates a barplot for a specific gene with Tukey HSD p-values from ANOVA.
#' @param df A long-format dataframe with 'Treatment' and 'Abundance' columns
#' @param pvalue_df A wide-format dataframe containing Tukey p-values
#' @param cohort_labels A vector of cohort labels
#' @param gene The gene name to plot
#' @return A ggplot2 object containing the barplot
#' @export
#' @examples
#' # df_long <- data.frame(Treatment=rep(c("A","B","C","D"),each=5), Abundance=rnorm(20))
#' # pvalue_df <- data.frame(Gene="ACTB", tukey_pvalue12=0.01, tukey_pvalue13=0.05, tukey_pvalue14=0.001)
#' # plot <- makeBarplot(df_long, pvalue_df, c("A","B","C","D"), "ACTB")
makeBarplot <- function(df, pvalue_df, cohort_labels, gene) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(tidyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggpubr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  ymax <- (ceiling(max(df$Abundance, na.rm = TRUE)/5) * 5) + 5
  ymin <- (floor(min(df$Abundance, na.rm = TRUE)/5) * 5)
  breaks <- seq(ymin, ymax, by = 5)
  
  pvalue_list <- list()
  
  # Create p-value dataframes for comparisons
  pvalue_df_1 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[2],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.05, "*", ""))),
    y.position = ymax - 5.8
  )
  if (pvalue_df_1$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_1))
  
  pvalue_df_2 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[3],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.05, "*", ""))),
    y.position = ymax - 3.5
  )
  if (pvalue_df_2$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_2))
  
  pvalue_df_3 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[4],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue14),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue14) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue14) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue14) < 0.05, "*", ""))),
    y.position = ymax - 1.5
  )
  if (pvalue_df_3$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_3))
  
  # Create the plot
  p <- ggplot(df, aes(x = Treatment, y = Abundance)) +
    geom_bar(aes(fill = Treatment), stat = "summary", fun = "mean", position = "dodge", width = 0.6, alpha = 0.6, color = 'black') +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.6), width = 0.2) +
    geom_jitter(aes(fill = Treatment), shape = 21, color = "black", size = 1.5, width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    ggprism::theme_prism(base_size = 8) +
    labs(title = gene, x = NULL, y = "Log2Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(breaks = breaks, expand = c(0, 0)) +
    theme(legend.position = "none")
  
  # Add p-values to the plot
  if (length(pvalue_list) > 0) {
    for (pval_df in pvalue_list) {
      p <- p + ggprism::add_pvalue(pval_df, label = "p.adj.signif", y.position = pval_df$y.position,
                                   tip.length = 0.04, bracket.size = 0.3, label.size = 5)
    }
  }
  
  return(p)
}



#' makeBarplot3
#'
#' @param df long dataframe containing the values
#' @param pvalue_df the main wide dataframe containing the pvalues from the ANOVA result
#' @param cohort_labels vector containing the cohort labels
#' @param gene the gene to be plotted
#'
#' @return ggplot2 object containing the barplot
#' @export
#'
#' @examples
#' \dontrun{
#'   df_long <- data.frame(
#'     Treatment = rep(c("A", "B", "C"), each = 3),
#'     Abundance = rnorm(9)
#'   )
#'   tukey_df <- data.frame(
#'     Gene = "Gene1",
#'     tukey_pvalue12 = 0.03,
#'     tukey_pvalue13 = 0.2,
#'     tukey_pvalue23 = 0.04
#'   )
#'   makeBarplot3(df_long, tukey_df, c("A", "B", "C"), "Gene1")
#' }
makeBarplot3 <- function(df, pvalue_df, cohort_labels, gene) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(tidyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggpubr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  ymax <- (ceiling(max(df$Abundance, na.rm = TRUE)/5) * 5) + 5
  ymin <- (floor(min(df$Abundance, na.rm = TRUE)/5) * 5)
  breaks <- seq(ymin, ymax, by = 5)
  
  pvalue_list <- list()
  
  # Create p-value dataframes for comparisons (for 3 groups)
  pvalue_df_1 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[2],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.05, "*", ""))),
    y.position = ymax - 5.8
  )
  if (pvalue_df_1$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_1))
  
  pvalue_df_2 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[3],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.05, "*", ""))),
    y.position = ymax - 3.5
  )
  if (pvalue_df_2$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_2))
  
  pvalue_df_3 <- data.frame(
    group1 = cohort_labels[2],
    group2 = cohort_labels[3],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.05, "*", ""))),
    y.position = ymax - 1.5
  )
  if (pvalue_df_3$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_3))
  
  # Create the plot
  p <- ggplot(df, aes(x = Treatment, y = Abundance)) +
    geom_bar(aes(fill = Treatment), stat = "summary", fun = "mean", position = "dodge", width = 0.6, alpha = 0.6, color = 'black') +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.6), width = 0.2) +
    geom_jitter(aes(fill = Treatment), shape = 21, color = "black", size = 1.5, width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    ggprism::theme_prism(base_size = 8) +
    labs(title = gene, x = NULL, y = "Log2Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(breaks = breaks, expand = c(0, 0)) +
    theme(legend.position = "none")
  
  # Add p-values to the plot
  if (length(pvalue_list) > 0) {
    for (pval_df in pvalue_list) {
      p <- p + ggprism::add_pvalue(pval_df, label = "p.adj.signif", y.position = pval_df$y.position,
                                   tip.length = 0.04, bracket.size = 0.3, label.size = 5)
    }
  }
  
  return(p)
}

  

  
#' Create Simple Barplot
#' @title Make Simple Barplot
#' @description
#' Creates a simple barplot without statistical annotations.
#' @param df A long-format dataframe with 'Treatment' and 'Abundance' columns
#' @param cohort_labels A vector of cohort labels
#' @param gene The gene name to plot
#' @param num Additional space to add to y-axis maximum (default 0)
#' @return A ggplot2 object containing the barplot
#' @export
#' @examples
#' # df_long <- data.frame(Treatment=rep(c("A","B"),each=5), Abundance=rnorm(10))
#' # plot <- makeBarplot_simp(df_long, c("A","B"), "ACTB")
makeBarplot_simp <- function(df, cohort_labels, gene, num = 0) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  ymax <- (ceiling(max(df$Abundance, na.rm = TRUE)/4) * 4)  + num
  ymin <- floor(min(df$Abundance, na.rm = TRUE))
  breaks <- seq(ymin, ymax, by = 5)
  
  # Create the plot
  p <- ggplot(df, aes(x = Treatment, y = Abundance)) +
    geom_bar(aes(fill = Treatment), stat = "summary", fun = "mean", position = "dodge", width = 0.6, alpha = 0.6, color = 'black') +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.6), width = 0.2) +
    geom_jitter(aes(fill = Treatment), shape = 21, color = "black", size = 1.5, width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    ggprism::theme_prism(base_size = 8) +
    labs(title = gene, x = NULL, y = "Log2Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none")
  
  return(p)
}

#' calculateROC
#'
#' @description
#' Given a dataframe of predictors and a target vector,
#' this function calculates ROC curves, AUC values, confidence intervals,
#' and associated metrics for each specified predictor.
#'
#' @param df A dataframe containing predictor columns.
#' @param identifiers A character vector of column names in \code{df}
#'   for which the ROC analysis should be performed.
#' @param target_vector A numeric or factor vector containing the target values
#'   (commonly 0/1 for binary classification).
#' @param ci_level A numeric value between 0 and 1 specifying the confidence
#'   level for the AUC CI (default is 0.95).
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item \strong{identifier}: Name of the predictor column
#'     \item \strong{auc}: The AUC for that predictor
#'     \item \strong{ci_lower}: Lower bound of the AUC confidence interval
#'     \item \strong{ci_upper}: Upper bound of the AUC confidence interval
#'     \item \strong{sensitivity}: Sensitivity at various thresholds
#'     \item \strong{specificity}: 1 - Specificity at those thresholds
#'   }
#'
#' @importFrom pROC roc ci.auc
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#' library(pROC)
#' library(dplyr)
#'
#' # Sample data
#' set.seed(123)
#' df_example <- data.frame(
#'   Predictor1 = runif(50, min = 0, max = 1),
#'   Predictor2 = runif(50, min = 0, max = 1)
#' )
#'
#' # Binary outcome (0 or 1)
#' target_vec <- rbinom(50, size = 1, prob = 0.3)
#'
#' # Calculate ROC for both columns
#' roc_results <- calculateROC(
#'   df = df_example,
#'   identifiers = c("Predictor1", "Predictor2"),
#'   target_vector = target_vec
#' )
#'
#' head(roc_results)
#' }
#'
#' @export
calculateROC <- function(df, identifiers, target_vector, ci_level = 0.95) {
  
  # Ensure pROC is installed
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Please install it with install.packages('pROC').")
  }
  
  # Check for missing columns
  if (!all(identifiers %in% colnames(df))) {
    missing_cols <- setdiff(identifiers, colnames(df))
    warning(
      "The following identifiers are not in 'df': ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Initialize an empty data frame for the combined results
  roc_df <- data.frame()
  
  for (identifier in identifiers) {
    
    # Skip if column does not exist
    if (!identifier %in% colnames(df)) {
      warning(paste("Column", identifier, "not found in dataframe. Skipping."))
      next
    }
    
    predictor <- df[[identifier]]
    
    # Check if predictor is numeric
    if (!is.numeric(predictor)) {
      warning(paste("Column", identifier, "is not numeric. Skipping ROC calculation."))
      next
    }
    
    # Calculate ROC curve
    roc_obj <- tryCatch(
      pROC::roc(response = target_vector, predictor = predictor),
      error = function(e) {
        warning(paste("Error calculating ROC for", identifier, ":", e$message))
        return(NULL)
      }
    )
    
    # If ROC calculation succeeded, extract metrics
    if (!is.null(roc_obj)) {
      auc_value <- roc_obj$auc
      ci_vals   <- pROC::ci.auc(roc_obj, conf.level = ci_level)
      
      # pROC's ci.auc typically returns c(lower, median, upper).
      ci_lower  <- ci_vals[1]
      ci_upper  <- ci_vals[3]
      std_error <- (ci_upper - ci_lower) / 3.92
      
      # Build a data frame of results for each threshold
      roc_data <- data.frame(
        identifier   = identifier,
        auc          = as.numeric(auc_value),
        ci_lower     = as.numeric(ci_lower),
        ci_upper     = as.numeric(ci_upper),
        std_error    = as.numeric(std_error),
        sensitivity  = roc_obj$sensitivities,
        specificity  = 1 - roc_obj$specificities
      )
      
      # Append to the master results
      roc_df <- dplyr::bind_rows(roc_df, roc_data)
    }
  }
  
  return(roc_df)
}


#' calculateROC_SVM_Poly
#'
#' @description
#' Fits a polynomial-kernel SVM model (via \code{e1071::svm()}) to a binary
#' target using all columns in \code{df} as predictors, then calculates ROC
#' curves, AUC, confidence intervals, and metrics similar to \code{calculateROC}.
#'
#' @param df A dataframe containing predictor columns (no target column).
#' @param target_vector A numeric or factor vector containing the binary target
#'   values (commonly 0/1).
#' @param ci_level A numeric value between 0 and 1 specifying the confidence
#'   level for the AUC CI (default is 0.95).
#'
#' @return A data frame similar in structure to \code{calculateROC}, containing:
#'   \itemize{
#'     \item \strong{identifier}: A character label (e.g., "SVM_Poly")
#'     \item \strong{auc}: The AUC for that SVM model
#'     \item \strong{ci_lower}: Lower bound of the AUC confidence interval
#'     \item \strong{ci_upper}: Upper bound of the AUC confidence interval
#'     \item \strong{std_error}: An approximate standard error of the AUC
#'     \item \strong{sensitivity}: Sensitivity at various thresholds
#'     \item \strong{specificity}: 1 - Specificity at those thresholds
#'   }
#'
#' @examples
#' \dontrun{
#' library(e1071)
#' library(pROC)
#'
#' # Example data
#' set.seed(123)
#' df_example <- data.frame(
#'   x1 = runif(50),
#'   x2 = runif(50)
#' )
#'
#' # Binary outcome (0 or 1)
#' target_vec <- rbinom(50, size = 1, prob = 0.3)
#'
#' # Fit polynomial-kernel SVM & compute ROC
#' svm_roc <- calculateROC_SVM_Poly(
#'   df = df_example,
#'   target_vector = target_vec
#' )
#'
#' head(svm_roc)
#' }
#'
#' @importFrom pROC roc ci.auc
#' @importFrom e1071 svm
#' @export
calculateSVM <- function(df, target_vector, ci_level = 0.95) {
  # Check for e1071
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required. Install with install.packages('e1071').")
  }
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Install with install.packages('pROC').")
  }
  
  # Ensure binary target
  target_factor <- as.factor(target_vector)
  if (length(levels(target_factor)) != 2) {
    stop("Target vector must have exactly 2 levels for binary classification.")
  }
  positive_class <- levels(target_factor)[2]  # We'll treat the 2nd level as "positive"
  
  # Fit polynomial-kernel SVM (all columns of df as predictors)
  svm_model <- e1071::svm(
    x = df,
    y = target_factor,
    kernel = "polynomial",
    probability = TRUE
  )
  
  # Predict probabilities
  pred_class <- predict(svm_model, df, probability = TRUE)
  # Extract the probability of the positive class
  prob_matrix <- attr(pred_class, "probabilities")
  prob_positive <- prob_matrix[, positive_class]
  
  # Calculate ROC
  roc_obj <- tryCatch(
    pROC::roc(response = target_vector, predictor = prob_positive),
    error = function(e) {
      stop("Error calculating ROC: ", e$message)
    }
  )
  
  # Extract AUC and CI
  auc_value <- roc_obj$auc
  ci_vals   <- pROC::ci.auc(roc_obj, conf.level = ci_level)
  ci_lower  <- ci_vals[1]
  ci_upper  <- ci_vals[3]
  
  # Approximate SE from the 95% CI width
  # 3.92 ~ 2 * 1.96 for a 95% CI normal approximation
  std_error <- (ci_upper - ci_lower) / 3.92
  
  # Build a data frame for each threshold in the ROC
  roc_df <- data.frame(
    identifier   = "Combined (SVM)",  # or customize as needed
    auc          = as.numeric(auc_value),
    ci_lower     = as.numeric(ci_lower),
    ci_upper     = as.numeric(ci_upper),
    std_error    = as.numeric(std_error),
    sensitivity  = roc_obj$sensitivities,
    specificity  = 1 - roc_obj$specificities
  )
  
  return(roc_df)
}

#' performANOVA6
#' @description
#' Takes in 6 data frames, each belonging to a cohort, and performs an ANOVA on each row.
#' The function returns a data frame with the p-value from the ANOVA, Tukey's HSD p-values
#' for pairwise comparisons, and fold changes for each pairwise comparison.
#'
#' @param df1 Data frame containing cohort 1
#' @param df2 Data frame containing cohort 2
#' @param df3 Data frame containing cohort 3
#' @param df4 Data frame containing cohort 4
#' @param df5 Data frame containing cohort 5
#' @param df6 Data frame containing cohort 6
#'
#' @return A data frame with the results of the ANOVA, Tukey’s HSD p-values, and fold changes
#' @export
#'
#' @examples
#' test <- performANOVA6(df1, df2, df3, df4, df5, df6)
performANOVA6 <- function(df1, df2, df3, df4, df5, df6) {
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(stats) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  # Number of rows (assuming df1,...,df6 have the same number of rows)
  n_rows <- nrow(df1)
  
  # Store overall ANOVA p-values
  pvalues <- numeric(n_rows)
  
  # Store Tukey's p-values for all 15 pairwise comparisons
  tukey_pvalue12 <- numeric(n_rows)
  tukey_pvalue13 <- numeric(n_rows)
  tukey_pvalue14 <- numeric(n_rows)
  tukey_pvalue15 <- numeric(n_rows)
  tukey_pvalue16 <- numeric(n_rows)
  tukey_pvalue23 <- numeric(n_rows)
  tukey_pvalue24 <- numeric(n_rows)
  tukey_pvalue25 <- numeric(n_rows)
  tukey_pvalue26 <- numeric(n_rows)
  tukey_pvalue34 <- numeric(n_rows)
  tukey_pvalue35 <- numeric(n_rows)
  tukey_pvalue36 <- numeric(n_rows)
  tukey_pvalue45 <- numeric(n_rows)
  tukey_pvalue46 <- numeric(n_rows)
  tukey_pvalue56 <- numeric(n_rows)
  
  # Store fold changes for all 15 pairwise comparisons
  foldchange12 <- numeric(n_rows)
  foldchange13 <- numeric(n_rows)
  foldchange14 <- numeric(n_rows)
  foldchange15 <- numeric(n_rows)
  foldchange16 <- numeric(n_rows)
  foldchange23 <- numeric(n_rows)
  foldchange24 <- numeric(n_rows)
  foldchange25 <- numeric(n_rows)
  foldchange26 <- numeric(n_rows)
  foldchange34 <- numeric(n_rows)
  foldchange35 <- numeric(n_rows)
  foldchange36 <- numeric(n_rows)
  foldchange45 <- numeric(n_rows)
  foldchange46 <- numeric(n_rows)
  foldchange56 <- numeric(n_rows)
  
  # Loop over each row
  tryCatch({
    for (i in seq_len(n_rows)) {
      values_df1 <- as.numeric(df1[i, ])
      values_df2 <- as.numeric(df2[i, ])
      values_df3 <- as.numeric(df3[i, ])
      values_df4 <- as.numeric(df4[i, ])
      values_df5 <- as.numeric(df5[i, ])
      values_df6 <- as.numeric(df6[i, ])
      
      # Combine all values
      values <- c(values_df1, values_df2, values_df3, values_df4, values_df5, values_df6)
      groups <- factor(c(
        rep("Group1", length(values_df1)),
        rep("Group2", length(values_df2)),
        rep("Group3", length(values_df3)),
        rep("Group4", length(values_df4)),
        rep("Group5", length(values_df5)),
        rep("Group6", length(values_df6))
      ))
      
      # Check for >=3 non-NA values in each group to run ANOVA
      enough_data <- (
        sum(!is.na(values_df1)) >= 3 &&
          sum(!is.na(values_df2)) >= 3 &&
          sum(!is.na(values_df3)) >= 3 &&
          sum(!is.na(values_df4)) >= 3 &&
          sum(!is.na(values_df5)) >= 3 &&
          sum(!is.na(values_df6)) >= 3
      )
      
      if (enough_data) {
        # Perform ANOVA
        anova_result <- aov(values ~ groups)
        pvalues[i] <- summary(anova_result)[[1]][["Pr(>F)"]][1]
        
        # Perform Tukey's HSD test
        # (Checking "if (pvalues[i] <= 1)" always succeeds unless p-value is NA;
        #  presumably you'd check p < 0.05 if you only want post-hoc for significant ANOVA)
        if (!is.na(pvalues[i]) && pvalues[i] <= 1) {
          tukey_result <- TukeyHSD(anova_result)
          # Access the matrix of pairwise comparisons
          comp <- tukey_result$groups
          
          # Extract Tukey p-values for pairwise comparisons
          tukey_pvalue12[i] <- comp["Group2-Group1", "p adj"]
          tukey_pvalue13[i] <- comp["Group3-Group1", "p adj"]
          tukey_pvalue14[i] <- comp["Group4-Group1", "p adj"]
          tukey_pvalue15[i] <- comp["Group5-Group1", "p adj"]
          tukey_pvalue16[i] <- comp["Group6-Group1", "p adj"]
          tukey_pvalue23[i] <- comp["Group3-Group2", "p adj"]
          tukey_pvalue24[i] <- comp["Group4-Group2", "p adj"]
          tukey_pvalue25[i] <- comp["Group5-Group2", "p adj"]
          tukey_pvalue26[i] <- comp["Group6-Group2", "p adj"]
          tukey_pvalue34[i] <- comp["Group4-Group3", "p adj"]
          tukey_pvalue35[i] <- comp["Group5-Group3", "p adj"]
          tukey_pvalue36[i] <- comp["Group6-Group3", "p adj"]
          tukey_pvalue45[i] <- comp["Group5-Group4", "p adj"]
          tukey_pvalue46[i] <- comp["Group6-Group4", "p adj"]
          tukey_pvalue56[i] <- comp["Group6-Group5", "p adj"]
          
          # Calculate fold changes = mean(GroupX) - mean(GroupY)
          foldchange12[i] <- mean(values_df2, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange13[i] <- mean(values_df3, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange14[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange15[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange16[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange23[i] <- mean(values_df3, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange24[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange25[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange26[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange34[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange35[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange36[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange45[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df4, na.rm = TRUE)
          foldchange46[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df4, na.rm = TRUE)
          foldchange56[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df5, na.rm = TRUE)
        }
      } else {
        # If not enough data for ANOVA, set p-values = 1, fold changes = 0
        pvalues[i]       <- 1
        tukey_pvalue12[i] <- 1
        tukey_pvalue13[i] <- 1
        tukey_pvalue14[i] <- 1
        tukey_pvalue15[i] <- 1
        tukey_pvalue16[i] <- 1
        tukey_pvalue23[i] <- 1
        tukey_pvalue24[i] <- 1
        tukey_pvalue25[i] <- 1
        tukey_pvalue26[i] <- 1
        tukey_pvalue34[i] <- 1
        tukey_pvalue35[i] <- 1
        tukey_pvalue36[i] <- 1
        tukey_pvalue45[i] <- 1
        tukey_pvalue46[i] <- 1
        tukey_pvalue56[i] <- 1
        
        foldchange12[i]  <- 0
        foldchange13[i]  <- 0
        foldchange14[i]  <- 0
        foldchange15[i]  <- 0
        foldchange16[i]  <- 0
        foldchange23[i]  <- 0
        foldchange24[i]  <- 0
        foldchange25[i]  <- 0
        foldchange26[i]  <- 0
        foldchange34[i]  <- 0
        foldchange35[i]  <- 0
        foldchange36[i]  <- 0
        foldchange45[i]  <- 0
        foldchange46[i]  <- 0
        foldchange56[i]  <- 0
      }
    }
  }, error = function(e) {
    # If an error occurs, fill the current row with safe defaults
    # (In a robust design, you'd handle or log the error more gracefully.)
    pvalues[i]       <- 1
    tukey_pvalue12[i] <- 1
    tukey_pvalue13[i] <- 1
    tukey_pvalue14[i] <- 1
    tukey_pvalue15[i] <- 1
    tukey_pvalue16[i] <- 1
    tukey_pvalue23[i] <- 1
    tukey_pvalue24[i] <- 1
    tukey_pvalue25[i] <- 1
    tukey_pvalue26[i] <- 1
    tukey_pvalue34[i] <- 1
    tukey_pvalue35[i] <- 1
    tukey_pvalue36[i] <- 1
    tukey_pvalue45[i] <- 1
    tukey_pvalue46[i] <- 1
    tukey_pvalue56[i] <- 1
    
    foldchange12[i]  <- 0
    foldchange13[i]  <- 0
    foldchange14[i]  <- 0
    foldchange15[i]  <- 0
    foldchange16[i]  <- 0
    foldchange23[i]  <- 0
    foldchange24[i]  <- 0
    foldchange25[i]  <- 0
    foldchange26[i]  <- 0
    foldchange34[i]  <- 0
    foldchange35[i]  <- 0
    foldchange36[i]  <- 0
    foldchange45[i]  <- 0
    foldchange46[i]  <- 0
    foldchange56[i]  <- 0
  })
  
  # Combine original data frames into one and append computed columns
  c_df <- dplyr::bind_cols(df1, df2, df3, df4, df5, df6) %>%
    dplyr::mutate(
      anova_pvalue   = pvalues,
      tukey_pvalue12 = tukey_pvalue12,  # Group1 vs Group2
      tukey_pvalue13 <- tukey_pvalue13,  # Group1 vs Group3
      tukey_pvalue14 <- tukey_pvalue14,  # Group1 vs Group4
      tukey_pvalue15 <- tukey_pvalue15,  # Group1 vs Group5
      tukey_pvalue16 <- tukey_pvalue16,  # Group1 vs Group6
      tukey_pvalue23 <- tukey_pvalue23,  # Group2 vs Group3
      tukey_pvalue24 <- tukey_pvalue24,  # Group2 vs Group4
      tukey_pvalue25 <- tukey_pvalue25,  # Group2 vs Group5
      tukey_pvalue26 <- tukey_pvalue26,  # Group2 vs Group6
      tukey_pvalue34 <- tukey_pvalue34,  # Group3 vs Group4
      tukey_pvalue35 <- tukey_pvalue35,  # Group3 vs Group5
      tukey_pvalue36 <- tukey_pvalue36,  # Group3 vs Group6
      tukey_pvalue45 <- tukey_pvalue45,  # Group4 vs Group5
      tukey_pvalue46 <- tukey_pvalue46,  # Group4 vs Group6
      tukey_pvalue56 <- tukey_pvalue56,  # Group5 vs Group6
      
      foldchange12 <- foldchange12,
      foldchange13   = foldchange13,
      foldchange14   = foldchange14,
      foldchange15   = foldchange15,
      foldchange16   = foldchange16,
      foldchange23   = foldchange23,
      foldchange24   = foldchange24,
      foldchange25   = foldchange25,
      foldchange26   = foldchange26,
      foldchange34   = foldchange34,
      foldchange35   = foldchange35,
      foldchange36   = foldchange36,
      foldchange45   = foldchange45,
      foldchange46   = foldchange46,
      foldchange56   = foldchange56
    )
  
  return(c_df)
}

#' Create Grouped Boxplot with Significance
#' @title Create Grouped Boxplot
#' @description
#' Creates a boxplot with beeswarm points and statistical significance annotations.
#' @param data A dataframe with 'group' and 'value' columns
#' @param alpha Significance threshold (default 0.05)
#' @return A ggplot2 object containing the boxplot
#' @export
#' @examples
#' # df <- data.frame(group=rep(c("A","B","C"),each=10), value=rnorm(30))
#' # plot <- create_grouped_boxplot(df, alpha=0.05)
create_grouped_boxplot <- function(data, alpha = 0.05) {
  # Check required packages
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("Please install ggbeeswarm: install.packages('ggbeeswarm')")
  }
  
  # Check required columns
  if (!all(c("group", "value") %in% colnames(data))) {
    stop("Data must contain 'group' and 'value' columns")
  }
  
  # Convert group to factor
  data$group <- factor(data$group)
  
  # Generate all possible group pairs
  all_pairs <- combn(levels(data$group), 2, simplify = FALSE)
  
  # Compute p-values for each pair
  pairwise_p <- purrr::map_dbl(all_pairs, ~{
    t.test(data$value[data$group == .x[1]], 
           data$value[data$group == .x[2]])$p.value
  })
  
  # Create pairwise data
  pairwise_data <- dplyr::tibble(
    pair = all_pairs,
    p_value = pairwise_p
  ) %>%
    dplyr::filter(p_value < alpha)
  
  # Return basic plot if no significant pairs
  if (nrow(pairwise_data) == 0) {
    warning("No significant pairs found at alpha = ", alpha)
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = group, y = value, fill = group)) +
        ggplot2::geom_boxplot() +
        ggbeeswarm::geom_beeswarm() +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       legend.position = "none")
    )
  }
  
  # Create annotations
  pairwise_data <- pairwise_data %>%
    dplyr::mutate(
      annotation = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*"
      )
    )
  
  # Calculate y positions
  y_range <- range(data$value, na.rm = TRUE)
  y_start <- max(data$value, na.rm = TRUE) + 0.2
  step <- (y_range[2] - y_range[1]) * 0.1
  y_positions <- y_start + (0:(nrow(pairwise_data) - 1)) * step
  
  # Create plot
  ggplot2::ggplot(data, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_boxplot(staplewidth = 0.5, outliers = FALSE) +
    ggbeeswarm::geom_quasirandom(aes(fill = group), shape = 21, color = "black", size = 2.2) +
    ggsignif::geom_signif(
      comparisons = pairwise_data$pair,
      annotations = pairwise_data$annotation,
      y_position = y_positions,
      tip_length = 0.01,
      map_signif_level = FALSE,
      textsize = 10,
      vjust = 0.6
    ) +
    ggplot2::ylim(c(min(data$value), max(y_positions) + step)) +
    ggprism::theme_prism(base_size = 18) +
    labs(x = NULL, y = "Log2Abundance") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      legend.position = "none"
    )
}

#' Perform Row-wise ANOVA, Extract Specific Tukey HSD Comparisons, and Calculate Log2 Fold Change
#'
#' This function iterates through rows of a combined data matrix, performs
#' ANOVA across groups defined by a grouping factor, runs Tukey's HSD post-hoc
#' test, calculates the log2 fold change (difference of means, assuming log-transformed input)
#' for user-specified comparisons, and extracts adjusted p-values.
#'
#' @param data_list A named list of numeric data frames or matrices. Each element
#'   represents a cohort (e.g., list(Q = Q_df, G = G_df, QZ = QZ_df, ...)).
#'   All data frames must have the same number of rows, corresponding to the
#'   same features (proteins) in the same order. Columns are replicates.
#'   Input data is assumed to be log-transformed (e.g., log2).
#' @param comparisons_list A character vector listing the desired comparisons
#'   in the format "Group1-Group2" (e.g., c("Q-G", "Q-QZ", "QZ-QAE")).
#'   The fold change will be calculated as mean(Group1) - mean(Group2).
#' @param feature_ids An optional vector of feature identifiers (e.g., protein IDs,
#'   gene names) corresponding to the rows of the data frames in `data_list`.
#'   If provided, these will be included in the output data frame. Must have
#'   the same length as the number of rows in the data frames.
#'
#' @return A data frame containing the results. If `feature_ids` were provided,
#'   the first column will be `FeatureID`. Subsequent columns correspond to the
#'   requested comparisons, with pairs of columns for log2 fold change (`_log2FC`)
#'   and adjusted p-value (`_p.adj`). Rows where analysis failed will have NA.
#' @export
#'
#' @examples
#' # --- Create Sample Data (simulate user's dataframes on log2 scale) ---
#' set.seed(42)
#' n_proteins <- 50
#' Q <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 10, sd = 1), nrow = n_proteins))
#' G <- as.data.frame(matrix(rnorm(n_proteins * 5, mean = 11, sd = 1), nrow = n_proteins)) # G higher mean (log2FC ~1)
#' QZ <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 10, sd = 1), nrow = n_proteins))
#' QAE <- as.data.frame(matrix(rnorm(n_proteins * 3, mean = 9, sd = 1), nrow = n_proteins)) # QAE lower mean than QZ (log2FC ~ -1)
#' GX <- as.data.frame(matrix(rnorm(n_proteins * 5, mean = 11, sd = 1), nrow = n_proteins)) # Same as G
#' GYM <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 12, sd = 1), nrow = n_proteins)) # GYM higher mean than GX (log2FC ~ -1 for GX-GYM)
#'
#' # Assign meaningful column names (optional but good practice)
#' colnames(Q) <- paste0("Q_", 1:ncol(Q))
#' colnames(G) <- paste0("G_", 1:ncol(G))
#' colnames(QZ) <- paste0("QZ_", 1:ncol(QZ))
#' colnames(QAE) <- paste0("QAE_", 1:ncol(QAE))
#' colnames(GX) <- paste0("GX_", 1:ncol(GX))
#' colnames(GYM) <- paste0("GYM_", 1:ncol(GYM))
#'
#' # Create the input list
#' data_for_analysis <- list(
#'   Q = Q, G = G, QZ = QZ, QAE = QAE, GX = GX, GYM = GYM
#' )
#'
#' # Define desired comparisons
#' my_comparisons <- c("Q-G", "Q-QZ", "QZ-QAE", "G-GX", "GX-GYM")
#'
#' # Generate dummy protein IDs
#' protein_ids <- paste0("Prot_", 1:n_proteins)
#'
#' # --- Run the function ---
#' analysis_results <- perform_rowwise_anova_tukey_fc(
#'   data_list = data_for_analysis,
#'   comparisons_list = my_comparisons,
#'   feature_ids = protein_ids
#' )
#'
#' # --- View Results ---
#' print(head(analysis_results))
#' # Note columns like Q-G_log2FC and Q-G_p.adj
#' summary(analysis_results)
#'
#' # Example: Filter for significant results with absolute log2FC > 1
#' # library(dplyr)
#' # significant_fc_Q_G <- analysis_results %>%
#' #   filter(`Q-G_p.adj` < 0.05 & abs(`Q-G_log2FC`) > 1)
#' # print(head(significant_fc_Q_G))

perform_rowwise_anova_tukey_fc <- function(data_list, comparisons_list, feature_ids = NULL) {
  
  # --- Input Validation ---
  if (!is.list(data_list) || is.data.frame(data_list)) {
    stop("Error: 'data_list' must be a named list of data frames or matrices.")
  }
  if (length(names(data_list)) != length(data_list) || any(names(data_list) == "")) {
    stop("Error: 'data_list' must be a *named* list.")
  }
  n_rows <- unique(sapply(data_list, nrow))
  if (length(n_rows) != 1) {
    stop("Error: All data frames in 'data_list' must have the same number of rows.")
  }
  # Check if data appears numeric (basic check on first data frame)
  if(!all(sapply(data_list[[1]], is.numeric))) {
    warning("Warning: Input data in the first element of 'data_list' does not appear fully numeric. Ensure all data frames contain only numeric values.")
  }
  if (!is.null(feature_ids) && length(feature_ids) != n_rows) {
    stop("Error: 'feature_ids' length must match the number of rows in the data frames.")
  }
  if(!is.character(comparisons_list)) {
    stop("Error: 'comparisons_list' must be a character vector.")
  }
  # Check comparison format
  if(!all(grepl("^[^-]+-[^-]+$", comparisons_list))) {
    stop("Error: 'comparisons_list' elements must be in the format 'Group1-Group2'.")
  }
  
  
  # --- Combine Data and Create Grouping Factor ---
  group_names <- names(data_list)
  combined_data <- do.call(cbind, data_list)
  # Ensure unique column names, necessary if original dfs had overlapping names
  colnames(combined_data) <- make.unique(unlist(sapply(data_list, colnames, USE.NAMES = FALSE)))
  
  group_vector <- factor(rep(group_names, times = sapply(data_list, ncol)))
  
  # --- Helper function to standardize comparison names (alphabetical order) ---
  # Used only for matching Tukey output, not for FC calculation direction
  standardize_comparison <- function(comp_str) {
    parts <- sort(strsplit(comp_str, "-", fixed = TRUE)[[1]])
    return(paste(parts, collapse = "-"))
  }
  
  # --- Function to process a single row ---
  process_row <- function(row_values, grouping_factor, target_comparisons_original) {
    
    # Create data frame for aov
    df_row <- data.frame(intensity = as.numeric(row_values), group = grouping_factor)
    df_row <- na.omit(df_row) # Remove NAs for this specific row
    
    # Check if enough data remains for ANOVA (at least 2 groups with >1 data point)
    group_counts <- table(df_row$group)
    valid_groups <- names(group_counts[group_counts >= 1]) # Need at least 1 point per group for mean calc
    valid_groups_anova <- names(group_counts[group_counts > 1]) # Need >1 point per group for ANOVA variance
    
    # Initialize result vector (interleaved p-val, fc)
    n_comp <- length(target_comparisons_original)
    result_vector <- rep(NA_real_, n_comp * 2)
    names(result_vector) <- character(n_comp * 2)
    
    # Calculate group means if possible
    group_means <- tapply(df_row$intensity, df_row$group, mean) # Will have NA for groups not present
    
    # --- Calculate Fold Changes ---
    fc_idx <- 2 # Start index for FC results in the interleaved vector
    for (comp in target_comparisons_original) {
      groups_to_compare <- strsplit(comp, "-", fixed = TRUE)[[1]]
      group1 <- groups_to_compare[1]
      group2 <- groups_to_compare[2]
      
      mean1 <- group_means[group1] # Returns NA if group1 not present in df_row
      mean2 <- group_means[group2] # Returns NA if group2 not present in df_row
      
      # Calculate FC only if both means are available
      if (!is.na(mean1) && !is.na(mean2)) {
        result_vector[fc_idx] <- mean1 - mean2
      } # else remains NA
      
      # Set names for FC column
      names(result_vector)[fc_idx] <- paste0(comp, "_log2FC")
      fc_idx <- fc_idx + 2 # Move to next FC slot
    }
    
    
    # --- Perform ANOVA and Tukey HSD (if possible) ---
    p_val_idx <- 1 # Start index for p-value results
    if (length(valid_groups_anova) < 2) {
      # Cannot perform ANOVA/Tukey, fill p-values with NA
      for (comp in target_comparisons_original) {
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        # result_vector[p_val_idx] is already NA
        p_val_idx <- p_val_idx + 2
      }
      return(result_vector) # Return vector with NAs for p-values and calculated FCs
    }
    
    # Try ANOVA/Tukey
    tukey_results <- tryCatch({
      aov_res <- aov(intensity ~ group, data = df_row)
      TukeyHSD(aov_res)
    }, error = function(e) {
      NULL # Return NULL if error occurs
    })
    
    # Extract Tukey p-values if successful
    if (!is.null(tukey_results) && "group" %in% names(tukey_results)) {
      tukey_table <- as.data.frame(tukey_results$group)
      tukey_comparisons_std <- sapply(rownames(tukey_table), standardize_comparison, USE.NAMES = FALSE)
      
      for (comp in target_comparisons_original) {
        req_std_comp <- standardize_comparison(comp) # Standardize requested comparison for matching
        match_idx <- which(tukey_comparisons_std == req_std_comp)
        
        if (length(match_idx) == 1) {
          result_vector[p_val_idx] <- tukey_table[match_idx, "p adj"] # Extract adjusted p-value
        } # else remains NA
        
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        p_val_idx <- p_val_idx + 2
      }
    } else {
      # ANOVA/Tukey failed, fill p-values with NA
      for (comp in target_comparisons_original) {
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        # result_vector[p_val_idx] is already NA
        p_val_idx <- p_val_idx + 2
      }
    }
    
    return(result_vector) # Return interleaved vector
    
  } # End process_row function
  
  # --- Apply process_row to each row ---
  results_matrix <- t(apply(combined_data, 1, process_row,
                            grouping_factor = group_vector,
                            target_comparisons_original = comparisons_list))
  results_df <- as.data.frame(results_matrix)
  
  # --- Add Feature IDs if provided ---
  if (!is.null(feature_ids)) {
    results_df <- cbind(FeatureID = feature_ids, results_df)
  }
  
  return(results_df)
  
} # End perform_rowwise_anova_tukey_fc function


