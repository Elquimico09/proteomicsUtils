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
#'   Set to NULL to disable p-value filtering.
#' @param log2fc_cutoff Numeric, absolute log2 fold change cutoff for filtering rows (default: 1).
#'   Set to NULL to disable log2 fold change filtering.
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
                           p_cutoff = NULL, log2fc_cutoff = NULL,
                           col = circlize::colorRamp2(c(-2, 0, 2), c("green", "black", "red")),
                           width = NULL, height = NULL) {
  required_cols <- c(
    if (!is.null(p_cutoff)) "pvalue",
    if (!is.null(log2fc_cutoff)) "log2fc"
  )
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste0("Input df is missing required column(s): ",
                paste(missing_cols, collapse = ", "),
                ". Run convertFormat() first."))
  }

  # Apply only the requested filters.
  if (!is.null(p_cutoff) && !is.null(log2fc_cutoff)) {
    df <- df %>%
      dplyr::filter(pvalue < p_cutoff & abs(log2fc) >= log2fc_cutoff)
  } else if (!is.null(p_cutoff)) {
    df <- df %>%
      dplyr::filter(pvalue < p_cutoff)
  } else if (!is.null(log2fc_cutoff)) {
    df <- df %>%
      dplyr::filter(abs(log2fc) >= log2fc_cutoff)
  }

  # if the df has a gene column, set it as the rownames
  if ("Gene" %in% colnames(df)) {
    rownames <- df$Gene
  }

  df <- df %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-dplyr::any_of(c("log2fc", "pvalue", "foldchange", "nlog10p"))) %>%
    dplyr::select(!tidyselect::contains("foldchange")) %>%
    dplyr::select(!tidyselect::contains("pvalue")) %>%
    dplyr::select(!tidyselect::contains("nlog10p"))
  require(ComplexHeatmap)
  colnames <- colnames(df)
  if (!is.null(rownames)) {
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
                                  width = width, height = height,, border_gp = gpar(col = "black"),
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
                                  width = width, height = height, border_gp = gpar(col = "black", lwd = 2),
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
