# Function: makeBoxplot
# Takes a long dataframe for plotting and a wide ANOVA results dataframe
# containing Tukey p-values from performANOVA/perform2WayANOVA.
#
# Args:
# @param df: dataframe containing the data
# @param pvalue_df: dataframe containing ANOVA/Tukey p-values
# @param cohort_labels: vector containing cohort labels (plot order)
# @param gene: gene name to plot
#
# Returns:
# @return: ggplot2 object containing the boxplot
#' Create Boxplot with Significant ANOVA Comparisons
#' @title Make Boxplot
#' @description
#' Creates a boxplot for a specific gene and adds significance brackets only for
#' significant Tukey HSD pairwise comparisons from `performANOVA()` or
#' `perform2WayANOVA()` output.
#' @param df A long-format dataframe with at least `Treatment` and `Abundance`
#'   columns. If a `Gene` column is present, rows are filtered by `gene`.
#' @param pvalue_df Output dataframe from `performANOVA()` or
#'   `perform2WayANOVA()` containing `tukey_pvalue_<group1>_<group2>` columns.
#' @param cohort_labels A character vector of cohort labels in plotting order.
#' @param gene The gene name to plot.
#' @param tip_length Numeric tip length for significance brackets (default 0.04).
#' @param label_size Numeric label size for significance annotations (default 5).
#' @return A ggplot2 object containing the boxplot.
#' @export
#' @examples
#' # p <- makeBoxplot(gene_long, anova_res, c("control", "treated", "placebo"), "ACTB", tip_length = 0.03, label_size = 4.5)
makeBoxplot <- function(df, pvalue_df, cohort_labels, gene, tip_length = 0.04, label_size = 5) {
  if (!all(c("Treatment", "Abundance") %in% names(df))) {
    stop("`df` must contain `Treatment` and `Abundance` columns.")
  }

  if ("Gene" %in% names(df)) {
    df <- df[df$Gene == gene, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop("No rows available to plot for the requested gene.")
  }

  df <- df[df$Treatment %in% cohort_labels, , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No rows in `df` matched the provided `cohort_labels`.")
  }

  df$Treatment <- factor(df$Treatment, levels = cohort_labels)

  if ("Gene" %in% names(pvalue_df)) {
    p_rows <- pvalue_df[pvalue_df$Gene == gene, , drop = FALSE]
  } else if (!is.null(rownames(pvalue_df))) {
    p_rows <- pvalue_df[rownames(pvalue_df) == gene, , drop = FALSE]
  } else {
    stop("Could not identify gene rows in `pvalue_df`. Add a `Gene` column or rownames.")
  }

  if (nrow(p_rows) == 0) {
    stop("Requested gene was not found in `pvalue_df`.")
  }

  p_row <- p_rows[1, , drop = FALSE]

  p_to_signif <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    ""
  }

  pair_mat <- utils::combn(cohort_labels, 2)
  sig_rows <- list()

  for (j in seq_len(ncol(pair_mat))) {
    g1 <- pair_mat[1, j]
    g2 <- pair_mat[2, j]

    col12 <- paste0("tukey_pvalue_", g1, "_", g2)
    col21 <- paste0("tukey_pvalue_", g2, "_", g1)
    col12_safe <- make.names(col12)
    col21_safe <- make.names(col21)

    pval <- NA_real_
    if (col12 %in% names(p_row)) {
      pval <- as.numeric(p_row[[col12]][1])
    } else if (col21 %in% names(p_row)) {
      pval <- as.numeric(p_row[[col21]][1])
    } else if (col12_safe %in% names(p_row)) {
      pval <- as.numeric(p_row[[col12_safe]][1])
    } else if (col21_safe %in% names(p_row)) {
      pval <- as.numeric(p_row[[col21_safe]][1])
    }

    if (!is.na(pval) && pval < 0.05) {
      sig_rows[[length(sig_rows) + 1]] <- data.frame(
        group1 = g1,
        group2 = g2,
        p.adj = pval,
        p.adj.signif = p_to_signif(pval),
        stringsAsFactors = FALSE
      )
    }
  }

  min_y <- min(df$Abundance, na.rm = TRUE)
  max_y <- max(df$Abundance, na.rm = TRUE)
  y_range <- max_y - min_y
  if (is.na(y_range) || y_range == 0) y_range <- 1

  y_step <- max(0.08 * y_range, 0.25)
  y_base <- max_y + (0.08 * y_range)

  if (length(sig_rows) > 0) {
    pvalue_ann <- do.call(rbind, sig_rows)
    pvalue_ann$y.position <- y_base + (seq_len(nrow(pvalue_ann)) - 1) * y_step
    ymax <- max(pvalue_ann$y.position) + y_step
  } else {
    pvalue_ann <- NULL
    ymax <- max_y + y_step
  }

  ymin <- min_y - (0.05 * y_range)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Treatment, y = Abundance)) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = Treatment),
      width = 0.6,
      alpha = 0.6,
      color = "black",
      outlier.shape = NA
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(fill = Treatment),
      shape = 21,
      color = "black",
      size = 1.5,
      width = 0.2,
      alpha = 0.75
    ) +
    ggprism::theme_prism(base_size = 8) +
    ggplot2::labs(title = gene, x = NULL, y = "Log2Abundance") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggplot2::coord_cartesian(ylim = c(ymin, ymax)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(legend.position = "none")

  if (exists("colors", inherits = TRUE)) {
    color_values <- get("colors", inherits = TRUE)
    if (length(color_values) >= length(cohort_labels)) {
      p <- p + ggplot2::scale_fill_manual(values = color_values[seq_along(cohort_labels)])
    }
  }

  if (!is.null(pvalue_ann) && nrow(pvalue_ann) > 0) {
    p <- p + ggprism::add_pvalue(
      pvalue_ann,
      label = "p.adj.signif",
      y.position = "y.position",
      tip.length = tip_length,
      bracket.size = 0.3,
      label.size = label_size
    )
  }

  return(p)
}
