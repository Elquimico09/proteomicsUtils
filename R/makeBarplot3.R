#' makeBarplot3
#'
#' Backward-compatible wrapper around `makeBarplot()`. Works with outputs from
#' `performANOVA()` and `perform2WayANOVA()` because it detects dynamic
#' `tukey_pvalue_<group1>_<group2>` columns.
#'
#' @param df long dataframe containing the values
#' @param pvalue_df the wide dataframe containing p-values from ANOVA results
#' @param cohort_labels vector containing the cohort labels
#' @param gene the gene to be plotted
#' @param tip_length Numeric tip length for significance brackets (default 0.04).
#' @param label_size Numeric label size for significance annotations (default 5).
#'
#' @return ggplot2 object containing the barplot
#' @export
#'
#' @examples
#' \dontrun{
#'   makeBarplot3(df_long, anova_results, c("A", "B", "C"), "Gene1", tip_length = 0.03, label_size = 4.5)
#' }
makeBarplot3 <- function(df, pvalue_df, cohort_labels, gene, tip_length = 0.04, label_size = 5) {
  makeBarplot(
    df = df,
    pvalue_df = pvalue_df,
    cohort_labels = cohort_labels,
    gene = gene,
    tip_length = tip_length,
    label_size = label_size
  )
}
