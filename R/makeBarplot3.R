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
