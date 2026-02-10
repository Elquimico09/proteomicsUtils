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