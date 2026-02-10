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
