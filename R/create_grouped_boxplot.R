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