##' Visualize Sample Intensity Distributions
##'
##' Creates a violin/boxplot visualization of per-sample intensity distributions.
##' Optionally applies a log2 transform and infers groups from sample name prefixes.
##'
##' @param data A data.frame containing sample intensity columns.
##' @param log_transform Logical; if TRUE (default) apply log2(x + 1).
##' @param sample_cols Character vector of sample column names. If NULL, uses
##'   numeric columns from column 4 onward.
##'
##' @return A ggplot object.
##'
##' @examples
##' # visualizeDist(df)
##' # visualizeDist(df, log_transform = FALSE, sample_cols = c("control1", "control2"))
##'
##' @importFrom dplyr %>% mutate case_when
##' @importFrom tidyr pivot_longer
##' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot scale_fill_manual
##' @importFrom ggplot2 labs theme_minimal theme element_text
##' @export
visualizeDist <- function(data, log_transform = TRUE, sample_cols = NULL) {
  # If sample_cols not specified, assume numeric columns after metadata (column 4 onwards)
  if (is.null(sample_cols)) {
    sample_cols <- colnames(data)[4:ncol(data)]
    # Keep only numeric columns
    sample_cols <- sample_cols[sapply(data[sample_cols], is.numeric)]
  }
  
  # Extract sample data
  sample_data <- data[, sample_cols, drop = FALSE]
  
  # Optional log2 transformation (common for proteomics intensity data)
  if (log_transform) {
    sample_data <- log2(sample_data + 1)
    y_label <- "Log2(Intensity + 1)"
  } else {
    y_label <- "Intensity"
  }
  
  # Pivot to long format for ggplot
  long_data <- sample_data %>%
    pivot_longer(cols = everything(), 
                 names_to = "Sample", 
                 values_to = "Intensity") %>%
    mutate(Sample = factor(Sample, levels = sample_cols))  # Preserve order
  
  # Extract group from sample name (control, chronic, acute)
  long_data <- long_data %>%
    mutate(Group = case_when(
      grepl("^control", Sample) ~ "Control",
      grepl("^chronic", Sample) ~ "Chronic",
      grepl("^acute", Sample) ~ "Acute",
      TRUE ~ "Other"
    ))
  
  # Create violin plot
  p <- ggplot(long_data, aes(x = Sample, y = Intensity, fill = Group)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    scale_fill_manual(values = c("Control" = "#4DAF4A", 
                                  "Chronic" = "#E41A1C", 
                                  "Acute" = "#377EB8",
                                  "Other" = "grey50")) +
    labs(x = "Sample", y = y_label, title = "Sample Intensity Distribution") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.position = "top")
  
  return(p)
}
