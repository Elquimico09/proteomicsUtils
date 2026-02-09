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
