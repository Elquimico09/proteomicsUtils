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