#' @title Initialize the GP object
#'
#' @description This function constructs an object for Gaussian process training.
#'
#' @param X a matrix where each row is an input data point and each column is an input dimension.
#' @param Y a matrix with only one column and each row being an input data point.
#' @param kernel a kernel object produced by the [kernel()] function.
#'
#' @return A GP object to be used by [train()] for GP training.
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
gp <- function(X, Y, kernel) {
  res <- pkg.env$dgpsi$gp(X, Y, kernel)
  return(res)
}
