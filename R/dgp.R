#' @title Initialize the DGP object
#'
#' @description This function constructs an object for deep Gaussian process training.
#'
#' @param X a matrix where each row is an input training data point and each column is an input dimension.
#' @param Y a matrix containing observed training output data..
#'     The matrix has it rows being output data points and columns being output dimensions
#'     (with the number of columns equals to the number of GP nodes in the final layer).
#' @param all_layer a list contains *L* (the number of layers) sub-lists, each of which contains
#'     the GPs defined by the [kernel()] function in that layer. The sub-lists are placed in the list
#'     in the same order of the specified DGP model. The final layer of DGP hierarchy can be set to a likelihood
#'     layer by putting an object created by a likelihood function (e.g., [Poisson()]) into the final sub-list of `all_layer`.
#'     Defaults to `NULL`. If a DGP structure is not provided, an input-connected two-layered DGP structure
#'     (for deterministic model emulation) with the number of GP nodes in the first layer equal to the dimension
#'     of `X` is automatically constructed.
#' @param check_rep whether to check the repetitions in the dataset, i.e., if one input
#'     position has multiple outputs. Defaults to `TRUE`.
#' @param rff whether to use random Fourier features to approximate the correlation matrices
#'     during the imputation in training. Defaults to `FALSE`.
#' @param M the number of features to be used by random Fourier approximation. It is only used
#'     when `rff` is set to `TRUE`. Defaults to `NULL`. If it is not specified, `M` is set to
#'     `max(100, ceil(sqrt(Data Size)*log(Data Size))))`.
#'
#' @return A DGP object to be used by [train()] for DGP training.
#' @details See examples in tutorials at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
dgp <- function(X, Y, all_layer = NULL, check_rep = TRUE, rff = FALSE, M = NULL) {
  if(!is.null(M)){
    M <- as.integer(M)
  }
  res <- pkg.env$dgpsi$dgp(X, Y, all_layer, check_rep, rff, M)
  return(res)
}
