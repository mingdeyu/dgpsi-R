#' @title Initialize the Poisson likelihood
#'
#' @description This function constructs an object for Poisson likelihood.
#'
#' @param input_dim a vector of length one that contains the indices of one GP in the feeding
#'     layer whose outputs feed into the likelihood node. When set to `NULL`,
#'     all outputs from GPs of the feeding layer feed into the likelihood node, and in this case
#'     one needs to ensure there is only one GP node specified in the feeding layer..
#'     Defaults to `NULL`.
#'
#' @return An object to represent the Poissonlikelihood node.
#' @note The Poisson likelihood node only needs one feeding GP node.
#' @details See examples in tutorials at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
Poisson <- function(input_dim = NULL) {
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  res <- pkg.env$dgpsi$Poisson(input_dim)
  return(res)
}


#' @title Initialize the heteroskedastic Gaussian likelihood
#'
#' @description This function constructs an object for heteroskedastic Gaussian likelihood.
#'
#' @param input_dim a vector of length two that contains the indices of two GPs in the feeding
#'     layer whose outputs feed into the likelihood node. When set to `NULL`,
#'     all outputs from GPs of feeding layer feed into the likelihood node, and in this case
#'     one needs to ensure there are only two GP nodes specified in the feeding layer..
#'     Defaults to `NULL`.
#'
#' @return An object to represent the heteroskedastic Gaussian likelihood node.
#' @note The heteroskedastic Gaussian likelihood node only needs two feeding GP nodes.
#' @details See examples in tutorials at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
Hetero <- function(input_dim = NULL) {
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  res <- pkg.env$dgpsi$Hetero(input_dim)
  return(res)
}

#' @title Initialize the negative Binomial likelihood
#'
#' @description This function constructs an object for negative Binomial likelihood.
#'
#' @param input_dim a vector of length two that contains the indices of two GPs in the feeding
#'     layer whose outputs feed into the likelihood node. When set to `NULL`,
#'     all outputs from GPs of feeding layer feed into the likelihood node, and in this case
#'     one needs to ensure there are only two GP nodes specified in the feeding layer.
#'     Defaults to `NULL`.
#'
#' @return An object to represent the negative Binomial likelihood node.
#' @note The negative Binomial likelihood node only needs two feeding GP nodes.
#' @details See examples in tutorials at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
NegBin <- function(input_dim = NULL) {
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  res <- pkg.env$dgpsi$NegBin(input_dim)
  return(res)
}


#' @title Calculate negative predicted log-likelihood
#'
#' @description This function compute the negative predicted log-likelihood from a
#'     trained DGP with likelihood layer.
#'
#' @param obj a DGP object produced by [emulator()]
#' @param x a matrix where each row is an input testing data point and each column is an input dimension.
#' @param y a matrix with only one column where each row is a scalar-valued testing output data point.
#'
#' @return A named list with two components. The first one, named `meanNLL`, is a scalar that gives the average negative
#'     predicted log-likelihood across all testing data points. The second one, named `NLL`, is a vector that gives
#'     the negative predicted log-likelihood for each testing data point.
#'
#' @details See examples in tutorials at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export
nllik <- function(obj, x, y) {
  res <- obj$nllik(x, y)
  named_res <- list("meanNLL" = res[[1]], "NLL" = res[[2]])
  return(named_res)
}
