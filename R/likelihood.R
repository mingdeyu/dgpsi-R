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
#' @return An object to represent the Poisson likelihood node.
#' @note The Poisson likelihood node only needs one feeding GP node.
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
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
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
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
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export
NegBin <- function(input_dim = NULL) {
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  res <- pkg.env$dgpsi$NegBin(input_dim)
  return(res)
}
