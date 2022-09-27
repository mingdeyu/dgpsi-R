#' @title Initialize a Poisson likelihood node
#'
#' @description This function constructs a likelihood object to represent a Poisson likelihood node.
#'
#' @param input_dim a vector of length one that contains the indices of one GP node in the feeding
#'     layer whose outputs feed into this likelihood node. When set to `NULL`,
#'     all outputs from GP nodes in the feeding layer feed into this likelihood node, and in such a case
#'     one needs to ensure that only one GP node is specified in the feeding layer.
#'     Defaults to `NULL`.
#'
#' @return A 'python' object to represent a Poisson likelihood node.
#' @note The Poisson likelihood node can only be linked to one feeding GP node.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # Check https://mingdeyu.github.io/dgpsi-R/ for examples
#' # on how to customize DGP structures using Poisson().
#' }
#' @md
#' @export
Poisson <- function(input_dim = NULL) {
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  res <- pkg.env$dgpsi$Poisson(input_dim)
  return(res)
}


#' @title Initialize a heteroskedastic Gaussian likelihood node
#'
#' @description This function constructs a likelihood object to represent a heteroskedastic Gaussian likelihood node.
#'
#' @param input_dim a vector of length two that contains the indices of two GP nodes in the feeding
#'     layer whose outputs feed into this likelihood node. When set to `NULL`,
#'     all outputs from GP nodes in the feeding layer feed into this likelihood node, and in such a case
#'     one needs to ensure that only two GP nodes are specified in the feeding layer.
#'     Defaults to `NULL`.
#'
#' @return A 'python' object to represent a heteroskedastic Gaussian likelihood node.
#' @note The heteroskedastic Gaussian likelihood node can only be linked to two feeding GP nodes.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # Check https://mingdeyu.github.io/dgpsi-R/ for examples
#' # on how to customize DGP structures using Hetero().
#' }
#' @md
#' @export
Hetero <- function(input_dim = NULL) {
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  res <- pkg.env$dgpsi$Hetero(input_dim)
  return(res)
}

#' @title Initialize a negative Binomial likelihood node
#'
#' @description This function constructs a likelihood object to represent a negative Binomial likelihood node.
#'
#' @param input_dim a vector of length two that contains the indices of two GP nodes in the feeding
#'     layer whose outputs feed into this likelihood node. When set to `NULL`,
#'     all outputs from GP nodes in the feeding layer feed into this likelihood node, and in such a case
#'     one needs to ensure that only two GP nodes are specified in the feeding layer.
#'     Defaults to `NULL`.
#'
#' @return A 'python' object to represent a negative Binomial likelihood node.
#' @note The negative Binomial likelihood node can only be linked to two feeding GP nodes.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # Check https://mingdeyu.github.io/dgpsi-R/ for examples
#' # on how to customize DGP structures using NegBin().
#' }
#' @md
#' @export
NegBin <- function(input_dim = NULL) {
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  res <- pkg.env$dgpsi$NegBin(input_dim)
  return(res)
}
