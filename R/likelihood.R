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
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
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
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
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
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  res <- pkg.env$dgpsi$NegBin(input_dim)
  return(res)
}

#' @title Initialize a Categorical likelihood node
#'
#' @description
#'
#' `r new_badge("new")`
#'
#' This function constructs a likelihood object to represent a categorical likelihood node.
#'
#' @param K an integer specifying the number of classes in the training data. If `K` is set to `NULL` (the default),
#' its value will be inferred from the training data when running [dgp()]. Defaults to `NULL`.
#' @param input_dim a vector that specifies which GP node outputs feed into this likelihood node. It should be either:
#' * A vector of length `1` (for binary output), or
#' * A vector of length `K` (for `K`-class output), containing the indices of one or `K` GP nodes in the feeding layer.
#'
#' If `input_dim` is set to `NULL` (the default), all GP node outputs from the feeding layer feed into the likelihood node.
#' In this case, you must ensure that the feeding layer contains exactly one or `K` GP nodes. Defaults to `NULL`.
#'
#' @return A 'python' object to represent a categorical likelihood node.
#' @note The categorical likelihood node can only be linked to one or `K` feeding GP nodes.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # Check https://mingdeyu.github.io/dgpsi-R/ for examples
#' # on how to customize DGP structures using Categorical().
#' }
#' @md
#' @export
Categorical <- function(K = NULL, input_dim = NULL) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }
  if(!is.null(K)){
    K <- reticulate::np_array(as.integer(K))
  }
  res <- pkg.env$dgpsi$Categorical(K, input_dim)
  return(res)
}
