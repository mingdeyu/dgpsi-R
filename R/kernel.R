#' @title Initialize a Gaussian process node
#'
#' @description This function constructs a kernel object to represent properties of a
#'     Gaussian process node.
#'
#' @param length a vector of lengthscales. The length of the vector equals to:
#' 1. either one if the lengthscales in the kernel function are assumed same across input dimensions; or
#' 2. the total number of input dimensions, which is the sum of the number of feeding GP nodes
#'    in the last layer (defined by the argument `input_dim`) and the number of connected global
#'    input dimensions (defined by the argument `connect`), if the lengthscales in the kernel function
#'    are assumed different across input dimensions.
#' @param scale the variance of a GP node. Defaults to `1`.
#' @param nugget the nugget term of a GP node. Defaults to `1e-6`.
#' @param name kernel function to be used. Either `"sexp"` for squared exponential kernel or
#'     `"matern2.5"` for Matérn-2.5 kernel. Defaults to `"sexp"`.
#' @param prior_name  prior options for the lengthscales and nugget term: gamma prior (`"ga"`), inverse gamma prior (`"inv_ga"`),
#'     or jointly robust prior (`"ref"`) for the lengthscales and nugget term. Set `NULL` to disable the prior. Defaults to `"ga"`.
#' @param prior_coef a vector that contains the coefficients for different priors:
#' * for the gamma prior, it is a vector of two values specifying the shape and rate parameters of the gamma distribution. Set to `NULL` for the
#'   default value `c(1.6,0.3)`.
#' * for the inverse gamma prior, it is a vector of two values specifying the shape and scale parameters of the inverse gamma distribution. Set
#'   to `NULL` for the default value `c(1.6,0.3)`.
#' * for the jointly robust prior, it is a vector of a single value specifying the `a` parameter in the prior. Set to `NULL` for the
#'   default value `c(0.2)`. See the reference below for the jointly robust prior.
#'
#' Defaults to `NULL`.
#' @param bounds a vector of length two that gives the lower bound (the first element of the vector) and the upper bound (the second element of the
#'     vector) of all lengthscales of the GP node. Defaults to `NULL` where no bounds are specified for the lengthscales.
#' @param nugget_est set to `TRUE` to estimate the nugget term or to `FALSE` to fix the nugget term as specified
#'     by the argument `nugget`. If set to `TRUE`, the value set to the argument `nugget` is used as the initial
#'     value. Defaults to `FALSE`.
#' @param scale_est set to `TRUE` to estimate the variance (i.e., scale) or to `FALSE` to fix the variance (i.e., scale) as specified
#'     by the argument `scale`. Defaults to `FALSE`.
#' @param input_dim a vector that contains either
#' 1. the indices of GP nodes in the feeding layer whose outputs feed into this GP node; or
#' 2. the indices of global input dimensions that are linked to the outputs of some feeding emulators,
#'    if this GP node is in the first layer of a GP or DGP, which will be used for the linked emulation.
#'
#' When set to `NULL`,
#' 1. all outputs from the GP nodes in the feeding layer feed into this GP node; or
#' 2. all global input dimensions feed into this GP node.
#'
#' Defaults to `NULL`.
#' @param connect a vector that contains the indices of dimensions in the global
#'      input connecting to this GP node as additional input dimensions. When set to `NULL`, no global input
#'      connection is implemented. Defaults to `NULL`. When this GP node is in the first layer of a GP or DGP emulator,
#'      which will consequently be used for linked emulation, `connect` gives the indices of global input dimensions
#'      that are not connected to some feeding emulators. In such a case, set `input_dim` to a vector of indices of
#'      the remaining input dimensions that are connected to the feeding emulators.
#'
#' @return A 'python' object to represent a GP node.
#' @references
#' Gu, M. (2019). Jointly robust prior for Gaussian stochastic process in emulation, calibration and variable selection. *Bayesian Analysis*, **14(3)**, 857-885.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # Check https://mingdeyu.github.io/dgpsi-R/ for examples
#' # on how to customize DGP structures using kernel().
#' }
#' @md
#' @export
kernel <- function(length, scale = 1., nugget = 1e-6, name = 'sexp', prior_name = 'ga', prior_coef = NULL, bounds = NULL, nugget_est = FALSE, scale_est = FALSE, input_dim = NULL, connect = NULL) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( name!='sexp' & name!='matern2.5' ) stop("'name' can only be either 'sexp' or 'matern2.5'.", call. = FALSE)
  if ( !is.null(prior_name) & prior_name!='ga' & prior_name!='inv_ga' ) stop("The provided 'prior_name' is not supported.", call. = FALSE)

  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }

  if(!is.null(connect)){
    connect <- reticulate::np_array(as.integer(connect - 1))
  }

  if(!is.null(bounds)){
    bounds <- reticulate::np_array(bounds)
  }

  if(!is.null(prior_coef)){
    prior_coef <- reticulate::np_array(prior_coef)
  }

  res <- pkg.env$dgpsi$kernel(reticulate::np_array(length), scale, nugget, name, prior_name, prior_coef, bounds, nugget_est, scale_est, input_dim, connect)
  return(res)
}
