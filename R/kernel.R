#' @title Initialize the kernel object
#'
#' @description This function constructs an object to represent properties of a
#'     Gaussian process node.
#'
#' @param length a vector of lengthscales. The length of the vector equals to:
#' 1. either one if the lengthscales in the kernel function are assumed same across input dimensions; or
#' 2. the total number of input dimensions, which is the sum of the number of feeding GPs
#'    in the last layer (defined by the argument `input_dim`) and the number of connected global
#'    input dimensions (defined by the argument `connect`), if the lengthscales in the kernel function
#'    are assumed different across input dimensions.
#' @param scale the variance of a GP. Defaults to `1`.
#' @param nugget the nugget term of a GP. Defaults to `1e-6`.
#' @param name kernel function to be used. Either `"sexp"` for squared exponential kernel or
#'     `"matern2.5"` for Matern2.5 kernel. Defaults to `"sexp"`.
#' @param prior_name  prior options for the lengthscales and nugget term. Either gamma (`"ga"`) or inverse gamma (`"inv_ga"`) distribution for
#'      the lengthscales and nugget term. Set `NULL` to disable the prior. Defaults to `"ga"`.
#' @param prior_coef a vector that contains two values specifying the shape and rate
#'      parameters of gamma prior, shape and scale parameters of inverse gamma prior. Defaults to ``c(1.6,0.3)``.
#' @param nugget_est set to `TRUE` to estimate nugget term or to `FALSE` to fix the nugget term as specified
#'      by the argument `nugget`. If set to `TRUE`, the value set to the argument `nugget` is used as the initial
#'      value. Defaults to `FALSE`.
#' @param scale_est set to `TRUE` to estimate the variance or to `FALSE` to fix the variance as specified
#'      by the argument `scale`. Defaults to `FALSE`.
#' @param input_dim a vector that contains either
#' 1. the indices of GPs in the feeding layer whose outputs feed into the GP; or
#' 2. the indices of dimensions in the global input if the GP is in the first layer.
#'
#' When set to `NULL`,
#' 1. all outputs from GPs in the feeding layer; or
#' 2. all global input dimensions feed into the GP.
#'
#' Defaults to `NULL`.
#' @param connect a vector that contains the indices of dimensions in the global
#'      input connecting to the GP as additional input dimensions to the input obtained from the output of
#'      GPs in the feeding layer (as determined by the argument `input_dim`). When set to `NULL`, no global input
#'      connection is implemented. Defaults to `NULL`. When the [kernel()] function is used in GP/DGP emulators for linked
#'      emulation and some input dimensions to the computer models are not connected to some feeding computer models,
#'      set `connect` to a vector of indices of these external global input dimensions, and accordingly, set
#'      `input_dim` to a vector of indices of the remaining input dimensions that are connected to the feeding
#'      computer models.
#'
#' @return A kernel object to represent a GP node.
#' @md
#' @export
kernel <- function(length, scale = 1., nugget = 1e-6, name = 'sexp', prior_name = 'ga', prior_coef = c(1.6,0.3), nugget_est = FALSE, scale_est = FALSE, input_dim = NULL, connect = NULL) {

  if(!is.null(input_dim)){
    input_dim <- reticulate::np_array(as.integer(input_dim - 1))
  }

  if(!is.null(connect)){
    connect <- reticulate::np_array(as.integer(connect - 1))
  }

  res <- pkg.env$dgpsi$kernel(reticulate::np_array(length), scale, nugget, name, prior_name, reticulate::np_array(prior_coef), nugget_est, scale_est, input_dim, connect)
  return(res)
}
