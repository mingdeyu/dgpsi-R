#' @title Training of GP and DGP
#'
#' @description This function trains a GP or DGP model constructed by [gp()] or [dgp()] function.
#'
#' @param obj a GP or DGP model produced by [gp()] or [dgp()] function.
#' @param N number of iterations for stochastic EM. Defaults to `500`. Only used if `obj` is a DGP object.
#' @param ess_burn number of burnin steps for the ESS-within-Gibbs
#'     at each I-step of the SEM. Defaults to `10`. Only used if `obj` is a DGP object.
#' @param disable (bool, optional): whether to disable the training progress bar.
#'     Defaults to `FALSE`. Only used if `obj` is a DGP object.
#'
#' @return Updated GP or DGP object to be used for GP or DGP predictions.
#' @note The function can be re-applied to the returned trained DGP object to continue
#'     the training of the DGP model with additional number of SEM iterations.
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export

train <- function(obj, N = 500, ess_burn = 10, disable = FALSE) {
  N <- as.integer(N)
  ess_burn <- as.integer(ess_burn)

  if (pkg.env$py_buildin$type(obj)$'__name__' == 'gp'){
    obj$train()
  } else if (pkg.env$py_buildin$type(obj)$'__name__' == 'dgp'){
    obj$train(N, ess_burn, disable)
  }
  return(obj)
}

#' @title Estimating the final DGP model
#'
#' @description This function estimates the final DGP model from the model trained using [train()] function.
#'
#' @param obj the DGP object produced by [train()] function.
#' @param burnin the number of SEM iterations to be discarded for
#'     point estimate calculation. Must be smaller than the SEM iterations
#'     implemented. If this is not specified, only the last 25% of iterations
#'     are used. Defaults to `NULL`.
#'
#' @return a DGP object to be used to construct the emulator by [emulator()] function.
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export

estimate <- function(obj, burnin = NULL) {
  if(!is.null(burnin)){
    burnin <- as.integer(burnin)
  }
  res <- obj$estimate(burnin)
  return(res)
}

#' @title Export the trained GP model
#'
#' @description This function exports the GP model trained using [train()] to be used for
#'    linked Gaussian process constructions.
#'
#' @param obj the GP model produced by [train()] function.
#'
#' @return A GP object representing the trained GP.
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export

export <- function(obj) {
  res <- obj$export()
  return(res)
}

#' @title Plot of DGP model parameter traces
#'
#' @description This function plots the traces of model parameters of a particular GP node
#'     in the trained DGP hierarchy.
#'
#' @param obj The DGP object produced by [train()] function.
#' @param layer_no the index of the interested layer.
#' @param ker_no the index of the interested GP in the layer specified by `layer_no`.
#' @param width the overall plot width. Defaults to `4`.
#' @param height the overall plot height. Defaults to `1`.
#' @param ticksize the size of sub-plot ticks. Defaults to `5`.
#' @param labelsize the font size of y labels. Defaults to `8`.
#' @param hspace the space between sub-plots. Defaults to `0.1`.
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R>.
#' @md
#' @export

trace_plot <- function(obj, layer_no, ker_no, width = 4., height = 1., ticksize = 5., labelsize = 8., hspace = 0.1) {
  obj$plot(as.integer(layer_no), as.integer(ker_no), width, height, ticksize, labelsize, hspace)
}


