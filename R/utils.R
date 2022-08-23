#' @title Combine layers
#'
#' @description This function combines customized layers into a DGP or linked (D)GP structure.
#'
#' @param ... a sequence of lists:
#' 1. For DGP emulations, each list represents a DGP layer and contains GP nodes (produced by [kernel()]), or
#'    likelihood nodes (produced by [Poisson()], [Hetero()], or [NegBin()]).
#' 2. For linked (D)GP emulations, each list represents a system layer and contains emulators (produced by [gp()] or
#'    [dgp()]) in that layer.
#'
#' @return A list defining a DGP structure (for `struc` of [dgp()]) or a linked (D)GP structure
#'     (for `struc` for [lgp()]).
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export
combine <- function(...) {
  res = list(...)
  return(res)
}


#' @title Save the constructed emulator
#'
#' @description This function saves the constructed emulator to a `.pkl` file.
#'
#' @param object an instance of the S3 class `gp`, `dgp`, or `lgp`.
#' @param pkl_file the path to and the name of the `.pkl` file to which
#'     the emulator `object` is saved.
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export
write <- function(object, pkl_file) {
  lst <- unclass(object)
  pkg.env$dgpsi$write(lst, pkl_file)
}


#' @title Load the stored emulator
#'
#' @description This function loads the `.pkl` file that stores the emulator.
#'
#' @param pkl_file the path to and the name of the `.pkl` file where the emulator is stored.
#'
#' @return A GP, DGP or linked (D)GP emulator S3 class.
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export
read <- function(pkl_file) {
  res <- pkg.env$dgpsi$read(pkl_file)
  type <- pkg.env$py_buildin$type(res$emulator_obj)$'__name__'
  if ( type=='emulator' ) {
    class(res) <- "dgp"
  } else if ( type=='gp' ) {
    class(res) <- "gp"
  } else if ( type=='lgp' ) {
    class(res) <- "lgp"
  }
  return(res)
}


#' @title Summary of a constructed GP, DGP, or linked (D)GP emulator
#'
#' @description This function summarizes key information of a GP, DGP or linked (D)GP emulator.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `lgp`.
#' @param ... N/A.
#'
#' @return A table summarizing key information contained in `object`.
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @name summary
NULL

#' @rdname summary
#' @method summary gp
#' @export
summary.gp <- function(object, ...) {
  pkg.env$dgpsi$summary(object$emulator_obj, 'pretty')
}

#' @rdname summary
#' @method summary dgp
#' @export
summary.dgp <- function(object, ...) {
  pkg.env$dgpsi$summary(object$emulator_obj, 'pretty')
}

#' @rdname summary
#' @method summary lgp
#' @export
summary.lgp <- function(object, ...) {
  pkg.env$dgpsi$summary(object$emulator_obj, 'pretty')
}


#' @title Set linked indices
#'
#' @description This function sets the linked indices of a GP or DGP emulator if the information is not provided
#'     when the emulator is constructed by [gp()] or [dgp()].
#'
#' @param object an instance of the S3 class `gp` or `dgp`.
#' @param idx indices of columns in the pooled output matrix (formed by column-combining outputs of all emulators
#'     in the feeding layer) that will feed into the GP or DGP emulator represented by `object`.
#'
#' @return An updated `object` with the information of `idx` incorporated.
#'
#' @note This function is useful when different models are emulated by different teams. Each team can create their (D)GP emulator
#'     even without knowing how different emulators are connected together. When this information is available and
#'     different emulators are collected, the connection information between emulators can then be assigned to
#'     individual emulators with this function.
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export
set_linked_idx <- function(object, idx) {
  idx <- reticulate::np_array(as.integer(idx - 1))
  object$container_obj$set_local_input(idx)
  return(object)
}


#' @title Calculate negative predicted log-likelihood
#'
#' @description This function computes the negative predicted log-likelihood from a
#'     DGP emulator with a likelihood layer.
#'
#' @param object an instance of the `dgp` class and when it should be produced by [dgp()] with one of the following two settings:
#' 1. if `struc = NULL`, `likelihood` is not `NULL`;
#' 2. if a customized structure is provided to `struc`, the final layer must be likelihood layer containing only one
#'    likelihood node produced by [Poisson()], [Hetero()], or [NegBin()].
#' @param x a matrix where each row is an input testing data point and each column is an input dimension.
#' @param y a matrix with only one column where each row is a scalar-valued testing output data point.
#'
#' @return An updated `object` is returned with an additional slot named `NLL` that contains two elements.
#'     The first one, named `meanNLL`, is a scalar that gives the average negative predicted log-likelihood
#'     across all testing data points. The second one, named `allNLL`, is a vector that gives the negative predicted
#'     log-likelihood for each testing data point.
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export
nllik <- function(object, x, y) {
  if ( class(object)!='dgp' ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( !is.matrix(x) ) stop("x must be a matrix", call. = FALSE)
  if ( !is.matrix(y) ) stop("y must be a matrix", call. = FALSE)
  if ( nrow(x)!=nrow(y) ) stop("x and y have different number of rows.", call. = FALSE)
  res <- object$emulator_obj$nllik(x, y)
  named_res <- list("meanNLL" = res[[1]], "allNLL" = res[[2]])
  object$NLL <- named_res
  return(object)
}


#' @title Plot of DGP model parameter traces
#'
#' @description This function plots the traces of model parameters of a particular GP node
#'     in a trained DGP emulator.
#'
#' @param object The DGP emulator produced by [dgp()] function.
#' @param layer_no the index of the interested layer.
#' @param ker_no the index of the interested GP in the layer specified by `layer_no`.
#' @param width the overall plot width. Defaults to `4`.
#' @param height the overall plot height. Defaults to `1`.
#' @param ticksize the size of sub-plot ticks. Defaults to `5`.
#' @param labelsize the font size of y labels. Defaults to `8`.
#' @param hspace the space between sub-plots. Defaults to `0.1`.
#'
#' @details See examples in Articles at <https://mingdeyu.github.io/dgpsi-R/>.
#' @md
#' @export

trace_plot <- function(object, layer_no, ker_no, width = 4., height = 1., ticksize = 5., labelsize = 8., hspace = 0.1) {
  if ( class(object)!='dgp' ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  object$trained_obj$plot(as.integer(layer_no), as.integer(ker_no), width, height, ticksize, labelsize, hspace)
}
