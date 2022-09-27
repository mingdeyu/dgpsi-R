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
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See lgp() for an example.
#' }
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
#' @return No return value. `object` will be save to a local `.pkl` file specified by `pkl_file`.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @note Since the constructed emulators are 'python' objects, [save()] from R will not work as it is only for R objects.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @export
write <- function(object, pkl_file) {
  pkl_file <- tools::file_path_sans_ext(pkl_file)
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
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @export
read <- function(pkl_file) {
  pkl_file <- tools::file_path_sans_ext(pkl_file)
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
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
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
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See lgp() for an example.
#' }
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
#' @param object an instance of the `dgp` class and it should be produced by [dgp()] with one of the following two settings:
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
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @note Any R vector detected in `x` and `y` will be treated as a column vector and automatically converted into a single-column
#'     R matrix.
#' @examples
#' \dontrun{
#'
#' # Check https://mingdeyu.github.io/dgpsi-R/ for examples
#' # on how to compute the negative predicted log-likelihood
#' # using nllik().
#' }
#' @md
#' @export
nllik <- function(object, x, y) {
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(y)&!is.vector(y) ) stop("'y' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x) ) x <- as.matrix(x)
  if ( is.vector(y) ) y <- as.matrix(y)

  if ( nrow(x)!=nrow(y) ) stop("'x' and 'y' have different number of data points.", call. = FALSE)
  n_dim_y <- ncol(y)
  if ( n_dim_y != 1 ) {
    stop("'y' must be a vector or a matrix with only one column.", call. = FALSE)
  }
  res <- object$emulator_obj$nllik(x, y)
  named_res <- list("meanNLL" = res[[1]], "allNLL" = res[[2]])
  object$NLL <- named_res
  return(object)
}


#' @title Plot of DGP model parameter traces
#'
#' @description This function plots the traces of model parameters of a chosen GP node
#'     in a DGP emulator.
#'
#' @param object an instance of the `dgp` class.
#' @param layer the index of a layer. Defaults to `NULL` for the final layer.
#' @param node the index of a GP node in the layer specified by `layer`. Defaults to `1` for the first GP node in the
#'     corresponding layer.
#'
#' @return A `ggplot` object.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See dgp() for an example.
#' }
#' @md
#' @export

trace_plot <- function(object, layer = NULL, node = 1) {
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)

  all_layer <- object$constructor_obj$all_layer

  layer_no <- length(all_layer)

  if ( is.null(layer) ){
    layer = layer_no
  }

  layer_l <- all_layer[[layer]]
  kernel_no <- length(layer_l)
  kernel <- layer_l[[node]]
  if (kernel$type == 'gp'){
    path <- kernel$para_path
    n_para <- ncol(path)
    dat <- as.data.frame(path)
    for ( i in 1:n_para){
      if ( i==1 ){
        colnames(dat)[i] <- 'Variance'
      } else if ( i==n_para ){
        colnames(dat)[i] <- 'Nugget'
      } else {
        colnames(dat)[i] <- paste('Lengthscale', i-1, sep = ' ')
      }
    }
    dat$Iteration <-  seq(nrow(path))
    mm <- reshape2::melt(dat, id.var="Iteration")

    color <- c("#E69F00", rep("#56B4E9",n_para-2), "#009E73")
    ggplot2::ggplot(mm, ggplot2::aes_(x = ~Iteration, y = ~value)) + ggplot2::geom_line(ggplot2::aes_(color = ~variable)) +
      ggplot2::facet_grid(variable ~ ., scales = "free_y") + ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title=paste("Node ", node, " (of ", kernel_no, ")", " in Layer ", layer, " (of ", layer_no, ")", sep = ""), x ="Iteration", y = "Parameter value") +
      ggplot2::scale_color_manual(values = color)
  } else {
    message('There is nothing to plot for a likelihood node, please choose a GP node instead.')
  }

}
