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

#' @title Pack GP and DGP emulators into a bundle
#'
#' @description This function packs GP emulators and DGP emulators (without likelihood layers) into a `bundle` class for
#'     sequential designs if each emulator emulates one output dimension of the underlying simulator.
#'
#' @param ... a sequence of emulators produced by [gp()] or [dgp()].
#'
#' @return An S3 class named `bundle` to be used by [design()] for sequential designs. It has:
#' - *N* slots named `emulator1,...,emulatorN`, each of which contains a GP or DGP emulator, where *N* is the number of emulators
#'   that are provided to the function.
#' - a slot called `data` which contains two elements `X` and `Y`. `X` contains *N* matrices named `emulator1,...,emulatorN` that are
#'   training input data for different emulators. `Y` contains *N* single-column matrices named `emulator1,...,emulatorN` that are
#'   training output data for different emulators.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # load packages and the Python env
#' library(lhs)
#' library(dgpsi)
#' init_py()
#'
#' # construct a function with a two-dimensional output
#' f <- function(x) {
#'  y1 = sin(30*((2*x-1)/2-0.4)^5)*cos(20*((2*x-1)/2-0.4))
#'  y2 = 1/3*sin(2*(2*x - 1))+2/3*exp(-30*(2*(2*x-1))^2)+1/3
#'  return(cbind(y1,y2))
#' }
#'
#' # generate the initial design
#' X <- maximinLHS(10,1)
#' Y <- f(X)
#'
#' # generate the validation data
#' validate_x <- maximinLHS(30,1)
#' validate_y <- f(validate_x)
#'
#' # training a 2-layered DGP emulator with respect to each output with the global connection off
#' m1 <- dgp(X, Y[,1], connect=F)
#' m2 <- dgp(X, Y[,2], connect=F)
#'
#' # specify the range of the input dimension
#' lim <- c(0, 1)
#'
#' # pack emulators to form an emulator bundle
#' m <- pack(m1, m2)
#'
#' # 1st wave of the sequential design with 10 steps with target RMSE 0.01
#' m <- design(m, N=10, limits = lim, f = f, x_test = validate_x, y_test = validate_y, target = 0.01)
#'
#' # 2rd wave of the sequential design with 10 steps, the same target, and the aggregation
#' # function that takes the average of the criterion scores across the two outputs
#' g <- function(x){
#'   return(rowMeans(x))
#' }
#' m <- design(m, N=10, limits = lim, f = f, x_test = validate_x,
#'                     y_test = validate_y, aggregate = g, target = 0.01)
#'
#' # draw sequential designs of the two packed emulators
#' draw(m, emulator = 1, type = 'design')
#' draw(m, emulator = 2, type = 'design')
#'
#' # inspect the traces of RMSEs of the two packed emulators
#' draw(m, emulator = 1, type = 'rmse')
#' draw(m, emulator = 2, type = 'rmse')
#'
#' # write and read the constructed emulator bundle
#' write(m, 'bundle_dgp')
#' m <- read('bundle_dgp')
#'
#' # unpack the bundle into individual emulators
#' m_unpacked <- unpack(m)
#'
#' # plot OOS validations of individual emulators
#' plot(m_unpacked[[1]], x_test = validate_x, y_test = validate_y[,1])
#' plot(m_unpacked[[2]], x_test = validate_x, y_test = validate_y[,2])
#' }
#' @md
#' @export
pack <- function(...) {
  res <- list(...)
  X_all <- list()
  Y_all <- list()
  if ( length(res)==1 ) stop("The function needs at least two emulators to pack.", call. = FALSE)
  training_input <- res[[1]]$data$X
  training_output <- c()
  for ( i in 1:length(res) ){
    if ( !inherits(res[[i]],"gp") & !inherits(res[[i]],"dgp") ) stop("The function only accepts GP or DGP emulators as inputs.", call. = FALSE)
    if ( inherits(res[[i]],"dgp") ){
      if ( res[[i]]$constructor_obj$all_layer[[res[[i]]$constructor_obj$n_layer]][[1]]$type == 'likelihood' ){
        stop("The function can only pack DGP emulators without likelihood layers.", call. = FALSE)
      }
    }
    if ( !identical(res[[i]]$data$X, training_input) ) stop("The function can only pack emulators with common training input data.", call. = FALSE)
    Y_dim <- ncol(res[[i]]$data$Y)
    if ( Y_dim!=1 ) stop(sprintf("The function is only applicable to emulators with 1D output. Your emulator %i has %i output dimensions.", i, Y_dim), call. = FALSE)
    X_all[[paste('emulator', i ,sep="")]] <- unname(training_input)
    Y_all[[paste('emulator', i ,sep="")]] <- unname(res[[i]]$data$Y)
    names(res)[i] <- paste('emulator', i, sep="")
  }
  res[['data']][['X']] <- X_all
  res[['data']][['Y']] <- Y_all
  class(res) <- "bundle"
  return(res)
}


#' @title Unpack a bundle of (D)GP emulators
#'
#' @description This function unpacks a bundle of (D)GP emulators safely so any further manipulations of unpacked individual emulators
#'     will not impact the ones in the bundle.
#'
#' @param object an instance of the class `bundle`.
#'
#' @return A named list that contains individual emulators (named `emulator1,...,emulatorS`) packed in `object`,
#'    where `S` is the number of emulators in `object`.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See pack() for an example.
#' }
#' @md
#' @export
unpack <- function(object) {
  if ( !inherits(object,"bundle") ){
    stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)
  }
  n_emulators <- length(object) - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1
  res <- list()
  for ( i in 1:n_emulators ){
    res[[paste('emulator', i, sep="")]] <- object[[paste('emulator', i, sep="")]]
    res[[paste('emulator', i, sep="")]]$constructor_obj <- pkg.env$copy$deepcopy(res[[paste('emulator', i, sep="")]]$constructor_obj)
    res[[paste('emulator', i, sep="")]]$container_obj <- pkg.env$copy$deepcopy(res[[paste('emulator', i, sep="")]]$container_obj)
    res[[paste('emulator', i, sep="")]]$emulator_obj <- pkg.env$copy$deepcopy(res[[paste('emulator', i, sep="")]]$emulator_obj)
  }
  return(res)
}


#' @title Save the constructed emulator
#'
#' @description This function saves the constructed emulator to a `.pkl` file.
#'
#' @param object an instance of the S3 class `gp`, `dgp`, `lgp`, or `bundle`.
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
#' # See gp(), dgp(), lgp(), or pack() for an example.
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
#' @return The S3 class of a GP emulator, a DGP emulator, a linked (D)GP emulator, or a bundle of (D)GP emulators.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), lgp(), or pack() for an example.
#' }
#' @md
#' @export
read <- function(pkl_file) {
  pkl_file <- tools::file_path_sans_ext(pkl_file)
  res <- pkg.env$dgpsi$read(pkl_file)
  if ('emulator_obj' %in% names(res)){
    type <- pkg.env$py_buildin$type(res$emulator_obj)$'__name__'
    if ( type=='emulator' ) {
      class(res) <- "dgp"
    } else if ( type=='gp' ) {
      class(res) <- "gp"
    } else if ( type=='lgp' ) {
      class(res) <- "lgp"
    }
  } else {
    N <- length(res) - 1
    if ( "design" %in% names(res) ) N <- N - 1
    class(res) <- "bundle"
    for ( i in 1:N ){
      type <- pkg.env$py_buildin$type(res[[paste('emulator',i, sep='')]]$emulator_obj)$'__name__'
      if ( type=='emulator' ) {
        class(res[[paste('emulator',i, sep='')]]) <- "dgp"
      } else if ( type=='gp' ) {
        class(res[[paste('emulator',i, sep='')]]) <- "gp"
      }
    }
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

#' @title Reset number of imputations for a DGP emulator
#'
#' @description This function resets the number of imputations for predictions from a DGP emulator.
#'
#' @param object an instance of the S3 class `dgp`.
#' @param B the number of imputations to produce predictions from `object`. Increase the value to account for
#'     more imputation uncertainties with slower predictions. Decrease the value for lower imputation uncertainties
#'     but faster predictions. Defaults to `10`.
#'
#' @return An updated `object` with the information of `B` incorporated.
#'
#' @note
#' * This function is useful when a DGP emulator has been trained and one wants to make faster predictions by decreasing
#'    the number of imputations without rebuilding the emulator.
#' * The following slots:
#'   - `loo` and `oos` created by [validate()]; and
#'   - `results` created by [predict()]
#'   in `object` will be removed and not contained in the returned object.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See design() for an example.
#' }
#' @md
#' @export
set_imp <- function(object, B = 10) {
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }
  B <- as.integer(B)

  linked_idx <- object$container_obj$local_input_idx
  constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
  burnin <- constructor_obj_cp$burnin
  est_obj <- constructor_obj_cp$estimate(burnin)

  new_object <- list()
  new_object[['data']][['X']] <- object$data$X
  new_object[['data']][['Y']] <- object$data$Y
  new_object[['constructor_obj']] <- constructor_obj_cp
  new_object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B)
  new_object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx)
  if ( "design" %in% names(object) ) new_object[['design']] <- object$design
  class(new_object) <- "dgp"
  return(new_object)
}


#' @title Trim the sequences of model parameters of a DGP emulator
#'
#' @description This function trim the sequences of model parameters of a DGP emulator
#'     that are generated during the training.
#'
#' @param object an instance of the S3 class `dgp`.
#' @param start the first iteration before which all iterations are trimmed from the sequences.
#' @param end the last iteration after which all iterations are trimmed from the sequences.
#'     Set to `NULL` to keep all iterations after (including) `start`. Defaults to `NULL`.
#' @param thin the interval between the `start` and `end` iterations to thin out the sequences.
#'     Defaults to 1.
#'
#' @return An updated `object` with trimmed sequences of model parameters.
#'
#' @note
#' * This function is useful when a DGP emulator has been trained and one wants to trim
#'   the sequences of model parameters and use the trimmed sequences to generate the point estimates
#'   of DGP model parameters for predictions.
#' * The following slots:
#'   - `loo` and `oos` created by [validate()]; and
#'   - `results` created by [predict()]
#'   in `object` will be removed and not contained in the returned object.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See dgp() for an example.
#' }
#' @md
#' @export
window <- function(object, start, end = NULL, thin = 1) {
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }
  if (length(start) != 1) stop("'start' must be an integer.", call. = FALSE)
  if ( !is.null(end) ){
    if (length(end) != 1) stop("'end' must be an integer.", call. = FALSE)
  }
  if (length(thin) != 1) stop("'thin' must be an integer.", call. = FALSE)
  if ( !is.null(end) ){
    if (start > end) stop("'start' must be before 'end'", call. = FALSE)
  }

  linked_idx <- object$container_obj$local_input_idx
  constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
  niter <- constructor_obj_cp$N + 1
  if ( is.null(end) ) end <- niter
  if (end > niter) end <- niter
  idx <- start:end
  idx <- idx[idx %% thin == 0]
  constructor_obj_cp$N <- length(idx) - 1
  for ( l in 1:constructor_obj_cp$n_layer ){
    n_kernel <- length(constructor_obj_cp$all_layer[[l]])
    for (k in 1:n_kernel){
      if (constructor_obj_cp$all_layer[[l]][[k]]$type == 'gp'){
        final_est <- constructor_obj_cp$all_layer[[l]][[k]]$para_path[idx[length(idx)],]
        constructor_obj_cp$all_layer[[l]][[k]]$scale <- reticulate::np_array(final_est[1])
        constructor_obj_cp$all_layer[[l]][[k]]$length <- reticulate::np_array(final_est[2:(length(final_est)-1)])
        constructor_obj_cp$all_layer[[l]][[k]]$nugget <- reticulate::np_array(final_est[length(final_est)])
        constructor_obj_cp$all_layer[[l]][[k]]$para_path <- constructor_obj_cp$all_layer[[l]][[k]]$para_path[idx,,drop=FALSE]
      }
    }
  }

  burnin <- 0L
  est_obj <- constructor_obj_cp$estimate(burnin)
  B <- length(object$emulator_obj$all_layer_set)

  new_object <- list()
  new_object[['data']][['X']] <- object$data$X
  new_object[['data']][['Y']] <- object$data$Y
  new_object[['constructor_obj']] <- constructor_obj_cp
  new_object[['emulator_obj']] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B)
  new_object[['container_obj']] <- pkg.env$dgpsi$container(est_obj, linked_idx)
  if ( "design" %in% names(object) ) new_object[['design']] <- object$design
  class(new_object) <- "dgp"
  return(new_object)
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
#'     R matrix. Thus, if `x` is a single testing data point with multiple dimensions, it must be given as a matrix.
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
