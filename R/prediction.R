#' @title Predictions from GP, DGP, or linked (D)GP emulators
#'
#' @description This function implements single-core or multi-core predictions (with or without multi-threading)
#'     from GP, DGP, or linked (D)GP emulators.
#'
#' @param object an instance of the `gp`, `dgp`, or `lgp` class.
#' @param x the testing input data:
#' * if `object` is an instance of the `gp` or `dgp` class, `x` is a matrix where each row is an input testing data point and each column is an input dimension.
#' * if `object` is an instance of the `lgp` class, `x` can be a matrix or a list:
#'    - if `x` is a matrix, it is the global testing input data that feed into the emulators in the first layer of a system.
#'      The rows of `x` represent different input data points and the columns represent input dimensions across all emulators in
#'      the first layer of the system. In this case, it is assumed that the only global input to the system is the input to the
#'      emulators in the first layer and there is no global input to emulators in other layers.
#'    - if `x` is a list, it should have *L* (the number of layers in an emulator system) elements. The first element
#'      is a matrix that represents the global testing input data that feed into the emulators in the first layer of the system. The
#'      remaining *L-1* elements are *L-1* sub-lists, each of which contains a number (the same number of emulators in
#'      the corresponding layer) of matrices (rows being testing input data points and columns being input dimensions) that represent the
#'      global testing input data to the emulators in the corresponding layer. The matrices must be placed in the sub-lists based on how
#'      their corresponding emulators are placed in `struc` argument of [lgp()]. If there is no global input data to a certain emulator,
#'      set `NULL` in the corresponding sub-list of `x`.
#' @param method the prediction approach: mean-variance (`"mean_var"`) or sampling (`"sampling"`) approach. Defaults to `"mean_var"`.
#' @param full_layer a bool indicating whether to output the predictions of all layers. Defaults to `FALSE`. Only used when `object` is a DGP and linked (D)GP emulator.
#' @param sample_size the number of samples to draw for each given imputation if `method = "sampling"`. Defaults to `50`.
#' @param cores the number of cores/workers to be used. If set to `NULL`,
#'     the number of cores is set to `(max physical cores available - 1)`. Defaults to `1`.
#' @param chunks the number of chunks that the testing input matrix `x` will be divided into for multi-cores to work on.
#'     Only used when `cores` is not `1`. If not specified (i.e., `chunks = NULL`), the number of chunks is set to the value of `cores`.
#'     Defaults to `NULL`.
#' @param threading a bool indicating whether to use the multi-threading to accelerate the predictions of DGP or linked (D)GP emulators. Turn this option on
#'     when you have a moderately large number of training data points as in such a case you could gain faster predictions. Defaults to `FALSE`.
#' @param ... N/A.
#' @return
#' * If `object` is an instance of the `gp` class:
#'   1. if `method = "mean_var"`: an updated `object` is returned with an additional slot called `results` that contains two matrices named `mean`
#'      for the predictive means and `var` for the predictive variances. Each matrix has only one column with its rows
#'      corresponding to testing positions (i.e., rows of `x`).
#'   2. if `method = "sampling"`: an updated `object` is returned with an additional slot called `results` that contains a matrix whose rows correspond
#'      to testing positions and columns correspond to `sample_size` number of samples drawn from the predictive distribution of GP.
#' * If `object` is an instance of the `dgp` class:
#'   1. if `method = "mean_var"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that contains two
#'      matrices named `mean` for the predictive means and `var` for the predictive variances respectively. Each matrix has its rows corresponding to testing
#'      positions and columns corresponding to DGP global output dimensions (i.e., the number of GP/likelihood nodes in the final layer).
#'   2. if `method = "mean_var"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains two
#'      sub-lists named `mean` for the predictive means and `var` for the predictive variances respectively. Each sub-list contains *L* (i.e., the number of layers)
#'      matrices named `layer1, layer2,..., layerL`. Each matrix has its rows corresponding to testing positions and columns corresponding to
#'      output dimensions (i.e., the number of GP/likelihood nodes from the associated layer).
#'   3. if `method = "sampling"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that contains *D* (i.e., the number
#'      of GP/likelihood nodes in the final layer) matrices named `output1, output2,..., outputD`. Each matrix has its rows corresponding to testing positions and
#'      columns corresponding to samples of size: `B * sample_size`, where `B` is the number of imputations specified in [dgp()].
#'   4. if `method = "sampling"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains *L* (i.e., the number
#'      of layers) sub-lists named `layer1, layer2,..., layerL`. Each sub-list represents samples drawn from the GP/likelihood nodes in the corresponding layer,
#'      and contains *D* (i.e., the number of GP/likelihood nodes in the corresponding layer) matrices named `output1, output2,..., outputD`. Each matrix gives samples
#'      of the output from one of *D* GP/likelihood nodes, and has its rows corresponding to testing positions and columns corresponding to samples
#'      of size: `B * sample_size`, where `B` is the number of imputations specified in [dgp()].
#' * If `object` is an instance of the `lgp` class:
#'   1. if `method = "mean_var"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that
#'      contains two sub-lists named `mean` for the predictive means and `var` for the predictive variances respectively. Each sub-list
#'      contains *M* number (same number of emulators in the final layer of the system) of matrices named `emulator1, emulator2,..., emulatorM`.
#'      Each matrix has its rows corresponding to global testing positions and columns corresponding to output dimensions of the associated emulator
#'      in the final layer.
#'   2. if `method = "mean_var"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains
#'      two sub-lists named `mean` for the predictive means and `var` for the predictive variances respectively. Each sub-list contains *L*
#'      (i.e., the number of layers in the emulated system) components named `layer1, layer2,..., layerL`. Each component represents a layer
#'      and contains *M* number (same number of emulators in the corresponding layer of the system) of matrices named `emulator1, emulator2,..., emulatorM`.
#'      Each matrix has its rows corresponding to global testing positions and columns corresponding to output dimensions of the associated
#'      GP/DGP emulator in the corresponding layer.
#'   3. if `method = "sampling"` and  `full_layer = FALSE`: an updated `object` is returned with an additional slot called `results` that contains
#'      *M* number (same number of emulators in the final layer of the system) of sub-lists named `emulator1, emulator2,..., emulatorM`. Each
#'      sub-list corresponds to an emulator in the final layer, and contains *D* matrices, named `output1, output2,..., outputD`, that correspond to the output
#'      dimensions of the GP/DGP emulator. Each matrix has its rows corresponding to testing positions and columns corresponding to samples
#'      of size: `B * sample_size`, where `B` is the number of imputations specified in [lgp()].
#'   4. if `method = "sampling"` and  `full_layer = TRUE`: an updated `object` is returned with an additional slot called `results` that contains
#'      *L* (i.e., the number of layers of the emulated system) sub-lists named `layer1, layer2,..., layerL`. Each sub-list represents a layer
#'      and contains *M* number (same number of emulators in the corresponding layer of the system) of components named `emulator1, emulator2,..., emulatorM`.
#'      Each component corresponds to an emulator in the associated layer, and contains *D* matrices, named `output1, output2,..., outputD`, that correspond to
#'      the output dimensions of the GP/DGP emulator. Each matrix has its rows corresponding to testing positions and columns corresponding to
#'      samples of size: `B * sample_size`, where `B` is the number of imputations specified in [lgp()].
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @note Any R vector detected in `x` will be treated as a column vector and automatically converted into a single-column R matrix.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @name predict
NULL


#' @rdname predict
#' @method predict dgp
#' @export
predict.dgp <- function(object, x, method = 'mean_var', full_layer = FALSE, sample_size = 50, cores = 1, chunks = NULL, threading = FALSE, ...) {
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x) ) x <- as.matrix(x)
  sample_size <- as.integer(sample_size)
  if( !is.null(chunks) ) {
    chunks <- as.integer(chunks)
    if ( chunks < 1 ) stop("The chunk number must be >= 1.", call. = FALSE)
  }
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  object$emulator_obj$set_nb_parallel(threading)
  if ( identical(cores,as.integer(1)) ){
    res <- object$emulator_obj$predict(x, method, full_layer, sample_size)
  } else {
    res <- object$emulator_obj$ppredict(x, method, full_layer, sample_size, chunks, cores)
  }

  if (method == 'mean_var'){
    if (full_layer) {
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
      for (l in 1:length(named_res$mean)) {
        names(named_res$mean)[l] <- paste('layer', l, sep="")
        names(named_res$var)[l] <- paste('layer', l, sep="")
      }
    } else {
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
    }
  } else if (method == 'sampling') {
    if (full_layer) {
      named_res <- res
      for (l in 1:length(named_res)) {
        names(named_res)[l] <- paste('layer', l, sep="")
        for (k in 1:length(named_res[[l]])) {
          names(named_res[[l]])[k] <- paste('output', k, sep="")
        }
      }
    } else {
      named_res <- res
      for (l in 1:length(named_res)) {
        names(named_res)[l] <- paste('output', l, sep="")
      }
    }
  }

  object$results <- named_res
  return(object)
}


#' @rdname predict
#' @method predict lgp
#' @export
predict.lgp <- function(object, x, method = 'mean_var', full_layer = FALSE, sample_size = 50, cores = 1, chunks = NULL, threading = FALSE, ...) {
  if ( !inherits(object,"lgp") ) stop("'object' must be an instance of the 'lgp' class.", call. = FALSE)
  if ( !is.list(x) ) {
    if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
    x <- unname(x)
    if ( is.vector(x) ) x <- as.matrix(x)
  } else {
    for ( l in 1:length(x) ){
      if ( l==1 ){
        if ( !is.matrix(x[[l]])&!is.vector(x[[l]]) ) {
          stop("The first element of 'x' must be a vector or a matrix.", call. = FALSE)
        } else {
          x[[l]] <- unname(x[[l]])
          if ( is.vector(x[[l]]) ) x[[l]] <- as.matrix(x[[l]])
          nrow_x <- nrow(x[[l]])
        }
      } else {
        for ( k in 1:length(x[[l]]) ){
          if ( !is.matrix(x[[l]][[k]])&!is.null(x[[l]][[k]])&!is.vector(x[[l]][[k]]) ) stop(sprintf("The element %i in the sublist %i of 'x' must be a vector, a matrix, or 'NULL'.", k, l), call. = FALSE)
          if ( is.matrix(x[[l]][[k]])|is.vector(x[[l]][[k]]) ){
            x[[l]][[k]] <- unname(x[[l]][[k]])
            if (is.vector(x[[l]][[k]])) x[[l]][[k]] <- as.matrix(x[[l]][[k]])
            if ( nrow(x[[l]][[k]])!=nrow_x ) {
              stop(sprintf("The element %i in the sublist %i of 'x' has inconsistent number of data points with the first element of 'x'.", k, l), call. = FALSE)
            }
          }
        }
      }
    }
  }

  sample_size <- as.integer(sample_size)
  if( !is.null(chunks) ) {
    chunks <- as.integer(chunks)
    if ( chunks < 1 ) stop("The chunk number must be >= 1.", call. = FALSE)
  }
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  object$emulator_obj$set_nb_parallel(threading)
  if ( identical(cores,as.integer(1)) ){
    res <- object$emulator_obj$predict(x, method, full_layer, sample_size)
  } else {
    res <- object$emulator_obj$ppredict(x, method, full_layer, sample_size, chunks, cores)
  }

  if (method == 'mean_var'){
    if (full_layer) {
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
      for (l in 1:length(named_res$mean)) {
        names(named_res$mean)[l] <- paste('layer', l, sep="")
        names(named_res$var)[l] <- paste('layer', l, sep="")
        for (k in 1:length(named_res$mean[[l]])) {
          names(named_res$mean[[l]])[k] <- paste('emulator', k, sep="")
          names(named_res$var[[l]])[k] <- paste('emulator', k, sep="")
        }
      }
    } else {
      named_res <- list("mean" = res[[1]], "var" = res[[2]])
      for (k in 1:length(named_res$mean)) {
        names(named_res$mean)[k] <- paste('emulator', k, sep="")
        names(named_res$var)[k] <- paste('emulator', k, sep="")
      }
    }
  } else if (method == 'sampling') {
    if (full_layer) {
      named_res <- list()
      for (l in 1:length(res)) {
        name1 <- paste('layer', l, sep="")
        for (k in 1:length(res[[l]])){
          name2 <- paste('emulator', k, sep="")
          for (n in 1:dim(res[[l]][[k]])[1]) {
            name3 <- paste('output', n, sep="")
            named_res[[name1]][[name2]][[name3]] <- res[[l]][[k]][n,,]
          }
        }
      }
    } else {
      named_res <- list()
      for (k in 1:length(res)) {
        name1 <- paste('emulator', k, sep="")
        for (n in 1:dim(res[[k]])[1]) {
          name2 <- paste('output', n, sep="")
          named_res[[name1]][[name2]] <- res[[k]][n,,]
        }
      }
    }
  }

  object$results <- named_res
  return(object)
}


#' @rdname predict
#' @method predict gp
#' @export
predict.gp <- function(object, x, method = 'mean_var', sample_size = 50, cores = 1, chunks = NULL, ...) {
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  if ( !is.matrix(x)&!is.vector(x) ) stop("'x' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x) ) x <- as.matrix(x)
  sample_size <- as.integer(sample_size)
  if( !is.null(chunks) ) {
    chunks <- as.integer(chunks)
    if ( chunks < 1 ) stop("The chunk number must be >= 1.", call. = FALSE)
  }
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  if ( identical(cores,as.integer(1)) ){
    res <- object$emulator_obj$predict(x, method, sample_size)
  } else {
    res <- object$emulator_obj$ppredict(x, method, sample_size, chunks, cores)
  }

  if (method == 'mean_var'){
    object$results <- list("mean" = res[[1]], "var" = res[[2]])
  } else if (method == 'sampling'){
    object$results <- res
  }

  return(object)
}

