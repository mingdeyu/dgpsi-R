#' @title Locate the next design point for a (D)GP emulator or a bundle of (D)GP emulators
#'
#' @description This function searches from a candidate set to locate the next design point to be added to a (D)GP or a bundle of
#'     (D)GP emulators.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param x_cand a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
#'     in which the next design point is determined. If `x_cand = NULL`, the candidate set will be generated using `n_cand` and
#'     `limits`. Defaults to `NULL`.
#' @param n_cand an integer that gives the size of the candidate set in which the next design point is determined. This argument
#'     is used when `x_cand = NULL`. Defaults to `200`
#' @param limits a two-column matrix that gives the ranges of each input dimension, or a vector of length two if there is only one
#'     input dimension. If a vector is provided, it will be converted to a two-column row matrix. The rows of the matrix correspond
#'     to input dimensions, and its first and second columns correspond to the minimum and maximum values of the input dimensions.
#'     If `limits = NULL`, the ranges of input dimensions will be determined from the training data contained in `object`. This
#'     argument is used when `x_cand = NULL`. Defaults to `NULL`.
#' @param batch_size an integer that gives the number of design points (for each simulator output dimension if `aggregate = NULL`) to be chosen.
#'     Defaults to `1`.
#' @param method the criterion used to locate the next design point: ALM (`"ALM"`) or MICE (`"MICE"`). See references below.
#'     Defaults to `"ALM"`.
#' @param nugget_s the value of the smoothing nugget term used when `method = "MICE"`. Defaults to `1.0`.
#' @param verb a bool indicating if the trace information will be printed during the function execution.
#'     Defaults to `TRUE`.
#' @param cores the number of cores/workers to be used for the criterion calculation. If set to `NULL`,
#'     the number of cores is set to `(max physical cores available - 1)`. Defaults to `1`.
#' @param threading a bool indicating whether to use the multi-threading to accelerate the criterion calculation.
#'     Turning this option on could improve the speed of criterion calculations when the emulator is built with a moderately large number of
#'     training data points and the Matérn-2.5 kernel.
#' @param aggregate an R function with only one argument which is a vector with the length equal to:
#' * the emulator output dimension if `object` is an instance of the `dgp` class; or
#' * the number of emulators contained in `object` if `object` is an instance of the `bundle` class.
#'
#' The vector represents values of the ALM or MICE criterion across different output dimensions at a design point in the candidate set. The
#'     output of `aggregate` should be a single value that gives an aggregation of the ALM or MICE criterion across different output dimensions.
#'     Set to `NULL` to disable the aggregation. Defaults to `NULL`.
#' @param ... N/A.
#'
#' @return
#' * If `object` is an instance of the `gp` class, a list is returned:
#'   - when `x_cand = NULL`, the list has one slot called `location` that contains
#'     - a single-row matrix that gives the position of determined design point, if `batch_size = 1`; or
#'     - a number of slots named `position1,...,positionS` where `S = batch_size`, each of which is a single-row matrix giving the position of
#'       a determined design point, if `batch_size` is greater than `1`.
#'   - when `x_cand` is not `NULL`, the list has an additional slot called `index` that contains
#'     - the row index of the next design point in `x_cand`, if `batch_size = 1`; or
#'     - a number of slots named `position1,...,positionS` giving indices of next design points in `x_cand`, if `batch_size` is greater than `1`.
#
#' * If `object` is an instance of the `dgp` class, a list is returned:
#'   - when `x_cand = NULL` and `aggregate = NULL`, the list has one slot called `location` that contains
#'     - a matrix whose rows give the positions of determined design points with respect to different output dimensions, if `batch_size = 1`; or
#'     - a number of slots named `position1,...,positionS` where `S = batch_size`, each of which is a matrix giving the position of determined
#'       design points with respect to different output dimensions, if `batch_size` is greater than `1`.
#'
#'     when `x_cand` is not `NULL`, the list has an additional slot called `index` that contains
#'     - a vector giving the row indices of determined design points in `x_cand` with respect to different output dimensions, if `batch_size = 1`; or
#'     - a number of vectors named `position1,...,positionS`, each of which gives indices of determined design points in `x_cand` with respect to
#'       different output dimensions, if `batch_size` is greater than `1`.
#'   - when `x_cand = NULL` and `aggregate` is provided, the list has one slot called `location` that contains
#'     - a single-row matrix that gives the position of the determined design point, if `batch_size = 1`; or
#'     - a number of slots named `position1,...,positionS` where `S = batch_size`, each of which is a single-row matrix giving the position of
#'       the determined design point, if `batch_size` is greater than `1`.
#'
#'     when `x_cand` is not `NULL`, the list has an additional slot called `index` that contains
#'     - the row index of the next design point in `x_cand`, if `batch_size = 1`; or
#'     - a number of slots named `position1,...,positionS`, each of which gives the index of a determined design point in `x_cand`, if `batch_size`
#'       is greater than `1`.
#'
#' * If `object` is an instance of the `bundle` class, a list is returned:
#'   - when `x_cand = NULL` and `aggregate = NULL`, the list has one slot called `location` that contains
#'     - a matrix whose rows give the positions of determined design points with respect to different emulators in the bundle, if `batch_size = 1`;
#'     - a number of slots named `position1,...,positionS` where `S = batch_size`, each of which is a matrix giving the position of determined
#'       design points with respect to different emulators in the bundle, if `batch_size` is greater than `1`.
#'
#'     when `x_cand` is not `NULL`, the list has an additional slot called `index` that contains
#'     - a vector giving the row indices of determined design points in `x_cand` with respect to different emulators in the bundle, if `batch_size = 1`; or
#'     - a number of vectors named `position1,...,positionS`, each of which gives indices of determined design points in `x_cand` with respect to
#'       different emulators in the bundle, if `batch_size` is greater than `1`.
#'   - when `x_cand = NULL` and `aggregate` is provided, the list has one slot called `location` that contains
#'     - a single-row matrix that gives the position of the determined design point, if `batch_size = 1`; or
#'     - a number of slots named `position1,...,positionS` where `S = batch_size`, each of which is a single-row matrix giving the position of
#'       the determined design point, if `batch_size` is greater than `1`.
#'
#'     when `x_cand` is not `NULL`, the list has an additional slot called `index` that contains
#'     - the row index of the next design point in `x_cand`, if `batch_size = 1`; or
#'     - a number of slots named `position1,...,positionS`, each of which gives the index of a determined design point in `x_cand`, if `batch_size`
#'       is greater than `1`.
#' @note
#' * The function is only applicable to GP emulators, DGP emulators without likelihood layers, or bundles of (D)GP emulators created by [pack()].
#' * Any R vector detected in `x_cand` will be treated as a column vector and automatically converted into a single-column
#'   R matrix.
#' @references
#' MacKay, D. J. (1992). Information-based objective functions for active data selection. *Neural Computation*, **4(4)**, 590-604.
#'
#' Beck, J., & Guillas, S. (2016). Sequential design with mutual information for computer experiments (MICE): emulation of a tsunami model.
#' *SIAM/ASA Journal on Uncertainty Quantification*, **4(1)**, 739-766.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # load packages and the Python env
#' library(lhs)
#' library(dgpsi)
#' init_py()
#'
#' # construct a 1D non-stationary function
#' f <- function(x) {
#'  sin(30*((2*x-1)/2-0.4)^5)*cos(20*((2*x-1)/2-0.4))
#' }
#'
#' # generate the initial design
#' X <- maximinLHS(10,1)
#' Y <- f(X)
#'
#' # training a 2-layered DGP emulator with the global connection off
#' m <- dgp(X, Y, connect = F)
#'
#' # specify the range of the input dimension
#' lim <- c(0,1)
#'
#' # locate the next design point
#' next_point <- locate(m, limits = lim)
#' X_new <- next_point$location
#'
#' # obtain the corresponding output at the located design point
#' Y_new <- f(X_new)
#'
#' # combine the new input-output pair to the existing data
#' X <- rbind(X, X_new)
#' Y <- rbind(Y, Y_new)
#'
#' # update the DGP emulator with the new input and output data and refit with 500 training iterations
#' m <- update(m, X, Y, refit = TRUE, N = 500)
#'
#' # plot the LOO validation
#' plot(m)
#' }
#' @md
#' @name locate
#' @export
locate <- function(object, x_cand, n_cand, limits, batch_size, method, nugget_s, verb, cores, ...){
  UseMethod("locate")
}

#' @rdname locate
#' @method locate gp
#' @export
locate.gp <- function(object, x_cand = NULL, n_cand = 200, limits = NULL, batch_size = 1, method = 'ALM', nugget_s = 1., verb = TRUE, cores = 1, ...) {
  #check class
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)

  #extract input dimensions
  training_input <- object$data$X
  n_dim_X <- ncol(training_input)

  if ( isTRUE(verb) ) message("Locating the next design point(s) ...", appendLF = FALSE)
  #locate next design point
  if ( is.null(x_cand) ){
    #check limits
    if ( !is.null(limits) ){
      if ( !is.matrix(limits)&!is.vector(limits) ) stop("'limits' must be a vector or a matrix.", call. = FALSE)
      if ( is.matrix(limits) ){
        if ( ncol(limits)!=2 ) stop("'limits' must be a matrix with two columns.", call. = FALSE)
      } else if ( is.vector(limits) ){
        if ( length(limits)!=2 ) stop("'limits' must be a vector of length two or a matrix with two columns.", call. = FALSE)
        limits <- matrix(limits, ncol=2, byrow=TRUE)
      }
      if ( nrow(limits)!=n_dim_X ) stop("You must provide ranges for all input dimensions through 'limits'.", call. = FALSE)
    }
  #  scale_training_input <- matrix(0,n_points,n_dim_X)
  #  if ( is.null(limits) ){
  #    for (i in 1:n_dim_X){
  #      scale_training_input[,i] <- (training_input[,i]-min(training_input[,i]))/(max(training_input[,i])-min(training_input[,i]))
  #    }
  #  } else {
  #    for (i in 1:n_dim_X){
  #      scale_training_input[,i] <- (training_input[,i]-limits[i,1])/(limits[i,2]-limits[i,1])
  #    }
  #  }
  #  lhd_full <- lhs::augmentLHS(scale_training_input, n_cand)
  #  lhd <- lhd_full[-seq(1,n_points),,drop=FALSE]
    x_cand <- lhs::maximinLHS(n_cand,n_dim_X)
    if ( is.null(limits) ){
      for (i in 1:n_dim_X){
        x_cand[,i] <- x_cand[,i]*(max(training_input[,i])-min(training_input[,i])) + min(training_input[,i])
      }
    } else {
      for (i in 1:n_dim_X){
        x_cand[,i] <- x_cand[,i]*(limits[i,2]-limits[i,1]) + limits[i,1]
      }
    }
    if ( identical(cores,as.integer(1)) ){
      res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size, core_num = cores)
    }
    dat <- list()
    if ( batch_size==1 ){
      dat[['location']] <- x_cand[res[[1]]+1,,drop=FALSE]
    } else {
      for (j in 1:batch_size){
        dat[['location']][[paste('position', j, sep="")]] <- x_cand[res[[1]][j]+1,,drop=FALSE]
      }
    }
  } else {
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
    if ( identical(cores,as.integer(1)) ){
      res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size, core_num = cores)
    }
    dat <- list()
    if ( batch_size==1 ){
      dat[['index']] <- res[[1]]+1
      dat[['location']] <- x_cand[res[[1]]+1,,drop=FALSE]
    } else {
      for (j in 1:batch_size){
        dat[['index']][[paste('position', j, sep="")]] <- res[[1]][j]+1
        dat[['location']][[paste('position', j, sep="")]] <- x_cand[res[[1]][j]+1,,drop=FALSE]
      }
    }
  }
  if ( isTRUE(verb) ) message(" done")
  if ( isTRUE(verb) ) {
    if ( batch_size==1 ){
      message(paste(c("The next design point(s):\n", sprintf("%.06f", dat[['location']])), collapse=" "))
    } else {
      message(paste(c("The next design points(s):")))
      for (i in 1:length(dat[['location']]) ){
        message(paste(" Position", i, ": ", paste(c(sprintf("%.06f", dat[['location']][[paste('position',i,sep="")]])), collapse=" "), sep=""))
      }
    }
  }
  return(dat)
}


#' @rdname locate
#' @method locate dgp
#' @export
locate.dgp <- function(object, x_cand = NULL, n_cand = 200, limits = NULL, batch_size = 1, method = 'ALM', nugget_s = 1., verb = TRUE, cores = 1, threading = FALSE, aggregate = NULL, ...) {
  #check class
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$type == 'likelihood' ){
    stop("The function is only applicable to DGP emulators without likelihood layers.", call. = FALSE)
  }
  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)

  object$emulator_obj$set_nb_parallel(threading)
  #extract input dimensions
  training_input <- object$data$X
  n_dim_X <- ncol(training_input)
  #extract output dimensions
  n_dim_Y <- ncol(object$data$Y)
  if ( n_dim_Y==1 & !is.null(aggregate) ){
    if ( isTRUE(verb) ) message("The emulator output is one-dimensional, the provided 'aggregate' function will be ignored.")
    aggregate <- NULL
  }
  if ( isTRUE(verb) ) message("Locating the next design point(s) ...", appendLF = FALSE)
  #locate next design point
  if ( is.null(x_cand) ){
    #check limits
    if ( !is.null(limits) ){
      if ( !is.matrix(limits)&!is.vector(limits) ) stop("'limits' must be a vector or a matrix.", call. = FALSE)
      if ( is.matrix(limits) ){
        if ( ncol(limits)!=2 ) stop("'limits' must be a matrix with two columns.", call. = FALSE)
      } else if ( is.vector(limits) ){
        if ( length(limits)!=2 ) stop("'limits' must be a vector of length two or a matrix with two columns.", call. = FALSE)
        limits <- matrix(limits, ncol=2, byrow=TRUE)
      }
      if ( nrow(limits)!=n_dim_X ) stop("You must provide ranges for all input dimensions through 'limits'.", call. = FALSE)
    }
    #scale_training_input <- matrix(0,n_points,n_dim_X)
    #if ( is.null(limits) ){
    #  for (i in 1:n_dim_X){
    #    scale_training_input[,i] <- (training_input[,i]-min(training_input[,i]))/(max(training_input[,i])-min(training_input[,i]))
    #  }
    #} else {
    #  for (i in 1:n_dim_X){
    #    scale_training_input[,i] <- (training_input[,i]-limits[i,1])/(limits[i,2]-limits[i,1])
    #  }
    #}
    #lhd_full <- lhs::augmentLHS(scale_training_input, n_cand)
    #lhd <- lhd_full[-seq(1,n_points),,drop=FALSE]
    #x_cand <- matrix(0,n_cand,n_dim_X)
    #lhd <- lhs::maximinLHS(n_cand,n_dim_X)
    x_cand <- lhs::maximinLHS(n_cand,n_dim_X)
    if ( is.null(limits) ){
      for (i in 1:n_dim_X){
        x_cand[,i] <- x_cand[,i]*(max(training_input[,i])-min(training_input[,i])) + min(training_input[,i])
      }
    } else {
      for (i in 1:n_dim_X){
        x_cand[,i] <- x_cand[,i]*(limits[i,2]-limits[i,1]) + limits[i,1]
      }
    }

    if ( is.null(aggregate) ){
      if ( identical(cores,as.integer(1)) ){
        res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size, core_num = cores)
      }
      dat <- list()
      if ( batch_size==1 ){
        dat[['location']] <- x_cand[res[[1]]+1,,drop=FALSE]
      } else {
        for (j in 1:batch_size){
          dat[['location']][[paste('position', j, sep="")]] <- x_cand[res[[1]][j,]+1,,drop=FALSE]
        }
      }
    } else {
      if ( identical(cores,as.integer(1)) ){
        res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE, core_num = cores)
      }
      agg_scores <- apply(res, 1, aggregate)
      dat <- list()
      if ( batch_size==1 ){
        idx <- which.max(agg_scores)
        dat[['location']] <- x_cand[idx,,drop=FALSE]
      } else {
        idx <- order(agg_scores, decreasing=TRUE)[1:batch_size]
        for (j in 1:batch_size){
          dat[['location']][[paste('position', j, sep="")]] <- x_cand[idx[j],,drop=FALSE]
        }
      }
    }
  } else {
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
    if ( is.null(aggregate) ){
      if ( identical(cores,as.integer(1)) ){
        res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size, core_num = cores)
      }
      dat <- list()
      if ( batch_size==1 ){
        dat[['index']] <- res[[1]]+1
        dat[['location']] <- x_cand[res[[1]]+1,,drop=FALSE]
      } else {
        for (j in 1:batch_size){
          dat[['index']][[paste('position', j, sep="")]] <- res[[1]][j,]+1
          dat[['location']][[paste('position', j, sep="")]] <- x_cand[res[[1]][j,]+1,,drop=FALSE]
        }
      }
    } else {
      if ( identical(cores,as.integer(1)) ){
        res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE, core_num = cores)
      }
      agg_scores <- apply(res, 1, aggregate)
      dat <- list()
      if ( batch_size==1 ){
        idx <- which.max(agg_scores)
        dat[['index']] <- idx
        dat[['location']] <- x_cand[idx,,drop=FALSE]
      } else {
        idx <- order(agg_scores, decreasing=TRUE)[1:batch_size]
        for (j in 1:batch_size){
          dat[['index']][[paste('position', j, sep="")]] <- idx[j]
          dat[['location']][[paste('position', j, sep="")]] <- x_cand[idx[j],,drop=FALSE]
        }
      }
    }
  }
  if ( isTRUE(verb) ) message(" done")
  if ( isTRUE(verb) ) {
    if ( batch_size==1 ){
      if ( nrow(dat[['location']])==1 ){
        message(paste(c("The next design point(s):\n", sprintf("%.06f", dat[['location']][1,])), collapse=" "))
      } else {
        message(paste(c("The next design point(s):")))
        for (i in 1:nrow(dat[['location']]) ){
          message(paste(" Output", i, ": ", paste(c(sprintf("%.06f", dat[['location']][i,])), collapse=" "), sep=""))
        }
      }
    } else {
      message(paste(c("The next design point(s):")))
      for (i in 1:batch_size ){
        pos_i <- dat[['location']][[paste('position',i,sep="")]]
        n_pos <- nrow(pos_i)
        if ( n_pos==1 ){
          message(paste(" Position", i, ": ", paste(c(sprintf("%.06f", pos_i[1,])), collapse=" "), sep=""))
        } else {
          for (j in 1:n_pos){
            message(paste(" Position", i, " for Output", j, ": ", paste(c(sprintf("%.06f", pos_i[j,])), collapse=" "), sep=""))
          }
        }
      }
    }
  }
  return(dat)
}


#' @rdname locate
#' @method locate bundle
#' @export
locate.bundle <- function(object, x_cand = NULL, n_cand = 200, limits = NULL, batch_size = 1, method = 'ALM', nugget_s = 1., verb = TRUE, cores = 1, threading = FALSE, aggregate = NULL, ...) {
  #check class
  if ( !inherits(object,"bundle") ) stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)

  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  batch_size <- as.integer(batch_size)

  n_dim_X <- ncol(object$data$X[[1]])

  x_cand_null <- is.null(x_cand)

  #check/create x_cand
  if ( x_cand_null ){
    #check limits
    if ( !is.null(limits) ){
      if ( !is.matrix(limits)&!is.vector(limits) ) stop("'limits' must be a vector or a matrix.", call. = FALSE)
      if ( is.matrix(limits) ){
        if ( ncol(limits)!=2 ) stop("'limits' must be a matrix with two columns.", call. = FALSE)
      } else if ( is.vector(limits) ){
        if ( length(limits)!=2 ) stop("'limits' must be a vector of length two or a matrix with two columns.", call. = FALSE)
        limits <- matrix(limits, ncol=2, byrow=TRUE)
      }
      if ( nrow(limits)!=n_dim_X ) stop("You must provide ranges for all input dimensions through 'limits'.", call. = FALSE)
    } else {
      all_training_input <- c()
      for ( j in 1:length(object$data$X) ){
        all_training_input <- rbind(all_training_input, object$data$X[[j]])
      }
    }
    x_cand <- lhs::maximinLHS(n_cand,n_dim_X)
    if ( is.null(limits) ){
      for (i in 1:n_dim_X){
        x_cand[,i] <- x_cand[,i]*(max(all_training_input[,i])-min(all_training_input[,i])) + min(all_training_input[,i])
      }
    } else {
      for (i in 1:n_dim_X){
        x_cand[,i] <- x_cand[,i]*(limits[i,2]-limits[i,1]) + limits[i,1]
      }
    }
  } else {
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
  }

  if ( is.null(aggregate) ) {
    if ( batch_size==1 ){
      location <- c()
      if ( !x_cand_null ) indice <- c()
    } else {
      location <- list()
      if ( !x_cand_null ) indice <- list()
    }
  } else {
    score <- c()
  }

  if ( isTRUE(verb) ) message("Locating the next design point(s) ...", appendLF = FALSE)
  n_emulators <- length(object)
  if ( "data" %in% names(object) ) n_emulators <- n_emulators - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1
  for ( i in 1:n_emulators ){
    obj_i <- object[[paste('emulator',i,sep='')]]
    if ( inherits(obj_i,"gp") ){
      if ( is.null(aggregate) ){
        res = obj_i$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size)
        if ( batch_size==1 ){
          location <- rbind(location, x_cand[res[[1]]+1,])
          if ( !x_cand_null ) indice <- c(indice, res[[1]]+1)
        } else {
          for ( j in 1:batch_size ){
            location[[paste('position',j,sep="")]] <- rbind(location[[paste('position',j,sep="")]], x_cand[res[[1]][j]+1,])
            if ( !x_cand_null ) indice[[paste('position',j,sep="")]] <- c(indice[[paste('position',j,sep="")]], res[[1]][j]+1)
          }
        }
      } else {
        res = obj_i$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE)
        score <- cbind(score, res)
      }
    } else if ( inherits(obj_i,"dgp") ){
      obj_i$emulator_obj$set_nb_parallel(threading)
      #locate next design point
      if ( is.null(aggregate) ){
        if ( identical(cores,as.integer(1)) ){
          res = obj_i$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size)
        } else {
          res = obj_i$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, batch_size = batch_size, core_num = cores)
        }
        if ( batch_size==1 ){
          location <- rbind(location, x_cand[res[[1]]+1,,drop=FALSE])
          if ( !x_cand_null ) indice <- c(indice, res[[1]]+1)
        } else {
          for ( j in 1:batch_size ){
            location[[paste('position',j,sep="")]] <- rbind(location[[paste('position',j,sep="")]], x_cand[res[[1]][j]+1,,drop=FALSE])
            if ( !x_cand_null ) indice[[paste('position',j,sep="")]] <- c(indice[[paste('position',j,sep="")]], res[[1]][j]+1)
          }
        }
      } else {
        if ( identical(cores,as.integer(1)) ){
          res = obj_i$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE)
        } else {
          res = obj_i$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE, core_num = cores)
        }
        score <- cbind(score, res)
      }
    }
  }
  if ( isTRUE(verb) ) message(" done")

  if ( is.null(aggregate) ) {
    if ( batch_size==1 ){
      if ( isTRUE(verb) ) {
        message(paste(c("The next design point(s):")))
        for (i in 1:nrow(location) ){
          message(paste(" Emulator", i, ": ", paste(c(sprintf("%.06f", location[i,])), collapse=" "), sep=""))
        }
      }
    } else {
      if ( isTRUE(verb) ) {
        message(paste(c("The next design point(s):")))
        for (i in 1:batch_size ){
          pos_i <- location[[paste('position',i,sep="")]]
          n_pos <- nrow(pos_i)
          for (j in 1:n_pos){
            message(paste(" Position", i, " for Emulator", j, ": ", paste(c(sprintf("%.06f", pos_i[j,])), collapse=" "), sep=""))
          }
        }
      }
    }
    if ( x_cand_null ) {
      dat <- list('location' = location)
    } else {
      dat <- list('idx' = indice, 'location' = location)
    }
  } else {
    agg_scores <- apply(score, 1, aggregate)
    dat <- list()
    if ( batch_size==1 ){
      idx <- which.max(agg_scores)
      if ( !x_cand_null ) dat[['index']] <- idx
      dat[['location']] <- x_cand[idx,,drop=FALSE]
    } else {
      idx <- order(agg_scores, decreasing=TRUE)[1:batch_size]
      for (j in 1:batch_size){
        if ( !x_cand_null ) dat[['index']][[paste('position', j, sep="")]] <- idx[j]
        dat[['location']][[paste('position', j, sep="")]] <- x_cand[idx[j],,drop=FALSE]
      }
    }

    if ( isTRUE(verb) ) {
      if ( batch_size==1 ){
        message(paste(c("The next design point(s):\n", sprintf("%.06f", dat[['location']][1,])), collapse=" "))
      } else {
        message(paste(c("The next design point(s):")))
        for (i in 1:batch_size ){
          pos_i <- dat[['location']][[paste('position',i,sep="")]]
          message(paste(" Position", i, ": ", paste(c(sprintf("%.06f", pos_i[1,])), collapse=" "), sep=""))
        }
      }
    }
  }
  return(dat)
}






