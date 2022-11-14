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
#'   - when `x_cand = NULL`, the list has one slot called `location` that contains a vector that gives the position of determined
#'     design point.
#'   - when `x_cand` is not `NULL`, the list has two slots. The first slot is called `index` that gives the row index of the next design point in
#'     `x_cand`. The second slot is called `location` that gives the position of the determined design point.
#
#' * If `object` is an instance of the `dgp` class, a list is returned:
#'   - when `x_cand = NULL` and `aggregate = NULL`, the list has one slot called `location` that contains a matrix whose rows give the positions of
#'     determined design points with respect to different output dimensions;
#'   - when `x_cand = NULL` and `aggregate` is provided, the list has one slot called `location` that contains a single-row matrix that gives the
#'     position of the determined design point.
#'   - when `x_cand` is not `NULL` and `aggregate = NULL`, the list has two slots:
#'     - the first slot is called `index` that contains a vector giving the row indices of determined design points in `x_cand` with respect to
#'       different output dimensions.
#'     - the second slot is called `location` that contains a matrix whose rows give the positions of determined design points with respect to
#'       different output dimensions.
#'   - when `x_cand` is not `NULL` and `aggregate` is provided, the list has two slots:
#'     - the first slot is called `index` that contains the row index of the next design point in `x_cand`.
#'     - the second slot is called `location` that contains a single-row matrix that gives the position of the determined design point.
#'
#' * If `object` is an instance of the `bundle` class, a list is returned:
#'   - when `x_cand = NULL` and `aggregate = NULL`, the list has one slot called `location` that contains a matrix whose rows give the positions of
#'     determined design points with respect to different emulators in the bundle;
#'   - when `x_cand = NULL` and `aggregate` is provided, the list has one slot called `location` that contains a single-row matrix that gives the
#'     position of the determined design point.
#'   - when `x_cand` is not `NULL` and `aggregate = NULL`, the list has two slots:
#'     - the first slot is called `index` that contains a vector giving the row indices of determined design points in `x_cand` with respect to
#'       different emulators in the bundle.
#'     - the second slot is called `location` that contains a matrix whose rows give the positions of determined design points with respect to
#'       different emulators in the bundle.
#'   - when `x_cand` is not `NULL` and `aggregate` is provided, the list has two slots:
#'     - the first slot is called `index` that contains the row index of the next design point in `x_cand`.
#'     - the second slot is called `location` that contains a single-row matrix that gives the position of the determined design point.
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
#' Y <- as.matrix(sapply(X, f))
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
locate <- function(object, x_cand, n_cand, limits, method, nugget_s, verb, cores, ...){
  UseMethod("locate")
}

#' @rdname locate
#' @method locate gp
#' @export
locate.gp <- function(object, x_cand = NULL, n_cand = 200, limits = NULL, method = 'ALM', nugget_s = 1., verb = TRUE, cores = 1, ...) {
  #check class
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  #extract input dimensions
  training_input <- object$data$X
  n_dim_X <- ncol(training_input)
  n_points <- nrow(training_input)
  if ( isTRUE(verb) ) message("Locating the next design point ...", appendLF = FALSE)
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
      res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, core_num = cores)
    }
    dat <- list('location' = x_cand[res[[1]]+1,])
  } else {
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
    if ( identical(cores,as.integer(1)) ){
      res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, core_num = cores)
    }
    dat <- list('index' = res[[1]]+1, 'location' = x_cand[res[[1]]+1,])
  }
  if ( isTRUE(verb) ) message(" done")
  if ( isTRUE(verb) ) message(paste(c("The next design point is:\n", sprintf("%.06f", dat[['location']])), collapse=" "))
  return(dat)
}


#' @rdname locate
#' @method locate dgp
#' @export
locate.dgp <- function(object, x_cand = NULL, n_cand = 200, limits = NULL, method = 'ALM', nugget_s = 1., verb = TRUE, cores = 1, threading = FALSE, aggregate = NULL, ...) {
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

  object$emulator_obj$set_nb_parallel(threading)
  #extract input dimensions
  training_input <- object$data$X
  n_dim_X <- ncol(training_input)
  n_points <- nrow(training_input)
  #extract output dimensions
  n_dim_Y <- ncol(object$data$Y)
  if ( n_dim_Y==1 & !is.null(aggregate) ){
    if ( isTRUE(verb) ) message("The emulator output is one-dimensional, the provided 'aggregate' function will be ignored.")
    aggregate <- NULL
  }
  if ( isTRUE(verb) ) message("Locating the next design point ...", appendLF = FALSE)
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
        res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, core_num = cores)
      }
      dat <- list('location' = x_cand[res[[1]]+1,,drop=FALSE])
    } else {
      if ( identical(cores,as.integer(1)) ){
        res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE, core_num = cores)
      }
      agg_scores <- apply(res, 1, aggregate)
      idx <- which.max(agg_scores)
      dat <- list('location' = x_cand[idx,,drop=FALSE])
    }
  } else {
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
    if ( is.null(aggregate) ){
      if ( identical(cores,as.integer(1)) ){
        res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, core_num = cores)
      }
      dat <- list('index' = res[[1]]+1, 'location' = x_cand[res[[1]]+1,,drop=FALSE])
    } else {
      if ( identical(cores,as.integer(1)) ){
        res = object$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE, core_num = cores)
      }
      agg_scores <- apply(res, 1, aggregate)
      idx <- which.max(agg_scores)
      dat <- list('index' = idx, 'location' = x_cand[idx,,drop=FALSE])
    }
  }
  if ( isTRUE(verb) ) message(" done")
  if ( isTRUE(verb) ) {
    if ( nrow(dat[['location']])==1 ){
       message(paste(c("The next design point is:\n", sprintf("%.06f", dat[['location']][1,])), collapse=" "))
    } else {
      message(paste(c("The next design point is:")))
      for (i in 1:nrow(dat[['location']]) ){
        message(paste(" Output", i, ": ", paste(c(sprintf("%.06f", dat[['location']][i,])), collapse=" "), sep=""))
      }
    }
  }
  return(dat)
}


#' @rdname locate
#' @method locate bundle
#' @export
locate.bundle <- function(object, x_cand = NULL, n_cand = 200, limits = NULL, method = 'ALM', nugget_s = 1., verb = TRUE, cores = 1, threading = FALSE, aggregate = NULL, ...) {
  #check class
  if ( !inherits(object,"bundle") ) stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)

  #check core number
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  training_input <- object$data$X
  n_dim_X <- ncol(training_input)
  n_points <- nrow(training_input)

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
    }
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
  } else {
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
  }

  if ( is.null(aggregate) ) {
    location <- c()
    if ( !x_cand_null ) indice <- c()
  } else {
    score <- c()
  }

  if ( isTRUE(verb) ) message("Locating the next design point ...", appendLF = FALSE)
  n_emulators <- length(object)
  if ( "data" %in% names(object) ) n_emulators <- n_emulators - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1
  for ( i in 1:n_emulators ){
    obj_i <- object[[paste('emulator',i,sep='')]]
    if ( inherits(obj_i,"gp") ){
      if ( is.null(aggregate) ){
        res = obj_i$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s)
        location <- rbind(location, x_cand[res[[1]]+1,])
        if ( !x_cand_null ) indice <- c(indice, res[[1]]+1)
      } else {
        res = obj_i$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s, score_only = TRUE)
        score <- cbind(score, res)
      }
    } else if ( inherits(obj_i,"dgp") ){
      obj_i$emulator_obj$set_nb_parallel(threading)
      #locate next design point
      if ( is.null(aggregate) ){
        if ( identical(cores,as.integer(1)) ){
          res = obj_i$emulator_obj$metric(x_cand = x_cand, method = method, nugget_s = nugget_s)
        } else {
          res = obj_i$emulator_obj$pmetric(x_cand = x_cand, method = method, nugget_s = nugget_s, core_num = cores)
        }
        location <- rbind(location, x_cand[res[[1]]+1,,drop=FALSE])
        if ( !x_cand_null ) indice <- c(indice, res[[1]]+1)
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
    if ( isTRUE(verb) ) {
      message(paste(c("The next design point is:")))
      for (i in 1:nrow(location) ){
        message(paste(" Emulator", i, ": ", paste(c(sprintf("%.06f", location[i,])), collapse=" "), sep=""))
      }
    }
    if ( x_cand_null ) {
      dat <- list('location' = location)
    } else {
      dat <- list('idx' = indice, 'location' = location)
    }
  } else {
    agg_scores <- apply(score, 1, aggregate)
    idx <- which.max(agg_scores)
    if ( x_cand_null ) {
      dat <- list('location' = x_cand[idx,,drop=FALSE])
    } else {
      dat <- list('index' = idx, 'location' = x_cand[idx,,drop=FALSE])
    }
    if ( isTRUE(verb) ) message(paste(c("The next design point is:\n", sprintf("%.06f", dat[['location']][1,])), collapse=" "))
  }
  return(dat)
}






