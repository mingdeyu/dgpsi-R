#' @title Locate the next design point for a (D)GP emulator or a bundle of (D)GP emulators using MICE
#'
#' @description This function searches from a candidate set to locate the next design point(s) to be added to a (D)GP emulator
#'     or a bundle of (D)GP emulators using the Mutual Information for Computer Experiments (MICE), see the reference below.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param x_cand a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
#'     from which the next design point(s) are determined.
#' @param batch_size an integer that gives the number of design points to be chosen.
#'     Defaults to `1`.
#' @param nugget_s the value of the smoothing nugget term used by MICE. Defaults to `1e-6`.
#' @param workers the number of workers/cores to be used for the criterion calculation. If set to `NULL`,
#'     the number of workers is set to `(max physical cores available - 1)`. Defaults to `1`.
#' @param threading a bool indicating whether to use the multi-threading to accelerate the criterion calculation for a DGP emulator.
#'     Turning this option on could improve the speed of criterion calculations when the DGP emulator is built with a moderately large number of
#'     training data points and the Matérn-2.5 kernel.
#' @param aggregate an R function that aggregates scores of the MICE across different output dimensions (if `object` is an instance
#'     of the `dgp` class) or across different emulators (if `object` is an instance of the `bundle` class). The function should be specified in the
#'     following basic form:
#' * the first argument is a matrix representing scores. The rows of the matrix correspond to different design points. The number of columns
#'   of the matrix equals to:
#'   - the emulator output dimension if `object` is an instance of the `dgp` class; or
#'   - the number of emulators contained in `object` if `object` is an instance of the `bundle` class.
#' * the output should be a vector that gives aggregations of scores at different design points.
#'
#' Set to `NULL` to disable the aggregation. Defaults to `NULL`.
#' @param ... any arguments (with names different from those of arguments used in [mice()]) that are used by `aggregate`
#'     can be passed here.
#'
#' @return
#' * If `object` is an instance of the `gp` class, a vector is returned with the length equal to `batch_size`, giving the positions (i.e., row numbers)
#'   of next design points from `x_cand`.
#' * If `object` is an instance of the `dgp` class, a matrix is returned with row number equal to `batch_size` and column number equal to one (if `aggregate`
#'   is not `NULL`) or the output dimension (if `aggregate` is `NULL`), giving positions (i.e., row numbers) of next design points from `x_cand` to be added
#'   to the DGP emulator across different outputs.
#' * If `object` is an instance of the `bundle` class, a matrix is returned with row number equal to `batch_size` and column number equal to the number of
#'   emulators in the bundle, giving positions (i.e., row numbers) of next design points from `x_cand` to be added to individual emulators.
#'
#' @note
#' * The column order of the first argument of `aggregate` must be consistent with the order of emulator output dimensions (if `object` is an instance of the
#'     `dgp` class), or the order of emulators placed in `object` if `object` is an instance of the `bundle` class;
#' * The function is only applicable to DGP emulators without likelihood layers.
#' * Any R vector detected in `x_cand` will be treated as a column vector and automatically converted into a single-column
#'   R matrix.
#' @references
#' Beck, J., & Guillas, S. (2016). Sequential design with mutual information for computer experiments (MICE): emulation of a tsunami model.
#' *SIAM/ASA Journal on Uncertainty Quantification*, **4(1)**, 739-766.
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
#' # generate a candidate set
#' x_cand <- maximinLHS(200,1)
#'
#' # locate the next design point using MICE
#' next_point <- mice(m, x_cand = x_cand)
#' X_new <- x_cand[next_point,,drop = F]
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
#' @name mice
#' @export
mice <- function(object, x_cand, ...){
  UseMethod("mice")
}

#' @rdname mice
#' @method mice gp
#' @export
mice.gp <- function(object, x_cand, batch_size = 1, nugget_s = 1e-6, workers = 1, ...) {
  #check class
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  training_input <- object$data$X
  training_output <- object$data$Y
  n_dim_X <- ncol(training_input)
  #check x_cand
  if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
  if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
  #check core number
  if( !is.null(workers) ) {
    workers <- as.integer(workers)
    if ( workers < 1 ) stop("The worker number must be >= 1.", call. = FALSE)
  }
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)
  #locate
  if ( batch_size==1 ){
    if ( identical(workers,as.integer(1)) ){
      res = object$emulator_obj$metric(x_cand = x_cand, method = 'MICE', nugget_s = nugget_s)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = 'MICE', nugget_s = nugget_s, core_num = workers)
    }
    idx <- res[[1]]+1
  } else {
    idx <- c()
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
    for (i in 1:batch_size){
      if ( identical(workers,as.integer(1)) ){
        res = constructor_obj_cp$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', nugget_s = nugget_s)
      } else {
        res = constructor_obj_cp$pmetric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', nugget_s = nugget_s, core_num = workers)
      }
      idx_i <- res[[1]]+1
      X_new <- x_cand[idx_x_cand,,drop=F][idx_i,,drop=F]
      Y_new <- constructor_obj_cp$predict(X_new)[[1]]
      training_input <- rbind(training_input, X_new)
      training_output <- rbind(training_output, Y_new)
      constructor_obj_cp$update_xy(training_input, training_output)
      idx <- c(idx,  idx_x_cand[idx_i])
      idx_x_cand <- idx_x_cand0[-idx]
    }
  }
  pkg.env$py_gc$collect()
  gc(full=T)
  return(idx)
}

#' @rdname mice
#' @method mice dgp
#' @export
mice.dgp <- function(object, x_cand, batch_size = 1, nugget_s = 1e-6, workers = 1, threading = FALSE, aggregate = NULL, ...) {
  #check class
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$type == 'likelihood' ){
    stop("The function is only applicable to DGP emulators without likelihood layers.", call. = FALSE)
  }
  object$emulator_obj$set_nb_parallel(threading)
  training_input <- object$data$X
  training_output <- object$data$Y
  n_dim_X <- ncol(training_input)
  n_dim_Y <- ncol(training_output)
  #check x_cand
  if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
  if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
  #check core number
  if( !is.null(workers) ) {
    workers <- as.integer(workers)
    if ( workers < 1 ) stop("The worker number must be >= 1.", call. = FALSE)
  }
  #check aggregate
  if ( !is.null(aggregate) ){
    add_arg <- list(...)
    gnames <- methods::formalArgs(aggregate)
  }
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)
  #locate
  if ( batch_size==1 ){
    if ( identical(workers,as.integer(1)) ){
      res = object$emulator_obj$metric(x_cand = x_cand, method = 'MICE', nugget_s = nugget_s, score_only = TRUE)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = 'MICE', nugget_s = nugget_s, score_only = TRUE, core_num = workers)
    }
    if ( is.null(aggregate) ){
      idx <- pkg.env$np$argmax(res, axis=0L) + 1
    } else {
      if ( ncol(res)==1 ){
        idx <- pkg.env$np$argmax(res, axis=0L) + 1
      } else {
        if ("..." %in% gnames){
          agg_res <- do.call(aggregate, c(list(res), add_arg))
        } else {
          gidx <- gnames %in% names(add_arg)
          gparam <- add_arg[gnames[gidx]]
          agg_res <- do.call(aggregate, c(list(res), gparam))
        }
        idx <- which.max(agg_res)
      }
    }
    idx <- matrix(idx, nrow = 1, byrow = T)
  } else {
    idx <- c()
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
    emulator_obj_cp <- pkg.env$copy$deepcopy(object$emulator_obj)
    B <- as.integer(length(emulator_obj_cp$all_layer_set))
    burnin <- constructor_obj_cp$burnin
    isblock <- constructor_obj_cp$block
    for (i in 1:batch_size){
      if ( identical(workers,as.integer(1)) ){
        res = emulator_obj_cp$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', nugget_s = nugget_s, score_only = TRUE)
      } else {
        res = emulator_obj_cp$pmetric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', nugget_s = nugget_s, score_only = TRUE, core_num = workers)
      }
      emulator_obj_cp$set_nb_parallel(FALSE)
      if ( is.null(aggregate) ){
        idx_i <- pkg.env$np$argmax(res, axis=0L) + 1
      } else {
        if ( ncol(res)==1 ){
          idx_i <- pkg.env$np$argmax(res, axis=0L) + 1
        } else {
          if ("..." %in% gnames){
            agg_res <- do.call(aggregate, c(list(res), add_arg))
          } else {
            gidx <- gnames %in% names(add_arg)
            gparam <- add_arg[gnames[gidx]]
            agg_res <- do.call(aggregate, c(list(res), gparam))
          }
          idx_i <- which.max(agg_res)
        }
      }
      X_new <- x_cand[idx_x_cand,,drop=F][unique(idx_i),,drop=F]
      Y_new <- emulator_obj_cp$predict(X_new)[[1]]

      training_input <- rbind(training_input, X_new)
      training_output <- rbind(training_output, Y_new)
      constructor_obj_cp$update_xy(training_input, training_output)

      est_obj <- constructor_obj_cp$estimate(burnin)
      emulator_obj_cp <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)

      idx <- c(idx,  idx_x_cand[idx_i])
      idx_x_cand <- idx_x_cand0[-unique(idx)]
    }
    idx <- matrix(idx, nrow = batch_size, byrow = T)
  }
  pkg.env$py_gc$collect()
  gc(full=T)
  return(idx)
}


#' @rdname mice
#' @method mice bundle
#' @export
mice.bundle <- function(object, x_cand, batch_size = 1, nugget_s = 1e-6, workers = 1, threading = FALSE, aggregate = NULL, ...) {
  #check class
  if ( !inherits(object,"bundle") ) stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)
  #check no of emulators
  n_emulators <- length(object)
  if ( "data" %in% names(object) ) n_emulators <- n_emulators - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1
  training_input <- object$data$X
  training_output <- object$data$Y
  n_dim_X <- ncol(training_input[[1]])
  #check x_cand
  if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
  if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
  #check core number
  if( !is.null(workers) ) {
    workers <- as.integer(workers)
    if ( workers < 1 ) stop("The worker number must be >= 1.", call. = FALSE)
  }
  #check aggregate
  if ( !is.null(aggregate) ){
    add_arg <- list(...)
    gnames <- methods::formalArgs(aggregate)
  }
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)

  #locate
  if ( batch_size==1 ){
    scores <- c()
    for ( i in 1:n_emulators ){
      obj_i <- object[[paste('emulator',i,sep='')]]
      if ( inherits(obj_i,"gp") ){
        res = obj_i$emulator_obj$metric(x_cand = x_cand, method = 'MICE', nugget_s = nugget_s, score_only = TRUE)
        scores <- cbind(scores, res)
      } else {
        obj_i$emulator_obj$set_nb_parallel(threading)
        if ( identical(workers,as.integer(1)) ){
          res = obj_i$emulator_obj$metric(x_cand = x_cand, method = 'MICE', nugget_s = nugget_s, score_only = TRUE)
        } else {
          res = obj_i$emulator_obj$pmetric(x_cand = x_cand, method = 'MICE', nugget_s = nugget_s, score_only = TRUE, core_num = workers)
        }
        scores <- cbind(scores, res)
      }
    }

    if ( is.null(aggregate) ){
      idx <- pkg.env$np$argmax(scores, axis=0L) + 1
    } else {
      if ("..." %in% gnames){
        agg_scores <- do.call(aggregate, c(list(scores), add_arg))
      } else {
        gidx <- gnames %in% names(add_arg)
        gparam <- add_arg[gnames[gidx]]
        agg_scores <- do.call(aggregate, c(list(scores), gparam))
      }
      idx <- which.max(agg_scores)
      idx <- rep(idx, n_emulators)
    }

    idx <- matrix(idx, nrow = 1, byrow = T)
  } else {
    idx <- c()
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    constructor_obj_list <- list()
    emulator_obj_list <- list()
    for ( i in 1:n_emulators){
      obj_i <- object[[paste('emulator',i,sep='')]]
      constructor_obj_list[[i]] <- pkg.env$copy$deepcopy(obj_i$constructor_obj)
      emulator_obj_list[[i]] <- pkg.env$copy$deepcopy(obj_i$emulator_obj)
    }

    for (i in 1:batch_size){

      scores <- c()
      for ( j in 1:n_emulators ){
        obj_j <- object[[paste('emulator',j,sep='')]]
        if ( inherits(obj_j,"gp") ){
          res = emulator_obj_list[[j]]$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', nugget_s = nugget_s, score_only = TRUE)
          scores <- cbind(scores, res)
        } else {
          emulator_obj_list[[j]]$set_nb_parallel(threading)
          if ( identical(workers,as.integer(1)) ){
            res = emulator_obj_list[[j]]$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', nugget_s = nugget_s, score_only = TRUE)
          } else {
            res = emulator_obj_list[[j]]$pmetric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', nugget_s = nugget_s, score_only = TRUE, core_num = workers)
          }
          emulator_obj_list[[j]]$set_nb_parallel(FALSE)
          scores <- cbind(scores, res)
        }
      }

      if ( is.null(aggregate) ){
        idx_i <- pkg.env$np$argmax(scores, axis=0L) + 1
      } else {
        if ("..." %in% gnames){
          agg_scores <- do.call(aggregate, c(list(scores), add_arg))
        } else {
          gidx <- gnames %in% names(add_arg)
          gparam <- add_arg[gnames[gidx]]
          agg_scores <- do.call(aggregate, c(list(scores), gparam))
        }
        idx_i <- which.max(agg_scores)
        idx_i <- rep(idx_i, n_emulators)
      }

      X_new <- x_cand[idx_x_cand,,drop=F][idx_i,,drop=F]

      for ( j in 1:n_emulators ){
        obj_j <- object[[paste('emulator',j,sep='')]]
        Y_new_j <- emulator_obj_list[[j]]$predict(X_new[j,,drop=F])[[1]]
        if ( inherits(obj_j,"gp") ){
          emulator_obj_list[[j]]$update_xy(rbind(training_input[[j]], X_new[j,]), rbind(training_output[[j]], Y_new_j))
        } else {
          B <- as.integer(length(obj_j$emulator_obj$all_layer_set))
          burnin <- obj_j$constructor_obj$burnin
          isblock <- obj_j$constructor_obj$block
          constructor_obj_list[[j]]$update_xy(rbind(training_input[[j]], X_new[j,]), rbind(training_output[[j]], Y_new_j))
          est_obj <- constructor_obj_list[[j]]$estimate(burnin)
          emulator_obj_list[[j]] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
        }
      }
      idx <- c(idx,  idx_x_cand[idx_i])
      idx_x_cand <- idx_x_cand0[-unique(idx)]
    }
    idx <- matrix(idx, nrow = batch_size, byrow = T)
  }
  pkg.env$py_gc$collect()
  gc(full=T)
  return(idx)
}

