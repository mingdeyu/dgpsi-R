#' @title Locate the next design point for a (D)GP emulator or a bundle of (D)GP emulators using ALM
#'
#' @description This function searches from a candidate set to locate the next design point(s) to be added to a (D)GP emulator
#'     or a bundle of (D)GP emulators using the Active Learning MacKay (ALM), see the reference below.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param x_cand a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
#'     from which the next design point(s) are determined. If `object` is an instance of the `bundle` class, `x_cand` could also
#'     be a list with the length equal to the number of emulators contained in the `object`. Each slot in `x_cand` is a matrix
#'     that gives a candidate set for each emulator included in the bundle. See *Note* section below for further information.
#' @param batch_size an integer that gives the number of design points to be chosen.
#'     Defaults to `1`.
#' @param M the size of the conditioning set for the Vecchia approximation in the criterion calculation. This argument is only used if the emulator `object`
#'     was constructed under the Vecchia approximation. Defaults to `50`.
#' @param workers the number of processes to be used for the criterion calculation. If set to `NULL`,
#'     the number of processes is set to `max physical cores available %/% 2`. Defaults to `1`.
#' @param aggregate an R function that aggregates scores of the ALM across different output dimensions (if `object` is an instance
#'     of the `dgp` class) or across different emulators (if `object` is an instance of the `bundle` class). The function should be specified in the
#'     following basic form:
#' * the first argument is a matrix representing scores. The rows of the matrix correspond to different design points. The number of columns
#'   of the matrix equals to:
#'   - the emulator output dimension if `object` is an instance of the `dgp` class; or
#'   - the number of emulators contained in `object` if `object` is an instance of the `bundle` class.
#' * the output should be a vector that gives aggregations of scores at different design points.
#'
#' Set to `NULL` to disable the aggregation. Defaults to `NULL`.
#' @param ... any arguments (with names different from those of arguments used in [alm()]) that are used by `aggregate`
#'     can be passed here.
#'
#' @return
#' * If `object` is an instance of the `gp` class, a vector is returned with the length equal to `batch_size`, giving the positions (i.e., row numbers)
#'   of next design points from `x_cand`.
#' * If `object` is an instance of the `dgp` class, a matrix is returned with row number equal to `batch_size` and column number equal to one (if `aggregate`
#'   is not `NULL`) or the output dimension (if `aggregate` is `NULL`), giving positions (i.e., row numbers) of next design points from `x_cand` to be added
#'   to the DGP emulator across different outputs. If `object` is a DGP emulator with either `Hetero` or `NegBin` likelihood layer, the returned matrix has
#'   two columns with the first column giving positions of next design points from `x_cand` that correspond to the mean parameter of the normal or negative Binomial
#'   distribution, and the second column giving positions of next design points from `x_cand` that correspond to the variance parameter of the normal distribution or
#'   the dispersion parameter of the negative Binomial distribution.
#' * If `object` is an instance of the `bundle` class, a matrix is returned with row number equal to `batch_size` and column number equal to the number of
#'   emulators in the bundle, giving positions (i.e., row numbers) of next design points from `x_cand` to be added to individual emulators.
#'
#' @note
#' * The column order of the first argument of `aggregate` must be consistent with the order of emulator output dimensions (if `object` is an instance of the
#'     `dgp` class), or the order of emulators placed in `object` if `object` is an instance of the `bundle` class;
#' * If `x_cand` is supplied as a list when `object` is an instance of `bundle` class and a `aggregate` function is provided, the matrices in `x_cand` must have
#'   common rows (i.e., the candidate sets of emulators in the bundle have common input locations) so the `aggregate` function can be applied.
#' * Any R vector detected in `x_cand` will be treated as a column vector and automatically converted into a single-column
#'   R matrix.
#' @references
#' MacKay, D. J. (1992). Information-based objective functions for active data selection. *Neural Computation*, **4(4)**, 590-604.
#'
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # load packages and the Python env
#' library(lhs)
#' library(dgpsi)
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
#' # locate the next design point using ALM
#' next_point <- alm(m, x_cand = x_cand)
#' X_new <- x_cand[next_point,,drop = F]
#'
#' # obtain the corresponding output at the located design point
#' Y_new <- f(X_new)
#'
#' # combine the new input-output pair to the existing data
#' X <- rbind(X, X_new)
#' Y <- rbind(Y, Y_new)
#'
#' # update the DGP emulator with the new input and output data and refit
#' m <- update(m, X, Y, refit = TRUE)
#'
#' # plot the LOO validation
#' plot(m)
#' }
#' @md
#' @name alm
#' @export
alm <- function(object, x_cand, ...){
  UseMethod("alm")
}

#' @rdname alm
#' @method alm gp
#' @export
alm.gp <- function(object, x_cand, batch_size = 1, M = 50, workers = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
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
  M <- as.integer(M)
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)
  #locate
  if ( batch_size==1 ){
    if ( identical(workers,as.integer(1)) ){
      res = object$emulator_obj$metric(x_cand = x_cand, method = 'ALM', m = M)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = 'ALM', m = M, core_num = workers)
    }
    idx <- res[[1]]+1
  } else {
    idx <- c()
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
    for (i in 1:batch_size){
      if ( identical(workers,as.integer(1)) ){
        res = constructor_obj_cp$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'ALM', m = M)
      } else {
        res = constructor_obj_cp$pmetric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'ALM', m = M, core_num = workers)
      }
      idx_i <- res[[1]]+1
      X_new <- x_cand[idx_x_cand,,drop=F][idx_i,,drop=F]
      Y_new <- constructor_obj_cp$predict(X_new, m = M)[[1]]
      training_input <- rbind(training_input, X_new)
      training_output <- rbind(training_output, Y_new)
      constructor_obj_cp$update_xy(training_input, training_output)
      idx <- c(idx,  idx_x_cand[idx_i])
      idx_x_cand <- idx_x_cand0[-idx]
    }
  }
  if ( batch_size!=1 ){
    pkg.env$py_gc$collect()
    gc(full=T)
  }
  return(idx)
}


#' @rdname alm
#' @method alm dgp
#' @export
alm.dgp <- function(object, x_cand, batch_size = 1, M = 50, workers = 1, aggregate = NULL, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  #if ( object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$type == 'likelihood' ){
  #  stop("The function is only applicable to DGP emulators without likelihood layers.", call. = FALSE)
  #}
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
  M <- as.integer(M)
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
      res = object$emulator_obj$metric(x_cand = x_cand, method = 'ALM', m = M, score_only = TRUE)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = 'ALM', m = M, score_only = TRUE, core_num = workers)
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
        res = emulator_obj_cp$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'ALM', m = M, score_only = TRUE)
      } else {
        res = emulator_obj_cp$pmetric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'ALM', m = M, score_only = TRUE, core_num = workers)
      }
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
      Y_new <- emulator_obj_cp$predict(X_new, m = M)[[1]]

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
  if ( batch_size!=1 ){
    pkg.env$py_gc$collect()
    gc(full=T)
  }
  return(idx)
}


#' @rdname alm
#' @method alm bundle
#' @export
alm.bundle <- function(object, x_cand, batch_size = 1, M = 50, workers = 1, aggregate = NULL, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"bundle") ) stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)
  #check no of emulators
  n_emulators <- length(object)
  if ( "data" %in% names(object) ) n_emulators <- n_emulators - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1
  if ( "id" %in% names(object) ) n_emulators <- n_emulators - 1
  training_input <- object$data$X
  training_output <- object$data$Y
  n_dim_X <- ncol(training_input[[1]])
  #check x_cand
  if ( !is.list(x_cand) ){
    if ( !is.matrix(x_cand) ) {
      if ( !is.vector(x_cand) ) {
        stop("'x_cand' must be a vector, a matrix, or a list.", call. = FALSE)
      }
    }
  }
  if ( is.list(x_cand) ){
    if (length(x_cand) != n_emulators) stop("When 'x_cand' is a list, the number of elements in it should match the number of emulators in the bundle.", call. = FALSE)
    for ( i in 1:n_emulators ){
      if ( is.vector(x_cand[[i]]) ) x_cand[[i]] <- as.matrix(x_cand[[i]])
      if ( ncol(x_cand[[i]])!=n_dim_X ) stop("Elements in 'x_cand' have different number of dimensions with the training input.", call. = FALSE)
    }
    islist <- TRUE
  } else {
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
    islist <- FALSE
  }
  #check core number
  if( !is.null(workers) ) {
    workers <- as.integer(workers)
    if ( workers < 1 ) stop("The worker number must be >= 1.", call. = FALSE)
  }
  M <- as.integer(M)
  #check aggregate
  if ( !is.null(aggregate) ){
    add_arg <- list(...)
    gnames <- methods::formalArgs(aggregate)
    if (islist){
      x_cand_dfs <- lapply(x_cand, as.data.frame)
      x_cand <- as.matrix(Reduce(function(x, y) merge(x, y, all = FALSE), x_cand_dfs))
      if (length(x_cand)==0) stop("Elements in 'x_cand' must have common positions when 'aggregate' is used.", call. = FALSE)
    }
  }
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)

  #locate
  if ( batch_size==1 ){
    scores <- vector('list', n_emulators)
    for ( i in 1:n_emulators ){
      obj_i <- object[[paste('emulator',i,sep='')]]
      if ( inherits(obj_i,"gp") ){
        res = obj_i$emulator_obj$metric(x_cand = if (is.list(x_cand)) {x_cand[[i]]} else {x_cand}, method = 'ALM', m = M, score_only = TRUE)
        scores[[i]] <- res
      } else {
        if ( identical(workers,as.integer(1)) ){
          res = obj_i$emulator_obj$metric(x_cand = if (is.list(x_cand)) {x_cand[[i]]} else {x_cand}, method = 'ALM', m = M, score_only = TRUE)
        } else {
          res = obj_i$emulator_obj$pmetric(x_cand = if (is.list(x_cand)) {x_cand[[i]]} else {x_cand}, method = 'ALM',  m = M, score_only = TRUE, core_num = workers)
        }
        scores[[i]] <- if(ncol(res) == 1) res else rowMeans(res)
      }
    }

    if ( is.null(aggregate) ){
      idx <- sapply(scores, function(x){pkg.env$np$argmax(x, axis=0L) + 1})
    } else {
      scores <- do.call(cbind, scores)
      if ("..." %in% gnames){
        agg_scores <- do.call(aggregate, c(list(scores), add_arg))
      } else {
        gidx <- gnames %in% names(add_arg)
        gparam <- add_arg[gnames[gidx]]
        agg_scores <- do.call(aggregate, c(list(scores), gparam))
      }
      idx <- which.max(agg_scores)
      if (islist){
        chosen_row <- x_cand[idx, ]
        idx <- unlist(lapply(x_cand_dfs, function(df) {
          which(colSums(t(as.matrix(df)) == chosen_row, na.rm = TRUE) == n_dim_X)
        }))
      } else {
        idx <- rep(idx, n_emulators)
      }
    }

    idx <- matrix(idx, nrow = 1, byrow = T)
  } else {
    idx <- c()
    idx_x_cand0 <- if (is.list(x_cand)) {lapply(x_cand, function(x){c(1:nrow(x))})} else replicate(n_emulators, c(1:nrow(x_cand)), simplify = FALSE)
    idx_x_cand <- idx_x_cand0
    constructor_obj_list <- list()
    emulator_obj_list <- list()
    for ( i in 1:n_emulators){
      obj_i <- object[[paste('emulator',i,sep='')]]
      constructor_obj_list[[i]] <- pkg.env$copy$deepcopy(obj_i$constructor_obj)
      emulator_obj_list[[i]] <- pkg.env$copy$deepcopy(obj_i$emulator_obj)
    }

    for (i in 1:batch_size){

      scores <- vector('list', n_emulators)
      for ( j in 1:n_emulators ){
        obj_j <- object[[paste('emulator',j,sep='')]]
        if ( inherits(obj_j,"gp") ){
          res = emulator_obj_list[[j]]$metric(x_cand = if (is.list(x_cand)) {x_cand[[j]][idx_x_cand[[j]],,drop=F]} else {x_cand[idx_x_cand[[j]],,drop=F]}, method = 'ALM', m = M, score_only = TRUE)
          scores[[j]] <- res
        } else {
          if ( identical(workers,as.integer(1)) ){
            res = emulator_obj_list[[j]]$metric(x_cand = if (is.list(x_cand)) {x_cand[[j]][idx_x_cand[[j]],,drop=F]} else {x_cand[idx_x_cand[[j]],,drop=F]}, method = 'ALM', m = M, score_only = TRUE)
          } else {
            res = emulator_obj_list[[j]]$pmetric(x_cand = if (is.list(x_cand)) {x_cand[[j]][idx_x_cand[[j]],,drop=F]} else {x_cand[idx_x_cand[[j]],,drop=F]}, method = 'ALM', m = M, score_only = TRUE, core_num = workers)
          }
          scores[[j]] <- if(ncol(res) == 1) res else rowMeans(res)
        }
      }

      if ( is.null(aggregate) ){
        idx_i <- sapply(scores, function(x){pkg.env$np$argmax(x, axis=0L) + 1})
      } else {
        scores <- do.call(cbind, scores)
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

      if (is.list(x_cand)){
        X_new <- lapply(1:n_emulators, function(k) {x_cand[[k]][idx_x_cand[[k]],,drop=F][idx_i[k],,drop=F]})
      } else {
        X_new <- lapply(1:n_emulators, function(k) {x_cand[idx_x_cand[[k]],,drop=F][idx_i[k],,drop=F]})
      }

      for ( j in 1:n_emulators ){
        obj_j <- object[[paste('emulator',j,sep='')]]
        Y_new_j <- emulator_obj_list[[j]]$predict(X_new[[j]], m = M)[[1]]
        if ( inherits(obj_j,"gp") ){
          emulator_obj_list[[j]]$update_xy(rbind(training_input[[j]], X_new[[j]]), rbind(training_output[[j]], Y_new_j))
        } else {
          B <- as.integer(length(obj_j$emulator_obj$all_layer_set))
          burnin <- obj_j$constructor_obj$burnin
          isblock <- obj_j$constructor_obj$block
          constructor_obj_list[[j]]$update_xy(rbind(training_input[[j]], X_new[[j]]), rbind(training_output[[j]], Y_new_j))
          est_obj <- constructor_obj_list[[j]]$estimate(burnin)
          emulator_obj_list[[j]] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
        }
      }
      if ( !is.null(aggregate) & islist ){
        chosen_row <- x_cand[idx_x_cand[[1]],,drop=F][idx_i[1], ]
        idx <- rbind( idx, unlist(lapply(x_cand_dfs, function(df) {
          which(colSums(t(as.matrix(df)) == chosen_row, na.rm = TRUE) == n_dim_X)
        })))
        idx_x_cand <- replicate(n_emulators, idx_x_cand[[1]][-idx_i[1]], simplify = FALSE)
      } else {
        idx <- rbind( idx,  sapply(1:n_emulators, function(k){idx_x_cand[[k]][idx_i[k]]}) )
        idx_x_cand <- lapply(1:n_emulators, function(k){idx_x_cand0[[k]][-idx[,k]]})
      }
    }
    #idx <- matrix(idx, nrow = batch_size, byrow = T)
  }
  if ( batch_size!=1 ){
    pkg.env$py_gc$collect()
    gc(full=T)
  }
  return(idx)
}
