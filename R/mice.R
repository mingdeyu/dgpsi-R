#' @title Locate the next design point for a (D)GP emulator or a bundle of (D)GP emulators using MICE
#'
#' @description
#'
#' `r new_badge("updated")`
#'
#' This function searches from a candidate set to locate the next design point(s) to be added to a (D)GP emulator
#'     or a bundle of (D)GP emulators using the Mutual Information for Computer Experiments (MICE), see the reference below.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param x_cand a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
#'     from which the next design point(s) are determined. If `object` is an instance of the `bundle` class and `aggregate` is not supplied, `x_cand` can also be a list.
#'     The list must have a length equal to the number of emulators in `object`, with each element being a matrix representing the candidate set for a corresponding
#'     emulator in the bundle. Defaults to `NULL`.
#' @param n_cand an integer specifying the size of the candidate set to be generated for selecting the next design point(s).
#'     This argument is used only when `x_cand` is `NULL`. Defaults to `200`.
#' @param batch_size an integer that gives the number of design points to be chosen.
#'     Defaults to `1`.
#' @param M `r new_badge("new")` the size of the conditioning set for the Vecchia approximation in the criterion calculation. This argument is only used if the emulator `object`
#'     was constructed under the Vecchia approximation. Defaults to `50`.
#' @param nugget_s the value of the smoothing nugget term used by MICE. Defaults to `1e-6`.
#' @param workers  the number of processes to be used for the criterion calculation. If set to `NULL`,
#'     the number of processes is set to `max physical cores available %/% 2`. Defaults to `1`.
#' @param limits `r new_badge("new")` a two-column matrix that gives the ranges of each input dimension, or a vector of length two if there is only one input dimension.
#'     If a vector is provided, it will be converted to a two-column row matrix. The rows of the matrix correspond to input dimensions, and its
#'     first and second columns correspond to the minimum and maximum values of the input dimensions. This
#'     argument is only used when `x_cand = NULL`. Defaults to `NULL`.
#' @param int `r new_badge("new")` a bool or a vector of bools that indicates if an input dimension is an integer type. If a single bool is given, it will be applied to
#'     all input dimensions. If a vector is provided, it should have a length equal to the input dimensions and will be applied to individual
#'     input dimensions. This argument is only used when `x_cand = NULL`. Defaults to `FALSE`.
#' @param aggregate an R function that aggregates scores of the MICE across different output dimensions (if `object` is an instance
#'     of the `dgp` class) or across different emulators (if `object` is an instance of the `bundle` class). The function should be specified in the
#'     following basic form:
#' * the first argument is a matrix representing scores. The rows of the matrix correspond to different design points. The number of columns
#'   of the matrix equals to:
#'   - the emulator output dimension if `object` is an instance of the `dgp` class; or
#'   - the number of emulators contained in `object` if `object` is an instance of the `bundle` class.
#' * the output should be a vector that gives aggregate scores at different design points.
#'
#' Set to `NULL` to disable aggregation. Defaults to `NULL`.
#' @param ... any arguments (with names different from those of arguments used in [mice()]) that are used by `aggregate`
#'     can be passed here.
#'
#' @return
#' 1. If `x_cand` is not `NULL`:
#'    - When `object` is an instance of the `gp` class, a vector of length `batch_size` is returned, containing the positions
#'      (row numbers) of the next design points from `x_cand`.
#'    - When `object` is an instance of the `dgp` class, a vector of length `batch_size * D` is returned, containing the positions
#'      (row numbers) of the next design points from `x_cand` to be added to the DGP emulator.
#'      * `D` is the number of output dimensions of the DGP emulator if no likelihood layer is included.
#'      * For a DGP emulator with a `Hetero` or `NegBin` likelihood layer, `D = 2`.
#'      * For a DGP emulator with a `Categorical` likelihood layer, `D = 1` for binary output or `D = K` for multi-class output with `K` classes.
#'    - When `object` is an instance of the `bundle` class, a matrix is returned with `batch_size` rows and a column for each emulator in
#'      the bundle, containing the positions (row numbers) of the next design points from `x_cand` for individual emulators.
#' 2. If `x_cand` is `NULL`:
#'    - When `object` is an instance of the `gp` class, a matrix with `batch_size` rows is returned, giving the next design points to be evaluated.
#'    - When `object` is an instance of the `dgp` class, a matrix with `batch_size * D` rows is returned, where:
#'      - `D` is the number of output dimensions of the DGP emulator if no likelihood layer is included.
#'      - For a DGP emulator with a `Hetero` or `NegBin` likelihood layer, `D = 2`.
#'      - For a DGP emulator with a `Categorical` likelihood layer, `D = 1` for binary output or `D = K` for multi-class output with `K` classes.
#'    - When `object` is an instance of the `bundle` class, a list is returned with a length equal to the number of emulators in the bundle. Each
#'      element of the list is a matrix with `batch_size` rows, where each row represents a design point to be added to the corresponding emulator.
#'
#' @note
#' The first column of the matrix supplied to the first argument of `aggregate` must correspond to the first output dimension of the DGP emulator
#'     if `object` is an instance of the `dgp` class, and so on for subsequent columns and dimensions. If `object` is an instance of the `bundle` class,
#'     the first column must correspond to the first emulator in the bundle, and so on for subsequent columns and emulators.
#' @references
#' Beck, J., & Guillas, S. (2016). Sequential design with mutual information for computer experiments (MICE): emulation of a tsunami model.
#' *SIAM/ASA Journal on Uncertainty Quantification*, **4(1)**, 739-766.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
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
#' # update the DGP emulator with the new input and output data and refit
#' m <- update(m, X, Y, refit = TRUE)
#'
#' # plot the LOO validation
#' plot(m)
#' }
#' @md
#' @name mice
#' @export
mice <- function(object, ...){
  UseMethod("mice")
}

#' @rdname mice
#' @method mice gp
#' @export
mice.gp <- function(object, x_cand = NULL, n_cand = 200, batch_size = 1, M = 50, nugget_s = 1e-6, workers = 1, limits = NULL, int = FALSE, ...) {
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
  if ( is.null(x_cand) ){
    is_cand <- FALSE
    limits <- check_limits(limits, n_dim_X)
    #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input, limits), n_cand), n_cand), limits)
    #x_cand <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
    x_cand <- lhs::maximinLHS(n_cand,n_dim_X)
    for (j in 1:n_dim_X){
      if ( int[j] ){
        if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
        x_cand[,j] <- floor( x_cand[,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
      } else {
        x_cand[,j] <- x_cand[,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
      }
    }
  } else {
    is_cand <- TRUE
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) {
      if ( ncol(object$data$X)!=1 ){
        x_cand <- matrix(x_cand, nrow = 1)
      } else {
        x_cand <- as.matrix(x_cand)
      }
    }
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
  }
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
      res = object$emulator_obj$metric(x_cand = x_cand, method = 'MICE', m = M, nugget_s = nugget_s)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = 'MICE', m = M, nugget_s = nugget_s, core_num = workers)
    }
    idx <- res[[1]]+1
  } else {
    idx <- c()
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
    for (i in 1:batch_size){
      if ( identical(workers,as.integer(1)) ){
        res = constructor_obj_cp$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', m = M, nugget_s = nugget_s)
      } else {
        res = constructor_obj_cp$pmetric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', m = M, nugget_s = nugget_s, core_num = workers)
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
  if (is_cand){
    return(idx)
  } else {
    final_res <- x_cand[idx,,drop=F]
    rownames(final_res) <- NULL
    return(final_res)
  }
}

#' @rdname mice
#' @method mice dgp
#' @export
mice.dgp <- function(object, x_cand = NULL, n_cand = 200, batch_size = 1, M = 50, nugget_s = 1e-6, workers = 1, limits = NULL, int = FALSE, aggregate = NULL, ...) {
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
  if ( is.null(x_cand) ){
    is_cand <- FALSE
    limits <- check_limits(limits, n_dim_X)
    int <- check_int(int, n_dim_X)
    #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input, limits), n_cand), n_cand), limits)
    x_cand <- lhs::maximinLHS(n_cand,n_dim_X)
    for (j in 1:n_dim_X){
      if ( int[j] ){
        if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
        x_cand[,j] <- floor( x_cand[,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
      } else {
        x_cand[,j] <- x_cand[,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
      }
    }
  } else {
    is_cand <- TRUE
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) {
      if ( ncol(object$data$X)!=1 ){
        x_cand <- matrix(x_cand, nrow = 1)
      } else {
        x_cand <- as.matrix(x_cand)
      }
    }
    if ( ncol(x_cand)!=n_dim_X ) stop("'x_cand' and the training input have different number of dimensions.", call. = FALSE)
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
  }
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)
  #locate
  if ( batch_size==1 ){
    if ( identical(workers,as.integer(1)) ){
      res = object$emulator_obj$metric(x_cand = x_cand, method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE)
    } else {
      res = object$emulator_obj$pmetric(x_cand = x_cand, method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE, core_num = workers)
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
    #idx <- matrix(idx, nrow = 1, byrow = T)
    idx <- as.vector(idx)
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
        res = emulator_obj_cp$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE)
      } else {
        res = emulator_obj_cp$pmetric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE, core_num = workers)
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
    #idx <- matrix(idx, nrow = batch_size, byrow = T)
  }
  if ( batch_size!=1 ){
    pkg.env$py_gc$collect()
    gc(full=T)
  }
  if (is_cand) {
    return(idx)
  } else {
    final_res <- x_cand[idx,,drop=F]
    rownames(final_res) <- NULL
    return(final_res)
  }
}


#' @rdname mice
#' @method mice bundle
#' @export
mice.bundle <- function(object, x_cand = NULL, n_cand = 200, batch_size = 1, M = 50, nugget_s = 1e-6, workers = 1, limits = NULL, int = FALSE, aggregate = NULL, ...) {
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
  if ( is.null(x_cand) ){
    is_cand <- FALSE
    limits <- check_limits(limits, n_dim_X)
    int <- check_int(int, n_dim_X)
    if ( is.null(aggregate) ) {
      x_cand <- vector('list', n_emulators)
      for (k in 1:n_emulators){
        #x_cand[[i]] <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input[[i]], limits), n_cand), n_cand), limits)
        #x_cand[[i]] <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
        x_cand[[k]] <- lhs::maximinLHS(n_cand,n_dim_X)
        for (j in 1:n_dim_X){
          if ( int[j] ){
            if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
            x_cand[[k]][,j] <- floor( x_cand[[k]][,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
          } else {
            x_cand[[k]][,j] <- x_cand[[k]][,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
          }
        }
      }
    } else {
      #total_input <- do.call(rbind, training_input)
      #total_input <- pkg.env$np$unique(total_input, axis=0L)
      #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(total_input, limits), n_cand), n_cand), limits)
      #x_cand <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
      x_cand <- lhs::maximinLHS(n_cand,n_dim_X)
      for (j in 1:n_dim_X){
        if ( int[j] ){
          if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
          x_cand[,j] <- floor( x_cand[,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
        } else {
          x_cand[,j] <- x_cand[,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
        }
      }
    }
  } else {
    is_cand <- TRUE
  }

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
      if ( is.vector(x_cand[[i]]) ) {
        if ( n_dim_X!=1 ){
          x_cand[[i]] <- matrix(x_cand[[i]], nrow = 1)
        } else {
          x_cand[[i]] <- as.matrix(x_cand[[i]])
        }
      }
      if ( ncol(x_cand[[i]])!=n_dim_X ) stop("Elements in 'x_cand' have different number of dimensions with the training input.", call. = FALSE)
    }
    islist <- TRUE
  } else {
    if ( is.vector(x_cand) ) {
      if ( n_dim_X!=1 ){
        x_cand <- matrix(x_cand, nrow = 1)
      } else {
        x_cand <- as.matrix(x_cand)
      }
    }
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
      if (length(x_cand)==0) stop("When using 'aggregate,' matrices in 'x_cand' must share at least some common rows.", call. = FALSE)
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
        res = obj_i$emulator_obj$metric(x_cand = if (is.list(x_cand)) {x_cand[[i]]} else {x_cand}, method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE)
        scores[[i]] <- res
      } else {
        if ( identical(workers,as.integer(1)) ){
          res = obj_i$emulator_obj$metric(x_cand = if (is.list(x_cand)) {x_cand[[i]]} else {x_cand}, method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE)
        } else {
          res = obj_i$emulator_obj$pmetric(x_cand = if (is.list(x_cand)) {x_cand[[i]]} else {x_cand}, method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE, core_num = workers)
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
          res = emulator_obj_list[[j]]$metric(x_cand = if (is.list(x_cand)) {x_cand[[j]][idx_x_cand[[j]],,drop=F]} else {x_cand[idx_x_cand[[j]],,drop=F]}, method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE)
          scores[[j]] <- res
        } else {
          if ( identical(workers,as.integer(1)) ){
            res = emulator_obj_list[[j]]$metric(x_cand = if (is.list(x_cand)) {x_cand[[j]][idx_x_cand[[j]],,drop=F]} else {x_cand[idx_x_cand[[j]],,drop=F]}, method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE)
          } else {
            res = emulator_obj_list[[j]]$pmetric(x_cand = if (is.list(x_cand)) {x_cand[[j]][idx_x_cand[[j]],,drop=F]} else {x_cand[idx_x_cand[[j]],,drop=F]}, method = 'MICE', m = M, nugget_s = nugget_s, score_only = TRUE, core_num = workers)
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
  if (is_cand){
    return(idx)
  } else {
    final_res <- vector('list', n_emulators)
    for (i in 1:n_emulators){
      final_res[[i]] <- if (is.null(aggregate)) {x_cand[[i]][idx[,i],,drop=F]} else {x_cand[idx[,i],,drop=F]}
      rownames(final_res[[i]]) <- NULL
    }
    return(unname(final_res))
  }
}

