#' @title Locate the next design point(s) for a (D)GP emulator or a bundle of (D)GP emulators using Active Learning MacKay (ALM)
#'
#' @description This function searches from a candidate set to locate the next design point(s) to be added to a (D)GP emulator
#'     or a bundle of (D)GP emulators using the Active Learning MacKay (ALM) criterion (see the reference below).
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param x_cand a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
#'     from which the next design point(s) are determined. If `object` is an instance of the `bundle` class and `aggregate` is not supplied, `x_cand` can also be a list.
#'     The list must have a length equal to the number of emulators in `object`, with each element being a matrix representing the candidate set for a corresponding
#'     emulator in the bundle. Defaults to `NULL`.
#' @param n_start the maximum number of local refinement starts to determine next design point(s)
#'     when `x_cand = NULL`. Defaults to `5`.
#' @param n_screen the number of points generated for screening when `x_cand = NULL`.
#'     Defaults to `100` times the number of input dimensions. Must be an integer
#'     at least `n_start`.
#' @param opt_control a named list controlling local refinement when `x_cand = NULL`.
#'     Supported entries are `maxit` (default `100`), `ftol` (nonnegative
#'     function tolerance, default `1e7 * .Machine$double.eps`), `gtol` (nonnegative
#'     projected-gradient tolerance, default `0`), and `step` (positive scalar numerical
#'     derivative step in normalized coordinates, default `NULL` for the solver's native step).
#' @param batch_size an integer that gives the number of design points to be chosen. Defaults to `1`.
#' @param M the size of the conditioning set for the Vecchia approximation in the criterion calculation. This argument is only used if the emulator `object`
#'     was constructed under the Vecchia approximation. Defaults to `50`.
#' @param limits a two-column matrix that gives the ranges of each input dimension, or a vector of length two if there is only one input dimension.
#'     If a vector is provided, it will be converted to a two-column row matrix. The rows of the matrix correspond to input dimensions, and its
#'     first and second columns correspond to the minimum and maximum values of the input dimensions. This
#'     argument is only used when `x_cand = NULL`. Defaults to `NULL`.
#' @param int a bool or a vector of bools that indicates if an input dimension is an integer type. If a single bool is given, it will be applied to
#'     all input dimensions. If a vector is provided, it should have a length equal to the input dimensions and will be applied to individual
#'     input dimensions. This argument is only used when `x_cand = NULL`. Defaults to `FALSE`.
#' @param aggregate an R function that aggregates scores of the ALM across different output dimensions (if `object` is an instance
#'     of the `dgp` class) or across different emulators (if `object` is an instance of the `bundle` class). The function should be specified in the
#'     following basic form:
#' * the first argument is a matrix representing scores. The rows of the matrix correspond to different design points. The number of columns
#'   of the matrix is equal to:
#'   - the emulator output dimension if `object` is an instance of the `dgp` class; or
#'   - the number of emulators contained in `object` if `object` is an instance of the `bundle` class.
#' * the output should be a vector that gives aggregate scores at different design points.
#'
#' Set to `NULL` to disable aggregation. Defaults to `NULL`.
#' @param ... any arguments (with names different from those of arguments used in [alm()]) that are used by `aggregate`
#'     can be passed here.
#'
#' @return
#' 1. If `x_cand` is not `NULL`:
#'    - When `object` is an instance of the `gp` class, a vector of length `batch_size` is returned, containing the positions
#'      (row numbers) of the next design points from `x_cand`.
#'    - When `object` is an instance of the `dgp` class, a vector of length `batch_size * D` is returned, containing the positions
#'      (row numbers) of the next design points from `x_cand` to be added to the DGP emulator.
#'      * `D` is the number of output dimensions of the DGP emulator if no likelihood layer is included.
#'      * For a DGP emulator with a `Hetero` or `NegBin` or `ZIP` likelihood layer, `D = 2`.
#'      * For a DGP emulator with a `ZINB` likelihood layer, `D = 3`.
#'      * For a DGP emulator with a `Categorical` likelihood layer, `D = 1` for binary output or `D = K` for multi-class output with `K` classes.
#'    - When `object` is an instance of the `bundle` class, a matrix is returned with `batch_size` rows and a column for each emulator in
#'      the bundle, containing the positions (row numbers) of the next design points from `x_cand` for individual emulators.
#' 2. If `x_cand` is `NULL`:
#'    - When `object` is an instance of the `gp` class, a matrix with `batch_size` rows is returned, giving the next design points to be evaluated.
#'    - When `object` is an instance of the `dgp` class, a matrix with `batch_size * D` rows is returned, where:
#'      - `D` is the number of output dimensions of the DGP emulator if no likelihood layer is included.
#'      - For a DGP emulator with a `Hetero` or `NegBin` or `ZIP` likelihood layer, `D = 2`.
#'      - For a DGP emulator with a `ZINB` likelihood layer, `D = 3`.
#'      - For a DGP emulator with a `Categorical` likelihood layer, `D = 1` for binary output or `D = K` for multi-class output with `K` classes.
#'    - When `object` is an instance of the `bundle` class, a list is returned with a length equal to the number of emulators in the bundle. Each
#'      element of the list is a matrix with `batch_size` rows, where each row represents a design point to be added to the corresponding emulator.
#'
#' @note
#' The first column of the matrix supplied to the first argument of `aggregate` must correspond to the first output dimension of the DGP emulator
#'     if `object` is an instance of the `dgp` class, and so on for subsequent columns and dimensions. If `object` is an instance of the `bundle` class,
#'     the first column must correspond to the first emulator in the bundle, and so on for subsequent columns and emulators.
#' @references
#' MacKay, D. J. (1992). Information-based objective functions for active data selection. *Neural Computation*, **4(4)**, 590-604.
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
#' # specify the input range
#' lim <- c(0,1)
#'
#' # locate the next design point using ALM
#' X_new <- alm(m, limits = lim)
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
alm <- function(object, ...){
  UseMethod("alm")
}

#' @rdname alm
#' @method alm gp
#' @export
alm.gp <- function(object, x_cand = NULL, n_start = 5, batch_size = 1, M = 50, limits = NULL, int = FALSE, n_screen = 100 * ncol(object$data$X), opt_control = list(), ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  if ( reticulate::py_is_null_xptr(object$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
  training_input <- object$data$X
  training_output <- object$data$Y
  n_dim_X <- ncol(training_input)
  #check x_cand
  if ( is.null(x_cand) ){
    is_cand <- FALSE
    limits <- check_limits(limits, n_dim_X)
    int <- check_int(int, n_dim_X)
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
  M <- as.integer(M)
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)
  #locate
  if ( batch_size==1 ){
    if (!is_cand){
      idx <- alm_search(object$emulator_obj, n_screen, n_start, limits,
                        int, M, NULL, list(), opt_control, is_gp = TRUE)
    } else {
      res = object$emulator_obj$metric(x_cand = x_cand, method = 'ALM', m = M)

      idx <- res[[1]]+1
    }
  } else {
    if (!is_cand){
      idx <- c()
      constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
      for (i in 1:batch_size){
        X_new <- alm_search(constructor_obj_cp, n_screen, n_start, limits,
                            int, M, NULL, list(), opt_control, is_gp = TRUE)
        Y_new <- constructor_obj_cp$predict(matrix(X_new, nrow = 1), m = M)[[1]]
        training_input <- rbind(training_input, X_new)
        training_output <- rbind(training_output, Y_new)
        constructor_obj_cp$update_xy(training_input, training_output)
        idx <- rbind(idx, X_new)
      }
    } else {
      idx <- c()
      idx_x_cand0 <- c(1:nrow(x_cand))
      idx_x_cand <- idx_x_cand0
      constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
      for (i in 1:batch_size){

        res = constructor_obj_cp$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'ALM', m = M)

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
  }
  if ( batch_size!=1 ){
    pkg.env$py_gc$collect()
    gc(full=T)
  }
  return(unname(idx))
}


#' @rdname alm
#' @method alm dgp
#' @export
alm.dgp <- function(object, x_cand = NULL, n_start = 5, batch_size = 1, M = 50, limits = NULL, int = FALSE, aggregate = NULL, n_screen = 100 * ncol(object$data$X), opt_control = list(), ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  if ( reticulate::py_is_null_xptr(object$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulator or, if it was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
  #if ( object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$type == 'likelihood' ){
  #  stop("The function is only applicable to DGP emulators without likelihood layers.", call. = FALSE)
  #}
  training_input <- object$data$X
  training_output <- object$data$Y
  n_dim_X <- ncol(training_input)
  #check x_cand
  if ( is.null(x_cand) ){
    is_cand <- FALSE
    limits <- check_limits(limits, n_dim_X)
    int <- check_int(int, n_dim_X)
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
  M <- as.integer(M)
  extra_param <- list()
  #check aggregate
  if ( !is.null(aggregate) ){
    add_arg <- list(...)
    gnames <- methods::formalArgs(aggregate)
    extra_param <- if ("..." %in% gnames) add_arg else
      add_arg[intersect(names(add_arg), gnames)]
  }
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)
  #locate
  if ( batch_size==1 ){
    if (!is_cand){
      idx <- alm_search(object$emulator_obj, n_screen, n_start, limits,
                        int, M, aggregate, extra_param, opt_control)
    } else {

      res = object$emulator_obj$metric(x_cand = x_cand, method = 'ALM', m = M, score_only = TRUE)

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
    }
  } else {
    if (!is_cand){
      idx <- c()
      constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
      emulator_obj_cp <- pkg.env$copy$deepcopy(object$emulator_obj)
      B <- as.integer(length(emulator_obj_cp$all_layer_set))
      burnin <- constructor_obj_cp$burnin
      isblock <- constructor_obj_cp$block
      for (i in 1:batch_size){
        X_new <- alm_search(emulator_obj_cp, n_screen, n_start, limits,
                            int, M, aggregate, extra_param, opt_control)
        Y_new <- alm_fantasy_response(emulator_obj_cp, constructor_obj_cp, X_new, M)
        training_input <- rbind(training_input, X_new)
        training_output <- rbind(training_output, Y_new)
        constructor_obj_cp$update_xy(training_input, training_output)
        est_obj <- constructor_obj_cp$estimate(burnin)
        emulator_obj_cp <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
        idx <- rbind(idx, X_new)
      }
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

        res = emulator_obj_cp$metric(x_cand = x_cand[idx_x_cand,,drop=F], method = 'ALM', m = M, score_only = TRUE)

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
        Y_new <- alm_fantasy_response(emulator_obj_cp, constructor_obj_cp, X_new, M)

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
  }
  if ( batch_size!=1 ){
    pkg.env$py_gc$collect()
    gc(full=T)
  }
  return(unname(idx))
}


#' @rdname alm
#' @method alm bundle
#' @export
alm.bundle <- function(object, x_cand = NULL, n_start = 5, batch_size = 1, M = 50, limits = NULL, int = FALSE, aggregate = NULL, n_screen = 100 * ncol(object$data$X[[1L]]), opt_control = list(), ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"bundle") ) stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)
  if ( reticulate::py_is_null_xptr(object$emulator1$constructor_obj) ) stop("The Python session originally associated with 'object' is no longer active. Please rebuild the emulators in the bundle or, if the bundle was saved using dgpsi::write(), load it into the R session with dgpsi::read().", call. = FALSE)
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
  } else {
    is_cand <- TRUE
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
  }
  M <- as.integer(M)
  extra_param <- list()
  #check aggregate
  if ( !is.null(aggregate) ){
    add_arg <- list(...)
    gnames <- methods::formalArgs(aggregate)
    extra_param <- if ("..." %in% gnames) add_arg else
      add_arg[intersect(names(add_arg), gnames)]
    if (is_cand){
      if (islist){
        x_cand_dfs <- lapply(x_cand, as.data.frame)
        x_cand <- as.matrix(Reduce(function(x, y) merge(x, y, all = FALSE), x_cand_dfs))
        if (length(x_cand)==0) stop("When using 'aggregate,' matrices in 'x_cand' must share at least some common rows.", call. = FALSE)
      }
    }
  }
  #check batch size
  batch_size <- as.integer(batch_size)
  if ( batch_size < 1 ) stop("'batch_size' must be >= 1.", call. = FALSE)

  is_gp <- logical(n_emulators)
  for (i in seq_len(n_emulators)) {
    is_gp[i] <- inherits(object[[paste0('emulator', i)]], "gp")
  }
  #locate
  if ( batch_size==1 ){
    if (!is_cand){
      emulator_obj_list <- vector('list', n_emulators)
      for (i in seq_len(n_emulators)) {
        emulator_obj_list[[i]] <- object[[paste0('emulator', i)]]$emulator_obj
      }
      idx <- alm_bundle_search(emulator_obj_list, is_gp, n_screen, n_start,
                               limits, int, M, aggregate, extra_param, opt_control)
    } else {
      scores <- vector('list', n_emulators)
      for ( i in 1:n_emulators ){
        obj_i <- object[[paste('emulator',i,sep='')]]
        if ( inherits(obj_i,"gp") ){
          res = obj_i$emulator_obj$metric(x_cand = if (is.list(x_cand)) {x_cand[[i]]} else {x_cand}, method = 'ALM', m = M, score_only = TRUE)
          scores[[i]] <- res
        } else {

          res = obj_i$emulator_obj$metric(x_cand = if (is.list(x_cand)) {x_cand[[i]]} else {x_cand}, method = 'ALM', m = M, score_only = TRUE)

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
    }
  } else {
    if (!is_cand){
      idx <- vector('list', n_emulators)
      constructor_obj_list <- list()
      emulator_obj_list <- list()
      for ( i in 1:n_emulators){
        obj_i <- object[[paste('emulator',i,sep='')]]
        constructor_obj_list[[i]] <- pkg.env$copy$deepcopy(obj_i$constructor_obj)
        emulator_obj_list[[i]] <- pkg.env$copy$deepcopy(obj_i$emulator_obj)
      }
      for (i in 1:batch_size){
        X_new <- alm_bundle_search(emulator_obj_list, is_gp, n_screen, n_start,
                                   limits, int, M, aggregate, extra_param, opt_control)

        for ( j in 1:n_emulators ){
          obj_j <- object[[paste('emulator',j,sep='')]]
          Y_new_j <- alm_fantasy_response(emulator_obj_list[[j]], constructor_obj_list[[j]],
                                          X_new[[j]], M, is_gp = is_gp[j])
          training_input[[j]] <- rbind(training_input[[j]], X_new[[j]])
          training_output[[j]] <- rbind(training_output[[j]], Y_new_j)
          if ( inherits(obj_j,"gp") ){
            emulator_obj_list[[j]]$update_xy(training_input[[j]], training_output[[j]])
          } else {
            B <- as.integer(length(obj_j$emulator_obj$all_layer_set))
            burnin <- obj_j$constructor_obj$burnin
            isblock <- obj_j$constructor_obj$block
            constructor_obj_list[[j]]$update_xy(training_input[[j]], training_output[[j]])
            est_obj <- constructor_obj_list[[j]]$estimate(burnin)
            emulator_obj_list[[j]] <- pkg.env$dgpsi$emulator(all_layer = est_obj, N = B, block = isblock)
          }
          idx[[j]] <- rbind(idx[[j]], X_new[[j]])
        }
      }
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

          res = emulator_obj_list[[j]]$metric(x_cand = if (is.list(x_cand)) {x_cand[[j]][idx_x_cand[[j]],,drop=F]} else {x_cand[idx_x_cand[[j]],,drop=F]}, method = 'ALM', m = M, score_only = TRUE)

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
        Y_new_j <- alm_fantasy_response(emulator_obj_list[[j]], constructor_obj_list[[j]],
                                          X_new[[j]], M, is_gp = is_gp[j])
        training_input[[j]] <- rbind(training_input[[j]], X_new[[j]])
        training_output[[j]] <- rbind(training_output[[j]], Y_new_j)
        if ( inherits(obj_j,"gp") ){
          emulator_obj_list[[j]]$update_xy(training_input[[j]], training_output[[j]])
        } else {
          B <- as.integer(length(obj_j$emulator_obj$all_layer_set))
          burnin <- obj_j$constructor_obj$burnin
          isblock <- obj_j$constructor_obj$block
          constructor_obj_list[[j]]$update_xy(training_input[[j]], training_output[[j]])
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
  }
  if ( batch_size!=1 ){
    pkg.env$py_gc$collect()
    gc(full=T)
  }
  return(unname(idx))
}

alm_bundle_search <- function(emulator_objs, is_gp, n_screen, n_start, limits, int, M, aggregate, extra_param, opt_control) {
  if (!is.null(aggregate)) {
    point <- alm_search(emulator_objs, n_screen, n_start, limits, int, M,
                        aggregate, extra_param, opt_control, bundle = TRUE)
    return(rep(list(point), length(emulator_objs)))
  }
  points <- vector('list', length(emulator_objs))
  for (i in seq_along(emulator_objs)) {
    points[[i]] <- alm_search(emulator_objs[[i]], n_screen, n_start, limits,
                             int, M, NULL, list(), opt_control,
                             is_gp = is_gp[i], mean_output = TRUE)
  }
  points
}

alm_bundle_scores <- function(emulator_objs, points, M) {
  scores <- matrix(NA_real_, nrow = nrow(points), ncol = length(emulator_objs))
  for (i in seq_along(emulator_objs)) {
    scores[, i] <- rowMeans(emulator_objs[[i]]$metric(
      x_cand = points, method = 'ALM', m = M, score_only = TRUE
    ))
  }
  scores
}

alm_opt_control <- function(opt_control, n_free) {
  allowed <- c("maxit", "ftol", "gtol", "step")
  if (!is.list(opt_control) || (length(opt_control) &&
      (is.null(names(opt_control)) || anyNA(names(opt_control)) ||
       any(!nzchar(names(opt_control))) || anyDuplicated(names(opt_control)) ||
       any(!names(opt_control) %in% allowed)))) {
    stop("'opt_control' must be a named list containing only maxit, ftol, gtol, and step.", call. = FALSE)
  }
  values <- list(maxit = 100L, ftol = 1e7 * .Machine$double.eps, gtol = 0, step = NULL)
  values[names(opt_control)] <- opt_control
  if (length(values$maxit) != 1L || !is.numeric(values$maxit) ||
      !is.finite(values$maxit) || values$maxit < 1 ||
      values$maxit != floor(values$maxit) || values$maxit > .Machine$integer.max) {
    stop("'opt_control$maxit' must be a positive integer within R's integer range.", call. = FALSE)
  }
  for (name in c("ftol", "gtol")) {
    value <- values[[name]]
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value < 0) {
      stop(sprintf("'opt_control$%s' must be one finite nonnegative number.", name), call. = FALSE)
    }
  }
  if (!is.null(values$step) && (length(values$step) != 1L ||
      !is.numeric(values$step) || !is.finite(values$step) || values$step <= 0)) {
    stop("'opt_control$step' must be NULL or one finite positive number.", call. = FALSE)
  }
  r <- list(maxit = as.integer(values$maxit), factr = values$ftol / .Machine$double.eps,
            pgtol = values$gtol)
  python <- list(maxiter = as.integer(values$maxit), ftol = values$ftol, gtol = values$gtol)
  if (!is.null(values$step)) {
    r$ndeps <- rep(values$step, n_free)
    python$eps <- values$step
  }
  list(r = r, python = python)
}

alm_search <- function(emulator_obj, n_screen, n_start, limits, int, M, aggregate, extra_param, opt_control = list(), is_gp = FALSE, mean_output = FALSE, bundle = FALSE) {
  if (length(n_start) != 1L || !is.numeric(n_start) ||
      !is.finite(n_start) || n_start < 1 || n_start != floor(n_start)) {
    stop("'n_start' must be a positive integer.", call. = FALSE)
  }
  if (length(n_screen) != 1L || !is.numeric(n_screen) ||
      !is.finite(n_screen) || n_screen < n_start ||
      n_screen != floor(n_screen)) {
    stop("'n_screen' must be an integer >= 'n_start'.", call. = FALSE)
  }
  if (is.numeric(limits) && !is.complex(limits)) storage.mode(limits) <- "double"
  if (!is.matrix(limits) || ncol(limits) != 2L || !is.numeric(limits) ||
      is.complex(limits) || any(!is.finite(limits)) ||
      any(limits[, 1] > limits[, 2]) || any(!is.finite(limits[, 2] - limits[, 1]))) {
    stop("'limits' must have finite numeric, ordered bounds and finite widths.", call. = FALSE)
  }
  if (!is.logical(int) || length(int) != nrow(limits) || anyNA(int)) {
    stop("'int' must be a logical mask without missing values.", call. = FALSE)
  }
  lower <- limits[, 1]
  upper <- limits[, 2]
  lower[int] <- ceiling(lower[int])
  upper[int] <- floor(upper[int])
  if (any(lower[int] > upper[int])) {
    stop("No feasible integer ALM screening points.", call. = FALSE)
  }
  span <- upper - lower
  free <- span > 0 & !int
  solver_control <- alm_opt_control(opt_control, sum(free))
  sampling_limits <- limits
  sampling_limits[int, 1] <- 0
  sampling_limits[int, 2] <- upper[int] - lower[int] + 1
  x_screen <- reverse_minmax(lhs::maximinLHS(n_screen, length(lower)), sampling_limits)
  if (any(int)) {
    x_screen[, int] <- sweep(floor(x_screen[, int, drop = FALSE]), 2L, lower[int], "+")
  }
  x_screen <- unique(x_screen[
    rowSums(sweep(x_screen, 2L, lower, "<") |
            sweep(x_screen, 2L, upper, ">")) == 0, , drop = FALSE
  ])
  if (!nrow(x_screen)) stop("No feasible ALM screening points.", call. = FALSE)
  scores <- if (bundle) alm_bundle_scores(emulator_obj, x_screen, M) else
    emulator_obj$metric(x_cand = x_screen, method = 'ALM', m = M, score_only = TRUE)
  if (!is.null(aggregate)) {
    scores <- matrix(alm_aggregate_scores(scores, aggregate, extra_param), ncol = 1)
  } else if (mean_output) {
    scores <- matrix(rowMeans(scores), ncol = 1)
  }

  span[span == 0] <- 1
  x_scaled <- sweep(sweep(x_screen, 2L, lower, "-"), 2L, span, "/")
  norm2 <- rowSums(x_scaled * x_scaled)

  points <- matrix(NA_real_, nrow = ncol(scores), ncol = ncol(x_screen))
  for (j in seq_len(ncol(scores))) {
    starts <- alm_select_starts(x_scaled, scores[, j], n_start, norm2)
    points[j, ] <- x_screen[starts[1], ]
    if (!any(free)) next

    if (is.null(aggregate)) {
      res <- tryCatch(pkg.env$dgpsi$utils$multistart(
        func = emulator_obj$metric,
        initials = x_screen[starts, , drop = FALSE],
        lb = reticulate::np_array(lower),
        up = reticulate::np_array(upper),
        args = if (is_gp) reticulate::tuple('ALM', 1., M, TRUE) else
          reticulate::tuple('ALM', NULL, 1., M, TRUE),
        method = "L-BFGS-B",
        out_dim = if (mean_output) -1L else as.integer(j - 1L),
        int_mask = int,
        options = reticulate::r_to_py(solver_control$python)
      ), error = identity)
      if (inherits(res, "error")) {
        next
      }
      refined <- matrix(res, nrow = 1)
    } else {
      refined <- matrix(NA_real_, nrow = length(starts), ncol = ncol(x_screen))
      for (k in seq_along(starts)) {
        res <- tryCatch(stats::optim(
          par = x_scaled[starts[k], free], fn = alm_objective,
          method = 'L-BFGS-B', lower = rep(0, sum(free)), upper = rep(1, sum(free)),
          control = solver_control$r,
          emulator_obj = emulator_obj, M = M, aggregate = aggregate,
          extra_param = extra_param, point = x_screen[starts[k], ], free = free,
          lower_free = lower[free], upper_free = upper[free], span_free = span[free],
          bundle = bundle
        ), error = identity)
        if (inherits(res, "error") || any(!is.finite(res$par))) {
          next
        }
        refined[k, ] <- x_screen[starts[k], ]
        refined[k, free] <- pmin(pmax(lower[free] + res$par * span[free], lower[free]), upper[free])
      }
    }

    if (any(int)) refined[, int] <- round(refined[, int, drop = FALSE])
    for (k in seq_along(lower)) {
      refined[, k] <- pmin(pmax(refined[, k], lower[k]), upper[k])
    }
    valid <- rowSums(!is.finite(refined)) == 0 &
      rowSums(refined[, int, drop = FALSE] !=
              round(refined[, int, drop = FALSE])) == 0
    refined <- refined[valid, , drop = FALSE]
    if (!nrow(refined)) next
    refined_scores <- tryCatch(if (bundle) alm_bundle_scores(emulator_obj, refined, M) else
      emulator_obj$metric(x_cand = refined, method = 'ALM', m = M, score_only = TRUE), error = identity)
    if (inherits(refined_scores, "error")) {
      next
    }
    values <- tryCatch(if (!is.null(aggregate))
      alm_aggregate_scores(refined_scores, aggregate, extra_param) else if (mean_output)
      rowMeans(refined_scores) else refined_scores[, j], error = identity)
    if (inherits(values, "error")) {
      next
    }
    values[!is.finite(values)] <- -Inf
    candidates <- rbind(x_screen[starts[1], , drop = FALSE], refined)
    points[j, ] <- candidates[which.max(c(scores[starts[1], j], values)), ]
  }
  points
}

alm_fantasy_response <- function(emulator_obj, constructor_obj, points, M, is_gp = FALSE) {
  response <- emulator_obj$predict(points, m = M)[[1]]
  if (is_gp) return(response)
  last_node <- constructor_obj$all_layer[[constructor_obj$n_layer]][[1]]
  if (last_node$name == "Categorical") {
    probabilities <- as.matrix(response)
    index <- if (last_node$num_classes == 2L) as.integer(probabilities[, 1] > 0.5) + 1L else
      max.col(probabilities, ties.method = "first")
    response <- matrix(last_node$class_encoder$classes_[index], ncol = 1L)
  }
  response
}

alm_aggregate_scores <- function(scores, aggregate, extra_param) {
  values <- numeric(nrow(scores))
  for (i in seq_len(nrow(scores))) {
    value <- do.call(aggregate, c(list(scores[i, , drop = FALSE]), extra_param))
    if (!is.numeric(value) || length(value) != 1L) {
      stop("'aggregate' must return one numeric score for a candidate.", call. = FALSE)
    }
    values[i] <- value
  }
  values
}

alm_objective <- function(x, emulator_obj, M, aggregate, extra_param, point, free, lower_free, upper_free, span_free, bundle = FALSE) {
  point[free] <- pmin(pmax(lower_free + x * span_free, lower_free), upper_free)
  scores <- if (bundle) alm_bundle_scores(emulator_obj, matrix(point, nrow = 1), M) else
    emulator_obj$metric(x_cand = matrix(point, nrow = 1), method = 'ALM', m = M, score_only = TRUE)
  -alm_aggregate_scores(scores, aggregate, extra_param)
}

alm_select_starts <- function(x, scores, n_start, norm2) {
  eligible <- which(is.finite(scores))
  if (!length(eligible)) stop("No finite ALM screening scores.", call. = FALSE)
  n_start <- min(n_start, length(eligible))
  y <- scores[eligible]
  first <- which.max(y)
  if (n_start == 1L) return(eligible[first])
  if (n_start == length(eligible)) {
    return(eligible[c(first, seq_along(eligible)[-first])])
  }

  # x is bound-normalized; squared row norms are shared across outputs.
  if (length(eligible) != nrow(x)) {
    norm2 <- norm2[eligible]
    x <- x[eligible, , drop = FALSE]
  }
  # Rescale scores and compare weight * squared distance in log space.
  magnitude <- max(abs(y))
  if (magnitude > 0) y <- y / magnitude
  score_sd <- stats::sd(y)
  log_weight <- if (score_sd > 0) (y - max(y)) / score_sd else numeric(length(y))

  starts <- integer(n_start)
  starts[1L] <- first
  chosen <- logical(length(y))
  chosen[first] <- TRUE
  nearest2 <- rep(Inf, length(y))
  for (k in seq_len(n_start - 1L)) {
    last <- starts[k]
    distance2 <- pmax(norm2 + norm2[last] - 2 * drop(x %*% x[last, ]), 0)
    nearest2 <- pmin(nearest2, distance2)
    remaining <- which(!chosen)
    next_id <- remaining[which.max(log_weight[remaining] + log(nearest2[remaining]))]
    starts[k + 1L] <- next_id
    chosen[next_id] <- TRUE
  }
  eligible[starts]
}
