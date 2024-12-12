#' @title Locate the next design point(s) for a (D)GP emulator or a bundle of (D)GP emulators using Active Learning MacKay (ALM)
#'
#' @description
#'
#' `r new_badge("updated")`
#'
#' This function searches from a candidate set to locate the next design point(s) to be added to a (D)GP emulator
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
#' @param n_start `r new_badge("new")` an integer that gives the number of initial design points to be used to determine next design point(s). This argument
#'     is only used when `x_cand` is `NULL`. Defaults to `20`.
#' @param batch_size an integer that gives the number of design points to be chosen. Defaults to `1`.
#' @param M `r new_badge("new")` the size of the conditioning set for the Vecchia approximation in the criterion calculation. This argument is only used if the emulator `object`
#'     was constructed under the Vecchia approximation. Defaults to `50`.
#' @param workers the number of processes to be used for design point selection. If set to `NULL`,
#'     the number of processes is set to `max physical cores available %/% 2`. Defaults to `1`. The argument does not currently support Windows machines when the `aggregate`
#'     function is provided, due to the significant overhead caused by initializing the Python environment for each worker under spawning.
#' @param limits a two-column matrix that gives the ranges of each input dimension, or a vector of length two if there is only one input dimension.
#'     If a vector is provided, it will be converted to a two-column row matrix. The rows of the matrix correspond to input dimensions, and its
#'     first and second columns correspond to the minimum and maximum values of the input dimensions. This
#'     argument is only used when `x_cand = NULL`. Defaults to `NULL`.
#' @param int `r new_badge("new")` a bool or a vector of bools that indicates if an input dimension is an integer type. If a single bool is given, it will be applied to
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
alm.gp <- function(object, x_cand = NULL, n_start = 20, batch_size = 1, M = 50, workers = 1, limits = NULL, int = FALSE, ...) {
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
    lower <- limits[, 1]
    upper <- limits[, 2]
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
    if (!is_cand){
      #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input, limits), n_start), n_start), limits)
      x_cand <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
      res <- pkg.env$dgpsi$utils$multistart(
        func = object$emulator_obj$metric,
        initials = x_cand,
        lb = reticulate::np_array(lower),
        up = reticulate::np_array(upper),
        args = reticulate::tuple('ALM', 1., M, TRUE),
        method = "L-BFGS-B",
        core_num = workers,
        int_mask = int
      )
      idx <- matrix(res, nrow = 1)
    } else {
      if ( identical(workers,as.integer(1)) ){
        res = object$emulator_obj$metric(x_cand = x_cand, method = 'ALM', m = M)
      } else {
        res = object$emulator_obj$pmetric(x_cand = x_cand, method = 'ALM', m = M, core_num = workers)
      }
      idx <- res[[1]]+1
    }
  } else {
    if (!is_cand){
      idx <- c()
      constructor_obj_cp <- pkg.env$copy$deepcopy(object$constructor_obj)
      for (i in 1:batch_size){
        #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input, limits), n_start), n_start), limits)
        x_cand <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
        res <- pkg.env$dgpsi$utils$multistart(
          func = constructor_obj_cp$metric,
          initials = x_cand,
          lb = reticulate::np_array(lower),
          up = reticulate::np_array(upper),
          args = reticulate::tuple('ALM', 1., M, TRUE),
          method = "L-BFGS-B",
          core_num = workers,
          int_mask = int
        )
        X_new <- res
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
alm.dgp <- function(object, x_cand = NULL, n_start = 20, batch_size = 1, M = 50, workers = 1, limits = NULL, int = FALSE, aggregate = NULL, ...) {
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
    lower <- limits[, 1]
    upper <- limits[, 2]
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
    if (!is_cand){
      #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input, limits), n_start), n_start), limits)
      x_cand <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
      if ( is.null(aggregate) ){
        if ( object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$type == 'likelihood' ){
          lik_name <- object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$name
          if (lik_name=='Poisson'){
            n_output_dim <- 1
          } else if (lik_name=='Hetero' | lik_name=='NegBin'){
            n_output_dim <- 2
          } else if (lik_name=='Categorical') {
            num_class <- object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$num_classes
            if (num_class == 2){
              n_output_dim <- 1
            } else {
              n_output_dim <- num_class
            }
          }
        } else {
          n_output_dim <- n_dim_Y
        }
        idx <- c()
        for (i in 1:n_output_dim){
          res <- pkg.env$dgpsi$utils$multistart(
            func = object$emulator_obj$metric,
            initials = x_cand,
            lb = reticulate::np_array(lower),
            up = reticulate::np_array(upper),
            args = reticulate::tuple('ALM', NULL, 1., M, TRUE),
            method = "L-BFGS-B",
            core_num = workers,
            out_dim = as.integer(i-1),
            int_mask = int
          )
          res <- pmin(pmax(res, lower), upper)
          idx <- rbind(idx, res)
        }
        #idx <- pkg.env$np$unique(idx, axis=0L)
      } else {
        if (Sys.info()["sysname"] == "Windows") workers <- 1
        total_cores <- parallel::detectCores(logical = FALSE)
        if (is.null(workers)) {
          workers <- max(1, floor(total_cores / 2))  # Default to half of the physical cores
        }
        num_thread <- floor(total_cores / workers)
        set_thread_num(num_thread)
        if ("..." %in% gnames){
          extra_param <- add_arg
        } else {
          gidx <- gnames %in% names(add_arg)
          extra_param <- add_arg[gnames[gidx]]
        }
        if (workers == 1) {
          fn_R <- function(x, m, extra_param, int){
            x[int] <- round(x[int])
            scores <- object$emulator_obj$metric(x_cand = matrix(x, nrow = 1), method = 'ALM', m = m, score_only = TRUE)
            agg_res <- do.call(aggregate, c(list(scores), extra_param))
            return(-agg_res)
          }
          results <- lapply(1:n_start, function(x0) {
            stats::optim(
              par = x_cand[x0,],
              fn = fn_R,
              method = 'L-BFGS-B',
              lower = lower,
              upper = upper,
              control = list(maxit = 100),
              extra_param = extra_param,
              m = M,
              int = int
            )
          })
        } else {
          fn_R <- function(x, m, extra_param, int){
            x[int] <- round(x[int])
            scores <- object$emulator_obj$metric(x_cand = matrix(x, nrow = 1), method = 'ALM', m = m, score_only = TRUE)
            agg_res <- do.call(aggregate, c(list(scores), extra_param))
            return(-agg_res)
          }
          results <- parallel::mclapply(1:n_start, function(x0) {
            stats::optim(
              par = x_cand[x0,],
              fn = fn_R,
              method = 'L-BFGS-B',
              lower = lower,
              upper = upper,
              control = list(maxit = 100),
              extra_param = extra_param,
              m = M,
              int = int
            )
          }, mc.cores = workers)
        }
        values <- sapply(results, function(res) res$value)
        parameters <- lapply(results, function(res) res$par)

        # Find the best result
        best_idx <- which.min(values)
        idx <- parameters[[best_idx]]
        idx[int] <- round(idx[int])
        idx <- pmin(pmax(idx, lower), upper)
        idx <- matrix(idx, nrow = 1)
      }
    } else {
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
      if ( !is.null(aggregate) ){
        if (Sys.info()["sysname"] == "Windows") workers <- 1
        total_cores <- parallel::detectCores(logical = FALSE)
        if (is.null(workers)) {
          workers <- max(1, floor(total_cores / 2))  # Default to half of the physical cores
        }
        num_thread <- floor(total_cores / workers)
        set_thread_num(num_thread)
        if ("..." %in% gnames){
          extra_param <- add_arg
        } else {
          gidx <- gnames %in% names(add_arg)
          extra_param <- add_arg[gnames[gidx]]
        }
      } else {
        if ( object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$type == 'likelihood' ){
          lik_name <- object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$name
          if (lik_name=='Poisson'){
            n_output_dim <- 1
          } else if (lik_name=='Hetero' | lik_name=='NegBin'){
            n_output_dim <- 2
          } else if (lik_name=='Categorical') {
            num_class <- object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$num_classes
            if (num_class == 2){
              n_output_dim <- 1
            } else {
              n_output_dim <- num_class
            }
          }
        } else {
          n_output_dim <- n_dim_Y
        }
      }
      for (i in 1:batch_size){
        #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input, limits), n_start), n_start), limits)
        x_cand <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
        X_new <- c()
        if ( is.null(aggregate) ){
          for (j in 1:n_output_dim){
            res <- pkg.env$dgpsi$utils$multistart(
              func = object$emulator_obj$metric,
              initials = x_cand,
              lb = reticulate::np_array(lower),
              up = reticulate::np_array(upper),
              args = reticulate::tuple('ALM', NULL, 1., M, TRUE),
              method = "L-BFGS-B",
              core_num = workers,
              out_dim = as.integer(j-1),
              int_mask = int
            )
            res <- pmin(pmax(res, lower), upper)
            X_new <- rbind(X_new, res)
          }
          #X_new <- pkg.env$np$unique(X_new, axis=0L)
        } else {
          fn_R <- function(x, m, extra_param, int){
            x[int] <- round(x[int])
            scores <- object$emulator_obj$metric(x_cand = matrix(x, nrow = 1), method = 'ALM', m = m, score_only = TRUE)
            agg_res <- do.call(aggregate, c(list(scores), extra_param))
            return(-agg_res)
          }
          if (workers == 1) {
            results <- lapply(1:n_start, function(x0) {
              stats::optim(
                par = x_cand[x0,],
                fn = fn_R,
                method = 'L-BFGS-B',
                lower = lower,
                upper = upper,
                control = list(maxit = 100),
                extra_param = extra_param,
                m = M,
                int = int
              )
            })
          } else {
            results <- parallel::mclapply(1:n_start, function(x0) {
              stats::optim(
                par = x_cand[x0,],
                fn = fn_R,
                method = 'L-BFGS-B',
                lower = lower,
                upper = upper,
                control = list(maxit = 100),
                extra_param = extra_param,
                m = M,
                int = int
              )
            }, mc.cores = workers)
          }
          values <- sapply(results, function(res) res$value)
          parameters <- lapply(results, function(res) res$par)

          # Find the best result
          best_idx <- which.min(values)
          res <- parameters[[best_idx]]
          res[int] <- round(res[int])
          res <- pmin(pmax(res, lower), upper)
          X_new <- rbind(X_new, matrix(res,nrow=1))
        }
        Y_new <- emulator_obj_cp$predict(X_new, m = M)[[1]]
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
alm.bundle <- function(object, x_cand = NULL, n_start = 20, batch_size = 1, M = 50, workers = 1, limits = NULL, int = FALSE, aggregate = NULL, ...) {
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
    lower <- limits[, 1]
    upper <- limits[, 2]
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

  #locate
  if ( batch_size==1 ){
    if (!is_cand){
      if ( is.null(aggregate) ){
        x_cand <- vector('list', n_emulators)
        idx <- vector('list', n_emulators)
        for (i in 1:n_emulators){
          #x_cand[[i]] <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input[[i]], limits), n_start), n_start), limits)
          x_cand[[i]] <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
          obj_i <- object[[paste('emulator',i,sep='')]]
          res <- pkg.env$dgpsi$utils$multistart(
            func = obj_i$emulator_obj$metric,
            initials = x_cand[[i]],
            lb = reticulate::np_array(lower),
            up = reticulate::np_array(upper),
            args = reticulate::tuple('ALM', NULL, 1., M, TRUE),
            method = "L-BFGS-B",
            core_num = workers,
            out_dim = -1L,
            int_mask = int
          )
          res <- pmin(pmax(res, lower), upper)
          idx[[i]] <- matrix(res,nrow=1)
        }
        #idx <- pkg.env$np$unique(idx, axis=0L)
      } else {
        idx <- vector('list', n_emulators)
        if (Sys.info()["sysname"] == "Windows") workers <- 1
        total_cores <- parallel::detectCores(logical = FALSE)
        if (is.null(workers)) {
          workers <- max(1, floor(total_cores / 2))  # Default to half of the physical cores
        }
        num_thread <- floor(total_cores / workers)
        set_thread_num(num_thread)
        if ("..." %in% gnames){
          extra_param <- add_arg
        } else {
          gidx <- gnames %in% names(add_arg)
          extra_param <- add_arg[gnames[gidx]]
        }
        #total_input <- do.call(rbind, training_input)
        #total_input <- pkg.env$np$unique(total_input, axis=0L)
        #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(total_input, limits), n_start), n_start), limits)
        x_cand <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
        if (workers == 1) {
          fn_R <- function(x, m, extra_param, int){
            scores <- vector('list', n_emulators)
            x[int] <- round(x[int])
            for (i in 1:n_emulators){
              obj_i <- object[[paste('emulator',i,sep='')]]
              scores_i <- obj_i$emulator_obj$metric(x_cand = matrix(x, nrow = 1), method = 'ALM', m = m, score_only = TRUE)
              scores[[i]] <- if(ncol(scores_i) == 1) as.vector(scores_i) else rowMeans(scores_i)
            }
            scores <- do.call(cbind, scores)
            agg_res <- do.call(aggregate, c(list(matrix(scores, nrow = 1)), extra_param))
            return(-agg_res)
          }
          results <- lapply(1:n_start, function(x0) {
            stats::optim(
              par = x_cand[x0,],
              fn = fn_R,
              method = 'L-BFGS-B',
              lower = lower,
              upper = upper,
              control = list(maxit = 100),
              extra_param = extra_param,
              m = M,
              int = int
            )
          })
        } else {
          fn_R <- function(x, m, extra_param, int){
            scores <- vector('list', n_emulators)
            x[int] <- round(x[int])
            for (i in 1:n_emulators){
              obj_i <- object[[paste('emulator',i,sep='')]]
              scores_i <- obj_i$emulator_obj$metric(x_cand = matrix(x, nrow = 1), method = 'ALM', m = m, score_only = TRUE)
              scores[[i]] <- if(ncol(scores_i) == 1) as.vector(scores_i) else rowMeans(scores_i)
            }
            scores <- do.call(cbind, scores)
            agg_res <- do.call(aggregate, c(list(matrix(scores, nrow = 1)), extra_param))
            return(-agg_res)
          }
          results <- parallel::mclapply(1:n_start, function(x0) {
            stats::optim(
              par = x_cand[x0,],
              fn = fn_R,
              method = 'L-BFGS-B',
              lower = lower,
              upper = upper,
              control = list(maxit = 100),
              extra_param = extra_param,
              m = M,
              int = int
            )
          }, mc.cores = workers)
        }
        values <- sapply(results, function(res) res$value)
        parameters <- lapply(results, function(res) res$par)

        # Find the best result
        best_idx <- which.min(values)
        res <- parameters[[best_idx]]
        res[int] <- round(res[int])
        res <- pmin(pmax(res, lower), upper)
        for (j in 1:n_emulators){
          idx[[j]] <- matrix(res,nrow=1)
        }
      }
    } else {
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
      if ( !is.null(aggregate) ){
        if (Sys.info()["sysname"] == "Windows") workers <- 1
        total_cores <- parallel::detectCores(logical = FALSE)
        if (is.null(workers)) {
          workers <- max(1, floor(total_cores / 2))  # Default to half of the physical cores
        }
        num_thread <- floor(total_cores / workers)
        set_thread_num(num_thread)
        if ("..." %in% gnames){
          extra_param <- add_arg
        } else {
          gidx <- gnames %in% names(add_arg)
          extra_param <- add_arg[gnames[gidx]]
        }
      }

      for (i in 1:batch_size){
        X_new <- vector('list', n_emulators)
        if ( is.null(aggregate) ){
          x_cand <- vector('list', n_emulators)
          for (j in 1:n_emulators){
            #x_cand[[j]] <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(training_input[[j]], limits), n_start), n_start), limits)
            x_cand[[j]] <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
            res <- pkg.env$dgpsi$utils$multistart(
              func = emulator_obj_list[[j]]$metric,
              initials = x_cand[[j]],
              lb = reticulate::np_array(lower),
              up = reticulate::np_array(upper),
              args = reticulate::tuple('ALM', NULL, 1., M, TRUE),
              method = "L-BFGS-B",
              core_num = workers,
              out_dim = -1L,
              int_mask = int
            )
            res <- pmin(pmax(res, lower), upper)
            X_new[[j]] <- matrix(res,nrow=1)
          }
          #X_new <- pkg.env$np$unique(X_new, axis=0L)
        } else {
          #total_input <- do.call(rbind, training_input)
          #total_input <- pkg.env$np$unique(total_input, axis=0L)
          #x_cand <- reverse_minmax(utils::tail(lhs::augmentLHS(minmax(total_input, limits), n_start), n_start), limits)
          x_cand <- reverse_minmax(lhs::maximinLHS(n_start,n_dim_X), limits)
          fn_R <- function(x, m, extra_param, int){
            scores <- vector('list', n_emulators)
            x[int] <- round(x[int])
            for (k in 1:n_emulators){
              scores_k <- emulator_obj_list[[k]]$metric(x_cand = matrix(x, nrow = 1), method = 'ALM', m = m, score_only = TRUE)
              scores[[k]] <- if(ncol(scores_k) == 1) as.vector(scores_k) else rowMeans(scores_k)
            }
            scores <- do.call(cbind, scores)
            agg_res <- do.call(aggregate, c(list(matrix(scores, nrow = 1)), extra_param))
            return(-agg_res)
          }

          if (workers == 1) {
            results <- lapply(1:n_start, function(x0) {
              stats::optim(
                par = x_cand[x0,],
                fn = fn_R,
                method = 'L-BFGS-B',
                lower = lower,
                upper = upper,
                control = list(maxit = 100),
                extra_param = extra_param,
                m = M,
                int = int
              )
            })
          } else {
            results <- parallel::mclapply(1:n_start, function(x0) {
              stats::optim(
                par = x_cand[x0,],
                fn = fn_R,
                method = 'L-BFGS-B',
                lower = lower,
                upper = upper,
                control = list(maxit = 100),
                extra_param = extra_param,
                m = M,
                int = int
              )
            }, mc.cores = workers)
          }
          values <- sapply(results, function(res) res$value)
          parameters <- lapply(results, function(res) res$par)

          # Find the best result
          best_idx <- which.min(values)
          res <- parameters[[best_idx]]
          res[int] <- round(res[int])
          res <- pmin(pmax(res, lower), upper)
          for (j in 1:n_emulators){
            X_new[[j]] <- matrix(res,nrow=1)
          }
        }

        for ( j in 1:n_emulators ){
          obj_j <- object[[paste('emulator',j,sep='')]]
          Y_new_j <- emulator_obj_list[[j]]$predict(X_new[[j]], m = M)[[1]]
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
  }
  if ( batch_size!=1 ){
    pkg.env$py_gc$collect()
    gc(full=T)
  }
  return(unname(idx))
}
