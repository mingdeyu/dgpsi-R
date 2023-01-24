#' @title Sequential design of a (D)GP emulator or a bundle of (D)GP emulators
#'
#' @description This function implements the sequential design of a (D)GP emulator or a bundle of (D)GP emulators.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param N the number of steps for the sequential design.
#' @param x_cand a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
#'     in which the next design point is determined. If `x_cand = NULL`, the candidate set will be generated using `n_cand`,
#'     `limits`, and `int`. Defaults to `NULL`.
#' @param y_cand a matrix (with each row being a simulator evaluation and column being an output dimension) that gives the realizations
#'    from the simulator at input positions in `x_cand`. Defaults to `NULL`.
#' @param n_cand an integer that gives
#' * the size of the candidate set in which the next design point is determined, if `x_cand = NULL`;
#' * the size of a sub-set to be sampled from the candidate set `x_cand` at each step of the sequential design to determine the next
#'   design point, if `x_cand` is not `NULL`.
#'
#' Defaults to `200`.
#' @param limits a two-column matrix that gives the ranges of each input dimension, or a vector of length two if there is only one
#'     input dimension. If a vector is provided, it will be converted to a two-column row matrix. The rows of the matrix correspond to input
#'     dimensions, and its first and second columns correspond to the minimum and maximum values of the input dimensions. If
#'     `limits = NULL`, the ranges of input dimensions will be determined from the training data contained in `object`. This argument
#'     is used when `x_cand = NULL` and `y_cand = NULL`. Defaults to `NULL`.
#' @param int a bool or a vector of bools that indicates if an input dimension is an integer type. If a bool is given, it will be applied to
#'     all input dimensions. If a vector is provided, it should have a length equal to the input dimensions and will be applied to individual
#'     input dimensions. Defaults to `FALSE`.
#' @param f an R function that represents the simulator. `f` needs to be specified with the following basic rules:
#' * the first argument of the function should be a matrix with rows being different design points and columns being input dimensions.
#' * the output of the function can either
#'   - a matrix with rows being different outputs (corresponding to the input design points) and columns being output dimensions. If there is
#'     only one output dimension, the matrix still needs to be returned with a single column.
#'   - a list with the first element being the output matrix described above and, optionally, additional named elements which will update values
#'     of any arguments with the same names passed via `...`. The list output can be useful if some additional arguments of `f` and `aggregate`
#'     need to be updated after each step of the sequential design.
#' @param freq a vector of two integers with the first element giving the frequency (in number of steps) to re-fit the
#'     emulator, and the second element giving the frequency to implement the emulator validation (for RMSE). Defaults to `c(1, 1)`.
#' @param x_test a matrix (with each row being an input testing data point and each column being an input dimension) that gives the testing
#'     input data to evaluate the emulator after each step of the sequential design. Set to `NULL` for the LOO-based emulator validation.
#'     Defaults to `NULL`. This argument is only used if `eval = NULL`.
#' @param y_test the testing output data that correspond to `x_test` for the emulator validation after each step of the sequential design:
#' * if `object` is an instance of the `gp` class, `y_test` is a matrix with only one column and each row being an testing output data point.
#' * if `object` is an instance of the `dgp` class, `y_test` is a matrix with its rows being testing output data points and columns being
#'   output dimensions.
#'
#' Set to `NULL` for the LOO-based emulator validation. Defaults to `NULL`. This argument is only used if `eval = NULL`.
#' @param target a numeric or a vector that gives the target RMSEs at which the sequential design is terminated. Defaults to `NULL`, in which
#'     case the sequential design stops after `N` steps. See *Note* section below for further information about `target`.
#' @param method an R function that give indices of designs points in a candidate set. The function must satisfy the following basic rules:
#' * the first argument is an emulator object that can be either an instance of
#'   - the `gp` class (produced by [gp()]);
#'   - the `dgp` class (produced by [dgp()]);
#'   - the `bundle` class (produced by [pack()]).
#' * the second argument is a matrix with rows representing a set of different design points.
#' * the output of the function
#'   - is a vector of indices if the first argument is an instance of the `gp` class;
#'   - is a matrix of indices if the first argument is an instance of the `dgp` class. If there are different design points to be added with
#'     respect to different outputs of the DGP emulator, the column number of the matrix should equal to the number of the outputs. If design
#'     points are common to all outputs of the DGP emulator, the matrix should be single-columned. If more than one design points are determined
#'     for a given output or for all outputs, the indices of these design points are placed in the matrix with extra rows.
#'   - is a matrix of indices if the first argument is an instance of the `bundle` class. Each row of the matrix gives the indices of the design
#'     points to be added to individual emulators in the bundle.
#'
#' See [alm()], [mice()], and [pei()] for examples on customizing `method`. Defaults to [mice()].
#' @param eval an R function that calculates the customized evaluating metric of the emulator.The function must satisfy the following basic rules:
#' * the first argument is an emulator object that can be either an instance of
#'   - the `gp` class (produced by [gp()]);
#'   - the `dgp` class (produced by [dgp()]);
#'   - the `bundle` class (produced by [pack()]).
#' * the output of the function can be
#'   - a single metric value, if the first argument is an instance of the `gp` class;
#'   - a single metric value or a vector of metric values with the length equal to the number of output dimensions, if the first argument is an
#'     instance of the `dgp` class;
#'   - a single metric value metric or a vector of metric values with the length equal to the number of emulators in the bundle, if the first
#'     argument is an instance of the `bundle` class.
#'
#' If no customized function is provided, the built-in evaluation metric, RMSE, will be calculated. Defaults to `NULL`. See *Note* section below for further information.
#' @param verb a bool indicating if the trace information will be printed during the sequential design.
#'     Defaults to `TRUE`.
#' @param check_point a vector of integers that indicates at which steps the sequential design will pause and ask for the confirmation
#'     from the user if the sequential design should continue or be terminated. Set to `NULL` to suspend the manual intervention. Defaults
#'     to `NULL`.
#' @param cores an integers that gives the number of cores to be used for emulator validations. If set to `NULL`, the number of cores is
#'     set to `(max physical cores available - 1)`. Defaults to `1`. This argument is only used if `eval = NULL`.
#' @param train_N an integer or a vector of integers that gives the number of training iterations to be used to re-fit the DGP emulator at each step
#'     of the sequential design:
#' * If `train_N` is an integer, then at each step the DGP emulator will re-fitted (based on the frequency of re-fit specified in `freq`) with `train_N` iterations.
#' * If `train_N` is a vector, then its size must be `N` even the re-fit frequency specified in `freq` is not one.
#'
#' Defaults to `100`.
#' @param ... any arguments (with names different from those of arguments used in [design()]) that are used by `f`, `method`, and `eval`
#'     can be passed here. [design()] will pass relevant arguments to `f`, `method`, and `eval` based on the names of additional arguments provided.
#'
#' @return
#' An updated `object` is returned with a slot called `design` that contains:
#'   - *S* slots, named `wave1, wave2,..., waveS`, that contain information of *S* waves of sequential designs that have been applied to the emulator.
#'     Each slot contains the following elements:
#'     - `N`, an integer that gives the numbers of steps implemented in the corresponding wave;
#'     - `rmse`, a matrix that gives the RMSEs of emulators constructed during the corresponding wave, if `eval = NULL`;
#'     - `metric`, a matrix that gives the customized evaluating metric values of emulators constructed during the corresponding wave,
#'        if a customized function is supplied to `eval`;
#'     - `freq`, an integer that gives the frequency that the emulator validations are implemented during the corresponding wave.
#'     - `enrichment`, a vector of size `N` that gives the number of new design points added after each step of the sequential design (if `object` is
#'        an instance of the `gp` or `dgp` class), or a matrix that gives the number of new design points added to emulators in a bundle after each step of
#'        the sequential design (if `object` is an instance of the `bundle` class).
#'
#'     If `target` is not `NULL`, the following additional elements are also included:
#'     - `target`, the target RMSE(s) to stop the sequential design.
#'     - `reached`, a bool (if `object` is an instance of the `gp` or `dgp` class) or a vector of bools (if `object` is an instance of the `bundle`
#'       class) that indicate if the target RMSEs are reached at the end of the sequential design.
#'   - a slot called `type` that gives the type of validations:
#'     - either LOO (`loo`) or OOS (`oos`) if `eval = NULL`. See [validate()] for more information about LOO and OOS.
#'     - the customized R function provided to `eval`.
#'   - two slots called `x_test` and `y_test` that contain the data points for the OOS validation if the `type` slot is `oos`.
#'
#' See *Note* section below for further information.
#' @note
#' * The re-fitting and validation of an emulator are forced after the final step of a sequential design even `N` is not multiples of elements in `freq`.
#' * Any `loo` or `oos` slot that already exists in `object` will be cleaned, and a new slot called `loo` or `oos` will be created in the returned object
#'   depending on whether `x_test` and `y_test` are provided. The new slot gives the validation information of the emulator constructed in the final step of
#'   the sequential design. See [validate()] for more information about the slots `loo` and `oos`.
#' * If `object` has previously been used by [design()] for sequential designs, the information of the current wave of the sequential design will replace
#'   those of old waves and be contained in the returned object, unless the following conditions are met:
#'   - the validation type (`loo`, `oos`, or the customized function provided to `eval`) of the current wave of the sequential design is the same as the
#'     validation types in previous waves, and
#'   - if the validation type is OOS, `x_test` and `y_test` in the current wave of the sequential design are identical to those in the previous waves.
#'
#'   When the above conditions are met, the information of the current wave of the sequential design will be added to
#'       the `design` slot of the returned object under the name `waveS`.
#' * If `object` is an instance of the `gp` class and `eval = NULL`, the matrix in the `rmse` slot is single-columned. If `object` is an instance of
#'   the `dgp` or `bundle` class and `eval = NULL`, the matrix in the `rmse` slot can have multiple columns that correspond to different output dimensions
#'   or different emulators in the bundle.
#' * If `object` is an instance of the `gp` class and `eval = NULL`, `target` needs to be a single value giving the RMSE threshold. If `object` is an instance
#'   of the `dgp` or `bundle` class and `eval = NULL`, `target` can be a vector of values that gives the RMSE thresholds for different output dimensions or
#'   different emulators. If a single value is provided, it will be used as the RMSE threshold for all output dimensions (if `object` is an instance of the `dgp`) or all emulators
#'   (if `object` is an instance of the `bundle`). If a customized function is supplied to `eval`, the user needs to ensure that the length of `target` is equal
#'   to that of the output from `eval`.
#' * When defining `f`, it is important to ensure that:
#'   - the column order of the first argument of `f` is consistent with the training input used for the emulator;
#'   - the column order of the output matrix of `f` is consistent with the order of emulator output dimensions (if `object` is an instance of the `dgp` class),
#'     or the order of emulators placed in `object` (if `object` is an instance of the `bundle` class).
#' * When defining `eval`, the output metric needs to be positive if [draw()] is used with `log = T`. And one needs to ensure that a lower metric value indicates
#'   a better emulation performance if `target` is set.
#' * Any R vector detected in `x_test` and `y_test` will be treated as a column vector and automatically converted into a single-column
#'   R matrix. Thus, if `x_test` or `y_test` is a single testing data point with multiple dimensions, it must be given as a matrix.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#'
#' @examples
#' \dontrun{
#'
#' # load packages and the Python env
#' library(lhs)
#' library(dgpsi)
#' init_py()
#'
#' # construct a 2D non-stationary function that takes a matrix as the input
#' f <- function(x) {
#'   sin(1/((0.7*x[,1,drop=F]+0.3)*(0.7*x[,2,drop=F]+0.3)))
#' }
#'
#' # generate the initial design
#' X <- maximinLHS(5,2)
#' Y <- f(X)
#'
#' # generate the validation data
#' validate_x <- maximinLHS(30,2)
#' validate_y <- f(validate_x)
#'
#' # training a 2-layered DGP emulator with the initial design
#' m <- dgp(X, Y)
#'
#' # specify the ranges of the input dimensions
#' lim_1 <- c(0, 1)
#' lim_2 <- c(0, 1)
#' lim <- rbind(lim_1, lim_2)
#'
#' # 1st wave of the sequential design with 10 steps
#' m <- design(m, N=10, limits = lim, f = f, x_test = validate_x, y_test = validate_y)
#'
#' # 2nd wave of the sequential design with 10 steps
#' m <- design(m, N=10, limits = lim, f = f, x_test = validate_x, y_test = validate_y)
#'
#' # 3rd wave of the sequential design with 10 steps
#' m <- design(m, N=10, limits = lim, f = f, x_test = validate_x, y_test = validate_y)
#'
#' # draw the design created by the sequential design
#' draw(m,'design')
#'
#' # inspect the trace of RMSEs during the sequential design
#' draw(m,'rmse')
#'
#' # reduce the number of imputations for faster OOS
#' m_faster <- set_imp(m, 10)
#'
#' # plot the OOS validation with the faster DGP emulator
#' plot(m_faster, x_test = validate_x, y_test = validate_y)
#' }
#' @md
#' @name design
#' @export
design <- function(object, N, x_cand, y_cand, n_cand, limits, int, f, freq, x_test, y_test, target, method, eval, verb, check_point, cores, ...){
  UseMethod("design")
}

#' @rdname design
#' @method design gp
#' @export
design.gp <- function(object, N, x_cand = NULL, y_cand = NULL, n_cand = 200, limits = NULL, int = FALSE, f = NULL, freq = c(1, 1), x_test = NULL, y_test = NULL, target = NULL, method = mice, eval = NULL, verb = TRUE, check_point = NULL, cores = 1, ...) {
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)

  N <- check_N(N)
  freq <- check_freq(freq)
  n_cand <- check_n_cand(n_cand)
  if (!is.null(x_test) & !is.null(y_test)) {
    xy_test <- check_xy_test(x_test, y_test)
    x_test <- xy_test[[1]]
    y_test <- xy_test[[2]]
  }
  if ( "design" %in% names(object) ){
    design_info <- object$design
    time_and_wave <- check_design(object, x_test, y_test, eval)
    first_time <- time_and_wave[[1]]
    n_wave <- time_and_wave[[2]]
  } else {
    first_time <- TRUE
    n_wave <- 0
  }

  X <- object$data$X
  Y <- object$data$Y
  n_dim_X <- ncol(X)
  n_dim_Y <- ncol(Y)

  N_acq <- c()
  mnames <- methods::formalArgs(method)
  if ( !is.null(eval) ){
    vnames <- methods::formalArgs(eval)
  }

  if ( is.null(x_cand) ) {
    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)
    fnames <- methods::formalArgs(f)

    limits <- check_limits(limits, n_dim_X)
    if ( is.null(limits) ){
      limits <- matrix(0, n_dim_X, 2)
      limits[,1] <- pkg.env$np$amin(X, axis=0L)
      limits[,2] <- pkg.env$np$amax(X, axis=0L)
    }
    if (identical(method, pei)) {
      add_arg <- list(pseudo_points = pp(X, limits))
      add_arg <- utils::modifyList(add_arg, list(...))
    } else {
      add_arg <- list(...)
    }

    int <- check_int(int, n_dim_X)

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if ( is.null(eval) ){
      if (is.null(x_test) & is.null(y_test)){
        type <- 'loo'
        object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores)
        rmse <- object$loo$rmse
      } else {
        type <- 'oos'
        object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, cores = cores)
        rmse <- object$oos$rmse
      }
    } else {
      type <- eval
      if ("..." %in% vnames){
        rmse <- do.call(eval, c(list(object), add_arg))
      } else {
        vidx <- vnames %in% names(add_arg)
        vparam <- add_arg[vnames[vidx]]
        rmse <- do.call(eval, c(list(object), vparam))
      }
      if (length(rmse)!=1) stop("'eval' must return a vector of length one.", call. = FALSE)
    }
    if ( verb ) message(" done")
    if ( is.null(eval) ){
      if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    } else {
      if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
    }
    rmse_records <- rmse

    target <- check_target(target, length(rmse))

    if ( !is.null(target) ){
      istarget <- rmse<=target
      if ( istarget ){
        run <- FALSE
      } else {
        run <- TRUE
      }
    } else {
      run <- TRUE
    }

    if ( run ){
      for ( i in 1:N ){
        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        x_cand <- lhs::maximinLHS(n_cand,n_dim_X)

        for (j in 1:n_dim_X){
          if ( int[j] ){
            #if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
            x_cand[,j] <- floor( x_cand[,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
          } else {
            x_cand[,j] <- x_cand[,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
          }
        }

        arg_list <- list(object, x_cand)
        if ("..." %in% mnames){
          res <- do.call(method, c(arg_list, add_arg))
        } else {
          midx <- mnames %in% names(add_arg)
          mparam <- add_arg[mnames[midx]]
          res <- do.call(method, c(arg_list, mparam))
        }

        new_X <- x_cand[res,,drop=F]

        if ( verb ) {
          message(" done")
          if ( nrow(new_X)==1 ){
            message(paste(c(" * Next design point:", sprintf("%.06f", new_X[1,])), collapse=" "))
          } else {
            for (j in 1:nrow(new_X) ){
              message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", new_X[j,])), collapse=" "))
            }
          }
        }

        N_acq <- c(N_acq, nrow(new_X))

        X <- rbind(X, new_X)
        if ("..." %in% fnames){
          new_output <- do.call(f, c(list(new_X), add_arg))
        } else {
          fidx <- fnames %in% names(add_arg)
          fparam <- add_arg[fnames[fidx]]
          new_output <- do.call(f, c(list(new_X), fparam))
        }

        if ( is.list(new_output) ){
          Y <- rbind(Y, new_output[[1]])
          if ( length(new_output)!=1 ){
            updated_arg <- new_output[-1]
            add_arg <- utils::modifyList(add_arg, updated_arg)
          }
        } else {
          Y <- rbind(Y, new_output)
        }

        if ( i %% freq[1]==0 | i==N){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = TRUE, verb = FALSE)
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = FALSE, verb = FALSE)
          if ( verb ) Sys.sleep(0.5)
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N ){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores)
              if ( verb ) message(" done")
              rmse <- object$loo$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, cores = cores)
              if ( verb ) message(" done")
              rmse <- object$oos$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            }
          } else {
            if ( verb ) message(" - Validating ...", appendLF = FALSE)
            if ("..." %in% vnames){
              rmse <- do.call(eval, c(list(object), add_arg))
            } else {
              vidx <- vnames %in% names(add_arg)
              vparam <- add_arg[vnames[vidx]]
              rmse <- do.call(eval, c(list(object), vparam))
            }
            if ( verb ) message(" done")
            if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
          }
          rmse_records <- rbind(rmse_records, rmse)
        }

        if( !is.null(target) ) {
          if ( rmse<=target ) {
            message(sprintf("Target reached! The sequential design stops at step %i.", i))
            istarget <- TRUE
            N <- i
            break
          }
        }

        if ( !is.null(check_point) ){
          if (i %in% check_point ){
            ans <- readline(prompt="Do you want to continue? (Y/N) ")
            if ( tolower(ans)=='n'|tolower(ans)=='no' ) {
              message(sprintf("The sequential design is terminated at step %i.", i))
              N <- i
              break
            }
          }
        }
      }
    }
  } else {
    add_arg <- list(...)
    xy_cand_list <- check_xy_cand(x_cand, y_cand, n_dim_Y)
    x_cand <- xy_cand_list[[1]]
    y_cand <- xy_cand_list[[2]]
    x_cand <- remove_dup(x_cand, X)
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- c()

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if ( is.null(eval) ){
      if (is.null(x_test) & is.null(y_test)){
        type <- 'loo'
        object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores)
        rmse <- object$loo$rmse
      } else {
        type <- 'oos'
        object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, cores = cores)
        rmse <- object$oos$rmse
      }
    } else {
      type <- eval
      if ("..." %in% vnames){
        rmse <- do.call(eval, c(list(object), add_arg))
      } else {
        vidx <- vnames %in% names(add_arg)
        vparam <- add_arg[vnames[vidx]]
        rmse <- do.call(eval, c(list(object), vparam))
      }
      if (length(rmse)!=1) stop("'eval' must return a vector of length one.", call. = FALSE)
    }
    if ( verb ) message(" done")
    if ( is.null(eval) ){
      if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    } else {
      if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
    }
    rmse_records <- rmse

    target <- check_target(target, length(rmse))

    if ( !is.null(target) ){
      istarget <- rmse<=target
      if ( istarget ){
        run <- FALSE
      } else {
        run <- TRUE
      }
    } else {
      run <- TRUE
    }

    if ( run ){
      for ( i in 1:N ){
        if ( n_cand<=length(idx_x_cand) ){
          idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
        } else {
          idx_sub_cand <- idx_x_cand
        }
        sub_cand <- x_cand[idx_sub_cand,,drop = FALSE]

        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        arg_list <- list(object, sub_cand)
        if ("..." %in% mnames){
          res <- do.call(method, c(arg_list, add_arg))
        } else {
          midx <- mnames %in% names(add_arg)
          mparam <- add_arg[mnames[midx]]
          res <- do.call(method, c(arg_list, mparam))
        }

        new_X <- sub_cand[res,,drop=F]

        if ( verb ) {
          message(" done")
          if ( nrow(new_X)==1 ){
            message(paste(c(" * Next design point:", sprintf("%.06f", new_X[1,])), collapse=" "))
          } else {
            for (j in 1:nrow(new_X) ){
              message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", new_X[j,])), collapse=" "))
            }
          }
        }

        N_acq <- c(N_acq, nrow(new_X))

        X <- rbind(X, new_X)
        Y <- rbind(Y, y_cand[idx_sub_cand,,drop=FALSE][res,,drop=FALSE])
        idx_x_acq <- c(idx_x_acq, idx_sub_cand[res])
        idx_x_cand <- idx_x_cand0[-idx_x_acq]

        if ( i %% freq[1]==0 | i==N ){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = TRUE, verb = FALSE)
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = FALSE, verb = FALSE)
          if ( verb ) Sys.sleep(0.5)
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N ){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores)
              if ( verb ) message(" done")
              rmse <- object$loo$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, cores = cores)
              if ( verb ) message(" done")
              rmse <- object$oos$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            }
          } else {
            if ( verb ) message(" - Validating ...", appendLF = FALSE)
            if ("..." %in% vnames){
              rmse <- do.call(eval, c(list(object), add_arg))
            } else {
              vidx <- vnames %in% names(add_arg)
              vparam <- add_arg[vnames[vidx]]
              rmse <- do.call(eval, c(list(object), vparam))
            }
            if ( verb ) message(" done")
            if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
          }
          rmse_records <- rbind(rmse_records, rmse)
        }

        if( !is.null(target) ) {
          if ( rmse<=target ) {
            message(sprintf("Target reached! The sequential design stops at step %i.", i))
            istarget <- TRUE
            N <- i
            break
          }
        }

        if ( !is.null(check_point) ){
          if (i %in% check_point ){
            ans <- readline(prompt="Do you want to continue? (Y/N) ")
            if ( tolower(ans)=='n'|tolower(ans)=='no' ) {
              message(sprintf("The sequential design is terminated at step %i.", i))
              N <- i
              break
            }
          }
        }
      }
    }
  }

  if ( run ){
    if ( isTRUE(first_time) ){
      object$design <- list()
      object[["design"]][["wave1"]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][["wave1"]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][["wave1"]][["freq"]] <- freq[2]
      object[["design"]][["wave1"]][["enrichment"]] <- N_acq
      if( !is.null(target) ) {
        object[["design"]][["wave1"]][["target"]] <- target
        object[["design"]][["wave1"]][["reached"]] <- istarget
      }
      object[["design"]][["type"]] <- type
      if ( identical(type, 'oos') ){
        object[["design"]][["x_test"]] <- x_test
        object[["design"]][["y_test"]] <- y_test
      }
    } else {
      object$design <- design_info
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <-  N_acq
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["target"]] <- target
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["reached"]] <- istarget
      }
    }
    if( !is.null(target) ) {
      if ( !istarget ) message("The target is not reached at the end of the sequential design.")
    }
  } else {
    message("Target already reached. The sequential design is not performed.")
  }

  return(object)
}

#' @rdname design
#' @method design dgp
#' @export
design.dgp <- function(object, N, x_cand = NULL, y_cand = NULL, n_cand = 200, limits = NULL, int = FALSE, f = NULL, freq = c(1, 1), x_test = NULL, y_test = NULL, target = NULL, method = mice, eval = NULL, verb = TRUE, check_point = NULL, cores = 1, train_N = 100, ...) {
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)

  N <- check_N(N)
  freq <- check_freq(freq)
  train_N <- check_train_N(train_N, N)
  n_cand <- check_n_cand(n_cand)
  if (!is.null(x_test) & !is.null(y_test)) {
    xy_test <- check_xy_test(x_test, y_test)
    x_test <- xy_test[[1]]
    y_test <- xy_test[[2]]
  }

  if ( "design" %in% names(object) ){
    design_info <- object$design
    time_and_wave <- check_design(object, x_test, y_test, eval)
    first_time <- time_and_wave[[1]]
    n_wave <- time_and_wave[[2]]
  } else {
    first_time <- TRUE
    n_wave <- 0
  }

  X <- object$data$X
  Y <- object$data$Y
  n_dim_X <- ncol(X)
  n_dim_Y <- ncol(Y)

  N_acq <- c()
  mnames <- methods::formalArgs(method)
  if ( !is.null(eval) ){
    vnames <- methods::formalArgs(eval)
  }
  #if candidate set is given
  if ( is.null(x_cand) ) {

    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)
    fnames <- methods::formalArgs(f)

    limits <- check_limits(limits, n_dim_X)
    if ( is.null(limits) ){
      limits <- matrix(0, n_dim_X, 2)
      limits[,1] <- pkg.env$np$amin(X, axis=0L)
      limits[,2] <- pkg.env$np$amax(X, axis=0L)
    }
    if (identical(method, pei)) {
      add_arg <- list(pseudo_points = pp(X, limits))
      add_arg <- utils::modifyList(add_arg, list(...))
    } else {
      add_arg <- list(...)
    }
    int <- check_int(int, n_dim_X)

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if ( is.null(eval) ){
      if (is.null(x_test) & is.null(y_test)){
        type <- 'loo'
        object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores, ...)
        rmse <- object$loo$rmse
      } else {
        type <- 'oos'
        object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, cores = cores, ...)
        rmse <- object$oos$rmse
      }
    } else {
      type <- eval
      if ("..." %in% vnames){
        rmse <- do.call(eval, c(list(object), add_arg))
      } else {
        vidx <- vnames %in% names(add_arg)
        vparam <- add_arg[vnames[vidx]]
        rmse <- do.call(eval, c(list(object), vparam))
      }
      if (length(rmse)!=1 & length(rmse)!=n_dim_Y) stop(sprintf("'eval' must return a single or %i metric values.", n_dim_Y), call. = FALSE)
    }
    if ( verb ) message(" done")
    if ( is.null(eval) ){
      if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    } else {
      if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
    }
    rmse_records <- rmse

    target <- check_target(target, length(rmse))

    if ( !is.null(target) ){
      istarget <- all(rmse<=target)
      if ( istarget ){
        run <- FALSE
      } else {
        run <- TRUE
      }
    } else {
      run <- TRUE
    }

    if ( run ){
      for ( i in 1:N ){
        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        x_cand <- lhs::maximinLHS(n_cand,n_dim_X)

        for (j in 1:n_dim_X){
          if ( int[j] ){
            #if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
            x_cand[,j] <- floor( x_cand[,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
          } else {
            x_cand[,j] <- x_cand[,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
          }
        }

        arg_list <- list(object, x_cand)
        if ("..." %in% mnames){
          res <- do.call(method, c(arg_list, add_arg))
        } else {
          midx <- mnames %in% names(add_arg)
          mparam <- add_arg[mnames[midx]]
          res <- do.call(method, c(arg_list, mparam))
        }

        if ( verb ) {
          message(" done")
          if ( nrow(res)==1 ){
            if ( ncol(res)==1 ){
              message(paste(c(" * Next design point:", sprintf("%.06f", x_cand[res[1,1],])), collapse=" "))
            } else {
              for ( j in 1:ncol(res) ){
                message(paste(c(sprintf(" * Next design point (Output%i):", j), sprintf("%.06f", x_cand[res[1,j],])), collapse=" "))
              }
            }
          } else {
            for ( j in 1:nrow(res) ){
              if ( ncol(res)==1 ){
                message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", x_cand[res[j,1],])), collapse=" "))
              } else {
                for ( k in 1:nrow(res) ){
                  message(paste(c(sprintf(" * Next design point (Position%i for Output%i):", j, k), sprintf("%.06f", x_cand[res[j,k],])), collapse=" "))
                }
              }
            }
          }
        }

        new_X <- x_cand[unique(as.vector(t(res))),,drop=FALSE]

        N_acq <- c(N_acq, nrow(new_X))

        X <- rbind(X, new_X)
        if ("..." %in% fnames){
          new_output <- do.call(f, c(list(new_X), add_arg))
        } else {
          fidx <- fnames %in% names(add_arg)
          fparam <- add_arg[fnames[fidx]]
          new_output <- do.call(f, c(list(new_X), fparam))
        }

        if ( is.list(new_output) ){
          Y <- rbind(Y, new_output[[1]])
          if ( length(new_output)!=1 ){
            updated_coeff <- new_output[-1]
            add_arg <- utils::modifyList(add_arg, updated_coeff)
          }
        } else {
          Y <- rbind(Y, new_output)
        }

        if ( i %% freq[1]==0 | i==N){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = TRUE, verb = FALSE, N = train_N[i], B = 10)
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = FALSE, verb = FALSE, B = 10)
          if ( verb ) Sys.sleep(0.5)
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores, ...)
              if ( verb ) message(" done")
              rmse <- object$loo$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, cores = cores, ...)
              if ( verb ) message(" done")
              rmse <- object$oos$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            }
          } else {
            if ( verb ) message(" - Validating ...", appendLF = FALSE)
            if ("..." %in% vnames){
              rmse <- do.call(eval, c(list(object), add_arg))
            } else {
              vidx <- vnames %in% names(add_arg)
              vparam <- add_arg[vnames[vidx]]
              rmse <- do.call(eval, c(list(object), vparam))
            }
            if ( verb ) message(" done")
            if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
          }
          rmse_records <- rbind(rmse_records, rmse)
        }

        if( !is.null(target) ) {
          if ( all(rmse<=target) ) {
            message(sprintf("Target reached! The sequential design stops at step %i.", i))
            istarget <- TRUE
            N <- i
            break
          }
        }

        if ( !is.null(check_point) ){
          if (i %in% check_point ){
            ans <- readline(prompt="Do you want to continue? (Y/N) ")
            if ( tolower(ans)=='n'|tolower(ans)=='no' ) {
              message(sprintf("The sequential design is terminated at step %i.", i))
              N <- i
              break
            }
          }
        }
      }
    }
  } else {
    add_arg <- list(...)
    xy_cand_list <- check_xy_cand(x_cand, y_cand, n_dim_Y)
    x_cand <- xy_cand_list[[1]]
    y_cand <- xy_cand_list[[2]]
    x_cand <- remove_dup(x_cand, X)
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- c()

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if ( is.null(eval) ){
      if (is.null(x_test) & is.null(y_test)){
        type <- 'loo'
        object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores, ...)
        rmse <- object$loo$rmse
      } else {
        type <- 'oos'
        object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, cores = cores, ...)
        rmse <- object$oos$rmse
      }
    } else {
      type <- eval
      if ("..." %in% vnames){
        rmse <- do.call(eval, c(list(object), add_arg))
      } else {
        vidx <- vnames %in% names(add_arg)
        vparam <- add_arg[vnames[vidx]]
        rmse <- do.call(eval, c(list(object), vparam))
      }
      if (length(rmse)!=1 & length(rmse)!=n_dim_Y) stop(sprintf("'eval' must return a single or %i metric values.", n_dim_Y), call. = FALSE)
    }
    if ( verb ) message(" done")
    if ( is.null(eval) ){
      if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    } else {
      if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
    }
    rmse_records <- rmse

    target <- check_target(target, length(rmse))

    if ( !is.null(target) ){
      istarget <- all(rmse<=target)
      if ( istarget ){
        run <- FALSE
      } else {
        run <- TRUE
      }
    } else {
      run <- TRUE
    }

    if ( run ){
      for ( i in 1:N ){
        if ( n_cand<=length(idx_x_cand) ){
          idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
        } else {
          idx_sub_cand <- idx_x_cand
        }
        sub_cand <- x_cand[idx_sub_cand,,drop = FALSE]

        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        arg_list <- list(object, sub_cand)
        if ("..." %in% mnames){
          res <- do.call(method, c(arg_list, add_arg))
        } else {
          midx <- mnames %in% names(add_arg)
          mparam <- add_arg[mnames[midx]]
          res <- do.call(method, c(arg_list, mparam))
        }

        if ( verb ) {
          message(" done")
          if ( nrow(res)==1 ){
            if ( ncol(res)==1 ){
              message(paste(c(" * Next design point:", sprintf("%.06f", sub_cand[res[1,1],])), collapse=" "))
            } else {
              for ( j in 1:ncol(res) ){
                message(paste(c(sprintf(" * Next design point (Output%i):", j), sprintf("%.06f", sub_cand[res[1,j],])), collapse=" "))
              }
            }
          } else {
            for ( j in 1:nrow(res) ){
              if ( ncol(res)==1 ){
                message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", sub_cand[res[j,1],])), collapse=" "))
              } else {
                for ( k in 1:nrow(res) ){
                  message(paste(c(sprintf(" * Next design point (Position%i for Output%i):", j, k), sprintf("%.06f", sub_cand[res[j,k],])), collapse=" "))
                }
              }
            }
          }
        }

        new_X <- sub_cand[unique(as.vector(t(res))),,drop=FALSE]

        N_acq <- c(N_acq, nrow(new_X))

        X <- rbind(X, new_X)
        Y <- rbind(Y, y_cand[idx_sub_cand,,drop=FALSE][unique(as.vector(t(res))),,drop=FALSE])
        idx_x_acq <- c(idx_x_acq, idx_sub_cand[unique(as.vector(t(res)))])
        idx_x_cand <- idx_x_cand0[-idx_x_acq]

        if ( i %% freq[1]==0 | i==N ){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = TRUE, verb = FALSE, N = train_N[i], B = 10)
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = FALSE, verb = FALSE, B = 10)
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N ){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores, ...)
              if ( verb ) message(" done")
              rmse <- object$loo$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, cores = cores, ...)
              if ( verb ) message(" done")
              rmse <- object$oos$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            }
          } else {
            if ( verb ) message(" - Validating ...", appendLF = FALSE)
            if ("..." %in% vnames){
              rmse <- do.call(eval, c(list(object), add_arg))
            } else {
              vidx <- vnames %in% names(add_arg)
              vparam <- add_arg[vnames[vidx]]
              rmse <- do.call(eval, c(list(object), vparam))
            }
            if ( verb ) message(" done")
            if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
          }
          rmse_records <- rbind(rmse_records, rmse)
        }

        if( !is.null(target) ) {
          if ( all(rmse<=target) ) {
            message(sprintf("Target reached! The sequential design stops at step %i.", i))
            istarget <- TRUE
            N <- i
            break
          }
        }

        if ( !is.null(check_point) ){
          if (i %in% check_point ){
            ans <- readline(prompt="Do you want to continue? (Y/N) ")
            if ( tolower(ans)=='n'|tolower(ans)=='no' ) {
              message(sprintf("The sequential design is terminated at step %i.", i))
              N <- i
              break
            }
          }
        }
      }
    }
  }

  if ( run ){
    if ( isTRUE(first_time) ){
      object$design <- list()
      object[["design"]][["wave1"]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][["wave1"]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][["wave1"]][["freq"]] <- freq[2]
      object[["design"]][["wave1"]][["enrichment"]] <- N_acq
      if( !is.null(target) ) {
        object[["design"]][["wave1"]][["target"]] <- target
        object[["design"]][["wave1"]][["reached"]] <- istarget
      }
      object[["design"]][["type"]] <- type
      if ( identical(type, 'oos') ){
        object[["design"]][["x_test"]] <- x_test
        object[["design"]][["y_test"]] <- y_test
      }
    } else {
      object$design <- design_info
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <- N_acq
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["target"]] <- target
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["reached"]] <- istarget
      }
    }
    if( !is.null(target) ) {
      if ( !istarget ) message("The target is not reached at the end of the sequential design.")
    }
  } else {
    message("Target already reached. The sequential design is not performed.")
  }

  return(object)
}


#' @rdname design
#' @method design bundle
#' @export
design.bundle <- function(object, N, x_cand = NULL, y_cand = NULL, n_cand = 200, limits = NULL, int = FALSE, f = NULL, freq = c(1, 1), x_test = NULL, y_test = NULL, target = NULL, method = mice, eval = NULL, verb = TRUE, check_point = NULL, cores = 1, train_N = 100, ...) {
  if ( !inherits(object,"bundle") ) stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)

  N <- check_N(N)
  freq <- check_freq(freq)
  train_N <- check_train_N(train_N, N)
  n_cand <- check_n_cand(n_cand)
  if (!is.null(x_test) & !is.null(y_test)) {
    xy_test <- check_xy_test(x_test, y_test)
    x_test <- xy_test[[1]]
    y_test <- xy_test[[2]]
  }
  if ( "design" %in% names(object) ){
    design_info <- object$design
    time_and_wave <- check_design(object, x_test, y_test, eval)
    first_time <- time_and_wave[[1]]
    n_wave <- time_and_wave[[2]]
  } else {
    first_time <- TRUE
    n_wave <- 0
  }

  X <- object$data$X
  Y <- object$data$Y
  n_dim_X <- ncol(X[[1]])
  n_emulators <- length(object) - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1

  N_acq_ind <- c()
  mnames <- methods::formalArgs(method)
  if ( !is.null(eval) ){
    vnames <- methods::formalArgs(eval)
  }

  #if candidate set is given
  if ( is.null(x_cand) ) {

    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)
    fnames <- methods::formalArgs(f)

    limits <- check_limits(limits, n_dim_X)
    if ( is.null(limits) ){
      all_training_input <- c()
      for ( k in 1:n_emulators ){
        all_training_input <- rbind(all_training_input, X[[k]])
      }
      limits <- matrix(0, n_dim_X, 2)
      limits[,1] <- pkg.env$np$amin(all_training_input, axis=0L)
      limits[,2] <- pkg.env$np$amax(all_training_input, axis=0L)
    }
    if (identical(method, pei)) {
      ppoints <- list()
      for ( k in 1:n_emulators ){
        ppoints[[k]] <- pp(X[[k]], limits)
      }
      add_arg <- list(pseudo_points = ppoints)
      add_arg <- utils::modifyList(add_arg, list(...))
    } else {
      add_arg <- list(...)
    }
    int <- check_int(int, n_dim_X)

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if ( is.null(eval) ){
      rmse <- c()
      for ( k in 1:n_emulators ){
        obj_k <- object[[paste('emulator',k,sep='')]]
        if (is.null(x_test) & is.null(y_test)){
          type <- 'loo'
          if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE)
          if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores, ...)
          object[[paste('emulator',k,sep='')]] <- obj_k
          rmse <- c(rmse, obj_k$loo$rmse)
        } else {
          type <- 'oos'
          if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE)
          if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, cores = cores, ...)
          object[[paste('emulator',k,sep='')]] <- obj_k
          rmse <- c(rmse, obj_k$oos$rmse)
        }
      }
    } else {
      type <- eval
      if ("..." %in% vnames){
        rmse <- do.call(eval, c(list(object), add_arg))
      } else {
        vidx <- vnames %in% names(add_arg)
        vparam <- add_arg[vnames[vidx]]
        rmse <- do.call(eval, c(list(object), vparam))
      }
      if (length(rmse)!=1 & length(rmse)!=n_emulators) stop(sprintf("'eval' must return a single or %i metric values.", n_emulators), call. = FALSE)
    }
    if ( verb ) message(" done")
    if ( is.null(eval) ){
      if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    } else {
      if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
    }
    rmse_records <- rmse

    target <- check_target(target, length(rmse))

    if ( !is.null(target) ){
      istarget <- rmse<=target
      if ( length(istarget)==1 ) istarget <- rep(istarget, n_emulators)
      if ( all(istarget) ){
        run <- FALSE
      } else {
        run <- TRUE
      }
    } else {
      run <- TRUE
    }

    if ( run ){
      for ( i in 1:N ){
        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        x_cand <- lhs::maximinLHS(n_cand,n_dim_X)

        for (j in 1:n_dim_X){
          if ( int[j] ){
            #if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
            x_cand[,j] <- floor( x_cand[,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
          } else {
            x_cand[,j] <- x_cand[,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
          }
        }

        arg_list <- list(object, x_cand)
        if ("..." %in% mnames){
          res <- do.call(method, c(arg_list, add_arg))
        } else {
          midx <- mnames %in% names(add_arg)
          mparam <- add_arg[mnames[midx]]
          res <- do.call(method, c(arg_list, mparam))
        }

        if ( verb ) {
          message(" done")
          if ( nrow(res)==1 ){
            for ( j in 1:n_emulators ){
              if ( is.null(target) ){
                message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", x_cand[res[1,j],])), collapse=" "))
              } else {
                if ( !istarget[j] ){
                  message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", x_cand[res[1,j],])), collapse=" "))
                } else {
                  message(sprintf(" * Next design point (Emulator%i): None (target reached)", j))
                }
              }
            }
          } else {
            for ( j in 1:nrow(res) ){
              for ( k in 1:n_emulators ){
                if ( is.null(target) ){
                  message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", x_cand[res[j,k],])), collapse=" "))
                } else {
                  if ( !istarget[k] ){
                    message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", x_cand[res[j,k],])), collapse=" "))
                  } else {
                    message(sprintf(" * Next design point (Position%i for Emulator%i): None (target reached)", j, k))
                  }
                }
              }
            }
          }
        }

        if ( nrow(res)==1 ){
          new_X <- x_cand[res[1,],,drop=FALSE]
          if ( !is.null(target) ){
            new_X <- new_X[!istarget,,drop=FALSE]
          }
          rep_new_X <- pkg.env$np$unique(new_X, return_inverse=TRUE, axis=0L)
          new_X_unique <- rep_new_X[[1]]
          rep_idx <- rep_new_X[[2]] + 1
          if ("..." %in% fnames){
            new_output_unique <- do.call(f, c(list(new_X_unique), add_arg))
          } else {
            fidx <- fnames %in% names(add_arg)
            fparam <- add_arg[fnames[fidx]]
            new_output_unique <- do.call(f, c(list(new_X_unique), fparam))
          }
          if ( is.list(new_output_unique) ){
            new_Y_unique <- new_output_unique[[1]]
            if ( length(new_output_unique)!=1 ){
              updated_coeff <- new_output_unique[-1]
              add_arg <- utils::modifyList(add_arg, updated_coeff)
            }
          } else {
            new_Y_unique <- new_output_unique
          }
          new_Y <- new_Y_unique[rep_idx,,drop=F]
          temp <- rep(1, n_emulators)
          if ( !is.null(target) ) temp[istarget] <- 0
          N_acq_ind <- rbind(N_acq_ind, temp)

          if ( !is.null(target) ) ctr <- 1
          for ( j in 1:n_emulators ){
            if ( is.null(target) ){
              X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X[j,,drop=FALSE])
              Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[j,j,drop=FALSE])
            } else {
              if ( !istarget[j] ) {
                X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X[ctr,,drop=FALSE])
                Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[ctr,j,drop=FALSE])
                ctr <- ctr + 1
              }
            }
          }
        } else {
          new_X <- c()
          for (j in 1:nrow(res) ){
            if ( is.null(target) ) {
              new_X <- rbind(new_X, x_cand[res[j,],,drop=FALSE])
            } else {
              new_X <- rbind(new_X, x_cand[res[j,],,drop=FALSE][!istarget,,drop=FALSE])
            }
          }
          rep_new_X <- pkg.env$np$unique(new_X, return_inverse=TRUE, axis=0L)
          new_X_unique <- rep_new_X[[1]]
          rep_idx <- rep_new_X[[2]] + 1
          if ("..." %in% fnames){
            new_output_unique <- do.call(f, c(list(new_X_unique), add_arg))
          } else {
            fidx <- fnames %in% names(add_arg)
            fparam <- add_arg[fnames[fidx]]
            new_output_unique <- do.call(f, c(list(new_X_unique), fparam))
          }
          if ( is.list(new_output_unique) ){
            new_Y_unique <- new_output_unique[[1]]
            if ( length(new_output_unique)!=1 ){
              updated_coeff <- new_output_unique[-1]
              add_arg <- utils::modifyList(add_arg, updated_coeff)
            }
          } else {
            new_Y_unique <- new_output_unique
          }

          new_Y <- new_Y_unique[rep_idx,,drop=F]

          temp <- rep(nrow(res), n_emulators)
          if ( !is.null(target) ) temp[istarget] <- 0
          N_acq_ind <- rbind(N_acq_ind, temp)

          if ( !is.null(target) ) {
            ctr <- 1
            active_emu <- sum(!istarget)
          }
          for ( j in 1:n_emulators ){
            if ( is.null(target) ){
              extract <- j + seq(0, by = n_emulators, length = nrow(res))
              X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X[extract,,drop=FALSE])
              Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[extract,j,drop=FALSE])
            } else {
              if ( !istarget[j] ) {
                extract <- ctr + seq(0, by = active_emu, length = nrow(res))
                X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X[extract,,drop=FALSE])
                Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[extract,j,drop=FALSE])
                ctr <- ctr + 1
              }
            }
          }
        }

        if ( i %% freq[1]==0 | i==N){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            if ( is.null(target) ) {
              obj_k <- object[[paste('emulator',k,sep='')]]
              if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, verb = FALSE)
              if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, verb = FALSE, N = train_N[i], B = 10)
              object[[paste('emulator',k,sep='')]] <- obj_k
            } else {
              if ( !istarget[k] ){
                obj_k <- object[[paste('emulator',k,sep='')]]
                if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, verb = FALSE)
                if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, verb = FALSE, N = train_N[i], B = 10)
                object[[paste('emulator',k,sep='')]] <- obj_k
              }
            }
          }
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            if ( is.null(target) ) {
              obj_k <- object[[paste('emulator',k,sep='')]]
              if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = FALSE, verb = FALSE)
              if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = FALSE, verb = FALSE, B = 10)
              object[[paste('emulator',k,sep='')]] <- obj_k
            } else {
              if ( !istarget[k] ){
                obj_k <- object[[paste('emulator',k,sep='')]]
                if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = FALSE, verb = FALSE)
                if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = FALSE, verb = FALSE, B = 10)
                object[[paste('emulator',k,sep='')]] <- obj_k
              }
            }
          }
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              for ( k in 1:n_emulators ){
                if ( is.null(target) ) {
                  obj_k <- object[[paste('emulator',k,sep='')]]
                  if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE)
                  if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores, ...)
                  object[[paste('emulator',k,sep='')]] <- obj_k
                  rmse[k] <- obj_k$loo$rmse
                } else {
                  if ( !istarget[k] ){
                    obj_k <- object[[paste('emulator',k,sep='')]]
                    if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE)
                    if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores, ...)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    rmse[k] <- obj_k$loo$rmse
                  }
                }
              }
              if ( verb ) message(" done")
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              for ( k in 1:n_emulators ){
                if ( is.null(target) ) {
                  obj_k <- object[[paste('emulator',k,sep='')]]
                  if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE)
                  if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE, cores = cores, ...)
                  object[[paste('emulator',k,sep='')]] <- obj_k
                  rmse[k] <- obj_k$oos$rmse
                } else {
                  if ( !istarget[k] ){
                    obj_k <- object[[paste('emulator',k,sep='')]]
                    if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE)
                    if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE, cores = cores, ...)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    rmse[k] <- obj_k$oos$rmse
                  }
                }
              }
              if ( verb ) message(" done")
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            }
          } else {
            if ( verb ) message(" - Validating ...", appendLF = FALSE)
            if ("..." %in% vnames){
              rmse_tmp <- do.call(eval, c(list(object), add_arg))
            } else {
              vidx <- vnames %in% names(add_arg)
              vparam <- add_arg[vnames[vidx]]
              rmse_tmp <- do.call(eval, c(list(object), vparam))
            }
            if ( is.null(target) ) {
              rmse <- rmse_tmp
            } else {
              if( length(rmse_tmp)==1 ){
                rmse <- rmse_tmp
              } else {
                rmse[!istarget] <- rmse_tmp[!istarget]
              }
            }
            if ( verb ) message(" done")
            if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
          }
          rmse_records <- rbind(rmse_records, rmse)
        }

        if( !is.null(target) ) {
          if ( all(rmse<=target) ) {
            message(sprintf("Target reached! The sequential design stops at step %i.", i))
            istarget <- rep(TRUE, n_emulators)
            N <- i
            break
          } else {
            istarget <- rmse<=target
            if ( length(istarget)==1 ) istarget <- rep(istarget, n_emulators)
          }
        }

        if ( !is.null(check_point) ){
          if (i %in% check_point ){
            ans <- readline(prompt="Do you want to continue? (Y/N) ")
            if ( tolower(ans)=='n'|tolower(ans)=='no' ) {
              message(sprintf("The sequential design is terminated at step %i.", i))
              N <- i
              break
            }
          }
        }

      }
    }
  } else {
    add_arg <- list(...)
    xy_cand_list <- check_xy_cand(x_cand, y_cand, n_emulators)
    x_cand <- xy_cand_list[[1]]
    y_cand <- xy_cand_list[[2]]

    for (j in 1:n_emulators){
      x_cand <- remove_dup(x_cand, X[[paste('emulator',j,sep="")]])
    }
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- c()

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if ( is.null(eval) ){
      rmse <- c()
      for ( k in 1:n_emulators ){
        obj_k <- object[[paste('emulator',k,sep='')]]
        if (is.null(x_test) & is.null(y_test)){
          type <- 'loo'
          if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE)
          if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores, ...)
          object[[paste('emulator',k,sep='')]] <- obj_k
          rmse <- c(rmse, obj_k$loo$rmse)
        } else {
          type <- 'oos'
          if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE)
          if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, cores = cores, ...)
          object[[paste('emulator',k,sep='')]] <- obj_k
          rmse <- c(rmse, obj_k$oos$rmse)
        }
      }
    } else {
      type <- eval
      if ("..." %in% vnames){
        rmse <- do.call(eval, c(list(object), add_arg))
      } else {
        vidx <- vnames %in% names(add_arg)
        vparam <- add_arg[vnames[vidx]]
        rmse <- do.call(eval, c(list(object), vparam))
      }
      if (length(rmse)!=1 & length(rmse)!=n_emulators) stop(sprintf("'eval' must return a single or %i metric values.", n_emulators), call. = FALSE)
    }
    if ( verb ) message(" done")
    if ( is.null(eval) ){
      if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    } else {
      if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
    }
    rmse_records <- rmse

    target <- check_target(target, length(rmse))

    if ( !is.null(target) ){
      istarget <- rmse<=target
      if ( length(istarget)==1 ) istarget <- rep(istarget, n_emulators)
      if ( all(istarget) ){
        run <- FALSE
      } else {
        run <- TRUE
      }
    } else {
      run <- TRUE
    }

    if ( run ){
      for ( i in 1:N ){
        if ( n_cand<=length(idx_x_cand) ){
          idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
        } else {
          idx_sub_cand <- idx_x_cand
        }
        sub_cand <- x_cand[idx_sub_cand,,drop = FALSE]
        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        arg_list <- list(object, sub_cand)
        if ("..." %in% mnames){
          res <- do.call(method, c(arg_list, add_arg))
        } else {
          midx <- mnames %in% names(add_arg)
          mparam <- add_arg[mnames[midx]]
          res <- do.call(method, c(arg_list, mparam))
        }

        if ( verb ) {
          message(" done")
          if ( nrow(res)==1 ){
            for ( j in 1:n_emulators ){
              if ( is.null(target) ){
                message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", sub_cand[res[1,j],])), collapse=" "))
              } else {
                if ( !istarget[j] ){
                  message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", sub_cand[res[1,j],])), collapse=" "))
                } else {
                  message(sprintf(" * Next design point (Emulator%i): None (target reached)", j))
                }
              }
            }
          } else {
            for ( j in 1:nrow(res) ){
              for ( k in 1:n_emulators ){
                if ( is.null(target) ){
                  message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", sub_cand[res[j,k],])), collapse=" "))
                } else {
                  if ( !istarget[k] ){
                    message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", sub_cand[res[j,k],])), collapse=" "))
                  } else {
                    message(sprintf(" * Next design point (Position%i for Emulator%i): None (target reached)", j, k))
                  }
                }
              }
            }
          }
        }

        if ( nrow(res)==1 ){
          new_X <- sub_cand[res[1,],,drop=FALSE]
          new_idx <- res[1,]
          new_Y <- y_cand[idx_sub_cand,,drop=FALSE][new_idx,,drop=FALSE]

          temp <- rep(1, n_emulators)
          if ( !is.null(target) ) temp[istarget] <- 0
          N_acq_ind <- rbind(N_acq_ind, temp)

          for ( j in 1:n_emulators ){
            if ( is.null(target) ){
              X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X[j,,drop=FALSE])
              Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[j,j,drop=FALSE])
            } else {
              if ( !istarget[j] ) {
                X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X[j,,drop=FALSE])
                Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[j,j,drop=FALSE])
              }
            }
          }
          if ( is.null(target) ){
            idx_x_acq <- c(idx_x_acq, idx_sub_cand[unique(new_idx)])
          } else {
            idx_x_acq <- c(idx_x_acq, idx_sub_cand[unique(new_idx[!istarget])])
          }
          idx_x_cand <- idx_x_cand0[-idx_x_acq]
        } else {
            new_X <- c()
            new_idx <- c()
            for (j in 1:nrow(res) ){
              new_X <- rbind(new_X, sub_cand[res[j,],,drop=FALSE])
              new_idx <- c(new_idx, res[j,])
            }
            new_Y <- y_cand[idx_sub_cand,,drop=FALSE][new_idx,,drop=FALSE]

            temp <- rep(nrow(res), n_emulators)
            if ( !is.null(target) ) temp[istarget] <- 0
            N_acq_ind <- rbind(N_acq_ind, temp)

            for ( j in 1:n_emulators ){
              extract <- j + seq(0, by = n_emulators, length = nrow(res))
              if ( is.null(target) ){
                X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X[extract,,drop=FALSE])
                Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[extract,j,drop=FALSE])
              } else {
                if ( !istarget[j] ) {
                  X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X[extract,,drop=FALSE])
                  Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[extract,j,drop=FALSE])
                }
              }
            }
            if ( is.null(target) ){
              idx_x_acq <- c(idx_x_acq, idx_sub_cand[unique(new_idx)])
            } else {
              mask <- rep(!istarget, nrow(res))
              idx_x_acq <- c(idx_x_acq, idx_sub_cand[unique(new_idx[mask])])
            }
            idx_x_cand <- idx_x_cand0[-idx_x_acq]
        }

        if ( i %% freq[1]==0 | i==N){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            if ( is.null(target) ){
              obj_k <- object[[paste('emulator',k,sep='')]]
              if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, verb = FALSE)
              if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, verb = FALSE, N = train_N[i], B = 10)
              object[[paste('emulator',k,sep='')]] <- obj_k
            } else {
              if ( !istarget[k] ){
                obj_k <- object[[paste('emulator',k,sep='')]]
                if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, verb = FALSE)
                if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, verb = FALSE, N = train_N[i], B = 10)
                object[[paste('emulator',k,sep='')]] <- obj_k
              }
            }
          }
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            if ( is.null(target) ){
              obj_k <- object[[paste('emulator',k,sep='')]]
              if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, verb = FALSE)
              if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, verb = FALSE, B = 10)
              object[[paste('emulator',k,sep='')]] <- obj_k
            } else {
              if ( !istarget[k] ){
                obj_k <- object[[paste('emulator',k,sep='')]]
                if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, verb = FALSE)
                if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, verb = FALSE, B = 10)
                object[[paste('emulator',k,sep='')]] <- obj_k
              }
            }
          }
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N){
          #rmse <- c()
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              for ( k in 1:n_emulators ){
                if ( is.null(target) ) {
                  obj_k <- object[[paste('emulator',k,sep='')]]
                  if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE)
                  if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores, ...)
                  object[[paste('emulator',k,sep='')]] <- obj_k
                  rmse[k] <- obj_k$loo$rmse
                } else {
                  if ( !istarget[k] ){
                    obj_k <- object[[paste('emulator',k,sep='')]]
                    if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE)
                    if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores, ...)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    rmse[k] <- obj_k$loo$rmse
                  }
                }
              }
              if ( verb ) message(" done")
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              for ( k in 1:n_emulators ){
                if ( is.null(target) ) {
                  obj_k <- object[[paste('emulator',k,sep='')]]
                  if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE)
                  if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE, cores = cores, ...)
                  object[[paste('emulator',k,sep='')]] <- obj_k
                  rmse[k] <- obj_k$oos$rmse
                } else {
                  if ( !istarget[k] ){
                    obj_k <- object[[paste('emulator',k,sep='')]]
                    if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE)
                    if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE, cores = cores, ...)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    rmse[k] <- obj_k$oos$rmse
                  }
                }
              }
              if ( verb ) message(" done")
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            }
          } else {
            if ( verb ) message(" - Validating ...", appendLF = FALSE)
            if ("..." %in% vnames){
              rmse_tmp <- do.call(eval, c(list(object), add_arg))
            } else {
              vidx <- vnames %in% names(add_arg)
              vparam <- add_arg[vnames[vidx]]
              rmse_tmp <- do.call(eval, c(list(object), vparam))
            }
            if ( is.null(target) ) {
              rmse <- rmse_tmp
            } else {
              if( length(rmse_tmp)==1 ){
                rmse <- rmse_tmp
              } else {
                rmse[!istarget] <- rmse_tmp[!istarget]
              }
            }
            if ( verb ) message(" done")
            if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
          }
          rmse_records <- rbind(rmse_records, rmse)
        }

        if( !is.null(target) ) {
          if ( all(rmse<=target) ) {
            message(sprintf("Target reached! The sequential design stops at step %i.", i))
            istarget <- rep(TRUE, n_emulators)
            N <- i
            break
          } else {
            istarget <- rmse<=target
            if ( length(istarget)==1 ) istarget <- rep(istarget, n_emulators)
          }
        }

        if ( !is.null(check_point) ){
          if (i %in% check_point ){
            ans <- readline(prompt="Do you want to continue? (Y/N) ")
            if ( tolower(ans)=='n'|tolower(ans)=='no' ) {
              message(sprintf("The sequential design is terminated at step %i.", i))
              N <- i
              break
            }
          }
        }

      }
    }
  }

  if ( run ){
    if ( isTRUE(first_time) ){
      object$design <- list()
      object[["design"]][["wave1"]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][["wave1"]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][["wave1"]][["freq"]] <- freq[2]
      object[["design"]][["wave1"]][["enrichment"]] <- unname(N_acq_ind)
      if( !is.null(target) ) {
        object[["design"]][["wave1"]][["target"]] <- target
        object[["design"]][["wave1"]][["reached"]] <- ifelse(length(target)==1, all(istarget), istarget)
      }
      object[["design"]][["type"]] <- type
      if ( identical(type, 'oos') ){
        object[["design"]][["x_test"]] <- x_test
        object[["design"]][["y_test"]] <- y_test
      }
      object[['data']][['X']] <- X
      object[['data']][['Y']] <- Y
    } else {
      object$design <- design_info
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
      if ( is.null(eval) ){
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
      } else {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["metric"]] <- unname(rmse_records)
      }
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <- unname(N_acq_ind)
      if( !is.null(target) ) {
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["target"]] <- target
        object[["design"]][[paste('wave', n_wave+1, sep="")]][["reached"]] <- ifelse(length(target)==1, all(istarget), istarget)
      }
      object[['data']][['X']] <- X
      object[['data']][['Y']] <- Y
    }
    if( !is.null(target) ) {
      if ( !all(istarget) ) message("Targets are not reached for all emulators at the end of the sequential design.")
    }
  } else {
    message("Target already reached. The sequential design is not performed.")
  }

  return(object)
}

#check argument N
check_N <- function(N){
  N <- as.integer(N)
  if ( N < 1 ) stop("The number of steps 'N' for a sequential design must be >= 1.", call. = FALSE)
  return(N)
}

#check argument freq
check_freq <- function(freq){
  freq <- as.integer(freq)
  if ( any(freq < 1) ) stop("All elements in 'freq' must be greater than or equal to 1.", call. = FALSE)
  return(freq)
}

#check argument n_cand
check_n_cand <- function(n_cand){
  n_cand <- as.integer(n_cand)
  if ( n_cand < 1 ) stop("'n_cand' must be greater than or equal to 1.", call. = FALSE)
  return(n_cand)
}

#check argument x_cand and y_cand
check_xy_cand <- function(x_cand, y_cand, n_dim_Y){
  if ( is.null(y_cand) ) stop("'y_cand' must be provided if 'x_cand' is not NULL.", call. = FALSE)
  if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(y_cand)&!is.vector(y_cand) ) stop("'y_cand' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
  if ( is.vector(y_cand) ) y_cand <- as.matrix(y_cand)
  if ( nrow(x_cand)!=nrow(y_cand) ) stop("'x_cand' and 'y_cand' have different number of data points.", call. = FALSE)
  if ( ncol(y_cand)!=n_dim_Y ) stop(sprintf("The dimension of 'y_cand' must be %i.", n_dim_Y), call. = FALSE)
  return(list(x_cand, y_cand))
}

#remove duplicates between x_cand and X
remove_dup <- function(x_cand, X){
  X_s <- apply(X, 1, paste, collapse = ", ")
  x_cand_s <- apply(x_cand, 1, paste, collapse = ", ")
  X_x_cand <- intersect(X_s, x_cand_s)
  if ( length(X_x_cand)!=0 ){
    x_cand <- x_cand[!x_cand_s %in% X_x_cand,,drop=FALSE]
  }
  return(x_cand)
}

#check argument x_test and y_test
check_xy_test <- function(x_test, y_test){
  x_test <- unname(x_test)
  y_test <- unname(y_test)
  if ( !is.matrix(x_test)&!is.vector(x_test) ) stop("'x_test' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(y_test)&!is.vector(y_test) ) stop("'y_test' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x_test) ) x_test <- as.matrix(x_test)
  if ( is.vector(y_test) ) y_test <- as.matrix(y_test)
  if ( nrow(x_test)!=nrow(y_test) ) stop("'x_test' and 'y_test' have different number of data points.", call. = FALSE)
  return(list(x_test, y_test))
}

#check argument train_N
check_train_N <- function(train_N, N){
  if ( length(train_N)==1 ){
    train_N <- as.integer(rep(train_N, N))
  } else {
    if ( length(train_N)!=N ) stop("The length of 'train_N' must be 'N'.", call. = FALSE)
    train_N <- as.integer(train_N)
  }
  return(train_N)
}

#check argument target
check_target <- function(target, n_dim_Y){
  if ( !is.null(target) ){
    if ( length(target)==1 ) {
      target <- rep(target, n_dim_Y)
    } else {
      if ( length(target)!=n_dim_Y ) {
        stop(sprintf("length(target) should equal to %i.", n_dim_Y), call. = FALSE)
      }
    }
  }
  return(target)
}

#check argument limits
check_limits <- function(limits, n_dim_X){
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
  return(limits)
}

check_int <- function(int, n_dim_X){
  if ( length(int)==1 ) {
    if ( n_dim_X!=1 ) {
      int <- rep(int, n_dim_X)
    }
  } else {
    if ( length(int)!=n_dim_X ) stop("The length of 'int' should equal to the number of input dimensions.", call. = FALSE)
  }
  return(int)
}

#check argument object and determine if it is the first time being used by design()
check_design <- function(object, x_test, y_test, eval){
  if ( is.null(eval) ){
    if ( is.null(x_test) & is.null(y_test) ){
      if ( !identical(object$design$type, "loo") ) {
        first_time <- TRUE
        n_wave <- 0
      } else {
        first_time <- FALSE
        n_wave <- length(object$design)
        if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
        if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
      }
    } else {
      if ( !identical(object$design$type, "oos") ) {
        first_time <- TRUE
        n_wave <- 0
      } else {
        if ( identical(object$design$x_test, x_test) & identical(object$design$y_test, y_test) ){
          first_time <- FALSE
          n_wave <- length(object$design)
          if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
          if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
        } else {
          first_time <- TRUE
          n_wave <- 0
        }
      }
    }
  } else {
    if ( identical(object$design$type, eval) ) {
      first_time <- FALSE
      n_wave <- length(object$design)
      if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
    } else {
      first_time <- TRUE
      n_wave <- 0
    }
  }
  return(list(first_time, n_wave))
}
