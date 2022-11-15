#' @title Sequential design of a (D)GP emulator or a bundle of (D)GP emulators
#'
#' @description This function implements the sequential design of a (D)GP emulator or a bundle of (D)GP emulators.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param N the number of steps for the sequential design, i.e., the number of design points to be added to the emulator `object`.
#' @param x_cand a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
#'     in which the next design point is determined. If `x_cand = NULL`, the candidate set will be generated using `n_cand` and
#'     `limits`. Defaults to `NULL`.
#' @param y_cand a matrix (with each row being a simulator evaluation and column being an output dimension) that gives the realizations
#'    from the simulator at input positions in `x_cand`. Defaults to `NULL`.
#' @param n_cand an integer that gives
#' * the size of the candidate set in which the next design point is determined, if `x_cand = NULL`;
#' * the size of a sub-set to be sampled from the candidate set `x_cand` at each step of the sequential design to determine the next
#'   design point, if `x_cand` is not `NULL`.
#'
#' Defaults to `200`
#' @param limits a two-column matrix that gives the ranges of each input dimension, or a vector of length two if there is only one
#'     input dimension. If a vector is provided, it will be converted to a two-column row matrix. The rows of the matrix correspond to input
#'     dimensions, and its first and second columns correspond to the minimum and maximum values of the input dimensions. If
#'     `limits = NULL`, the ranges of input dimensions will be determined from the training data contained in `object`. This argument
#'     is used when `x_cand = NULL` and `y_cand = NULL`. Defaults to `NULL`.
#' @param batch_size an integer or a vector of integers that gives the number of design points (for each simulator output dimension if
#'     `aggregate = NULL`) to be chosen at each step of the sequential design. If an integer is given, it will be applied to all steps
#'     of the sequential design. If a vector is given, it must have a length of `N`. Defaults to `1`.
#' @param f an R function that represents the simulator. See *Note* section below for details on specifications of `f`.
#' @param freq a vector of two integers with the first element giving the frequency (in number of steps) to re-fit the
#'     emulator, and the second element giving the frequency to implement the emulator validation (for RMSE). Defaults to `c(1, 1)`.
#' @param x_test a matrix (with each row being an input testing data point and each column being an input dimension) that gives the testing
#'     input data to evaluate the emulator after each step of the sequential design. Set to `NULL` for the LOO-based emulator validation.
#'     Defaults to `NULL`.
#' @param y_test the testing output data that correspond to `x_test` for the emulator validation after each step of the sequential design:
#' * if `object` is an instance of the `gp` class, `y_test` is a matrix with only one column and each row being an testing output data point.
#' * if `object` is an instance of the `dgp` class, `y_test` is a matrix with its rows being testing output data points and columns being
#'   output dimensions.
#'
#' Set to `NULL` for the LOO-based emulator validation. Defaults to `NULL`.
#' @param method the criterion used to locate the design points: ALM (`"ALM"`) or MICE (`"MICE"`). See references below.
#'     Defaults to `"ALM"`.
#' @param nugget_s the value of the smoothing nugget term used when `method = "MICE"`. Defaults to `1.0`.
#' @param verb a bool indicating if the trace information will be printed during the sequential design.
#'     Defaults to `TRUE`.
#' @param check_point a vector of integers that indicates at which steps the sequential design will pause and ask for the confirmation
#'     from the user if the sequential design should continue or be terminated. Set to `NULL` to suspend the manual intervention. Defaults
#'     to `NULL`.
#' @param cores a vector of two integers that gives the number of cores/workers to be used for criterion calculations and emulator
#'     validations, respectively. If either element is set to `NULL`, the number of cores is set to `(max physical cores available - 1)`.
#'     Defaults to `c(1, 1)`.
#' @param threading a bool indicating whether to use the multi-threading to accelerate criterion calculations and emulator
#'     validations if `object` is an instance of the `dgp` class. Turning this option on could improve the speed of validations when the DGP emulator
#'     is built with a moderately large number of training data points and the Matérn-2.5 kernel.
#' @param train_N an integer or a vector of integers that gives the number of training iterations to be used to re-fit the DGP emulator at each step
#'     of the sequential design:
#' * If `train_N` is an integer, then at each step the DGP emulator will re-fitted (based on the frequency of re-fit specified in `freq`) with `train_N` iterations.
#' * If `train_N` is a vector, then its size must be `N` even the re-fit frequency specified in `freq` is not one.
#'
#' Defaults to `100`.
#' @param aggregate an R function that aggregates the ALM or MICE criterion across different output dimensions of a DGP emulator, or across outputs of different
#'     emulators contained in an emulator bundle. See *Note* section below for details on specifications of `aggregate`.
#' @param ... N/A.
#'
#' @return
#' An updated `object` is returned with a slot called `design` that contains:
#'   - *S* slots, named `wave1, wave2,..., waveS`, that contain information of *S* waves of sequential designs that have been applied to the emulator.
#'     Each slot contains the following elements:
#'     - `N`, an integer that gives the numbers of steps implemented in the corresponding wave;
#'     - `rmse`, a matrix that gives the RMSEs of emulators constructed during the corresponding wave;
#'     - `freq`, an integer that gives the frequency that the emulator validations are implemented during the corresponding wave.
#'     - `aggregation`, a bool that indicates if the aggregation is applied (i.e., if `aggregate` is supplied) to locate the design points for a DGP
#'       emulator or a bundle of (D)GP emulators during the corresponding wave.
#'     - `batch_size`, a vector of size `N` that contains the information given by the argument `batch_size`.
#'     - `enrichment`, a vector of size `N` that gives the number of new design points added after each step of the sequential design.
#'   - a slot called `type` that gives the type of validations, either LOO (`loo`) or OOS (`oos`) used to calculate the RMSEs of emulators constructed
#'     in the sequential design. See [validate()] for more information about LOO and OOS.
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
#'   - the validation type (LOO or OOS) of the current wave of the sequential design is the same as the validation types in previous waves, and
#'   - if the validation type is OOS, `x_test` and `y_test` in the current wave of the sequential design are identical to those in the previous waves.
#'
#'   When the above conditions are met, the information (`N` and `rmse`) of the current wave of the sequential design will be added to
#'       the `design` slot of the returned object under the name `waveS`.
#' * If `object` is an instance of the `gp` class, the matrix in the `rmse` slot is single-columned. If `object` is an instance of
#'   the `dgp` class, the matrix in the `rmse` slot can have multiple columns that correspond to different output dimensions.
#' * `aggregate` can be defined as an R function with either one or two arguments.
#'   - if `aggregate` has one argument, the argument should be a vector representing values of ALM or MICE criterion across
#'     - different output dimensions of a DGP emulator; or
#'     - outputs of different emulators contained an emulator bundle,
#'     at a design point.
#'   - if `aggregate` has two arguments, the additional argument should also be a vector that contains parameters that are involved in the criterion
#'     aggregation. The values of the aggregation parameters should be specified by setting them as the default values of the second argument of
#'     `aggregate`.
#' * `f` should be defined as a function with only one argument which is a matrix with each row giving values of an input position to the simulator.
#'   `f` can return a single or two outputs.
#'   - if `f` has one output, it should return a matrix with each row being the output from the simulator given the corresponding input.
#'   - if `f` has two outputs, it should return a list with two slots:
#'     - the first slot is a matrix that gives values of outputs generated from the simulator.
#'     - the second slot is a vector that gives values of aggregation parameters that will be used to update the criterion aggregation via `aggregate`
#'       after each step of the sequential design.
#'
#'     In the two-output case, `aggregate` must be provided and defined with two arguments.
#' * When defining `f` and `aggregate`, it is important to ensure that:
#'   - the positions of input dimensions (i.e, the order of columns of the input matrix) of `f` are consistent with those of the emulator;
#'   - the positions of output dimensions (i.e, the order of columns of the output matrix) of `f` are consistent with those of the emulator, or the order
#'     of emulators placed in `object` if `object` is an instance of the `bundle` class;
#'   - the positions of elements in the first argument of `aggregate` are consistent with the positions of output dimensions of the emulator, or the
#'     order of emulators placed in `object` if `object` is an instance of the `bundle` class;
#'   - the positions of elements in the second output (if defined) of `f` are consistent with the positions of elements in the second argument of `aggregate`.
#' * Any R vector detected in `x_test` and `y_test` will be treated as a column vector and automatically converted into a single-column
#'   R matrix. Thus, if `x_test` or `y_test` is a single testing data point with multiple dimensions, it must be given as a matrix.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @references
#' MacKay, D. J. (1992). Information-based objective functions for active data selection. *Neural Computation*, **4(4)**, 590-604.
#'
#' Beck, J., & Guillas, S. (2016). Sequential design with mutual information for computer experiments (MICE): emulation of a tsunami model.
#' *SIAM/ASA Journal on Uncertainty Quantification*, **4(1)**, 739-766.
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
design <- function(object, N, x_cand, y_cand, n_cand, limits, batch_size, f, freq, x_test, y_test, method, nugget_s, verb, check_point, cores, ...){
  UseMethod("design")
}

#' @rdname design
#' @method design gp
#' @export
design.gp <- function(object, N, x_cand = NULL, y_cand = NULL, n_cand = 200, limits = NULL, batch_size = 1, f = NULL, freq = c(1, 1), x_test = NULL, y_test = NULL, method = 'ALM', nugget_s = 1., verb = TRUE, check_point = NULL, cores = c(1, 1), ...) {
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  #if candidate set is given
  N <- as.integer(N)
  if ( N < 1 ) stop("The number of steps 'N' for a sequential design must be >= 1.", call. = FALSE)

  freq <- as.integer(freq)
  if ( any(freq < 1) ) stop("All elements in 'freq' must be greater than or equal to 1.", call. = FALSE)

  n_cand <- as.integer(n_cand)
  if ( n_cand < 1 ) stop("'n_cand' must be greater than or equal to 1.", call. = FALSE)

  if ( length(batch_size)==1 ) {
    batch_size <- rep(batch_size, N)
  } else {
    if ( length(batch_size)!=N ) {
      stop(sprintf("length(batch_size) should equal to %i.", N), call. = FALSE)
    }
  }
  batch_size <- as.integer(batch_size)

  if (!is.null(x_test) & !is.null(y_test)) {
    x_test <- unname(x_test)
    y_test <- unname(y_test)
    if ( !is.matrix(x_test)&!is.vector(x_test) ) stop("'x_test' must be a vector or a matrix.", call. = FALSE)
    if ( !is.matrix(y_test)&!is.vector(y_test) ) stop("'y_test' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_test) ) x_test <- as.matrix(x_test)
    if ( is.vector(y_test) ) y_test <- as.matrix(y_test)
    if ( nrow(x_test)!=nrow(y_test) ) stop("'x_test' and 'y_test' have different number of data points.", call. = FALSE)
  }

  if ( "design" %in% names(object) ){
    design_info <- object$design
    if ( is.null(x_test) & is.null(y_test) ){
      if ( object$design$type=="oos" ) {
        first_time <- TRUE
      } else {
        first_time <- FALSE
        n_wave <- length(object$design)
        if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
        if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
      }
    } else {
      if ( object$design$type=="loo" ) {
        first_time <- TRUE
      } else {
        if ( identical(object$design$x_test, x_test) & identical(object$design$y_test, y_test) ){
          first_time <- FALSE
          n_wave <- length(object$design)
          if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
          if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
        } else {
          first_time <- TRUE
        }
      }
    }
  } else {
    first_time <- TRUE
  }

  N_acq <- c()

  if ( is.null(x_cand) ) {
    #if candidate set is not given
    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)

    X <- object$data$X
    Y <- object$data$Y

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if (is.null(x_test) & is.null(y_test)){
      type <- 'loo'
      object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores[2])
      rmse <- object$loo$rmse
    } else {
      type <- 'oos'
      object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, cores = cores[2])
      rmse <- object$oos$rmse
    }
    if ( verb ) message(" done")
    message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    rmse_records <- rmse

    for ( i in 1:N ){
      if ( verb ) {
        message(sprintf("Iteration %i:", i))
        message(" - Locating ...", appendLF = FALSE)
      }
      res <- locate(object, n_cand = n_cand, limits = limits, batch_size = batch_size[i], method = method, nugget_s = nugget_s, verb = FALSE, cores = cores[1])
      if ( verb ) {
        message(" done")
        if ( batch_size[i]==1 ){
          message(paste(c(" * Next design point:", sprintf("%.06f", res$location)), collapse=" "))
        } else {
          for (j in 1:batch_size[i] ){
            message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", res[['location']][[paste('position',j,sep="")]])), collapse=" "))
          }
        }
      }

      if ( batch_size[i]==1 ){
        new_X <- res$location
      } else {
        new_X <- c()
        for (j in 1:batch_size[i] ){
          new_X <- rbind(new_X, res[['location']][[paste('position',j,sep="")]])
        }
      }

      N_acq <- c(N_acq, nrow(new_X))
      X <- rbind(X, new_X)
      Y <- rbind(Y, f(new_X))

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
        if (is.null(x_test) & is.null(y_test)){
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores[2])
          if ( verb ) message(" done")
          rmse <- object$loo$rmse
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        } else {
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, cores = cores[2])
          if ( verb ) message(" done")
          rmse <- object$oos$rmse
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        }
        rmse_records <- rbind(rmse_records, rmse)
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
  } else {
    if ( is.null(y_cand) ) stop("'y_cand' must be provided if 'x_cand' is not NULL.", call. = FALSE)
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( !is.matrix(y_cand)&!is.vector(y_cand) ) stop("'y_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( is.vector(y_cand) ) y_cand <- as.matrix(y_cand)
    if ( nrow(x_cand)!=nrow(y_cand) ) stop("'x_cand' and 'y_cand' have different number of data points.", call. = FALSE)

    X <- object$data$X
    Y <- object$data$Y

    if ( ncol(y_cand)!=ncol(Y) ) stop("The dimension of 'y_cand' must be one.", call. = FALSE)

    X_s <- apply(X, 1, paste, collapse = ", ")
    x_cand_s <- apply(x_cand, 1, paste, collapse = ", ")
    X_x_cand <- intersect(X_s, x_cand_s)
    if ( length(X_x_cand)!=0 ){
      x_cand <- x_cand[!x_cand_s %in% X_x_cand,,drop=FALSE]
    }
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- c()

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if (is.null(x_test) & is.null(y_test)){
      type <- 'loo'
      object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores[2])
      rmse <- object$loo$rmse
    } else {
      type <- 'oos'
      object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, cores = cores[2])
      rmse <- object$oos$rmse
    }
    if ( verb ) message(" done")
    message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    rmse_records <- rmse

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
      res <- locate(object, x_cand = sub_cand, batch_size = batch_size[i], method = method, nugget_s = nugget_s, verb = FALSE, cores = cores[1])
      if ( verb ) {
        message(" done")
        if ( batch_size[i]==1 ){
          message(paste(c(" * Next design point:", sprintf("%.06f", res$location)), collapse=" "))
        } else {
          for (j in 1:batch_size[i] ){
            message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", res[['location']][[paste('position',j,sep="")]])), collapse=" "))
          }
        }
      }

      if ( batch_size[i]==1 ){
        new_X <- res$location
        new_idx <- res$index
      } else {
        new_X <- c()
        new_idx <- c()
        for (j in 1:batch_size[i] ){
          new_X <- rbind(new_X, res[['location']][[paste('position',j,sep="")]])
          new_idx <- c(new_idx, res[['index']][[paste('position',j,sep="")]])
        }
      }

      N_acq <- c(N_acq, nrow(new_X))

      X <- rbind(X, new_X)
      Y <- rbind(Y, y_cand[idx_sub_cand,,drop=FALSE][new_idx,,drop=FALSE])
      idx_x_acq <- c(idx_x_acq, idx_sub_cand[new_idx])
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
        if (is.null(x_test) & is.null(y_test)){
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores[2])
          if ( verb ) message(" done")
          rmse <- object$loo$rmse
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        } else {
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, cores = cores[2])
          if ( verb ) message(" done")
          rmse <- object$oos$rmse
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        }
        rmse_records <- rbind(rmse_records, rmse)
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

  if ( isTRUE(first_time) ){
    object$design <- list()
    object[["design"]][["wave1"]][["N"]] <- N
    object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
    object[["design"]][["wave1"]][["freq"]] <- freq[2]
    object[["design"]][["wave1"]][["batch_size"]] <- batch_size
    object[["design"]][["wave1"]][["enrichement"]] <- N_acq
    object[["design"]][["type"]] <- type
    if ( type == 'oos' ){
      object[["design"]][["x_test"]] <- x_test
      object[["design"]][["y_test"]] <- y_test
    }
  } else {
    object$design <- design_info
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["batch_size"]] <- batch_size
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichement"]] <-  N_acq
  }

  return(object)
}

#' @rdname design
#' @method design dgp
#' @export
design.dgp <- function(object, N, x_cand = NULL, y_cand = NULL, n_cand = 200, limits = NULL, batch_size = 1, f = NULL, freq = c(1, 1), x_test = NULL, y_test = NULL, method = 'ALM', nugget_s = 1., verb = TRUE, check_point = NULL, cores = c(1, 1), threading = FALSE, train_N = 100, aggregate = NULL, ...) {
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)

  N <- as.integer(N)
  if ( N < 1 ) stop("The number of steps 'N' for a sequential design must be >= 1.", call. = FALSE)

  freq <- as.integer(freq)
  if ( any(freq < 1) ) stop("All elements in 'freq' must be greater than or equal to 1.", call. = FALSE)

  if ( length(batch_size)==1 ) {
    batch_size <- rep(batch_size, N)
  } else {
    if ( length(batch_size)!=N ) {
      stop(sprintf("length(batch_size) should equal to %i.", N), call. = FALSE)
    }
  }
  batch_size <- as.integer(batch_size)

  if ( length(train_N)==1 ){
    train_N <- as.integer(rep(train_N, N))
  } else {
    if ( length(train_N)!=N ) stop("The length of 'train_N' must be 'N'.", call. = FALSE)
    train_N <- as.integer(train_N)
  }

  n_cand <- as.integer(n_cand)
  if ( n_cand < 1 ) stop("'n_cand' must be greater than or equal to 1.", call. = FALSE)

  if (!is.null(x_test) & !is.null(y_test)) {
    x_test <- unname(x_test)
    y_test <- unname(y_test)
    if ( !is.matrix(x_test)&!is.vector(x_test) ) stop("'x_test' must be a vector or a matrix.", call. = FALSE)
    if ( !is.matrix(y_test)&!is.vector(y_test) ) stop("'y_test' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_test) ) x_test <- as.matrix(x_test)
    if ( is.vector(y_test) ) y_test <- as.matrix(y_test)
    if ( nrow(x_test)!=nrow(y_test) ) stop("'x_test' and 'y_test' have different number of data points.", call. = FALSE)
  }

  if ( "design" %in% names(object) ){
    design_info <- object$design
    if ( is.null(x_test) & is.null(y_test) ){
      if ( object$design$type=="oos" ) {
        first_time <- TRUE
      } else {
        first_time <- FALSE
        n_wave <- length(object$design)
        if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
        if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
      }
    } else {
      if ( object$design$type=="loo" ) {
        first_time <- TRUE
      } else {
        if ( identical(object$design$x_test, x_test) & identical(object$design$y_test, y_test) ){
          first_time <- FALSE
          n_wave <- length(object$design)
          if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
          if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
        } else {
          first_time <- TRUE
        }
      }
    }
  } else {
    first_time <- TRUE
  }

  N_acq <- c()

  #if candidate set is given
  if ( is.null(x_cand) ) {

    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)

    if ( !is.null(aggregate) ) {
      g <- function(x){
        res <- aggregate(x)
        return(res)
      }
    } else {
      g <- NULL
    }

    X <- object$data$X
    Y <- object$data$Y

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if (is.null(x_test) & is.null(y_test)){
      type <- 'loo'
      object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores[2], threading = threading)
      rmse <- object$loo$rmse
    } else {
      type <- 'oos'
      object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, cores = cores[2], threading = threading)
      rmse <- object$oos$rmse
    }
    if ( verb ) message(" done")
    message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    rmse_records <- rmse

    for ( i in 1:N ){
      if ( verb ) {
        message(sprintf("Iteration %i:", i))
        message(" - Locating ...", appendLF = FALSE)
      }
      res <- locate(object, n_cand = n_cand, limits = limits, batch_size = batch_size[i], method = method, nugget_s = nugget_s, verb = FALSE, cores = cores[1], threading = threading, aggregate = g)
      if ( verb ) {
        message(" done")
        if ( batch_size[i]==1 ){
          if ( nrow(res$location)==1 ){
            message(paste(c(" * Next design point:", sprintf("%.06f", res$location[1,])), collapse=" "))
          } else {
            for ( j in 1:nrow(res$location) ){
              message(paste(c(sprintf(" * Next design point (Output%i):", j), sprintf("%.06f", res$location[j,])), collapse=" "))
            }
          }
        } else {
          for ( j in 1:batch_size[i] ){
            pos_j <- res[['location']][[paste('position',j,sep="")]]
            n_row <- nrow(pos_j)
            if ( n_row==1 ){
              message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", pos_j[1,])), collapse=" "))
            } else {
              for ( k in 1:n_row ){
                message(paste(c(sprintf(" * Next design point (Position%i for Output%i):", j, k), sprintf("%.06f", pos_j[k,])), collapse=" "))
              }
            }
          }
        }
      }

      if ( batch_size[i]==1 ){
        new_X <- res$location
      } else {
        new_X <- c()
        for (j in 1:batch_size[i] ){
          new_X <- rbind(new_X, res[['location']][[paste('position',j,sep="")]])
        }
        new_X <- unique(new_X, MARGIN = 1)
      }

      N_acq <- c(N_acq, nrow(new_X))

      X <- rbind(X, new_X)
      new_output <- f(new_X)

      if ( is.list(new_output) ){
        Y <- rbind(Y, new_output[[1]])
        coeff <- new_output[[2]]
        g <- function(x){
          res <- aggregate(x,coeff)
          return(res)
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
        if (is.null(x_test) & is.null(y_test)){
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores[2], threading = threading)
          if ( verb ) message(" done")
          rmse <- object$loo$rmse
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        } else {
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, cores = cores[2], threading = threading)
          if ( verb ) message(" done")
          rmse <- object$oos$rmse
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        }
        rmse_records <- rbind(rmse_records, rmse)
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
  } else {
    if ( is.null(y_cand) ) stop("'y_cand' must be provided if 'x_cand' is not NULL.", call. = FALSE)
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( !is.matrix(y_cand)&!is.vector(y_cand) ) stop("'y_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( is.vector(y_cand) ) y_cand <- as.matrix(y_cand)
    if ( nrow(x_cand)!=nrow(y_cand) ) stop("'x_cand' and 'y_cand' have different number of data points.", call. = FALSE)
    if ( n_cand>nrow(x_cand) ) stop("'n_cand' must be smaller than or equal to the size of the candidate set 'x_cand'.", call. = FALSE)

    if ( !is.null(aggregate) ) {
      g <- function(x){
        res <- aggregate(x)
        return(res)
      }
    } else {
      g <- NULL
    }

    X <- object$data$X
    Y <- object$data$Y

    if ( ncol(y_cand)!=ncol(Y) ) stop("The dimension of 'y_cand' must equal to that of the training output data.", call. = FALSE)

    X_s <- apply(X, 1, paste, collapse = ", ")
    x_cand_s <- apply(x_cand, 1, paste, collapse = ", ")
    X_x_cand <- intersect(X_s, x_cand_s)
    if ( length(X_x_cand)!=0 ){
      x_cand <- x_cand[!x_cand_s %in% X_x_cand,,drop=FALSE]
    }
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- c()

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    if (is.null(x_test) & is.null(y_test)){
      type <- 'loo'
      object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores[2], threading = threading)
      rmse <- object$loo$rmse
    } else {
      type <- 'oos'
      object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, cores = cores[2], threading = threading)
      rmse <- object$oos$rmse
    }
    if ( verb ) message(" done")
    message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    rmse_records <- rmse

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
      res <- locate(object, x_cand = sub_cand, batch_size = batch_size[i], method = method, nugget_s = nugget_s, verb = FALSE, cores = cores[1], threading = threading, aggregate = g)
      if ( verb ) {
        message(" done")
        if ( batch_size[i]==1 ){
          if ( nrow(res$location)==1 ){
            message(paste(c(" * Next design point:", sprintf("%.06f", res$location[1,])), collapse=" "))
          } else {
            for ( j in 1:nrow(res$location) ){
              message(paste(c(sprintf(" * Next design point (Output%i):", j), sprintf("%.06f", res$location[j,])), collapse=" "))
            }
          }
        } else {
          for ( j in 1:batch_size[i] ){
            pos_j <- res[['location']][[paste('position',j,sep="")]]
            n_row <- nrow(pos_j)
            if ( n_row==1 ){
              message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", pos_j[1,])), collapse=" "))
            } else {
              for ( k in 1:n_row ){
                message(paste(c(sprintf(" * Next design point (Position%i for Output%i):", j, k), sprintf("%.06f", pos_j[k,])), collapse=" "))
              }
            }
          }
        }
      }

      if ( batch_size[i]==1 ){
        new_X <- res$location
        new_idx <- res$index
      } else {
        new_X <- c()
        new_idx <- c()
        for (j in 1:batch_size[i] ){
          new_X <- rbind(new_X, res[['location']][[paste('position',j,sep="")]])
          new_idx <- c(new_idx, res[['index']][[paste('position',j,sep="")]])
        }
        new_X <- unique(new_X, MARGIN = 1)
        new_idx <- unique(new_idx)
      }

      N_acq <- c(N_acq, nrow(new_X))

      X <- rbind(X, new_X)
      Y <- rbind(Y, y_cand[idx_sub_cand,,drop=FALSE][new_idx,,drop=FALSE])
      idx_x_acq <- c(idx_x_acq, idx_sub_cand[new_idx])
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
        if (is.null(x_test) & is.null(y_test)){
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores[2], threading = threading)
          if ( verb ) message(" done")
          rmse <- object$loo$rmse
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        } else {
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, cores = cores[2], threading = threading)
          if ( verb ) message(" done")
          rmse <- object$oos$rmse
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        }
        rmse_records <- rbind(rmse_records, rmse)
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

  if ( isTRUE(first_time) ){
    object$design <- list()
    object[["design"]][["wave1"]][["N"]] <- N
    object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
    object[["design"]][["wave1"]][["freq"]] <- freq[2]
    object[["design"]][["wave1"]][["batch_size"]] <- batch_size
    object[["design"]][["wave1"]][["enrichment"]] <- N_acq
    if ( is.null(aggregate) ){
      object[["design"]][["wave1"]][["aggregation"]] <- FALSE
    } else {
      object[["design"]][["wave1"]][["aggregation"]] <- TRUE
    }
    object[["design"]][["type"]] <- type
    if ( type == 'oos' ){
      object[["design"]][["x_test"]] <- x_test
      object[["design"]][["y_test"]] <- y_test
    }
  } else {
    object$design <- design_info
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["batch_size"]] <- batch_size
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <- N_acq
    if ( is.null(aggregate) ){
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["aggregation"]] <- FALSE
    } else {
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["aggregation"]] <- TRUE
    }
  }

  return(object)
}


#' @rdname design
#' @method design bundle
#' @export
design.bundle <- function(object, N, x_cand = NULL, y_cand = NULL, n_cand = 200, limits = NULL, batch_size = 1, f = NULL, freq = c(1, 1), x_test = NULL, y_test = NULL, method = 'ALM', nugget_s = 1., verb = TRUE, check_point = NULL, cores = c(1, 1), threading = FALSE, train_N = 100, aggregate = NULL, ...) {
  if ( !inherits(object,"bundle") ) stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)

  N <- as.integer(N)
  if ( N < 1 ) stop("The number of steps 'N' for a sequential design must be >= 1.", call. = FALSE)

  freq <- as.integer(freq)
  if ( any(freq < 1) ) stop("All elements in 'freq' must be greater than or equal to 1.", call. = FALSE)

  if ( length(batch_size)==1 ) {
    batch_size <- rep(batch_size, N)
  } else {
    if ( length(batch_size)!=N ) {
      stop(sprintf("length(batch_size) should equal to %i.", N), call. = FALSE)
    }
  }
  batch_size <- as.integer(batch_size)

  if ( length(train_N)==1 ){
    train_N <- as.integer(rep(train_N, N))
  } else {
    if ( length(train_N)!=N ) stop("The length of 'train_N' must be 'N'.", call. = FALSE)
    train_N <- as.integer(train_N)
  }

  n_cand <- as.integer(n_cand)
  if ( n_cand < 1 ) stop("'n_cand' must be greater than or equal to 1.", call. = FALSE)

  if (!is.null(x_test) & !is.null(y_test)) {
    x_test <- unname(x_test)
    y_test <- unname(y_test)
    if ( !is.matrix(x_test)&!is.vector(x_test) ) stop("'x_test' must be a vector or a matrix.", call. = FALSE)
    if ( !is.matrix(y_test)&!is.vector(y_test) ) stop("'y_test' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_test) ) x_test <- as.matrix(x_test)
    if ( is.vector(y_test) ) y_test <- as.matrix(y_test)
    if ( nrow(x_test)!=nrow(y_test) ) stop("'x_test' and 'y_test' have different number of data points.", call. = FALSE)
  }

  if ( "design" %in% names(object) ){
    design_info <- object$design
    if ( is.null(x_test) & is.null(y_test) ){
      if ( object$design$type=="oos" ) {
        first_time <- TRUE
      } else {
        first_time <- FALSE
        n_wave <- length(object$design)
        if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
        if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
      }
    } else {
      if ( object$design$type=="loo" ) {
        first_time <- TRUE
      } else {
        if ( identical(object$design$x_test, x_test) & identical(object$design$y_test, y_test) ){
          first_time <- FALSE
          n_wave <- length(object$design)
          if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
          if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
        } else {
          first_time <- TRUE
        }
      }
    }
  } else {
    first_time <- TRUE
  }

  N_acq <- c()

  n_emulators <- length(object) - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1

  #if candidate set is given
  if ( is.null(x_cand) ) {

    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)

    if ( !is.null(aggregate) ) {
      g <- function(x){
        res <- aggregate(x)
        return(res)
      }
    } else {
      g <- NULL
    }

    X <- object$data$X
    Y <- object$data$Y

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    rmse <- c()
    for ( k in 1:n_emulators ){
      obj_k <- object[[paste('emulator',k,sep='')]]
      if (is.null(x_test) & is.null(y_test)){
        type <- 'loo'
        if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE)
        if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores[2], threading = threading)
        object[[paste('emulator',k,sep='')]] <- obj_k
        rmse <- c(rmse, obj_k$loo$rmse)
      } else {
        type <- 'oos'
        if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE)
        if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, cores = cores[2], threading = threading)
        object[[paste('emulator',k,sep='')]] <- obj_k
        rmse <- c(rmse, obj_k$oos$rmse)
      }
    }
    if ( verb ) message(" done")
    message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    rmse_records <- rmse

    for ( i in 1:N ){
      if ( verb ) {
        message(sprintf("Iteration %i:", i))
        message(" - Locating ...", appendLF = FALSE)
      }
      res <- locate(object, n_cand = n_cand, limits = limits, batch_size = batch_size[i], method = method, nugget_s = nugget_s, verb = FALSE, cores = cores[1], threading = threading, aggregate = g)
      if ( verb ) {
        message(" done")
        if ( batch_size[i]==1 ){
          if ( nrow(res$location)==1 ){
            message(paste(c(" * Next design point:", sprintf("%.06f", res$location[1,])), collapse=" "))
          } else {
            for ( j in 1:nrow(res$location) ){
              message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", res$location[j,])), collapse=" "))
            }
          }
        } else {
          for ( j in 1:batch_size[i] ){
            pos_j <- res[['location']][[paste('position',j,sep="")]]
            n_row <- nrow(pos_j)
            if ( n_row==1 ){
              message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", pos_j[1,])), collapse=" "))
            } else {
              for ( k in 1:n_row ){
                message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", pos_j[k,])), collapse=" "))
              }
            }
          }
        }
      }

      if ( batch_size[i]==1 ){
        new_X <- res$location
      } else {
        new_X <- c()
        for (j in 1:batch_size[i] ){
          new_X <- rbind(new_X, res[['location']][[paste('position',j,sep="")]])
        }
        new_X <- unique(new_X, MARGIN = 1)
      }

      N_acq <- c(N_acq, nrow(new_X))

      X <- rbind(X, new_X)
      new_output <- f(new_X)

      if ( is.list(new_output) ){
        Y <- rbind(Y, new_output[[1]])
        coeff <- new_output[[2]]
        g <- function(x){
          res <- aggregate(x,coeff)
          return(res)
        }
      } else {
        Y <- rbind(Y, new_output)
      }

      if ( i %% freq[1]==0 | i==N){
        if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
        for ( k in 1:n_emulators ){
          obj_k <- object[[paste('emulator',k,sep='')]]
          if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X, Y[,k], refit = TRUE, verb = FALSE)
          if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X, Y[,k], refit = TRUE, verb = FALSE, N = train_N[i], B = 10)
          object[[paste('emulator',k,sep='')]] <- obj_k
        }
        if ( verb ) message(" done")
      } else {
        if ( verb ) message(" - Updating ...", appendLF = FALSE)
        for ( k in 1:n_emulators ){
          obj_k <- object[[paste('emulator',k,sep='')]]
          if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, verb = FALSE)
          if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, verb = FALSE, B = 10)
          object[[paste('emulator',k,sep='')]] <- obj_k
        }
        if ( verb ) message(" done")
      }

      if ( i %% freq[2]==0 | i==N){
        rmse <- c()
        if (is.null(x_test) & is.null(y_test)){
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            obj_k <- object[[paste('emulator',k,sep='')]]
            if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE)
            if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores[2], threading = threading)
            object[[paste('emulator',k,sep='')]] <- obj_k
            rmse <- c(rmse, obj_k$loo$rmse)
          }
          if ( verb ) message(" done")
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        } else {
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            obj_k <- object[[paste('emulator',k,sep='')]]
            if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE)
            if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE, cores = cores[2], threading = threading)
            object[[paste('emulator',k,sep='')]] <- obj_k
            rmse <- c(rmse, obj_k$oos$rmse)
          }
          if ( verb ) message(" done")
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        }
        rmse_records <- rbind(rmse_records, rmse)
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
  } else {
    if ( is.null(y_cand) ) stop("'y_cand' must be provided if 'x_cand' is not NULL.", call. = FALSE)
    if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
    if ( !is.matrix(y_cand)&!is.vector(y_cand) ) stop("'y_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(x_cand) ) x_cand <- as.matrix(x_cand)
    if ( is.vector(y_cand) ) y_cand <- as.matrix(y_cand)
    if ( nrow(x_cand)!=nrow(y_cand) ) stop("'x_cand' and 'y_cand' have different number of data points.", call. = FALSE)
    if ( n_cand>nrow(x_cand) ) stop("'n_cand' must be smaller than or equal to the size of the candidate set 'x_cand'.", call. = FALSE)

    if ( !is.null(aggregate) ) {
      g <- function(x){
        res <- aggregate(x)
        return(res)
      }
    } else {
      g <- NULL
    }

    X <- object$data$X
    Y <- object$data$Y

    if ( ncol(y_cand)!=ncol(Y) ) stop("The dimension of 'y_cand' must equal to the number of emulators in 'object'.", call. = FALSE)

    X_s <- apply(X, 1, paste, collapse = ", ")
    x_cand_s <- apply(x_cand, 1, paste, collapse = ", ")
    X_x_cand <- intersect(X_s, x_cand_s)
    if ( length(X_x_cand)!=0 ){
      x_cand <- x_cand[!x_cand_s %in% X_x_cand,,drop=FALSE]
    }
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- c()

    if ( verb ) message("Initializing ...", appendLF = FALSE)
    rmse <- c()
    for ( k in 1:n_emulators ){
      obj_k <- object[[paste('emulator',k,sep='')]]
      if (is.null(x_test) & is.null(y_test)){
        type <- 'loo'
        if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE)
        if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, cores = cores[2], threading = threading)
        object[[paste('emulator',k,sep='')]] <- obj_k
        rmse <- c(rmse, obj_k$loo$rmse)
      } else {
        type <- 'oos'
        if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE)
        if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, cores = cores[2], threading = threading)
        object[[paste('emulator',k,sep='')]] <- obj_k
        rmse <- c(rmse, obj_k$oos$rmse)
      }
    }
    if ( verb ) message(" done")
    message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
    rmse_records <- rmse

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
      res <- locate(object, x_cand = sub_cand, method = method, batch_size = batch_size[i], nugget_s = nugget_s, verb = FALSE, cores = cores[1], threading = threading, aggregate = g)
      if ( verb ) {
        message(" done")
        if ( batch_size[i]==1 ){
          if ( nrow(res$location)==1 ){
            message(paste(c(" * Next design point:", sprintf("%.06f", res$location[1,])), collapse=" "))
          } else {
            for ( j in 1:nrow(res$location) ){
              message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", res$location[j,])), collapse=" "))
            }
          }
        } else {
          for ( j in 1:batch_size[i] ){
            pos_j <- res[['location']][[paste('position',j,sep="")]]
            n_row <- nrow(pos_j)
            if ( n_row==1 ){
              message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", pos_j[1,])), collapse=" "))
            } else {
              for ( k in 1:n_row ){
                message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", pos_j[k,])), collapse=" "))
              }
            }
          }
        }
      }

      if ( batch_size[i]==1 ){
        new_X <- res$location
        new_idx <- res$index
      } else {
        new_X <- c()
        new_idx <- c()
        for (j in 1:batch_size[i] ){
          new_X <- rbind(new_X, res[['location']][[paste('position',j,sep="")]])
          new_idx <- c(new_idx, res[['index']][[paste('position',j,sep="")]])
        }
        new_X <- unique(new_X, MARGIN = 1)
        new_idx <- unique(new_idx)
      }

      N_acq <- c(N_acq, nrow(new_X))

      X <- rbind(X, new_X)
      Y <- rbind(Y, y_cand[idx_sub_cand,,drop=FALSE][new_idx,,drop=FALSE])
      idx_x_acq <- c(idx_x_acq, idx_sub_cand[new_idx])
      idx_x_cand <- idx_x_cand0[-idx_x_acq]

      if ( i %% freq[1]==0 | i==N){
        if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
        for ( k in 1:n_emulators ){
          obj_k <- object[[paste('emulator',k,sep='')]]
          if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X, Y[,k], refit = TRUE, verb = FALSE)
          if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X, Y[,k], refit = TRUE, verb = FALSE, N = train_N[i], B = 10)
          object[[paste('emulator',k,sep='')]] <- obj_k
        }
        if ( verb ) message(" done")
      } else {
        if ( verb ) message(" - Updating ...", appendLF = FALSE)
        for ( k in 1:n_emulators ){
          obj_k <- object[[paste('emulator',k,sep='')]]
          if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, verb = FALSE)
          if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, verb = FALSE, B = 10)
          object[[paste('emulator',k,sep='')]] <- obj_k
        }
        if ( verb ) message(" done")
      }

      if ( i %% freq[2]==0 | i==N){
        rmse <- c()
        if (is.null(x_test) & is.null(y_test)){
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            obj_k <- object[[paste('emulator',k,sep='')]]
            if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE)
            if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, cores = cores[2], threading = threading)
            object[[paste('emulator',k,sep='')]] <- obj_k
            rmse <- c(rmse, obj_k$loo$rmse)
          }
          if ( verb ) message(" done")
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        } else {
          if ( verb ) message(" - Validating ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            obj_k <- object[[paste('emulator',k,sep='')]]
            if ( inherits(obj_k,"gp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE)
            if ( inherits(obj_k,"dgp") ) obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, force = TRUE, cores = cores[2], threading = threading)
            object[[paste('emulator',k,sep='')]] <- obj_k
            rmse <- c(rmse, obj_k$oos$rmse)
          }
          if ( verb ) message(" done")
          message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        }
        rmse_records <- rbind(rmse_records, rmse)
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

  if ( isTRUE(first_time) ){
    object$design <- list()
    object[["design"]][["wave1"]][["N"]] <- N
    object[["design"]][["wave1"]][["rmse"]] <- unname(rmse_records)
    object[["design"]][["wave1"]][["freq"]] <- freq[2]
    object[["design"]][["wave1"]][["batch_size"]] <- batch_size
    object[["design"]][["wave1"]][["enrichment"]] <- N_acq
    if ( is.null(aggregate) ){
      object[["design"]][["wave1"]][["aggregation"]] <- FALSE
    } else {
      object[["design"]][["wave1"]][["aggregation"]] <- TRUE
    }
    object[["design"]][["type"]] <- type
    if ( type == 'oos' ){
      object[["design"]][["x_test"]] <- x_test
      object[["design"]][["y_test"]] <- y_test
    }
    object[['data']][['X']] <- X
    object[['data']][['Y']] <- Y
  } else {
    object$design <- design_info
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["N"]] <- N
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["rmse"]] <- unname(rmse_records)
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["freq"]] <- freq[2]
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["batch_size"]] <- batch_size
    object[["design"]][[paste('wave', n_wave+1, sep="")]][["enrichment"]] <- N_acq
    if ( is.null(aggregate) ){
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["aggregation"]] <- FALSE
    } else {
      object[["design"]][[paste('wave', n_wave+1, sep="")]][["aggregation"]] <- TRUE
    }
    object[['data']][['X']] <- X
    object[['data']][['Y']] <- Y
  }

  return(object)
}


