#' @title Sequential design of a (D)GP emulator or a bundle of (D)GP emulators
#'
#' @description
#'
#' `r new_badge("updated")`
#'
#' This function implements sequential design and active learning for a (D)GP emulator or
#'     a bundle of (D)GP emulators, supporting an array of popular methods as well as user-specified approaches.
#'     It can also be used as a wrapper for Bayesian optimization methods.
#'
#' @param object can be one of the following:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param N the number of iterations for the sequential design.
#' @param x_cand a matrix (with each row being a design point and column being an input dimension) that gives a candidate set
#'     from which the next design points are determined. Defaults to `NULL`.
#' @param y_cand a matrix (with each row being a simulator evaluation and column being an output dimension) that gives the realizations
#'    from the simulator at input positions in `x_cand`. Defaults to `NULL`.
#' @param n_sample `r new_badge("new")` an integer that gives the size of a sub-set to be sampled from the candidate set `x_cand` at each step of the sequential design to determine the next
#'   design point, if `x_cand` is not `NULL`.
#'
#' Defaults to `200`.
#' @param n_cand `r lifecycle::badge("deprecated")` this argument is deprecated. Use `n_sample` instead.
#' @param limits a two-column matrix that gives the ranges of each input dimension, or a vector of length two if there is only one
#'     input dimension. If a vector is provided, it will be converted to a two-column row matrix. The rows of the matrix correspond to input
#'     dimensions, and its first and second columns correspond to the minimum and maximum values of the input dimensions. Set
#'     `limits = NULL` if `x_cand` is supplied. This argument is only used when `x_cand` is not supplied, i.e., `x_cand = NULL`. Defaults to `NULL`. If you provide
#'     a custom `method` function with an argument called `limits`, the value of `limits` will be passed to your function.
#' @param f an R function representing the simulator. `f` must adhere to the following rules:
#' - **First argument**: a matrix where rows correspond to different design points, and columns represent input dimensions.
#' - **Function output**:
#'   - a matrix where rows correspond to different outputs (matching the input design points) and columns represent output dimensions.
#'     If there is only one output dimension, the function should return a matrix with a single column.
#'   - alternatively, a list where:
#'     - the first element is the output matrix as described above.
#'     - additional named elements can optionally update values of arguments with matching names passed via `...`. This list output is
#'       useful if additional arguments to `f`, `method`, or `eval` need to be updated after each sequential design iteration.
#'
#' See the *Note* section below for additional details. This argument is required and must be supplied when `y_cand = NULL`. Defaults to `NULL`.
#' @param reps an integer that gives the number of repetitions of the located design points to be created and used for evaluations of `f`. Set the
#'     argument to an integer greater than `1` only if `f` is a stochastic function that can generate different responses given for the same input and the
#'     supplied emulator `object` can deal with stochastic responses, e.g., a (D)GP emulator with `nugget_est = TRUE` or a DGP emulator with a
#'     likelihood layer. The argument is only used when `f` is supplied. Defaults to `1`.
#' @param freq a vector of two integers with the first element indicating the number of iterations taken between re-estimating
#'     the emulator hyperparameters, and the second element defining the number of iterations to take between re-calculation of evaluating metrics
#'     on the validation set (see `x_test` below) via the `eval` function. Defaults to `c(1, 1)`.
#' @param x_test a matrix (with each row being an input testing data point and each column being an input dimension) that gives the testing
#'     input data to evaluate the emulator after each `freq[2]` iterations of the sequential design. Set to `NULL` for LOO-based emulator validation.
#'     Defaults to `NULL`. This argument is only used if `eval = NULL`.
#' @param y_test the testing output data corresponding to `x_test` for emulator validation after each `freq[2]` iterations of the sequential design:
#' * if `object` is an instance of the `gp` class, `y_test` is a matrix with only one column and each row contains a testing output data point from the corresponding row of `x_test`.
#' * if `object` is an instance of the `dgp` class, `y_test` is a matrix with its rows containing testing output data points corresponding to the same rows of `x_test` and columns representing the
#'   output dimensions.
#' * if `object` is an instance of the `bundle` class, `y_test` is a matrix with each row representing the outputs for the corresponding row of `x_test` and each column representing the output of the different emulators in the bundle.
#'
#' Set to `NULL` for LOO-based emulator validation. Defaults to `NULL`. This argument is only used if `eval = NULL`.
#' @param reset A bool or a vector of bools indicating whether to reset the hyperparameters of the emulator(s) to their initial values (as set during initial construction) before re-fitting.
#'     The re-fitting occurs based on the frequency specified by `freq[1]`. This option is useful when hyperparameters are suspected to have converged to a local optimum affecting validation performance.
#' - If a single bool is provided, it applies to every iteration of the sequential design.
#' - If a vector is provided, its length must equal `N` (even if the re-fit frequency specified in `freq[1]` is not 1) and it will apply to the corresponding iterations of the sequential design.
#'
#' Defaults to `FALSE`.
#' @param target a number or vector specifying the target evaluation metric value(s) at which the sequential design should terminate.
#'     Defaults to `NULL`, in which case the sequential design stops after `N` steps. See the *Note* section below for further details about `target`.
#' @param method `r new_badge("updated")` an R function that determines the next design points to be evaluated by `f`. The function must adhere to the following rules:
#' - **First argument**: an emulator object, which can be one of the following:
#'   - an instance of the `gp` class (produced by [gp()]);
#'   - an instance of the `dgp` class (produced by [dgp()]);
#'   - an instance of the `bundle` class (produced by [pack()]).
#' - **Second argument** (if `x_cand` is not `NULL`): a *candidate matrix* representing a set of potential design points from which the `method` function selects the next points.
#' - **Function output**:
#'   - If `x_cand` is not `NULL`:
#'     - for `gp` or `dgp` objects, the output must be a vector of row indices corresponding to the selected design points from the *candidate matrix* (the second argument).
#'     - for `bundle` objects, the output must be a matrix containing the row indices of the selected design points from the *candidate matrix*. Each column corresponds to
#'       the indices for an individual emulator in the bundle.
#'   - If `x_cand` is `NULL`:
#'     - for `gp` or `dgp` objects, the output must be a matrix where each row represents a new design point to be added.
#'     - for `bundle` objects, the output must be a list with a length equal to the number of emulators in the bundle. Each element in the list is a matrix where rows
#'       represent the new design points for the corresponding emulator.
#'
#' See [alm()], [mice()], and [vigf()] for examples of built-in `method` functions. Defaults to [vigf()].
#' @param batch_size `r new_badge("new")` an integer specifying the number of design points to select in a single iteration. Defaults to `1`.
#'     This argument is used by the built-in `method` functions [alm()], [mice()], and [vigf()].
#'     If you provide a custom `method` function with an argument named `batch_size`, the value of `batch_size` will be passed to your function.
#' @param eval an R function that computes a customized metric for evaluating emulator performance. The function must adhere to the following rules:
#' - **First argument**: an emulator object, which can be one of the following:
#'   - an instance of the `gp` class (produced by [gp()]);
#'   - an instance of the `dgp` class (produced by [dgp()]);
#'   - an instance of the `bundle` class (produced by [pack()]).
#' - **Function output**:
#'   - for `gp` objects, the output must be a single metric value.
#'   - for `dgp` objects, the output can be a single metric value or a vector of metric values with a length equal to the number of output dimensions.
#'   - for `bundle` objects, the output can be a single metric value or a vector of metric values with a length equal to the number of emulators in the bundle.
#'
#' If no custom function is provided, a built-in evaluation metric (RMSE or log-loss, in the case of DGP emulators with categorical likelihoods) will be used.
#' Defaults to `NULL`. See the *Note* section below for additional details.
#' @param verb a bool indicating if trace information will be printed during the sequential design.
#'     Defaults to `TRUE`.
#' @param autosave a list that contains configuration settings for the automatic saving of the emulator:
#' * `switch`: a bool indicating whether to enable automatic saving of the emulator during sequential design. When set to `TRUE`,
#'   the emulator in the final iteration is always saved. Defaults to `FALSE`.
#' * `directory`: a string specifying the directory path where the emulators will be stored. Emulators will be stored in a sub-directory
#'   of `directory` named 'emulator-`id`'. Defaults to './check_points'.
#' * `fname`: a string representing the base name for the saved emulator files. Defaults to 'check_point'.
#' * `save_freq`: an integer indicating the frequency of automatic saves, measured in the number of iterations. Defaults to `5`.
#' * `overwrite`: a bool value controlling the file saving behavior. When set to `TRUE`, each new automatic save overwrites the previous one,
#'   keeping only the latest version. If `FALSE`, each automatic save creates a new file, preserving all previous versions. Defaults to `FALSE`.
#' @param new_wave a bool indicating whether the current call to [design()] will create a new wave of sequential designs or add the next sequence of designs to the most recent wave.
#'     This argument is relevant only if waves already exist in the emulator. Creating new waves can improve the visualization of sequential design performance across different calls
#'     to [design()] via [draw()], and allows for specifying a different evaluation frequency in `freq`. However, disabling this option can help limit the number of waves visualized
#'     in [draw()] to avoid issues such as running out of distinct colors for large numbers of waves. Defaults to `TRUE`.
#' @param M_val `r new_badge("new")` an integer that gives the size of the conditioning set for the Vecchia approximation in emulator validations. This argument is only used if the emulator `object`
#'     was constructed under the Vecchia approximation. Defaults to `50`.
#' @param cores an integer that gives the number of processes to be used for emulator validation. If set to `NULL`, the number of processes is set to
#'     `max physical cores available %/% 2`. Defaults to `1`. This argument is only used if `eval = NULL`.
#' @param train_N the number of training iterations to be used for re-fitting the DGP emulator at each step of the sequential design:
#' - If `train_N` is an integer, the DGP emulator will be re-fitted at each step (based on the re-fit frequency specified in `freq[1]`) using `train_N` iterations.
#' - If `train_N` is a vector, its length must be `N`, even if the re-fit frequency specified in `freq[1]` is not 1.
#' - If `train_N` is `NULL`, the DGP emulator will be re-fitted at each step (based on the re-fit frequency specified in `freq[1]`) using:
#'   - `100` iterations if the DGP emulator was constructed without the Vecchia approximation, or
#'   - `50` iterations if the Vecchia approximation was used.
#'
#' Defaults to `NULL`.
#' @param refit_cores the number of processes to be used to re-fit GP components (in the same layer of a DGP emulator)
#'     at each M-step during the re-fitting. If set to `NULL`, the number of processes is set to `(max physical cores available - 1)`
#'     if the DGP emulator was constructed without the Vecchia approximation. Otherwise, the number of processes is set to `max physical cores available %/% 2`.
#'     Only use multiple processes when there is a large number of GP components in different layers and optimization of GP components
#'     is computationally expensive. Defaults to `1`.
#' @param pruning a bool indicating if dynamic pruning of DGP structures will be implemented during the sequential design after the total number of
#'     design points exceeds `min_size` in `control`. The argument is only applicable to DGP emulators (i.e., `object` is an instance of `dgp` class)
#'     produced by `dgp()`. Defaults to `TRUE`.
#' @param control a list that can supply any of the following components to control the dynamic pruning of the DGP emulator:
#' * `min_size`, the minimum number of design points required to trigger dynamic pruning. Defaults to 10 times the number of input dimensions.
#' * `threshold`, the \eqn{R^2} value above which a GP node is considered redundant. Defaults to `0.97`.
#' * `nexceed`, the minimum number of consecutive iterations that the \eqn{R^2} value of a GP node must exceed `threshold` to trigger the removal of that node from
#'   the DGP structure. Defaults to `3`.
#'
#' The argument is only used when `pruning = TRUE`.
#' @param ... Any arguments with names that differ from those used in [design()] but are required by `f`, `method`, or `eval` can be passed here.
#'     [design()] will forward relevant arguments to `f`, `method`, and `eval` based on the names of the additional arguments provided.
#'
#' @return
#' An updated `object` is returned with a slot called `design` that contains:
#'   - *S* slots, named `wave1, wave2,..., waveS`, that contain information of *S* waves of sequential design that have been applied to the emulator.
#'     Each slot contains the following elements:
#'     - `N`, an integer that gives the numbers of iterations implemented in the corresponding wave;
#'     - `rmse`, a matrix providing the evaluation metric values for emulators constructed during the corresponding wave, when `eval = NULL`.
#'        Each row of the matrix represents an iteration.
#'        - for an `object` of class `gp`, the matrix contains a single column of RMSE values.
#'        - for an `object` of class `dgp` without a categorical likelihood, each row contains mean/median squared errors corresponding to different output dimensions.
#'        - for an `object` of class `dgp` with a categorical likelihood, the matrix contains a single column of log-loss values.
#'        - for an `object` of class `bundle`, each row contains either mean/median squared errors or log-loss values for the emulators in the bundle.
#'     - `metric`: a matrix providing the values of custom evaluation metrics, as computed by the user-supplied `eval` function, for emulators constructed during the corresponding wave.
#'     - `freq`, an integer that gives the frequency that the emulator validations are implemented during the corresponding wave.
#'     - `enrichment`, a vector of size `N` that gives the number of new design points added after each step of the sequential design (if `object` is
#'        an instance of the `gp` or `dgp` class), or a matrix that gives the number of new design points added to emulators in a bundle after each step of
#'        the sequential design (if `object` is an instance of the `bundle` class).
#'
#'     If `target` is not `NULL`, the following additional elements are also included:
#'     - `target`: the target evaluating metric computed by the `eval` or built-in function to stop the sequential design.
#'     - `reached`: indicates whether the `target` was reached at the end of the sequential design:
#'        - a bool if `object` is an instance of the `gp` or `dgp` class.
#'        - a vector of bools if `object` is an instance of the `bundle` class, with its length determined as follows:
#'          - equal to the number of emulators in the bundle when `eval = NULL`.
#'          - equal to the length of the output from `eval` when a custom `eval` function is provided.
#'   - a slot called `type` that gives the type of validation:
#'     - either LOO ('loo') or OOS ('oos') if `eval = NULL`. See [validate()] for more information about LOO and OOS.
#'     - 'customized' if a customized R function is provided to `eval`.
#'   - two slots called `x_test` and `y_test` that contain the data points for the OOS validation if the `type` slot is 'oos'.
#'   - If `y_cand = NULL` and `x_cand` is supplied, and there are `NA`s returned from the supplied `f` during the sequential design, a slot called `exclusion` is included
#'     that records the located design positions that produced `NA`s via `f`. The sequential design will use this information to
#'     avoid re-visiting the same locations in later runs of `design()`.
#'
#' See *Note* section below for further information.
#' @note
#' * Validation of an emulator is forced after the final step of a sequential design even if `N` is not a multiple of the second element in `freq`.
#' * Any `loo` or `oos` slot that already exists in `object` will be cleaned, and a new slot called `loo` or `oos` will be created in the returned object
#'   depending on whether `x_test` and `y_test` are provided. The new slot gives the validation information of the emulator constructed in the final step of
#'   the sequential design. See [validate()] for more information about the slots `loo` and `oos`.
#' * If `object` has previously been used by [design()] for sequential design, the information of the current wave of the sequential design will replace
#'   those of old waves and be contained in the returned object, unless
#'   - the validation type (LOO or OOS depending on whether `x_test` and `y_test` are supplied or not) of the current wave of the sequential design is the
#'     same as the validation types (shown in the `type` of the `design` slot of `object`) in previous waves, and if the validation type is OOS,
#'     `x_test` and `y_test` in the current wave must also be identical to those in the previous waves;
#'   - both the current and previous waves of the sequential design supply customized evaluation functions to `eval`. Users need to ensure the customized evaluation
#'     functions are consistent among different waves. Otherwise, the trace plot of RMSEs produced by [draw()] will show values of different evaluation metrics in
#'     different waves.
#'
#'   For the above two cases, the information of the current wave of the sequential design will be added to the `design` slot of the returned object under the name `waveS`.
#' * If `object` is an instance of the `gp` class and `eval = NULL`, the matrix in the `rmse` slot is single-columned. If `object` is an instance of
#'   the `dgp` or `bundle` class and `eval = NULL`, the matrix in the `rmse` slot can have multiple columns that correspond to different output dimensions
#'   or different emulators in the bundle.
#' * If `object` is an instance of the `gp` class and `eval = NULL`, `target` needs to be a single value giving the RMSE threshold. If `object` is an instance
#'   of the `dgp` or `bundle` class and `eval = NULL`, `target` can be a vector of values that gives the thresholds of evaluating metrics for different output dimensions or
#'   different emulators. If a single value is provided, it will be used as the threshold for all output dimensions (if `object` is an instance of the `dgp`) or all emulators
#'   (if `object` is an instance of the `bundle`). If a customized function is supplied to `eval` and `target` is given as a vector, the user needs to ensure that the length
#'   of `target` is equal to that of the output from `eval`.
#' * When defining `f`, it is important to ensure that:
#'   - the column order of the first argument of `f` is consistent with the training input used for the emulator;
#'   - the column order of the output matrix of `f` is consistent with the order of emulator output dimensions (if `object` is an instance of the `dgp` class),
#'     or the order of emulators placed in `object` (if `object` is an instance of the `bundle` class).
#' * The output matrix produced by `f` may include `NA`s. This is especially beneficial as it allows the sequential design process to continue without interruption,
#'   even if errors or `NA` outputs are encountered from `f` at certain input locations identified by the sequential design. Users should ensure that any errors
#'   within `f` are handled by appropriately returning `NA`s.
#' * When defining `eval`, the output metric needs to be positive if [draw()] is used with `log = T`. And one needs to ensure that a lower metric value indicates
#'   a better emulation performance if `target` is set.
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#'
#' @examples
#' \dontrun{
#'
#' # load packages and the Python env
#' library(lhs)
#' library(dgpsi)
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
#' m_faster <- set_imp(m, 5)
#'
#' # plot the OOS validation with the faster DGP emulator
#' plot(m_faster, x_test = validate_x, y_test = validate_y)
#' }
#' @md
#' @name design
#' @export
design <- function(object, N, x_cand, y_cand, n_sample, n_cand, limits, f, reps, freq, x_test, y_test, reset, target, method, batch_size, eval, verb, autosave, new_wave, M_val, cores, ...){
  UseMethod("design")
}

#' @rdname design
#' @method design gp
#' @export
design.gp <- function(object, N, x_cand = NULL, y_cand = NULL, n_sample = 200, n_cand = lifecycle::deprecated(), limits = NULL, f = NULL, reps = 1, freq = c(1, 1), x_test = NULL, y_test = NULL, reset = FALSE, target = NULL, method = vigf, batch_size = 1, eval = NULL, verb = TRUE, autosave = list(), new_wave = TRUE, M_val = 50, cores = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"gp") ) stop("'object' must be an instance of the 'gp' class.", call. = FALSE)

  object <- copy_in_design(object)

  if (lifecycle::is_present(n_cand)) {
    lifecycle::deprecate_warn("2.5.0", "design(n_cand)", "design(n_sample)")
  }
  n_cand <- n_sample

  N <- check_N(N)
  reps <- check_reps(reps)
  freq <- check_freq(freq)
  reset <- check_reset(reset, N)
  n_cand <- check_n_cand(n_cand)

  X <- object$data$X
  Y <- object$data$Y
  n_dim_X <- ncol(X)
  n_dim_Y <- ncol(Y)

  if (!is.null(x_test) & !is.null(y_test)) {
    xy_test <- check_xy_test(x_test, y_test, n_dim_X, n_dim_Y)
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

  autosave_config <- list(switch = FALSE, directory = './check_points', fname = 'check_point', save_freq = 5, overwrite = FALSE)
  nms_autosave_config <- names(autosave_config)
  autosave_config[(nam_autosave <- names(autosave))] <- autosave
  if (length(noNms_autosave <- nam_autosave[!nam_autosave %in% nms_autosave_config])){
    warning("unknown names in 'autosave': ", paste(noNms_autosave,collapse=", "), call. = FALSE, immediate. = TRUE)
  }

  if (autosave_config$switch) {
    autosave_config$directory <- paste0(autosave_config$directory, '/emulator-', object$id)
    if (dir.exists(autosave_config$directory)){
      if (!first_time & !new_wave){
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave)
      } else {
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave+1)
      }
      if (!dir.exists(autosave_config$directory)){
        dir.create(autosave_config$directory, recursive = TRUE)
      }
    } else {
      if (!first_time & !new_wave){
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave)
      } else {
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave+1)
      }
      dir.create(autosave_config$directory, recursive = TRUE)
    }
  }

  if (!first_time & !new_wave){
    start_point <- object[["design"]][[paste('wave', n_wave, sep="")]][["N"]]
    N <- N + start_point
    if ( object[["design"]][[paste('wave', n_wave, sep="")]][["freq"]] != freq[2] ) stop("The frequency of emulator validation must match that in the last wave if 'new_wave = FALSE'.", call. = FALSE)
  } else {
    start_point <- 0
  }

  N_acq <- c()
  mnames <- methods::formalArgs(method)
  if ( !is.null(eval) ){
    vnames <- methods::formalArgs(eval)
  }

  if ( is.null(y_cand) ) {
    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)
    fnames <- methods::formalArgs(f)

    if ( "design" %in% names(object) ){
      if ( 'exclusion' %in% names(object$design) ){
        target_points <- object$design$exclusion
      } else {
        target_points <- c()
      }
    } else {
      target_points <- c()
    }

    if ( is.null(x_cand) ) {
      limits <- check_limits(limits, n_dim_X)
      #if ( is.null(limits) ){
      #  limits <- matrix(0, n_dim_X, 2)
      #  limits[,1] <- pkg.env$np$amin(X, axis=0L)
      #  limits[,2] <- pkg.env$np$amax(X, axis=0L)
      #}
      #if (identical(method, pei)) {
      #  add_arg <- list(pseudo_points = pp(X, limits))
      #  add_arg <- utils::modifyList(add_arg, list(...))
      #} else {
        add_arg <- list(...)
        add_arg$batch_size <- batch_size
        add_arg$limits <- limits
      #}
      #int <- check_int(int, n_dim_X)
    } else {
      add_arg <- list(...)
      add_arg$batch_size <- batch_size
      xy_cand_list <- check_xy_cand(x_cand, y_cand, n_dim_X, n_dim_Y)
      xy_cand_list <- remove_dup(xy_cand_list, rbind(X, target_points))
      x_cand_origin <- xy_cand_list[[1]]
      x_cand_rep <- pkg.env$np$unique(x_cand_origin, return_inverse=TRUE, axis=0L)
      if ( nrow(x_cand_rep[[1]])==nrow(x_cand_origin) ){
        x_cand <- x_cand_origin
      } else {
        x_cand <- x_cand_rep[[1]]
      }
      idx_x_cand0 <- c(1:nrow(x_cand))
      idx_x_cand <- idx_x_cand0
      idx_x_acq <- c()
    }

    if ( start_point==0 ){
      if ( verb ) message("Initializing ...", appendLF = FALSE)
      if ( is.null(eval) ){
        if (is.null(x_test) & is.null(y_test)){
          type <- 'loo'
          object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, cores = cores)
          rmse <- object$loo$rmse
        } else {
          type <- 'oos'
          object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, cores = cores)
          rmse <- object$oos$rmse
        }
      } else {
        type <- 'customized'
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
    } else {
      if ( is.null(eval) ){
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]][n_rmse,]
      } else {
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]][n_rmse,]
      }
      rmse_records <- c()
    }

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
      for ( i in (1+start_point):N ){

        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        tryagain <- TRUE

        while(tryagain){

          if ( !is.null(x_cand) ){
            if ( n_cand<=length(idx_x_cand) ){
              if (n_dim_X == 1){
                idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
              } else {
                idx_idx <- suppressPackageStartupMessages(clhs::clhs(as.data.frame(x_cand[idx_x_cand,,drop = FALSE]), size = n_cand, progress = FALSE, simple = TRUE))
                idx_sub_cand <- idx_x_cand[idx_idx]
              }
              #idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
            } else {
              idx_sub_cand <- idx_x_cand
            }
            sub_cand <- x_cand[idx_sub_cand,,drop = FALSE]
          } #else {
            #sub_cand <- lhs::maximinLHS(n_cand,n_dim_X)
            #if( !is.null(target_points) ){
            #  eud_dist <- sqrt(t(t(rowSums(sub_cand^2)-2*sub_cand%*%t(target_points))+rowSums(target_points^2)))
            #  rows_to_keep <- rowSums(eud_dist > 0.01) == ncol(eud_dist)
            #  sub_cand <- sub_cand[rows_to_keep,,drop = F]
            #}
            #sub_cand_lhd <- sub_cand
            #for (j in 1:n_dim_X){
            #  if ( int[j] ){
            #    #if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
            #    sub_cand[,j] <- floor( sub_cand[,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
            #  } else {
            #    sub_cand[,j] <- sub_cand[,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
            #  }
            #}
          #}
          if ( !is.null(x_cand) ){
            arg_list <- list(object, sub_cand)
          } else {
            arg_list <- list(object)
          }

          if ("..." %in% mnames){
            res <- do.call(method, c(arg_list, add_arg))
          } else {
            midx <- mnames %in% names(add_arg)
            mparam <- add_arg[mnames[midx]]
            res <- do.call(method, c(arg_list, mparam))
          }

          if ( verb ) {

            message(" done")
            if ( !is.null(x_cand) ){
              if ( length(res)==1 ){
                message(paste(c(" * Next design point:", sprintf("%.06f", sub_cand[res,])), collapse=" "))
              } else {
                for (j in 1:length(res) ){
                  message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", sub_cand[res[j],])), collapse=" "))
                }
              }
            } else {
              if ( nrow(res)==1 ){
                message(paste(c(" * Next design point:", sprintf("%.06f", res[1,])), collapse=" "))
              } else {
                for (j in 1:nrow(res) ){
                  message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", res[j,])), collapse=" "))
                }
              }
            }
          }

          if ( !is.null(x_cand) ){
            new_X <- sub_cand[rep(res, reps),,drop=F]
          } else {
            new_X <- res[rep(1:nrow(res), reps),,drop=F]
          }

          if ("..." %in% fnames){
            new_output <- do.call(f, c(list(new_X), add_arg))
          } else {
            fidx <- fnames %in% names(add_arg)
            fparam <- add_arg[fnames[fidx]]
            new_output <- do.call(f, c(list(new_X), fparam))
          }

          if (ifelse(is.list(new_output), all(!is.na(new_output[[1]])), all(!is.na(new_output)))) {
            tryagain <- FALSE
          } else if ( ifelse(is.list(new_output), all(is.na(new_output[[1]])), all(is.na(new_output))) ){
            if ( !is.null(x_cand) ) {
              target_points <- rbind(target_points, sub_cand[res,,drop=FALSE])
              idx_x_acq <- c(idx_x_acq, idx_sub_cand[res])
              idx_x_cand <- idx_x_cand0[-idx_x_acq]
            }
            if ( verb ) {
              message(" - Re-Locating ...", appendLF = FALSE)
            }
          } else {
            if ( is.list(new_output) ){
              keep_idx <- !is.na(new_output[[1]])
              new_X <- new_X[keep_idx,,drop=F]
              new_output[[1]] <- new_output[[1]][keep_idx,,drop=F]
            } else {
              keep_idx <- !is.na(new_output)
              new_X <- new_X[keep_idx,,drop=F]
              new_output <- new_output[keep_idx,,drop=F]
            }

            if ( !is.null(x_cand) ){
              keep_dat_idx <- unique(rep(res, reps)[keep_idx])
              target_points <- rbind(target_points, sub_cand[setdiff(res, keep_dat_idx),,drop=FALSE])
            }
            tryagain <- FALSE
          }

        }

        N_acq <- c(N_acq, nrow(new_X))

        X <- rbind(X, new_X)
        if ( !is.null(x_cand) ){
          idx_x_acq <- c(idx_x_acq, idx_sub_cand[res])
          idx_x_cand <- idx_x_cand0[-idx_x_acq]
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

        if ( i %% freq[1]==0 ){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = TRUE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = FALSE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
          if ( verb ) Sys.sleep(0.5)
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N ){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, force = TRUE, M = M_val, cores = cores)
              if ( verb ) message(" done")
              rmse <- object$loo$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, force = TRUE, M = M_val, cores = cores)
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
            message(sprintf("Target reached! The sequential design stopped at step %i.", i))
            istarget <- TRUE
            N <- i
            break
          }
        }

        if (autosave_config$switch && i %% autosave_config$save_freq == 0) {
          if (i != N) {
            if ( verb ) message(" - Auto-saving ...", appendLF = FALSE)
            save_file_name <- if (autosave_config$overwrite) {
              autosave_config$fname
            } else {
              paste0(autosave_config$fname, "_", i)
            }
            object_temp <- pack_gp(object, first_time, i, eval, new_wave, rmse_records, freq, N_acq, target, type, y_cand,
                                   n_wave = if(isFALSE(first_time)) n_wave else NULL,
                                   design_info = if(isFALSE(first_time)) design_info else NULL,
                                   x_test = if(identical(type, 'oos')) x_test else NULL,
                                   y_test = if(identical(type, 'oos')) y_test else NULL,
                                   istarget = if(!is.null(target)) istarget else NULL,
                                   target_points = if(is.null(y_cand)) target_points else NULL)
            save_file_path <- file.path(autosave_config$directory, save_file_name)
            write(object_temp, save_file_path)
            if ( verb ) message(" done")
          }
        }
      }
    }
  } else {
    add_arg <- list(...)
    add_arg$batch_size <- batch_size
    xy_cand_list <- check_xy_cand(x_cand, y_cand, n_dim_X, n_dim_Y)
    xy_cand_list <- remove_dup(xy_cand_list, X)
    x_cand_origin <- xy_cand_list[[1]]
    x_cand_rep <- pkg.env$np$unique(x_cand_origin, return_inverse=TRUE, axis=0L)
    if ( nrow(x_cand_rep[[1]])==nrow(x_cand_origin) ){
      x_cand <- x_cand_origin
      is.dup <- FALSE
    } else {
      x_cand <- x_cand_rep[[1]]
      is.dup <- TRUE
    }
    y_cand <- xy_cand_list[[2]]
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- c()

    if ( start_point==0 ){
      if ( verb ) message("Initializing ...", appendLF = FALSE)
      if ( is.null(eval) ){
        if (is.null(x_test) & is.null(y_test)){
          type <- 'loo'
          object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, cores = cores)
          rmse <- object$loo$rmse
        } else {
          type <- 'oos'
          object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, cores = cores)
          rmse <- object$oos$rmse
        }
      } else {
        type <- 'customized'
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
    } else {
      if ( is.null(eval) ){
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]][n_rmse,]
      } else {
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]][n_rmse,]
      }
      rmse_records <- c()
    }

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
      for ( i in (1+start_point):N ){
        if ( n_cand<=length(idx_x_cand) ){
          if (n_dim_X == 1){
            idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
          } else {
            idx_idx <- suppressPackageStartupMessages(clhs::clhs(as.data.frame(x_cand[idx_x_cand,,drop = FALSE]), size = n_cand, progress = FALSE, simple = TRUE))
            idx_sub_cand <- idx_x_cand[idx_idx]
          }
          #idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
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
          if ( length(res)==1 ){
            message(paste(c(" * Next design point:", sprintf("%.06f", sub_cand[res,])), collapse=" "))
          } else {
            for (j in 1:length(res) ){
              message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", sub_cand[res[j],])), collapse=" "))
            }
          }
        }

        new_X_temp <- sub_cand[res,,drop=FALSE]
        if ( is.dup ){
          new_idx <- extract_all(new_X_temp, x_cand_origin)
          new_X <- x_cand_origin[new_idx,,drop=FALSE]
          N_acq <- c(N_acq, nrow(new_X))
          X <- rbind(X, new_X)
          Y <- rbind(Y, y_cand[new_idx,,drop=FALSE])
          idx_x_acq <- c(idx_x_acq, idx_sub_cand[res])
          idx_x_cand <- idx_x_cand0[-idx_x_acq]
        } else {
          new_idx <- extract_all(new_X_temp, sub_cand)
          new_X <- sub_cand[new_idx,,drop=FALSE]
          N_acq <- c(N_acq, nrow(new_X))
          X <- rbind(X, new_X)
          Y <- rbind(Y, y_cand[idx_sub_cand,,drop=FALSE][new_idx,,drop=FALSE])
          idx_x_acq <- c(idx_x_acq, idx_sub_cand[new_idx])
          idx_x_cand <- idx_x_cand0[-idx_x_acq]
        }

        if ( i %% freq[1]==0 ){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = TRUE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = FALSE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
          if ( verb ) Sys.sleep(0.5)
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N ){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE, cores = cores)
              if ( verb ) message(" done")
              rmse <- object$loo$rmse
              if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, force = TRUE, cores = cores)
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
            message(sprintf("Target reached! The sequential design stopped at step %i.", i))
            istarget <- TRUE
            N <- i
            break
          }
        }

        if (autosave_config$switch && i %% autosave_config$save_freq == 0) {
          if (i != N) {
            if ( verb ) message(" - Auto-saving ...", appendLF = FALSE)
            save_file_name <- if (autosave_config$overwrite) {
              autosave_config$fname
            } else {
              paste0(autosave_config$fname, "_", i)
            }
            object_temp <- pack_gp(object, first_time, i, eval, new_wave, rmse_records, freq, N_acq, target, type, y_cand,
                                   n_wave = if(isFALSE(first_time)) n_wave else NULL,
                                   design_info = if(isFALSE(first_time)) design_info else NULL,
                                   x_test = if(identical(type, 'oos')) x_test else NULL,
                                   y_test = if(identical(type, 'oos')) y_test else NULL,
                                   istarget = if(!is.null(target)) istarget else NULL,
                                   target_points = if(is.null(y_cand)) target_points else NULL)
            save_file_path <- file.path(autosave_config$directory, save_file_name)
            write(object_temp, save_file_path)
            if ( verb ) message(" done")
          }
        }
      }
    }
  }

  if ( run ){
    object <- pack_gp(object, first_time, N, eval, new_wave, rmse_records, freq, N_acq, target, type, y_cand,
                      n_wave = if(isFALSE(first_time)) n_wave else NULL,
                      design_info = if(isFALSE(first_time)) design_info else NULL,
                      x_test = if(identical(type, 'oos')) x_test else NULL,
                      y_test = if(identical(type, 'oos')) y_test else NULL,
                      istarget = if(!is.null(target)) istarget else NULL,
                      target_points = if(is.null(y_cand)) target_points else NULL)

    if (autosave_config$switch) {
      if ( verb ) message(" - Auto-saving the final iteration ...", appendLF = FALSE)
      save_file_name <- if (autosave_config$overwrite) {
        autosave_config$fname
      } else {
        paste0(autosave_config$fname, "_", N)
      }
      save_file_path <- file.path(autosave_config$directory, save_file_name)
      write(object, save_file_path)
      if ( verb ) message(" done")
    }

    if( !is.null(target) ) {
      if ( !istarget ) message("The target was not reached at the end of the sequential design.")
    }
  } else {
    message("Target already reached. Sequential design not performed.")
  }

  return(object)
}

#' @rdname design
#' @method design dgp
#' @export
design.dgp <- function(object, N, x_cand = NULL, y_cand = NULL, n_sample = 200, n_cand = lifecycle::deprecated(), limits = NULL, f = NULL, reps = 1, freq = c(1, 1), x_test = NULL, y_test = NULL, reset = FALSE, target = NULL, method = vigf, batch_size = 1, eval = NULL, verb = TRUE, autosave = list(), new_wave = TRUE, M_val = 50, cores = 1, train_N = NULL, refit_cores = 1, pruning = TRUE, control = list(), ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"dgp") ) stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)

  if (object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$name == "Categorical") {
    is.categorical <- TRUE
  } else {
    is.categorical <- FALSE
  }
  object <- copy_in_design(object)

  if (lifecycle::is_present(n_cand)) {
    lifecycle::deprecate_warn("2.5.0", "design(n_cand)", "design(n_sample)")
  }
  n_cand <- n_sample

  N <- check_N(N)
  freq <- check_freq(freq)
  reset <- check_reset(reset, N)
  if ( !is.null(train_N) ) train_N <- check_train_N(train_N, N)
  reps <- check_reps(reps)
  if( !is.null(refit_cores) ) {
    refit_cores <- as.integer(refit_cores)
    if ( refit_cores < 1 ) stop("'refit_cores' must be >= 1.", call. = FALSE)
  }
  n_cand <- check_n_cand(n_cand)

  X <- object$data$X
  Y <- object$data$Y
  n_dim_X <- ncol(X)
  n_dim_Y <- ncol(Y)

  if (!is.null(x_test) & !is.null(y_test)) {
    xy_test <- check_xy_test(x_test, y_test, n_dim_X, n_dim_Y)
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

  autosave_config <- list(switch = FALSE, directory = './check_points', fname = 'check_point', save_freq = 5, overwrite = FALSE)
  nms_autosave_config <- names(autosave_config)
  autosave_config[(nam_autosave <- names(autosave))] <- autosave
  if (length(noNms_autosave <- nam_autosave[!nam_autosave %in% nms_autosave_config])){
    warning("unknown names in 'autosave': ", paste(noNms_autosave,collapse=", "), call. = FALSE, immediate. = TRUE)
  }

  if (autosave_config$switch) {
    autosave_config$directory <- paste0(autosave_config$directory, '/emulator-', object$id)
    if (dir.exists(autosave_config$directory)){
      if (!first_time & !new_wave){
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave)
      } else {
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave+1)
      }
      if (!dir.exists(autosave_config$directory)){
        dir.create(autosave_config$directory, recursive = TRUE)
      }
    } else {
      if (!first_time & !new_wave){
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave)
      } else {
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave+1)
      }
      dir.create(autosave_config$directory, recursive = TRUE)
    }
  }

  if (!first_time & !new_wave){
    start_point <- object[["design"]][[paste('wave', n_wave, sep="")]][["N"]]
    N <- N + start_point
    if ( object[["design"]][[paste('wave', n_wave, sep="")]][["freq"]] != freq[2] ) stop("The frequency of emulator validation must match to that in the last wave if 'new_wave = FALSE'.", call. = FALSE)
  } else {
    start_point <- 0
  }

  if (pruning){
    pruning <- check_auto(object)
    if (pruning){
      drop_list <- create_drop_list(object)
      #con <- list(threshold = 0.97, nexceed = 3, warmup = 15)
      con <- list(min_size = 10*n_dim_X, threshold = 0.97, nexceed = 3)
      nmsC <- names(con)
      con[(namc <- names(control))] <- control
      if (length(noNms <- namc[!namc %in% nmsC])){
        warning("unknown names in 'control': ", paste(noNms,collapse=", "), call. = FALSE, immediate. = TRUE)
      }

      if (con$min_size<1) stop("'min_size' in 'control' must be at least 1.", call. = FALSE)
      if (con$threshold>1 || con$threshold<0) stop("'threshold' in 'control' must be between 0 and 1.", call. = FALSE)
      if (con$nexceed<1) stop("'nexceed' in 'control' must be at least 1.", call. = FALSE)

      required_step <- con$nexceed + ifelse((nrow(X)>con$min_size), 0, con$min_size-nrow(X) )
      if (required_step > N){
        warning("'N' needs to be greater than ", paste(required_step), " to trigger dynamic pruning of the DGP emulator.", call. = FALSE, immediate. = TRUE)
      }
    }
  }

  N_acq <- c()
  mnames <- methods::formalArgs(method)
  if ( !is.null(eval) ){
    vnames <- methods::formalArgs(eval)
  }
  #if candidate set is given
  if ( is.null(y_cand) ) {

    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)
    fnames <- methods::formalArgs(f)

    if ( "design" %in% names(object) ){
      if ( 'exclusion' %in% names(object$design) ){
        target_points <- object$design$exclusion
      } else {
        target_points <- c()
      }
    } else {
      target_points <- c()
    }

    if ( is.null(x_cand) ) {
      limits <- check_limits(limits, n_dim_X)
      #if ( is.null(limits) ){
      #  limits <- matrix(0, n_dim_X, 2)
      #  limits[,1] <- pkg.env$np$amin(X, axis=0L)
      #  limits[,2] <- pkg.env$np$amax(X, axis=0L)
      #}
      #if (identical(method, pei)) {
      #  add_arg <- list(pseudo_points = pp(X, limits))
      #  add_arg <- utils::modifyList(add_arg, list(...))
      #} else {
        add_arg <- list(...)
        add_arg$batch_size <- batch_size
        add_arg$limits <- limits
      #}
      #int <- check_int(int, n_dim_X)
    } else {
      add_arg <- list(...)
      add_arg$batch_size <- batch_size
      xy_cand_list <- check_xy_cand(x_cand, y_cand, n_dim_X, n_dim_Y)
      xy_cand_list <- remove_dup(xy_cand_list, rbind(X, target_points))
      x_cand_origin <- xy_cand_list[[1]]
      x_cand_rep <- pkg.env$np$unique(x_cand_origin, return_inverse=TRUE, axis=0L)
      if ( nrow(x_cand_rep[[1]])==nrow(x_cand_origin) ){
        x_cand <- x_cand_origin
      } else {
        x_cand <- x_cand_rep[[1]]
      }
      idx_x_cand0 <- c(1:nrow(x_cand))
      idx_x_cand <- idx_x_cand0
      idx_x_acq <- c()
    }

    if (start_point == 0){
      if ( verb ) message("Initializing ...", appendLF = FALSE)
      if ( is.null(eval) ){
        if (is.null(x_test) & is.null(y_test)){
          type <- 'loo'
          object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, cores = cores, ...)
          if (is.categorical){
            rmse <- object$loo$log_loss
          } else {
            rmse <- object$loo$rmse
          }
        } else {
          type <- 'oos'
          object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, cores = cores, ...)
          if (is.categorical){
            rmse <- object$oos$log_loss
          } else {
            rmse <- object$oos$rmse
          }
        }
      } else {
        type <- 'customized'
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
        if (is.categorical){
          if ( verb ) message(paste(c(" * Log-Loss:", sprintf("%.06f", rmse)), collapse=" "))
        } else {
          if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        }
      } else {
        if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
      }
      rmse_records <- rmse
    } else {
      if ( is.null(eval) ){
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]][n_rmse,]
      } else {
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]][n_rmse,]
      }
      rmse_records <- c()
    }

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
      for ( i in (1+start_point):N ){
        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        tryagain <- TRUE

        while(tryagain){
          if ( !is.null(x_cand) ){
            #sub_cand <- lhs::maximinLHS(n_cand,n_dim_X)
            #if( !is.null(target_points) ){
            #  eud_dist <- sqrt(t(t(rowSums(sub_cand^2)-2*sub_cand%*%t(target_points))+rowSums(target_points^2)))
            #  rows_to_keep <- rowSums(eud_dist > 0.01) == ncol(eud_dist)
            #  sub_cand <- sub_cand[rows_to_keep,,drop = F]
            #}
            #sub_cand_lhd <- sub_cand
            #for (j in 1:n_dim_X){
            #  if ( int[j] ){
            #    #if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
            #    sub_cand[,j] <- floor( sub_cand[,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
            #  } else {
            #    sub_cand[,j] <- sub_cand[,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
            #  }
            #}
          #} else {
            if ( n_cand<=length(idx_x_cand) ){
              if (n_dim_X == 1){
                idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
              } else {
                idx_idx <- suppressPackageStartupMessages(clhs::clhs(as.data.frame(x_cand[idx_x_cand,,drop = FALSE]), size = n_cand, progress = FALSE, simple = TRUE))
                idx_sub_cand <- idx_x_cand[idx_idx]
              }
            } else {
              idx_sub_cand <- idx_x_cand
            }
            sub_cand <- x_cand[idx_sub_cand,,drop = FALSE]
          }

          if ( !is.null(x_cand) ){
            arg_list <- list(object, sub_cand)
          } else {
            arg_list <- list(object)
          }

          if ("..." %in% mnames){
            res <- do.call(method, c(arg_list, add_arg))
          } else {
            midx <- mnames %in% names(add_arg)
            mparam <- add_arg[mnames[midx]]
            res <- do.call(method, c(arg_list, mparam))
          }

          if ( verb ) {
            message(" done")
            if ( !is.null(x_cand) ){
              if ( length(res)==1 ){
                message(paste(c(" * Next design point:", sprintf("%.06f", sub_cand[res,])), collapse=" "))
              } else {
                for ( j in 1:length(res) ){
                  message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", sub_cand[res[j],])), collapse=" "))
                }
              }
            } else {
              if ( nrow(res)==1 ){
                message(paste(c(" * Next design point:", sprintf("%.06f", res[1,])), collapse=" "))
              } else {
                for (j in 1:nrow(res) ){
                  message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", res[j,])), collapse=" "))
                }
              }
            }
          }

          if ( !is.null(x_cand) ){
            new_X <- sub_cand[rep(unique(res), reps),,drop=FALSE]
          } else {
            res_unique <- pkg.env$np$unique(res, axis=0L)
            new_X <- res_unique[rep(1:nrow(res_unique), reps),,drop=F]
          }

          if ("..." %in% fnames){
            new_output <- do.call(f, c(list(new_X), add_arg))
          } else {
            fidx <- fnames %in% names(add_arg)
            fparam <- add_arg[fnames[fidx]]
            new_output <- do.call(f, c(list(new_X), fparam))
          }

          if (ifelse(is.list(new_output), all(!is.na(new_output[[1]])), all(!is.na(new_output)))) {
            tryagain <- FALSE
          } else if ( ifelse(is.list(new_output), all(!pkg.env$np$prod(!is.na(new_output[[1]]), axis=as.integer(1))), all(!pkg.env$np$prod(!is.na(new_output), axis=as.integer(1)))) )  {
            if ( !is.null(x_cand) ){
              target_points <- rbind(target_points, sub_cand[unique(res),,drop=FALSE])
              idx_x_acq <- c(idx_x_acq, idx_sub_cand[unique(res)])
              idx_x_cand <- idx_x_cand0[-idx_x_acq]
            }
            if ( verb ) {
              message(" - Re-Locating ...", appendLF = FALSE)
            }
          } else {
            if ( is.list(new_output) ){
              keep_idx <- pkg.env$np$prod(!is.na(new_output[[1]]), axis=as.integer(1))
              new_X <- new_X[keep_idx,,drop=F]
              new_output[[1]] <- new_output[[1]][keep_idx,,drop=F]
            } else {
              keep_idx <- pkg.env$np$prod(!is.na(new_output), axis=as.integer(1))
              new_X <- new_X[keep_idx,,drop=F]
              new_output <- new_output[keep_idx,,drop=F]
            }

            if ( !is.null(x_cand) ) {
              keep_dat_idx <- unique(rep(unique(res), reps)[keep_idx])
              target_point_idx <- setdiff(unique(res), keep_dat_idx)
              if (length(target_point_idx)!=0) target_points <- rbind(target_points, sub_cand[target_point_idx,,drop=FALSE])
            }
            tryagain <- FALSE
          }

        }

        N_acq <- c(N_acq, nrow(new_X))

        X <- rbind(X, new_X)
        if ( !is.null(x_cand) ){
          idx_x_acq <- c(idx_x_acq, idx_sub_cand[unique(res)])
          idx_x_cand <- idx_x_cand0[-idx_x_acq]
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

        if ( i %% freq[1]==0 ){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = TRUE, reset = reset[i-start_point], verb = FALSE, N = if(!is.null(train_N)) train_N[i-start_point] else NULL, cores = refit_cores, B = ifelse(length(object$emulator_obj$all_layer_set)<10, length(object$emulator_obj$all_layer_set), 10), update_in_design = NULL)
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = FALSE, reset = reset[i-start_point], verb = FALSE, B = ifelse(length(object$emulator_obj$all_layer_set)<10, length(object$emulator_obj$all_layer_set), 10), update_in_design = NULL)
          if ( verb ) Sys.sleep(0.5)
          if ( verb ) message(" done")
        }

        if (pruning){
          if (nrow(X) > con$min_size){
            r2 <- object$constructor_obj$aggregate_r2()
            crop_id_list <- vector("list", length = length(drop_list))
            N_cropped <- 0
            for (l in length(drop_list):1){
              above_threshold <- (r2[[l+1]][[1]]>con$threshold)
              drop_list[[l]] <- drop_list[[l]] + above_threshold
              drop_list[[l]][!above_threshold] <- 0
              crop_id <- drop_list[[l]]==con$nexceed
              crop_id_list[[l]] <- crop_id
              N_cropped <- N_cropped + sum(crop_id)
            }
            if (N_cropped!=0) {
              object <- crop(object, crop_id_list, refit_cores, verb)
              if ( inherits(object,"dgp") ) {
                pruning <- check_auto(object)
                if (pruning) drop_list <- create_drop_list(object)
              }
            }
          }
        }

        if ( i %% freq[2]==0 | i==N){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              if ( inherits(object,"dgp") ) {
                object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                if (object$constructor_obj$all_layer[[length(object$constructor_obj$all_layer)]][[1]]$name == "Categorical"){
                  rmse <- object$loo$log_loss
                  if ( verb ) message(" done")
                  if ( verb ) message(paste(c(" * Log-Loss:", sprintf("%.06f", rmse)), collapse=" "))
                } else {
                  rmse <- object$loo$rmse
                  if ( verb ) message(" done")
                  if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
                }
              } else if ( inherits(object,"gp") ) {
                object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE)
                rmse <- object$loo$rmse
                if ( verb ) message(" done")
                if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
              } else if ( inherits(object,"bundle") ) {
                n_emulators <- length(object) - 1
                if ( "id" %in% names(object) ) n_emulators <- n_emulators - 1
                for ( k in 1:n_emulators ){
                  object[[paste('emulator',k,sep='')]] <- validate(object[[paste('emulator',k,sep='')]], x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE)
                  rmse[k] <- object[[paste('emulator',k,sep='')]]$loo$rmse
                }
                if ( verb ) message(" done")
                if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
              }
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              if ( inherits(object,"dgp") ) {
                object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                if (object$constructor_obj$all_layer[[length(object$constructor_obj$all_layer)]][[1]]$name == "Categorical"){
                  rmse <- object$oos$log_loss
                  if ( verb ) message(" done")
                  if ( verb ) message(paste(c(" * Log-Loss:", sprintf("%.06f", rmse)), collapse=" "))
                } else {
                  rmse <- object$oos$rmse
                  if ( verb ) message(" done")
                  if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
                }
              } else if ( inherits(object,"gp") ) {
                object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, force = TRUE)
                rmse <- object$oos$rmse
                if ( verb ) message(" done")
                if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
              } else if ( inherits(object,"bundle") ){
                n_emulators <- length(object) - 1
                if ( "id" %in% names(object) ) n_emulators <- n_emulators - 1
                for ( k in 1:n_emulators ){
                  object[[paste('emulator',k,sep='')]] <- validate(object[[paste('emulator',k,sep='')]], x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, force = TRUE)
                  rmse[k] <- object[[paste('emulator',k,sep='')]]$oos$rmse
                }
                if ( verb ) message(" done")
                if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
              }
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
            message(sprintf("Target reached! Sequential design stopped at step %i.", i))
            if ( inherits(object,"bundle") ) {
              istarget <- rep(TRUE, n_emulators)
            } else {
              istarget <- TRUE
            }
            N <- i
            break
          } else {
            if ( inherits(object,"bundle") ) {
              istarget <- rmse<=target
              if ( length(istarget)==1 ) istarget <- rep(istarget, n_emulators)
            }
          }
        }

        if (autosave_config$switch && i %% autosave_config$save_freq == 0) {
          if (i != N) {
            if ( verb ) message(" - Auto-saving ...", appendLF = FALSE)
            save_file_name <- if (autosave_config$overwrite) {
              autosave_config$fname
            } else {
              paste0(autosave_config$fname, "_", i)
            }
            object_temp <- pack_dgp(object, first_time, i, eval, new_wave, rmse_records, freq, N_acq, n_dim_Y, target, type, y_cand,
                                    n_wave = if(isFALSE(first_time)) n_wave else NULL,
                                    design_info = if(isFALSE(first_time)) design_info else NULL,
                                    x_test = if(identical(type, 'oos')) x_test else NULL,
                                    y_test = if(identical(type, 'oos')) y_test else NULL,
                                    istarget = if(!is.null(target)) istarget else NULL,
                                    target_points = if(is.null(y_cand)) target_points else NULL)
            save_file_path <- file.path(autosave_config$directory, save_file_name)
            write(object_temp, save_file_path)
            if ( verb ) message(" done")
          }
        }

        if (inherits(object,"bundle") | inherits(object,"gp")) {
          remaining_steps <- N-i
          if (remaining_steps > 0) {
            N <- i
            break
          }
        }

      }
    }
  } else {
    add_arg <- list(...)
    add_arg$batch_size <- batch_size
    xy_cand_list <- check_xy_cand(x_cand, y_cand, n_dim_X, n_dim_Y)
    xy_cand_list <- remove_dup(xy_cand_list, X)
    x_cand_origin <- xy_cand_list[[1]]
    x_cand_rep <- pkg.env$np$unique(x_cand_origin, return_inverse=TRUE, axis=0L)
    if ( nrow(x_cand_rep[[1]])==nrow(x_cand_origin) ){
      x_cand <- x_cand_origin
      is.dup <- FALSE
    } else {
      x_cand <- x_cand_rep[[1]]
      is.dup <- TRUE
    }
    y_cand <- xy_cand_list[[2]]
    idx_x_cand0 <- c(1:nrow(x_cand))
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- c()

    if (start_point == 0){
      if ( verb ) message("Initializing ...", appendLF = FALSE)
      if ( is.null(eval) ){
        if (is.null(x_test) & is.null(y_test)){
          type <- 'loo'
          object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, cores = cores, ...)
          if (is.categorical){
            rmse <- object$loo$log_loss
          } else {
            rmse <- object$loo$rmse
          }
        } else {
          type <- 'oos'
          object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, cores = cores, ...)
          if (is.categorical){
            rmse <- object$oos$log_loss
          } else {
            rmse <- object$oos$rmse
          }
        }
      } else {
        type <- 'customized'
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
        if (is.categorical){
          if ( verb ) message(paste(c(" * Log-Loss:", sprintf("%.06f", rmse)), collapse=" "))
        } else {
          if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
        }
      } else {
        if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
      }
      rmse_records <- rmse
    } else {
      if ( is.null(eval) ){
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]][n_rmse,]
        rmse_records <- c()
      } else {
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]][n_rmse,]
        rmse_records <- c()
      }
    }

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
      for ( i in (1+start_point):N ){
        if ( n_cand<=length(idx_x_cand) ){
          if (n_dim_X == 1){
            idx_sub_cand <- sample(idx_x_cand, n_cand, replace = FALSE)
          } else {
            idx_idx <- suppressPackageStartupMessages(clhs::clhs(as.data.frame(x_cand[idx_x_cand,,drop = FALSE]), size = n_cand, progress = FALSE, simple = TRUE))
            idx_sub_cand <- idx_x_cand[idx_idx]
          }
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
          if ( length(res)==1 ){
            message(paste(c(" * Next design point:", sprintf("%.06f", sub_cand[res,])), collapse=" "))
          } else {
            for ( j in 1:length(res) ){
              message(paste(c(sprintf(" * Next design point (Position%i):", j), sprintf("%.06f", sub_cand[res[j],])), collapse=" "))
            }
          }
        }

        new_X_temp <- sub_cand[unique(res),,drop=FALSE]
        if ( is.dup ){
          new_idx <- extract_all(new_X_temp, x_cand_origin)
          new_X <- x_cand_origin[new_idx,,drop=FALSE]
          N_acq <- c(N_acq, nrow(new_X))
          X <- rbind(X, new_X)
          Y <- rbind(Y, y_cand[new_idx,,drop=FALSE])
          idx_x_acq <- c(idx_x_acq, idx_sub_cand[unique(res)])
          idx_x_cand <- idx_x_cand0[-idx_x_acq]
        } else {
          new_idx <- extract_all(new_X_temp, sub_cand)
          new_X <- sub_cand[new_idx,,drop=FALSE]
          N_acq <- c(N_acq, nrow(new_X))
          X <- rbind(X, new_X)
          Y <- rbind(Y, y_cand[idx_sub_cand,,drop=FALSE][new_idx,,drop=FALSE])
          idx_x_acq <- c(idx_x_acq, idx_sub_cand[new_idx])
          idx_x_cand <- idx_x_cand0[-idx_x_acq]
        }

        if ( i %% freq[1]==0 ){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = TRUE, reset = reset[i-start_point], verb = FALSE, N = if(!is.null(train_N)) train_N[i-start_point] else NULL, cores = refit_cores, B = ifelse(length(object$emulator_obj$all_layer_set)<10, length(object$emulator_obj$all_layer_set), 10), update_in_design = NULL)
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          object <- update(object, X, Y, refit = FALSE, reset = reset[i-start_point], verb = FALSE, B = ifelse(length(object$emulator_obj$all_layer_set)<10, length(object$emulator_obj$all_layer_set), 10), update_in_design = NULL)
          if ( verb ) message(" done")
        }

        if (pruning){
          if (nrow(X) > con$min_size){
            r2 <- object$constructor_obj$aggregate_r2()
            crop_id_list <- vector("list", length = length(drop_list))
            N_cropped <- 0
            for (l in length(drop_list):1){
              above_threshold <- (r2[[l+1]][[1]]>con$threshold)
              drop_list[[l]] <- drop_list[[l]] + above_threshold
              drop_list[[l]][!above_threshold] <- 0
              crop_id <- drop_list[[l]]==con$nexceed
              crop_id_list[[l]] <- crop_id
              N_cropped <- N_cropped + sum(crop_id)
            }
            if (N_cropped!=0) {
              object <- crop(object, crop_id_list, refit_cores, verb)
              if ( inherits(object,"dgp") ) {
                pruning <- check_auto(object)
                if (pruning) drop_list <- create_drop_list(object)
              }
            }
          }
        }

        if ( i %% freq[2]==0 | i==N){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              if ( inherits(object,"dgp") ) {
                object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                if (object$constructor_obj$all_layer[[length(object$constructor_obj$all_layer)]][[1]]$name == "Categorical"){
                  rmse <- object$loo$log_loss
                  if ( verb ) message(" done")
                  if ( verb ) message(paste(c(" * Log-Loss:", sprintf("%.06f", rmse)), collapse=" "))
                } else{
                  rmse <- object$loo$rmse
                  if ( verb ) message(" done")
                  if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
                }
              } else if ( inherits(object,"gp") ) {
                object <- validate(object, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE)
                rmse <- object$loo$rmse
                if ( verb ) message(" done")
                if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
              } else if ( inherits(object,"bundle") ) {
                n_emulators <- length(object) - 1
                if ( "id" %in% names(object) ) n_emulators <- n_emulators - 1
                for ( k in 1:n_emulators ){
                  object[[paste('emulator',k,sep='')]] <- validate(object[[paste('emulator',k,sep='')]], x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE)
                  rmse[k] <- object[[paste('emulator',k,sep='')]]$loo$rmse
                }
                if ( verb ) message(" done")
                if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
              }
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              if ( inherits(object,"dgp") ) {
                object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                if (object$constructor_obj$all_layer[[length(object$constructor_obj$all_layer)]][[1]]$name == "Categorical"){
                  rmse <- object$oos$log_loss
                  if ( verb ) message(" done")
                  if ( verb ) message(paste(c(" * Log-Loss:", sprintf("%.06f", rmse)), collapse=" "))
                } else {
                  rmse <- object$oos$rmse
                  if ( verb ) message(" done")
                  if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
                }
              } else if ( inherits(object,"gp") ) {
                object <- validate(object, x_test = x_test, y_test = y_test, verb = FALSE, M = M_val, force = TRUE)
                rmse <- object$oos$rmse
                if ( verb ) message(" done")
                if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
              } else if ( inherits(object,"bundle") ){
                n_emulators <- length(object) - 1
                if ( "id" %in% names(object) ) n_emulators <- n_emulators - 1
                for ( k in 1:n_emulators ){
                  object[[paste('emulator',k,sep='')]] <- validate(object[[paste('emulator',k,sep='')]], x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, force = TRUE)
                  rmse[k] <- object[[paste('emulator',k,sep='')]]$oos$rmse
                }
                if ( verb ) message(" done")
                if ( verb ) message(paste(c(" * RMSE:", sprintf("%.06f", rmse)), collapse=" "))
              }
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
            message(sprintf("Target reached! Sequential design stopped at step %i.", i))
            if ( inherits(object,"bundle") ) {
              istarget <- rep(TRUE, n_emulators)
            } else {
              istarget <- TRUE
            }
            N <- i
            break
          } else {
            if ( inherits(object,"bundle") ) {
              istarget <- rmse<=target
              if ( length(istarget)==1 ) istarget <- rep(istarget, n_emulators)
            }
          }
        }

        if (autosave_config$switch && i %% autosave_config$save_freq == 0) {
          if (i != N) {
            if ( verb ) message(" - Auto-saving ...", appendLF = FALSE)
            save_file_name <- if (autosave_config$overwrite) {
              autosave_config$fname
            } else {
              paste0(autosave_config$fname, "_", i)
            }
            object_temp <- pack_dgp(object, first_time, i, eval, new_wave, rmse_records, freq, N_acq, n_dim_Y, target, type, y_cand,
                                    n_wave = if(isFALSE(first_time)) n_wave else NULL,
                                    design_info = if(isFALSE(first_time)) design_info else NULL,
                                    x_test = if(identical(type, 'oos')) x_test else NULL,
                                    y_test = if(identical(type, 'oos')) y_test else NULL,
                                    istarget = if(!is.null(target)) istarget else NULL,
                                    target_points = if(is.null(y_cand)) target_points else NULL)
            save_file_path <- file.path(autosave_config$directory, save_file_name)
            write(object_temp, save_file_path)
            if ( verb ) message(" done")
          }
        }

        if (inherits(object,"bundle") | inherits(object,"gp")) {
          remaining_steps <- N-i
          if (remaining_steps > 0) {
            N <- i
            break
          }
        }

      }
    }
  }

  if ( run ){
    object <- pack_dgp(object, first_time, N, eval, new_wave, rmse_records, freq, N_acq, n_dim_Y, target, type, y_cand,
                       n_wave = if(isFALSE(first_time)) n_wave else NULL,
                       design_info = if(isFALSE(first_time)) design_info else NULL,
                       x_test = if(identical(type, 'oos')) x_test else NULL,
                       y_test = if(identical(type, 'oos')) y_test else NULL,
                       istarget = if(!is.null(target)) istarget else NULL,
                       target_points = if(is.null(y_cand)) target_points else NULL)

    if (!inherits(object,"dgp")){
      if (remaining_steps > 0){
        object <- design(object, remaining_steps, x_cand = if (is.null(x_cand)) {NULL} else {xy_cand_list[[1]]}, y_cand = if (is.null(x_cand)) {NULL} else {xy_cand_list[[2]]},
                         n_sample = n_cand, limits = limits, f = f, reps = reps, freq = freq, x_test = x_test, y_test = y_test, reset = utils::tail(reset, remaining_steps), target = target,
                         method = method, batch_size = batch_size, eval = eval, verb = verb, autosave = autosave, new_wave = FALSE, M_val = M_val, cores = cores, train_N = if(!is.null(train_N)) utils::tail(train_N, remaining_steps) else NULL, refit_cores = refit_cores, ...)
        return(object)
      }
    }

    if (autosave_config$switch) {
      if ( verb ) message(" - Auto-saving the final iteration ...", appendLF = FALSE)
      save_file_name <- if (autosave_config$overwrite) {
        autosave_config$fname
      } else {
        paste0(autosave_config$fname, "_", N)
      }
      save_file_path <- file.path(autosave_config$directory, save_file_name)
      write(object, save_file_path)
      if ( verb ) message(" done")
    }

    if( !is.null(target) ) {
      if ( !all(istarget) ) message("Target not reached at the end of the sequential design.")
    }
  } else {
    message("Target already reached. Sequential design not performed.")
  }

  return(object)
}


#' @rdname design
#' @method design bundle
#' @export
design.bundle <- function(object, N, x_cand = NULL, y_cand = NULL, n_sample = 200, n_cand = lifecycle::deprecated(), limits = NULL, f = NULL, reps = 1, freq = c(1, 1), x_test = NULL, y_test = NULL, reset = FALSE, target = NULL, method = vigf, batch_size = 1, eval = NULL, verb = TRUE, autosave = list(), new_wave = TRUE, M_val = 50, cores = 1, train_N = NULL, refit_cores = 1,  ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( !inherits(object,"bundle") ) stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)

  if (lifecycle::is_present(n_cand)) {
    lifecycle::deprecate_warn("2.5.0", "design(n_cand)", "design(n_sample)")
  }
  n_cand <- n_sample

  N <- check_N(N)
  reset <- check_reset(reset, N)
  freq <- check_freq(freq)
  if( !is.null(train_N) ) train_N <- check_train_N(train_N, N)
  reps <- check_reps(reps)
  if( !is.null(refit_cores) ) {
    refit_cores <- as.integer(refit_cores)
    if ( refit_cores < 1 ) stop("'refit_cores' must be >= 1.", call. = FALSE)
  }
  n_cand <- check_n_cand(n_cand)

  X <- object$data$X
  Y <- object$data$Y
  n_dim_X <- ncol(X[[1]])
  n_emulators <- length(object) - 1
  if ( "id" %in% names(object) ) n_emulators <- n_emulators - 1
  if ( "design" %in% names(object) ) n_emulators <- n_emulators - 1

  is.categorical <- rep(FALSE, n_emulators)
  for ( k in 1:n_emulators ){
    obj_k <- object[[paste('emulator',k,sep='')]]
    if (inherits(obj_k,"dgp")){
      if (obj_k$constructor_obj$all_layer[[obj_k$constructor_obj$n_layer]][[1]]$name == "Categorical") {
        is.categorical[k] <- TRUE
      }
    }
  }

  if (!is.null(x_test) & !is.null(y_test)) {
    xy_test <- check_xy_test(x_test, y_test, n_dim_X, n_emulators)
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

  autosave_config <- list(switch = FALSE, directory = './check_points', fname = 'check_point', save_freq = 5, overwrite = FALSE)
  nms_autosave_config <- names(autosave_config)
  autosave_config[(nam_autosave <- names(autosave))] <- autosave
  if (length(noNms_autosave <- nam_autosave[!nam_autosave %in% nms_autosave_config])){
    warning("unknown names in 'autosave': ", paste(noNms_autosave,collapse=", "), call. = FALSE, immediate. = TRUE)
  }

  if (autosave_config$switch) {
    autosave_config$directory <- paste0(autosave_config$directory, '/emulator-', object$id)
    if (dir.exists(autosave_config$directory)){
      if (!first_time & !new_wave){
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave)
      } else {
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave+1)
      }
      if (!dir.exists(autosave_config$directory)){
        dir.create(autosave_config$directory, recursive = TRUE)
      }
    } else {
      if (!first_time & !new_wave){
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave)
      } else {
        autosave_config$directory <- paste0(autosave_config$directory,'/wave_',n_wave+1)
      }
      dir.create(autosave_config$directory, recursive = TRUE)
    }
  }

  if (!first_time & !new_wave){
    start_point <- object[["design"]][[paste('wave', n_wave, sep="")]][["N"]]
    N <- N + start_point
    if ( object[["design"]][[paste('wave', n_wave, sep="")]][["freq"]] != freq[2] ) stop("The frequency of emulator validation must match to that in the last wave if 'new_wave = FALSE'.", call. = FALSE)
  } else {
    start_point <- 0
  }

  for ( k in 1:n_emulators ){
    object[[paste('emulator',k,sep='')]] <- copy_in_design(object[[paste('emulator',k,sep='')]])
  }

  N_acq_ind <- c()
  mnames <- methods::formalArgs(method)
  if ( !is.null(eval) ){
    vnames <- methods::formalArgs(eval)
  }

  #if candidate set is given
  if ( is.null(y_cand) ) {

    if ( is.null(f) ) stop("'f' must be provided.", call. = FALSE)
    fnames <- methods::formalArgs(f)

    if ( "design" %in% names(object) ){
      if ( 'exclusion' %in% names(object$design) ){
        target_points <- object$design$exclusion
      } else {
        target_points <- vector('list', n_emulators)
        target_list_names <- sapply(1:n_emulators, function(i) paste('emulator',i,sep=''))
        names(target_points) <- target_list_names
      }
    } else {
      target_points <- vector('list', n_emulators)
      target_list_names <- sapply(1:n_emulators, function(i) paste('emulator',i,sep=''))
      names(target_points) <- target_list_names
    }

    if ( is.null(x_cand) ) {
      limits <- check_limits(limits, n_dim_X)
      #if ( is.null(limits) ){
      #  all_training_input <- c()
      #  for ( k in 1:n_emulators ){
      #    all_training_input <- rbind(all_training_input, X[[k]])
      #  }
      #  limits <- matrix(0, n_dim_X, 2)
      #  limits[,1] <- pkg.env$np$amin(all_training_input, axis=0L)
      #  limits[,2] <- pkg.env$np$amax(all_training_input, axis=0L)
      #}
      #if (identical(method, pei)) {
      #  ppoints <- list()
      #  for ( k in 1:n_emulators ){
      #    ppoints[[k]] <- pp(X[[k]], limits)
      #  }
      #  add_arg <- list(pseudo_points = ppoints)
      #  add_arg <- utils::modifyList(add_arg, list(...))
      #} else {
        add_arg <- list(...)
        add_arg$batch_size <- batch_size
        add_arg$limits <- limits
      #}
      #int <- check_int(int, n_dim_X)
    } else {
      add_arg <- list(...)
      add_arg$batch_size <- batch_size
      xy_cand_list <- check_xy_cand(x_cand, y_cand, n_dim_X, n_emulators)
      x_cand <- vector('list', n_emulators)
      for (j in 1:n_emulators){
        xy_cand_list_j <- remove_dup( xy_cand_list, rbind(X[[paste('emulator',j,sep="")]], target_points[[paste('emulator',j,sep="")]]) )
        x_cand_origin <- xy_cand_list_j[[1]]
        x_cand_rep <- pkg.env$np$unique(x_cand_origin, return_inverse=TRUE, axis=0L)
        if ( nrow(x_cand_rep[[1]])==nrow(x_cand_origin) ){
          x_cand[[j]] <- x_cand_origin
        } else {
          x_cand[[j]] <- x_cand_rep[[1]]
        }
      }
      idx_x_cand0 <- lapply(1:n_emulators, function(k){c(1:nrow(x_cand[[k]]))})
      idx_x_cand <- idx_x_cand0
      idx_x_acq <- vector('list', n_emulators)
    }

    if (start_point == 0){
      if ( verb ) message("Initializing ...", appendLF = FALSE)
      if ( is.null(eval) ){
        rmse <- c()
        for ( k in 1:n_emulators ){
          obj_k <- object[[paste('emulator',k,sep='')]]
          if (is.null(x_test) & is.null(y_test)){
            type <- 'loo'
            if ( inherits(obj_k,"gp") ) {
              obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val)
              object[[paste('emulator',k,sep='')]] <- obj_k
              rmse <- c(rmse, obj_k$loo$rmse)
            }
            if ( inherits(obj_k,"dgp") ) {
              obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, cores = cores, ...)
              object[[paste('emulator',k,sep='')]] <- obj_k
              if (is.categorical[k]) {
                rmse <- c(rmse, obj_k$loo$log_loss)
              } else {
                rmse <- c(rmse, obj_k$loo$rmse)
              }
            }
          } else {
            type <- 'oos'
            if ( inherits(obj_k,"gp") ) {
              obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val)
              object[[paste('emulator',k,sep='')]] <- obj_k
              rmse <- c(rmse, obj_k$oos$rmse)
            }
            if ( inherits(obj_k,"dgp") ) {
              obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, cores = cores, ...)
              object[[paste('emulator',k,sep='')]] <- obj_k
              if (is.categorical[k]) {
                rmse <- c(rmse, obj_k$oos$log_loss)
              } else {
                rmse <- c(rmse, obj_k$oos$rmse)
              }
            }
          }
        }
      } else {
        type <- 'customized'
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
        if ( verb ) {
          metric_labels <- ifelse(is.categorical, "Log-Loss", "RMSE")

          # Create combined strings for each label and corresponding rmse value
          formatted_pairs <- mapply(function(label, value) {
            paste(label, sprintf("%.06f", value), sep = ": ")
          }, metric_labels, rmse)

          # Combine all pairs into a single message
          message_text <- paste(" * ", paste(formatted_pairs, collapse = ", "))

          # Print the message
          message(message_text)
        }
      } else {
        if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
      }
      rmse_records <- rmse
    } else {
      if ( is.null(eval) ){
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]][n_rmse,]
        rmse_records <- c()
      } else {
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]][n_rmse,]
        rmse_records <- c()
      }
    }

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
      for ( i in (1+start_point):N ){
        if ( verb ) {
          message(sprintf("Iteration %i:", i))
          message(" - Locating ...", appendLF = FALSE)
        }

        tryagain <- TRUE

        while(tryagain){

          if ( !is.null(x_cand) ){
            #sub_cand <- lhs::maximinLHS(n_cand,n_dim_X)
            #sub_cand <- replicate(n_emulators, sub_cand, simplify = FALSE)

            #for (j in 1:n_emulators){
            #  if( !is.null(target_points[[j]]) ){
            #    eud_dist <- sqrt(t(t(rowSums(sub_cand[[j]]^2)-2*sub_cand[[j]]%*%t(target_points[[j]]))+rowSums(target_points[[j]]^2)))
            #    rows_to_keep <- rowSums(eud_dist > 0.01) == ncol(eud_dist)
            #    sub_cand[[j]] <- sub_cand[[j]][rows_to_keep,,drop = F]
            #  }
            #}
            #sub_cand_lhd <- sub_cand

            #for (k in 1:n_emulators){
            #  for (j in 1:n_dim_X){
            #    if ( int[j] ){
            #      #if ( !is.integer(limits[j,2])|!is.integer(limits[j,1]) ) stop(sprintf("The upper and lower limits specified for the input dimension %i should be intgers.", j), call. = FALSE)
            #      sub_cand[[k]][,j] <- floor( sub_cand[[k]][,j]*(limits[j,2]-limits[j,1]+1) ) + limits[j,1]
            #    } else {
            #      sub_cand[[k]][,j] <- sub_cand[[k]][,j]*(limits[j,2]-limits[j,1]) + limits[j,1]
            #    }
            #  }
            #}
          #} else {
            sub_cand <- vector('list', n_emulators)
            idx_sub_cand <- vector('list', n_emulators)
            for (j in 1:n_emulators){
              if ( n_cand<=length(idx_x_cand[[j]]) ){
                if (n_dim_X == 1){
                  idx_sub_cand[[j]] <- sample(idx_x_cand[[j]], n_cand, replace = FALSE)
                } else {
                  idx_idx <- suppressPackageStartupMessages(clhs::clhs(as.data.frame(x_cand[[j]][idx_x_cand[[j]],,drop = FALSE]), size = n_cand, progress = FALSE, simple = TRUE))
                  idx_sub_cand[[j]] <- idx_x_cand[[j]][idx_idx]
                }
              } else {
                idx_sub_cand[[j]] <- idx_x_cand[[j]]
              }
              sub_cand[[j]] <- x_cand[[j]][idx_sub_cand[[j]],,drop = FALSE]
            }
          }

          if ( !is.null(x_cand) ){
            arg_list <- list(object, sub_cand)
          } else {
            arg_list <- list(object)
          }

          if ("..." %in% mnames){
            res <- do.call(method, c(arg_list, add_arg))
          } else {
            midx <- mnames %in% names(add_arg)
            mparam <- add_arg[mnames[midx]]
            res <- do.call(method, c(arg_list, mparam))
          }

          if ( verb ) {
            message(" done")
            if ( !is.null(x_cand) ){
              if ( nrow(res)==1 ){
                for ( j in 1:n_emulators ){
                  if ( is.null(target) ){
                    message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", sub_cand[[j]][res[1,j],])), collapse=" "))
                  } else {
                    if ( !istarget[j] ){
                      message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", sub_cand[[j]][res[1,j],])), collapse=" "))
                    } else {
                      message(sprintf(" * Next design point (Emulator%i): None (target reached)", j))
                    }
                  }
                }
              } else {
                for ( j in 1:nrow(res) ){
                  for ( k in 1:n_emulators ){
                    if ( is.null(target) ){
                      message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", sub_cand[[k]][res[j,k],])), collapse=" "))
                    } else {
                      if ( !istarget[k] ){
                        message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", sub_cand[[k]][res[j,k],])), collapse=" "))
                      } else {
                        message(sprintf(" * Next design point (Position%i for Emulator%i): None (target reached)", j, k))
                      }
                    }
                  }
                }
              }
            } else {
              for ( j in 1:n_emulators ){
                if ( is.null(target) ){
                  if (nrow(res[[j]])==1){
                    message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", res[[j]][1,])), collapse=" "))
                  } else {
                    for ( k in 1:nrow(res[[j]]) ){
                      message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", k, j), sprintf("%.06f", res[[j]][k,])), collapse=" "))
                    }
                  }
                } else {
                  if ( !istarget[j] ){
                    if (nrow(res[[j]])==1){
                      message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", res[[j]][1,])), collapse=" "))
                    } else {
                      for ( k in 1:nrow(res[[j]]) ){
                        message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", k, j), sprintf("%.06f", res[[j]][k,])), collapse=" "))
                      }
                    }
                  } else {
                    if (nrow(res[[j]])==1){
                      message(sprintf(" * Next design point (Emulator%i): None (target reached)", j))
                    } else {
                      for ( k in 1:nrow(res[[j]]) ){
                        message(sprintf(" * Next design point (Position%i for Emulator%i): None (target reached)", k, j))
                      }
                    }
                  }
                }
              }
            }
          }

          if (!is.null(x_cand)){
            if ( nrow(res)==1 ){
              single_batch <-  TRUE
            } else {
              single_batch <-  FALSE
            }
          } else {
            if (nrow(res[[1]])==1) {
              single_batch <-  TRUE
            } else {
              single_batch <-  FALSE
            }
          }

          if ( single_batch ){
            if (!is.null(x_cand)){
              new_X <- do.call(rbind, lapply(1:n_emulators, function(i) {
                sub_cand[[i]][res[i],,drop = FALSE]
              }))
            } else {
              new_X <- do.call(rbind, res)
            }

            if ( !is.null(target) ){
              new_X <- new_X[!istarget,,drop=FALSE]
            }
            rep_new_X <- pkg.env$np$unique(new_X, return_inverse=TRUE, axis=0L)
            new_X_unique <- matrix( rep( t( rep_new_X[[1]] ) , reps ) , ncol = ncol(rep_new_X[[1]]) , byrow = TRUE )
            rep_idx0 <- rep_new_X[[2]] + 1
            rep_idx <- c()
            for (j in 1:reps){
              rep_idx <- c(rep_idx, rep_idx0+(j-1)*nrow(rep_new_X[[1]]))
            }
            if ("..." %in% fnames){
              new_output_unique <- do.call(f, c(list(new_X_unique), add_arg))
            } else {
              fidx <- fnames %in% names(add_arg)
              fparam <- add_arg[fnames[fidx]]
              new_output_unique <- do.call(f, c(list(new_X_unique), fparam))
            }
          } else {
            new_X <- c()
            if (!is.null(x_cand)){
              for (j in 1:nrow(res) ){
                if ( is.null(target) ) {
                  new_X <- rbind(new_X, do.call(rbind, lapply(1:n_emulators, function(i) {
                    sub_cand[[i]][res[j,i],,drop = FALSE]
                  })))
                } else {
                  new_X <- rbind(new_X, do.call(rbind, lapply(1:n_emulators, function(i) {
                    sub_cand[[i]][res[j,i],,drop = FALSE]
                  }))[!istarget,,drop=FALSE])
                }
              }
            } else {
              for (j in 1:nrow(res[[1]]) ){
                if ( is.null(target) ) {
                  new_X <- rbind(new_X, do.call(rbind, lapply(1:n_emulators, function(i) {
                    res[[i]][j,,drop = FALSE]
                  })))
                } else {
                  new_X <- rbind(new_X, do.call(rbind, lapply(1:n_emulators, function(i) {
                    res[[i]][j,,drop = FALSE]
                  }))[!istarget,,drop=FALSE])
                }
              }
            }
            rep_new_X <- pkg.env$np$unique(new_X, return_inverse=TRUE, axis=0L)
            new_X_unique <- matrix( rep( t( rep_new_X[[1]] ) , reps ) , ncol = ncol(rep_new_X[[1]]) , byrow = TRUE )
            rep_idx0 <- rep_new_X[[2]] + 1
            rep_idx <- c()
            for (j in 1:reps){
              rep_idx <- c(rep_idx, rep_idx0+(j-1)*nrow(rep_new_X[[1]]))
            }

            if ("..." %in% fnames){
              new_output_unique <- do.call(f, c(list(new_X_unique), add_arg))
            } else {
              fidx <- fnames %in% names(add_arg)
              fparam <- add_arg[fnames[fidx]]
              new_output_unique <- do.call(f, c(list(new_X_unique), fparam))
            }
          }

          if (ifelse(is.list(new_output_unique), all(!is.na(new_output_unique[[1]])), all(!is.na(new_output_unique)))) {
            tryagain <- FALSE
          } else if ( ifelse(is.list(new_output_unique), all(is.na(new_output_unique[[1]])), all(is.na(new_output_unique)) ))  {
            if ( !is.null(x_cand) ){
              #exclude_lhd_input <- c()
              #for (j in 1:nrow(res) ){
              #  if ( is.null(target) ) {
              #    exclude_lhd_input <- rbind(exclude_lhd_input, do.call(rbind, lapply(1:n_emulators, function(i) {
              #      sub_cand_lhd[[i]][res[j,i],,drop = FALSE]
              #    })))
              #  } else {
              #    exclude_lhd_input <- rbind(exclude_lhd_input, do.call(rbind, lapply(1:n_emulators, function(i) {
              #      sub_cand_lhd[[i]][res[j,i],,drop = FALSE]
              #    }))[!istarget,,drop=FALSE])
              #  }
              #}
              #unique_exclude_lhd_input <- pkg.env$np$unique(exclude_lhd_input, return_inverse=TRUE, axis=0L)
              #target_points <- lapply(target_points, function(x) rbind(x, unique_exclude_lhd_input[[1]]))
            #} else {
              target_points <- lapply(target_points, function(x) rbind(x, rep_new_X))

              if ( is.null(target) ){
                idx_x_acq <- Map(c, idx_x_acq, lapply(x_cand, find_matching_indices, mat2 = rep_new_X))
              } else {
                idx_x_acq[!istarget] <- Map(c, idx_x_acq[!istarget], lapply(x_cand[!istarget], find_matching_indices, mat2 = rep_new_X))
              }

              idx_x_cand[!istarget] <- lapply((1:n_emulators)[!istarget], function(k) idx_x_cand0[[k]][-idx_x_acq[[k]]])
            }
            if ( verb ) {
              message(" - Re-Locating ...", appendLF = FALSE)
            }
          } else {
            if ( is.list(new_output_unique) ){
              keep_idx <- !is.na(new_output_unique[[1]])
              keep_X <- lapply(1:n_emulators, function(k) if (all(!keep_idx[,k])) {NULL} else {new_X_unique[keep_idx[,k],,drop=F]})
            } else {
              keep_idx <- !is.na(new_output_unique)
              keep_X <- lapply(1:n_emulators, function(k) if (all(!keep_idx[,k])) {NULL} else {new_X_unique[keep_idx[,k],,drop=F]})
            }
            if ( !is.null(x_cand) ){
            #  exclude_points <- lapply(keep_X, function(mat) if (is.null(mat)) {new_X_unique} else {unname(as.matrix(dplyr::setdiff(as.data.frame(new_X_unique), as.data.frame(mat))))})
            #  exclude_points_lhd <- lapply(1:n_emulators, function(k) if (length(exclude_points[[k]])==0) {exclude_points[[k]]} else {sub_cand_lhd[[k]][find_matching_indices(sub_cand[[k]], exclude_points[[k]]),,drop=F]})
            #  target_points <- lapply(1:n_emulators, function(k) rbind(target_points[[k]], exclude_points_lhd[[k]]))
            #} else {
              exclude_points <- lapply(keep_X, function(mat) if (is.null(mat)) {new_X_unique} else {unname(as.matrix(dplyr::setdiff(as.data.frame(new_X_unique), as.data.frame(mat))))})
              target_points <- lapply(1:n_emulators, function(k) rbind(target_points[[k]], exclude_points[[k]]))
            }
            tryagain <- FALSE
          }

        }

        if ( single_batch ){
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

          temp_before <- unlist(lapply(1:n_emulators, function(k) nrow(X[[paste('emulator',k,sep="")]])))

          if ( is.null(target) ){
            new_X_list <- lapply(1:n_emulators, function(k) matrix(rep(new_X[k, ], reps), nrow = reps, byrow = TRUE))
            new_Y_list <- lapply(1:n_emulators, function(k) new_Y[seq(k, by = n_emulators, length = reps), k, drop = FALSE] )
            if ( !is.null(x_cand) ) idx_x_acq <- Map(c, idx_x_acq, lapply(Map(find_matching_indices, x_cand, new_X_list), unique))
          } else {
            new_X_list <- vector('list', n_emulators)
            new_Y_list <- vector('list', n_emulators)
            ct <- 1
            for (k in 1:n_emulators){
              if (!istarget[k]){
                new_X_list[[k]] <- matrix(rep(new_X[ct, ], reps), nrow = reps, byrow = TRUE)
                new_Y_list[[k]] <- new_Y[seq(ct, by = sum(!istarget), length = reps), k, drop = FALSE]
                ct <- ct + 1
              }
            }
            if ( !is.null(x_cand) ) idx_x_acq[!istarget] <- Map(c, idx_x_acq[!istarget], lapply(Map(find_matching_indices, x_cand[!istarget], new_X_list[!istarget]), unique))
          }
          if ( !is.null(x_cand) ) idx_x_cand[!istarget] <- lapply((1:n_emulators)[!istarget], function(k) idx_x_cand0[[k]][-idx_x_acq[[k]]])

          for ( j in 1:n_emulators ){
            if ( !is.null(new_Y_list[[j]]) ){
              na.idx <- is.na(new_Y_list[[j]])
              if ( !all(na.idx) ){
                X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X_list[[j]][!na.idx,,drop=FALSE])
                Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y_list[[j]][!na.idx,,drop=FALSE])
              }
            }
          }

          temp_after <- unlist(lapply(1:n_emulators, function(k) nrow(X[[paste('emulator',k,sep="")]])))
          N_acq_ind <- rbind(N_acq_ind, temp_after-temp_before)
        } else {
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

          temp_before <- unlist(lapply(1:n_emulators, function(k) nrow(X[[paste('emulator',k,sep="")]])))

          if ( is.null(target) ){
            new_X_list <- vector('list', n_emulators)
            new_Y_list <- vector('list', n_emulators)
            for (k in 1:n_emulators){
              extract0 <- k + seq(0, by = n_emulators, length = nrow(res))
              extract <- rep(extract0, reps) +  rep(seq(0, by = nrow(new_X), length.out = reps), each = length(extract0))
              extracted_X <- new_X[extract0,,drop=FALSE]
              new_X_list[[k]] <- matrix( rep( t( extracted_X ) , reps ) , ncol = ncol(extracted_X) , byrow = TRUE )
              new_Y_list[[k]] <- new_Y[extract,k,drop=FALSE]
            }
            if ( !is.null(x_cand) ) idx_x_acq <- Map(c, idx_x_acq, lapply(Map(find_matching_indices, x_cand, new_X_list), unique))
          } else {
            ct <- 1
            active_emu <- sum(!istarget)
            new_X_list <- vector('list', n_emulators)
            new_Y_list <- vector('list', n_emulators)
            for (k in 1:n_emulators){
              if (!istarget[k]){
                extract0 <- ct + seq(0, by = active_emu, length = nrow(res))
                extract <- rep(extract0, reps) +  rep(seq(0, by = nrow(new_X), length.out = reps), each = length(extract0))
                extracted_X <- new_X[extract0,,drop=FALSE]
                new_X_list[[k]] <- matrix( rep( t( extracted_X ) , reps ) , ncol = ncol(extracted_X) , byrow = TRUE )
                new_Y_list[[k]] <- new_Y[extract,k,drop=FALSE]
                ct <- ct + 1
              }
            }

            if ( !is.null(x_cand) ) idx_x_acq[!istarget] <- Map(c, idx_x_acq[!istarget], lapply(Map(find_matching_indices, x_cand[!istarget], new_X_list[!istarget]), unique))
          }
          if ( !is.null(x_cand) ) idx_x_cand[!istarget] <- lapply((1:n_emulators)[!istarget], function(k) idx_x_cand0[[k]][-idx_x_acq[[k]]])

          for ( j in 1:n_emulators ){
            if ( !is.null(new_Y_list[[j]]) ){
              na.idx <- is.na(new_Y_list[[j]])
              if ( !all(na.idx) ){
                X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X_list[[j]][!na.idx,,drop=FALSE])
                Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y_list[[j]][!na.idx,,drop=FALSE])
              }
            }
          }

          temp_after <- unlist(lapply(1:n_emulators, function(k) nrow(X[[paste('emulator',k,sep="")]])))
          N_acq_ind <- rbind(N_acq_ind, temp_after-temp_before)
        }

        if ( i %% freq[1]==0 ){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            if (N_acq_ind[nrow(N_acq_ind),k]!=0){
              obj_k <- object[[paste('emulator',k,sep='')]]
              if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
              if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, reset = reset[i-start_point], verb = FALSE, N = if(!is.null(train_N)) train_N[i-start_point] else NULL, cores = refit_cores, B = ifelse(length(obj_k$emulator_obj$all_layer_set)<10, length(obj_k$emulator_obj$all_layer_set), 10), update_in_design = NULL)
              object[[paste('emulator',k,sep='')]] <- obj_k
            }
          }
          if ( verb ) message(" done")
        } else {
          if ( verb ) message(" - Updating ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            if ( N_acq_ind[nrow(N_acq_ind),k]!=0 ) {
              obj_k <- object[[paste('emulator',k,sep='')]]
              if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = FALSE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
              if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = FALSE, reset = reset[i-start_point], verb = FALSE, B = ifelse(length(obj_k$emulator_obj$all_layer_set)<10, length(obj_k$emulator_obj$all_layer_set), 10), update_in_design = NULL)
              object[[paste('emulator',k,sep='')]] <- obj_k
            }
          }
          if ( verb ) message(" done")
        }

        if ( i %% freq[2]==0 | i==N){
          if ( is.null(eval) ){
            if (is.null(x_test) & is.null(y_test)){
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              for ( k in 1:n_emulators ){
                if ( N_acq_ind[nrow(N_acq_ind),k]!=0 ) {
                  obj_k <- object[[paste('emulator',k,sep='')]]
                  if ( inherits(obj_k,"gp") ) {
                    obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    rmse[k] <- obj_k$loo$rmse
                  }
                  if ( inherits(obj_k,"dgp") ) {
                    obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    if (is.categorical[k]){
                      rmse[k] <- obj_k$loo$log_loss
                    } else {
                      rmse[k] <- obj_k$loo$rmse
                    }
                  }
                }
              }
              if ( verb ) message(" done")
              if ( verb ) {
                metric_labels <- ifelse(is.categorical, "Log-Loss", "RMSE")

                # Create combined strings for each label and corresponding rmse value
                formatted_pairs <- mapply(function(label, value) {
                  paste(label, sprintf("%.06f", value), sep = ": ")
                }, metric_labels, rmse)

                # Combine all pairs into a single message
                message_text <- paste(" * ", paste(formatted_pairs, collapse = ", "))

                # Print the message
                message(message_text)
              }
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              for ( k in 1:n_emulators ){
                if ( N_acq_ind[nrow(N_acq_ind),k]!=0 ) {
                  obj_k <- object[[paste('emulator',k,sep='')]]
                  if ( inherits(obj_k,"gp") ) {
                    obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, force = TRUE)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    rmse[k] <- obj_k$oos$rmse
                  }
                  if ( inherits(obj_k,"dgp") ) {
                    obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    if (is.categorical[k]){
                      rmse[k] <- obj_k$oos$log_loss
                    } else {
                      rmse[k] <- obj_k$oos$rmse
                    }
                  }
                }
              }
              if ( verb ) message(" done")
              if ( verb ) {
                metric_labels <- ifelse(is.categorical, "Log-Loss", "RMSE")

                # Create combined strings for each label and corresponding rmse value
                formatted_pairs <- mapply(function(label, value) {
                  paste(label, sprintf("%.06f", value), sep = ": ")
                }, metric_labels, rmse)

                # Combine all pairs into a single message
                message_text <- paste(" * ", paste(formatted_pairs, collapse = ", "))

                # Print the message
                message(message_text)
              }
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
            if( length(rmse_tmp)==1 ){
              rmse <- rmse_tmp
            } else {
              idx_rmse <- N_acq_ind[nrow(N_acq_ind),]!=0
              rmse[idx_rmse] <- rmse_tmp[idx_rmse]
            }
            if ( verb ) message(" done")
            if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
          }
          rmse_records <- rbind(rmse_records, rmse)
        }

        if( !is.null(target) ) {
          if ( all(rmse<=target) ) {
            message(sprintf("Target reached! Sequential design stopped at step %i.", i))
            istarget <- rep(TRUE, n_emulators)
            N <- i
            break
          } else {
            istarget <- rmse<=target
            if ( length(istarget)==1 ) istarget <- rep(istarget, n_emulators)
          }
        }

        if (autosave_config$switch && i %% autosave_config$save_freq == 0) {
          if (i != N) {
            if ( verb ) message(" - Auto-saving ...", appendLF = FALSE)
            save_file_name <- if (autosave_config$overwrite) {
              autosave_config$fname
            } else {
              paste0(autosave_config$fname, "_", i)
            }
            object_temp <- pack_bundle(object, first_time, i, eval, new_wave, rmse_records, freq, N_acq_ind, target, type, X, Y, y_cand,
                                       n_wave = if(isFALSE(first_time)) n_wave else NULL,
                                       design_info = if(isFALSE(first_time)) design_info else NULL,
                                       x_test = if(identical(type, 'oos')) x_test else NULL,
                                       y_test = if(identical(type, 'oos')) y_test else NULL,
                                       istarget = if(!is.null(target)) istarget else NULL,
                                       target_points = if(is.null(y_cand)) target_points else NULL)
            save_file_path <- file.path(autosave_config$directory, save_file_name)
            write(object_temp, save_file_path)
            if ( verb ) message(" done")
          }
        }

      }
    }
  } else {
    add_arg <- list(...)
    add_arg$batch_size <- batch_size
    xy_cand_list <- check_xy_cand(x_cand, y_cand, n_dim_X, n_emulators)

    x_cand <- vector('list', n_emulators)
    x_cand_origin <- vector('list', n_emulators)
    y_cand <- vector('list', n_emulators)
    is.dup <- vector('list', n_emulators)
    for (j in 1:n_emulators){
      xy_cand_list_j <- remove_dup( xy_cand_list, X[[paste('emulator',j,sep="")]] )
      y_cand[[j]] <- xy_cand_list_j[[2]]
      x_cand_origin[[j]] <- xy_cand_list_j[[1]]
      x_cand_rep <- pkg.env$np$unique(x_cand_origin[[j]], return_inverse=TRUE, axis=0L)
      if ( nrow(x_cand_rep[[1]])==nrow(x_cand_origin[[j]]) ){
        x_cand[[j]] <- x_cand_origin[[j]]
        is.dup[[j]] <- FALSE
      } else {
        x_cand[[j]] <- x_cand_rep[[1]]
        is.dup[[j]] <- TRUE
      }
    }
    idx_x_cand0 <- lapply(1:n_emulators, function(k){c(1:nrow(x_cand[[k]]))})
    idx_x_cand <- idx_x_cand0
    idx_x_acq <- vector('list', n_emulators)

    if (start_point == 0){
      if ( verb ) message("Initializing ...", appendLF = FALSE)
      if ( is.null(eval) ){
        rmse <- c()
        for ( k in 1:n_emulators ){
          obj_k <- object[[paste('emulator',k,sep='')]]
          if (is.null(x_test) & is.null(y_test)){
            type <- 'loo'
            if ( inherits(obj_k,"gp") ) {
              obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val)
              object[[paste('emulator',k,sep='')]] <- obj_k
              rmse <- c(rmse, obj_k$loo$rmse)
            }
            if ( inherits(obj_k,"dgp") ) {
              obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, cores = cores, ...)
              object[[paste('emulator',k,sep='')]] <- obj_k
              if (is.categorical[k]){
                rmse <- c(rmse, obj_k$loo$log_loss)
              } else {
                rmse <- c(rmse, obj_k$loo$rmse)
              }
            }
          } else {
            type <- 'oos'
            if ( inherits(obj_k,"gp") ) {
              obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val)
              object[[paste('emulator',k,sep='')]] <- obj_k
              rmse <- c(rmse, obj_k$oos$rmse)
            }
            if ( inherits(obj_k,"dgp") ) {
              obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, cores = cores, ...)
              object[[paste('emulator',k,sep='')]] <- obj_k
              if (is.categorical[k]){
                rmse <- c(rmse, obj_k$oos$log_loss)
              } else {
                rmse <- c(rmse, obj_k$oos$rmse)
              }
            }
          }
        }
      } else {
        type <- 'customized'
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
        if ( verb ){
          metric_labels <- ifelse(is.categorical, "Log-Loss", "RMSE")

          # Create combined strings for each label and corresponding rmse value
          formatted_pairs <- mapply(function(label, value) {
            paste(label, sprintf("%.06f", value), sep = ": ")
          }, metric_labels, rmse)

          # Combine all pairs into a single message
          message_text <- paste(" * ", paste(formatted_pairs, collapse = ", "))

          # Print the message
          message(message_text)
        }
      } else {
        if ( verb ) message(paste(c(" * Metric:", sprintf("%.06f", rmse)), collapse=" "))
      }
      rmse_records <- rmse
    } else {
      if ( is.null(eval) ){
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["rmse"]][n_rmse,]
        rmse_records <- c()
      } else {
        n_rmse <- nrow(object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]])
        rmse <- object[["design"]][[paste('wave', n_wave, sep="")]][["metric"]][n_rmse,]
        rmse_records <- c()
      }
    }

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
      for ( i in (1+start_point):N ){
        sub_cand <- vector('list', n_emulators)
        idx_sub_cand <- vector('list', n_emulators)
        for (j in 1:n_emulators){
          if ( n_cand<=length(idx_x_cand[[j]]) ){
            if (n_dim_X == 1){
              idx_sub_cand[[j]] <- sample(idx_x_cand[[j]], n_cand, replace = FALSE)
            } else {
              idx_idx <- suppressPackageStartupMessages(clhs::clhs(as.data.frame(x_cand[[j]][idx_x_cand[[j]],,drop = FALSE]), size = n_cand, progress = FALSE, simple = TRUE))
              idx_sub_cand[[j]] <- idx_x_cand[[j]][idx_idx]
            }
          } else {
            idx_sub_cand[[j]] <- idx_x_cand[[j]]
          }
          sub_cand[[j]] <- x_cand[[j]][idx_sub_cand[[j]],,drop = FALSE]
        }


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
                message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", sub_cand[[j]][res[1,j],])), collapse=" "))
              } else {
                if ( !istarget[j] ){
                  message(paste(c(sprintf(" * Next design point (Emulator%i):", j), sprintf("%.06f", sub_cand[[j]][res[1,j],])), collapse=" "))
                } else {
                  message(sprintf(" * Next design point (Emulator%i): None (target reached)", j))
                }
              }
            }
          } else {
            for ( j in 1:nrow(res) ){
              for ( k in 1:n_emulators ){
                if ( is.null(target) ){
                  message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", sub_cand[[j]][res[j,k],])), collapse=" "))
                } else {
                  if ( !istarget[k] ){
                    message(paste(c(sprintf(" * Next design point (Position%i for Emulator%i):", j, k), sprintf("%.06f", sub_cand[[j]][res[j,k],])), collapse=" "))
                  } else {
                    message(sprintf(" * Next design point (Position%i for Emulator%i): None (target reached)", j, k))
                  }
                }
              }
            }
          }
        }

        if ( nrow(res)==1 ){
          temp_before <- unlist(lapply(1:n_emulators, function(k) nrow(X[[paste('emulator',k,sep="")]])))
          for ( j in 1:n_emulators ){
            if ( is.dup[[j]] ){
              new_X_temp <- sub_cand[[j]][res[1,j],,drop=FALSE]
              new_idx <- extract_all(new_X_temp, x_cand_origin[[j]])
              new_X <- x_cand_origin[[j]][new_idx,,drop=FALSE]
              new_Y <- y_cand[[j]][new_idx,,drop=FALSE]
            } else {
              new_idx <- res[1,j]
              new_X <- sub_cand[[j]][new_idx,,drop=FALSE]
              new_Y <- y_cand[[j]][idx_sub_cand[[j]],,drop=FALSE][new_idx,,drop=FALSE]
            }
            if ( is.null(target) ){
              X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X)
              Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[,j,drop=FALSE])
            } else {
              if ( !istarget[j] ) {
                X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X)
                Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[,j,drop=FALSE])
              }
            }

            if ( is.null(target) ){
              idx_x_acq[[j]] <- c(idx_x_acq[[j]], idx_sub_cand[[j]][res[1,j]])
              idx_x_cand[[j]] <- idx_x_cand0[[j]][-idx_x_acq[[j]]]
            } else {
              if ( !istarget[j] ) {
                idx_x_acq[[j]] <- c(idx_x_acq[[j]], idx_sub_cand[[j]][res[1,j]])
                idx_x_cand[[j]] <- idx_x_cand0[[j]][-idx_x_acq[[j]]]
              }
            }
          }
          temp_after <- unlist(lapply(1:n_emulators, function(k) nrow(X[[paste('emulator',k,sep="")]])))
          N_acq_ind <- rbind(N_acq_ind, temp_after-temp_before)
        } else {
          temp_before <- unlist(lapply(1:n_emulators, function(k) nrow(X[[paste('emulator',k,sep="")]])))
          for ( j in 1:n_emulators ){
            if ( is.dup[[j]] ){
              new_X_temp <- sub_cand[[j]][res[,j],,drop=FALSE]
              new_idx <- extract_all(new_X_temp, x_cand_origin[[j]])
              new_X <- x_cand_origin[[j]][new_idx,,drop=FALSE]
              new_Y <- y_cand[[j]][new_idx,,drop=FALSE]
            } else {
              new_idx <- res[,j]
              new_X <- sub_cand[[j]][new_idx,,drop=FALSE]
              new_Y <- y_cand[[j]][idx_sub_cand[[j]],,drop=FALSE][new_idx,,drop=FALSE]
            }
            if ( is.null(target) ){
              X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X)
              Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[,j,drop=FALSE])
            } else {
              if ( !istarget[j] ) {
                X[[paste('emulator',j,sep="")]] <- rbind(X[[paste('emulator',j,sep="")]], new_X)
                Y[[paste('emulator',j,sep="")]] <- rbind(Y[[paste('emulator',j,sep="")]], new_Y[,j,drop=FALSE])
              }
            }

            if ( is.null(target) ){
              idx_x_acq[[j]] <- c(idx_x_acq[[j]], idx_sub_cand[[j]][res[,j]])
              idx_x_cand[[j]] <- idx_x_cand0[[j]][-idx_x_acq[[j]]]
            } else {
              if ( !istarget[j] ) {
                idx_x_acq[[j]] <- c(idx_x_acq[[j]], idx_sub_cand[[j]][res[,j]])
                idx_x_cand[[j]] <- idx_x_cand0[[j]][-idx_x_acq[[j]]]
              }
            }
          }
          temp_after <- unlist(lapply(1:n_emulators, function(k) nrow(X[[paste('emulator',k,sep="")]])))
          N_acq_ind <- rbind(N_acq_ind, temp_after-temp_before)
        }

        if ( i %% freq[1]==0 ){
          if ( verb ) message(" - Updating and re-fitting ...", appendLF = FALSE)
          for ( k in 1:n_emulators ){
            if ( is.null(target) ){
              obj_k <- object[[paste('emulator',k,sep='')]]
              if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
              if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, reset = reset[i-start_point], verb = FALSE, N = if(!is.null(train_N)) train_N[i-start_point] else NULL, cores = refit_cores, B = ifelse(length(obj_k$emulator_obj$all_layer_set)<10, length(obj_k$emulator_obj$all_layer_set), 10), update_in_design = NULL)
              object[[paste('emulator',k,sep='')]] <- obj_k
            } else {
              if ( !istarget[k] ){
                obj_k <- object[[paste('emulator',k,sep='')]]
                if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
                if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X[[paste('emulator',k,sep="")]], Y[[paste('emulator',k,sep="")]], refit = TRUE, reset = reset[i-start_point], verb = FALSE, N = if(!is.null(train_N)) train_N[i-start_point] else NULL, cores = refit_cores, B = ifelse(length(obj_k$emulator_obj$all_layer_set)<10, length(obj_k$emulator_obj$all_layer_set), 10), update_in_design = NULL)
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
              if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
              if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, reset = reset[i-start_point], verb = FALSE, B = ifelse(length(obj_k$emulator_obj$all_layer_set)<10, length(obj_k$emulator_obj$all_layer_set), 10), update_in_design = NULL)
              object[[paste('emulator',k,sep='')]] <- obj_k
            } else {
              if ( !istarget[k] ){
                obj_k <- object[[paste('emulator',k,sep='')]]
                if ( inherits(obj_k,"gp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, reset = reset[i-start_point], verb = FALSE, update_in_design = NULL)
                if ( inherits(obj_k,"dgp") ) obj_k <- update(obj_k, X, Y[,k], refit = FALSE, reset = reset[i-start_point], verb = FALSE, B = ifelse(length(obj_k$emulator_obj$all_layer_set)<10, length(obj_k$emulator_obj$all_layer_set), 10), update_in_design = NULL)
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
                  if ( inherits(obj_k,"gp") ) {
                    obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    rmse[k] <- obj_k$loo$rmse
                  }
                  if ( inherits(obj_k,"dgp") ) {
                    obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    if (is.categorical[k]){
                      rmse[k] <- obj_k$loo$log_loss
                    } else {
                      rmse[k] <- obj_k$loo$rmse
                    }
                  }
                } else {
                  if ( !istarget[k] ){
                    obj_k <- object[[paste('emulator',k,sep='')]]
                    if ( inherits(obj_k,"gp") ) {
                      obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE)
                      object[[paste('emulator',k,sep='')]] <- obj_k
                      rmse[k] <- obj_k$loo$rmse
                    }
                    if ( inherits(obj_k,"dgp") ) {
                      obj_k <- validate(obj_k, x_test = NULL, y_test = NULL, verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                      object[[paste('emulator',k,sep='')]] <- obj_k
                      if (is.categorical[k]){
                        rmse[k] <- obj_k$loo$log_loss
                      } else {
                        rmse[k] <- obj_k$loo$rmse
                      }
                    }
                  }
                }
              }
              if ( verb ) message(" done")
              if ( verb ) {
                metric_labels <- ifelse(is.categorical, "Log-Loss", "RMSE")

                # Create combined strings for each label and corresponding rmse value
                formatted_pairs <- mapply(function(label, value) {
                  paste(label, sprintf("%.06f", value), sep = ": ")
                }, metric_labels, rmse)

                # Combine all pairs into a single message
                message_text <- paste(" * ", paste(formatted_pairs, collapse = ", "))

                # Print the message
                message(message_text)
              }
            } else {
              if ( verb ) message(" - Validating ...", appendLF = FALSE)
              for ( k in 1:n_emulators ){
                if ( is.null(target) ) {
                  obj_k <- object[[paste('emulator',k,sep='')]]
                  if ( inherits(obj_k,"gp") ) {
                    obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, force = TRUE)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    rmse[k] <- obj_k$oos$rmse
                  }
                  if ( inherits(obj_k,"dgp") ) {
                    obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                    object[[paste('emulator',k,sep='')]] <- obj_k
                    if (is.categorical[k]){
                      rmse[k] <- obj_k$oos$log_loss
                    } else {
                      rmse[k] <- obj_k$oos$rmse
                    }
                  }
                } else {
                  if ( !istarget[k] ){
                    obj_k <- object[[paste('emulator',k,sep='')]]
                    if ( inherits(obj_k,"gp") ) {
                      obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, force = TRUE)
                      object[[paste('emulator',k,sep='')]] <- obj_k
                      rmse[k] <- obj_k$oos$rmse
                    }
                    if ( inherits(obj_k,"dgp") ) {
                      obj_k <- validate(obj_k, x_test = x_test, y_test = y_test[,k], verb = FALSE, M = M_val, force = TRUE, cores = cores, ...)
                      object[[paste('emulator',k,sep='')]] <- obj_k
                      if (is.categorical[k]){
                        rmse[k] <- obj_k$oos$log_loss
                      } else {
                        rmse[k] <- obj_k$oos$rmse
                      }
                    }
                  }
                }
              }
              if ( verb ) message(" done")
              if ( verb ) {
                metric_labels <- ifelse(is.categorical, "Log-Loss", "RMSE")

                # Create combined strings for each label and corresponding rmse value
                formatted_pairs <- mapply(function(label, value) {
                  paste(label, sprintf("%.06f", value), sep = ": ")
                }, metric_labels, rmse)

                # Combine all pairs into a single message
                message_text <- paste(" * ", paste(formatted_pairs, collapse = ", "))

                # Print the message
                message(message_text)
              }
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
            message(sprintf("Target reached! Sequential design stopped at step %i.", i))
            istarget <- rep(TRUE, n_emulators)
            N <- i
            break
          } else {
            istarget <- rmse<=target
            if ( length(istarget)==1 ) istarget <- rep(istarget, n_emulators)
          }
        }

        if (autosave_config$switch && i %% autosave_config$save_freq == 0) {
          if (i != N) {
            if ( verb ) message(" - Auto-saving ...", appendLF = FALSE)
            save_file_name <- if (autosave_config$overwrite) {
              autosave_config$fname
            } else {
              paste0(autosave_config$fname, "_", i)
            }
            object_temp <- pack_bundle(object, first_time, i, eval, new_wave, rmse_records, freq, N_acq_ind, target, type, X, Y, y_cand,
                                       n_wave = if(isFALSE(first_time)) n_wave else NULL,
                                       design_info = if(isFALSE(first_time)) design_info else NULL,
                                       x_test = if(identical(type, 'oos')) x_test else NULL,
                                       y_test = if(identical(type, 'oos')) y_test else NULL,
                                       istarget = if(!is.null(target)) istarget else NULL,
                                       target_points = if(is.null(y_cand)) target_points else NULL)
            save_file_path <- file.path(autosave_config$directory, save_file_name)
            write(object_temp, save_file_path)
            if ( verb ) message(" done")
          }
        }

      }
    }
  }

  if ( run ){
    object <- pack_bundle(object, first_time, N, eval, new_wave, rmse_records, freq, N_acq_ind, target, type, X, Y, y_cand,
                          n_wave = if(isFALSE(first_time)) n_wave else NULL,
                          design_info = if(isFALSE(first_time)) design_info else NULL,
                          x_test = if(identical(type, 'oos')) x_test else NULL,
                          y_test = if(identical(type, 'oos')) y_test else NULL,
                          istarget = if(!is.null(target)) istarget else NULL,
                          target_points = if(is.null(y_cand)) target_points else NULL)

    if (autosave_config$switch) {
      if ( verb ) message(" - Auto-saving the final iteration ...", appendLF = FALSE)
      save_file_name <- if (autosave_config$overwrite) {
        autosave_config$fname
      } else {
        paste0(autosave_config$fname, "_", N)
      }
      save_file_path <- file.path(autosave_config$directory, save_file_name)
      write(object, save_file_path)
      if ( verb ) message(" done")
    }

    if( !is.null(target) ) {
      if ( !all(istarget) ) message("Targets not reached for all emulators at the end of the sequential design.")
    }
  } else {
    message("Target already reached. Sequential design not performed.")
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

#check argument reps
check_reps <- function(reps){
  reps <- as.integer(reps)
  if ( reps < 1 ) stop("'reps' must be greater than or equal to 1.", call. = FALSE)
  return(reps)
}

# find match indices between two matrices
find_matching_indices <- function(mat1, mat2) {
  mat1_string <- apply(mat1, 1, paste, collapse = "_")
  mat2_string <- apply(mat2, 1, paste, collapse = "_")
  idx <- match(mat2_string, mat1_string)
  idx <- idx[!is.na(idx)]
  return(idx)
}

#check argument x_cand and y_cand
check_xy_cand <- function(x_cand, y_cand, n_dim_X, n_dim_Y){
  # if ( is.null(y_cand) ) stop("'y_cand' must be provided if 'x_cand' is not NULL.", call. = FALSE)
  if ( !is.matrix(x_cand)&!is.vector(x_cand) ) stop("'x_cand' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x_cand) ) {
    if ( n_dim_X!=1 ){
      x_cand <- matrix(x_cand, nrow = 1)
    } else {
      x_cand <- as.matrix(x_cand)
    }
  }
  if ( !is.null(y_cand) ) {
    if ( !is.matrix(y_cand)&!is.vector(y_cand) ) stop("'y_cand' must be a vector or a matrix.", call. = FALSE)
    if ( is.vector(y_cand) ) {
      if ( n_dim_Y!=1 ){
        y_cand <- matrix(y_cand, nrow = 1)
      } else {
        y_cand <- as.matrix(y_cand)
      }
    }
    if ( nrow(x_cand)!=nrow(y_cand) ) stop("'x_cand' and 'y_cand' have different number of data points.", call. = FALSE)
    if ( ncol(y_cand)!=n_dim_Y ) stop(sprintf("The dimension of 'y_cand' must be %i.", n_dim_Y), call. = FALSE)
  }
  return(list(x_cand, y_cand))
}

#remove duplicates between x_cand and X
remove_dup <- function(xy_cand_list, X){
  x_cand <- xy_cand_list[[1]]
  y_cand <- xy_cand_list[[2]]
  X_s <- apply(X, 1, paste, collapse = ", ")
  x_cand_s <- apply(x_cand, 1, paste, collapse = ", ")
  X_x_cand <- intersect(X_s, x_cand_s)
  if ( length(X_x_cand)!=0 ){
    idx <- !x_cand_s %in% X_x_cand
    x_cand <- x_cand[idx,,drop=FALSE]
    if ( !is.null(y_cand) ) {
      y_cand <- y_cand[idx,,drop=FALSE]
    }
  }
  return(list(x_cand, y_cand))
}

#find all data points in the candidate dataset that have the same input locations as the selected design points
extract_all <- function(X ,x_cand){
  x_cand_s <- apply(x_cand, 1, paste, collapse = ", ")
  X_s <- apply(X, 1, paste, collapse = ", ")
  idx <- x_cand_s %in% X_s
  return(idx)
}

#check argument x_test and y_test
check_xy_test <- function(x_test, y_test, n_dim_X, n_dim_Y){
  x_test <- unname(x_test)
  y_test <- unname(y_test)
  if ( !is.matrix(x_test)&!is.vector(x_test) ) stop("'x_test' must be a vector or a matrix.", call. = FALSE)
  if ( !is.matrix(y_test)&!is.vector(y_test) ) stop("'y_test' must be a vector or a matrix.", call. = FALSE)
  if ( is.vector(x_test) ) {
    if ( n_dim_X!=1 ){
      x_test <- matrix(x_test, nrow = 1)
    } else {
      x_test <- as.matrix(x_test)
    }
  }
  if ( is.vector(y_test) ) {
    if ( n_dim_Y!=1 ){
      y_test <- matrix(y_test, nrow = 1)
    } else {
      y_test <- as.matrix(y_test)
    }
  }
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
  } else {
    stop("'limits' must be provided.", call. = FALSE)
  }
  return(limits)
}

check_int <- function(int, n_dim_X){
  if ( length(int)==1 ) {
    if ( n_dim_X!=1 ) {
      int <- rep(int, n_dim_X)
    }
  } else {
    if ( length(int)!=n_dim_X ) stop("The length of 'int' should be equal to the number of input dimensions.", call. = FALSE)
  }
  return(int)
}

check_reset <- function(reset, N){
  if ( length(reset)==1 ) {
    reset <- rep(reset, N)
  } else {
    if ( length(reset)!=N ) stop("The length of 'reset' should be equal to the number of steps of the sequential design.", call. = FALSE)
  }
  return(reset)
}

check_auto <- function(object){
  auto_pruning <- T
  # exclude user-defined structure
    n_layer <- object$constructor_obj$n_layer
    if (object$constructor_obj$all_layer[[n_layer]][[1]]$type!='gp') {
      n_layer <- n_layer - 1
      if (n_layer == 1) {
        auto_pruning <- F
        return(auto_pruning)
      }
    }
    for (i in 2:n_layer){
      for (ker in object$constructor_obj$all_layer[[i]]){
        if (is.null(ker$global_input)) {
          auto_pruning <- F
          return(auto_pruning)
        }
      }
    }

  return(auto_pruning)
}

create_drop_list <- function(object){
  n_layer <- object$constructor_obj$n_layer
  if (object$constructor_obj$all_layer[[n_layer]][[1]]$type!='gp') {
    n_layer <- n_layer - 1
  }
  drop_list <- vector("list", length = n_layer-1)
  for (i in 1:(n_layer-1)){
    drop_list[[i]] <- rep( 0, length(object$constructor_obj$all_layer[[i]]) )
  }
  return(drop_list)
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
        if ( "exclusion" %in% names(object$design) ) n_wave <- n_wave - 1
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
          if ( "exclusion" %in% names(object$design) ) n_wave <- n_wave - 1
          if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) n_wave <- n_wave - 2
        } else {
          first_time <- TRUE
          n_wave <- 0
        }
      }
    }
  } else {
    if ( identical(object$design$type, 'customized') ) {
      first_time <- FALSE
      n_wave <- length(object$design)
      if ( "type" %in% names(object$design) ) n_wave <- n_wave - 1
      if ( "exclusion" %in% names(object$design) ) n_wave <- n_wave - 1
    } else {
      first_time <- TRUE
      n_wave <- 0
    }
  }
  return(list(first_time, n_wave))
}

minmax <- function(training_input, limits) {
  # Extract min and max values from the limits matrix
  mins <- limits[, 1]
  maxs <- limits[, 2]

  # Apply Min-Max normalization
  normalized <- sweep(training_input, 2, mins, FUN = "-")  # Subtract min values
  normalized <- sweep(normalized, 2, maxs - mins, FUN = "/")  # Divide by (max - min)

  return(normalized)
}

reverse_minmax <- function(normalized_data, limits) {
  # Extract min and max values from the limits matrix
  mins <- limits[, 1]
  maxs <- limits[, 2]

  # Revert Min-Max normalization
  original_data <- sweep(normalized_data, 2, maxs - mins, FUN = "*")  # Multiply by (max - min)
  original_data <- sweep(original_data, 2, mins, FUN = "+")  # Add min values

  return(original_data)
}

#generic_wrapper <- function(r_func) {
#  function(...) {
#    # Capture the arguments
#    args <- list(...)
#
#    # Convert Python-native arguments to R-native if necessary
#    args <- lapply(args, function(arg) {
#      if (inherits(arg, "python.builtin.object")) {
#        reticulate::py_to_r(arg)
#      } else {
#        arg
#      }
#    })
#
#    # Call the user-provided R function with converted arguments
#    result <- do.call(r_func, args)
#
#    # Convert the result back to Python-native types
#    reticulate::r_to_py(result)
#  }
#}
