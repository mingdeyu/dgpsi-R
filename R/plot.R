#' @title Validation plots of a constructed GP, DGP, or linked (D)GP emulator
#'
#' @description This function draws validation plots of a GP, DGP, or linked (D)GP emulator.
#'
#' @param x can be one of the following emulator classes:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `lgp`.
#' @param x_test same as that of [validate()].
#' @param y_test same as that of [validate()].
#' @param dim if `dim = NULL`, the index of an emulator's input will be shown on the x-axis in validation plots. Otherwise, `dim` indicates
#'     which dimension of an emulator's input will be shown on the x-axis in validation plots:
#' * If `x` is an instance of the `gp` of `dgp` class, `dim` is an integer.
#' * If `x` is an instance of the `lgp` class, `dim` can be
#'   1. an integer referring to the dimension of the global input to emulators in the first layer of a linked emulator system; or
#'   2. a vector of three integers referring to the dimension (specified by the third integer) of the global input to an emulator
#'   (specified by the second integer) in a layer (specified by the first integer) that is not the first layer of a linked emulator
#'   system.
#'
#' This argument is only used when `style = 1` and the emulator input is at least two-dimensional. Defaults to `NULL`.
#' @param method same as that of [validate()].
#' @param style either `1` or `2`, indicating two different types of validation plots.
#' @param min_max a bool indicating if min-max normalization will be used to scale the testing output, RMSE, predictive mean and std from the
#'     emulator. Defaults to `TRUE`.
#' @param color a character string indicating the color map to use when `style = 2`:
#' * `'magma'` (or `'A'`)
#' * `'inferno'` (or `'B'`)
#' * `'plasma'` (or '`C`')
#' * `'viridis'` (or `'D'`)
#' * `'cividis'` (or `'E'`)
#' * `'rocket'` (or `'F'`)
#' * `'mako'` (or `'G'`)
#' * `'turbo'` (or `'H'`)
#'
#' Defaults to `'turbo'` (or `'H'`).
#' @param type either `'line'` or `'points`, indicating whether to draw testing data in the OOS validation plot as a line or
#'     individual points when the input of the emulator is one-dimensional and `style = 1`. Defaults to `'points'`
#' @param verb a bool indicating if the trace information on plotting will be printed during the function execution.
#'     Defaults to `TRUE`.
#' @param force same as that of [validate()].
#' @param cores same as that of [validate()].
#' @param threading same as that of [validate()].
#' @param ... N/A.
#'
#' @return A `patchwork` object.
#'
#' @note
#' * [plot()] calls [validate()] internally to obtain validation results for plotting. However, [plot()] will not export the
#'   emulator object with validation results. Instead, it only returns the plotting object. For small-scale validations (i.e., small
#'   training or testing data points), direct execution of [plot()] is fine. However, for moderate- to large-scale validations,
#'   it is recommended to first run [validate()] to obtain and store validation results in the emulator object, and then supply the
#'   object to [plot()]. This is because if an emulator object has the validation results stored, each time when [plot()]
#'   is invoked, unnecessary evaluations of repetitive LOO or OOS validation will not be implemented.
#' * [plot()] uses information provided in `x_test` and `y_test` to produce the OOS validation plots. Therefore, if validation results
#'   are already stored in `x`, unless `x_test` and `y_test` are identical to those used by [validate()], [plot()] will re-evaluate OOS
#'   validations before plotting.
#' * Any R vector detected in `x_test` and `y_test` will be treated as a column vector and automatically converted into a single-column
#'   R matrix. Thus, if `x_test` or `y_test` is a single testing data point with multiple dimensions, it must be given as a matrix.
#' * The returned `patchwork` object contains the `ggplot2` objects. One can modify the included individual ggplots
#'   by accessing them with double-bracket indexing. See <https://patchwork.data-imaginist.com/> for further information.
#' @details See further examples and tutorials at <https://mingdeyu.github.io/dgpsi-R/>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @name plot
NULL


#' @rdname plot
#' @method plot dgp
#' @export
plot.dgp <- function(x, x_test = NULL, y_test = NULL, dim = NULL, method = 'mean_var', style = 1, min_max = TRUE, color = 'turbo', type = 'points', verb = TRUE, force = FALSE, cores = 1, threading = FALSE, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( style!=1&style!=2 ) stop("'style' must be either 1 or 2.", call. = FALSE)
  if ( type!='points'&type!='line' ) stop("'type' must be either 'points' or 'line'.", call. = FALSE)
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  if ( isTRUE(verb) ) message("Initializing ...", appendLF = FALSE)
  results <- validate(object = x, x_test = x_test, y_test = y_test, method = method, verb = FALSE, force = force, cores = cores, threading = threading)
  if ( isTRUE(verb) ) message(" done")

  # For LOO
  if ( is.null(x_test) & is.null(y_test) ){
    if ( isTRUE(verb) ) message("Post-processing LOO results ...", appendLF = FALSE)
    loo_res <- results$loo
    p_list <- list()
    if ( style==1 ){
      # create indices
      if ( ncol(loo_res$x_train)==1 ){
        idx <- loo_res$x_train[,1]
        isdup <- TRUE
      } else {
        if ( is.null(dim) ) {
          rep1 <- pkg.env$np$unique(loo_res$x_train, return_index=TRUE, axis=0L)
          rep2 <- pkg.env$np$unique(loo_res$x_train, return_inverse=TRUE, axis=0L)
          idx <- seq(1,length(rep1[[2]]))[pkg.env$np$argsort(pkg.env$np$argsort(rep1[[2]]))+1][rep2[[2]]+1]
          isdup <- TRUE
        } else {
          idx <- loo_res$x_train[,dim]
          isdup <- FALSE
        }
      }
      for (l in 1:ncol(loo_res$y_train) ) {
        dat <- list()
        dat[["idx"]] <- idx
        # extract mean or median
        if ( "mean" %in% names(loo_res) ){
          dat[["mean"]] <- loo_res$mean[,l]
          method <- "mean_var"
        } else if ( "median" %in% names(loo_res) ){
          dat[["median"]] <- loo_res$median[,l]
          method <- "sampling"
        }
        # extract other attributes
        dat[["lower"]] <- loo_res$lower[,l]
        dat[["upper"]] <- loo_res$upper[,l]
        dat[["y_validate"]] <- loo_res$y_train[,l]
        dat[["coverage"]] <- (dat[["y_validate"]]<=dat[["upper"]]) & (dat[["y_validate"]]>=dat[["lower"]])
        if ( min_max ){
          p_list[[l]] <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
            ggplot2::ggtitle(sprintf("O%i: NRMSE = %.2f%%", l, loo_res$nrmse[l]*100)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        } else {
          p_list[[l]] <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
            ggplot2::ggtitle(sprintf("O%i: RMSE = %.6f", l, loo_res$rmse[l])) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        }
      }
    } else if ( style==2 ) {
      for (l in 1:ncol(loo_res$y_train) ) {
        dat <- list()
        # extract mean or median
        if ( "mean" %in% names(loo_res) ){
          dat[["mean"]] <- loo_res$mean[,l]
          method <- "mean_var"
        } else if ( "median" %in% names(loo_res) ){
          dat[["median"]] <- loo_res$median[,l]
          method <- "sampling"
        }
        dat[["y_validate"]] <- loo_res$y_train[,l]
        dat[["std"]] <- loo_res$std[,l]
        if ( min_max ){
          p_list[[l]] <- plot_style_2(as.data.frame(dat), method, min_max, color) +
            ggplot2::ggtitle(sprintf("O%i: NRMSE = %.2f%%", l, loo_res$nrmse[l]*100)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        } else {
          p_list[[l]] <- plot_style_2(as.data.frame(dat), method, min_max, color) +
            ggplot2::ggtitle(sprintf("O%i: RMSE = %.6f", l, loo_res$rmse[l])) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        }
      }
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Plotting ...", appendLF = FALSE)
    if ( method == "mean_var" ){
      p_patch <- patchwork::wrap_plots(p_list) +
        patchwork::plot_annotation(
          title = 'Leave-One-Out (LOO) Cross Validation',
          caption = if ( style== 1) {
            if ( min_max ) {
              'Oi = Output i of the DGP emulator
               NRMSE = Normalized Root Mean Squared Error
               CI = Credible Interval'
            } else {
              'Oi = Output i of the DGP emulator
               RMSE = Root Mean Squared Error
               CI = Credible Interval'
            }
          } else {
            if ( min_max ) {
              'Oi = Output i of the DGP emulator
               NRMSE = Normalized Root Mean Squared Error'
            } else {
              'Oi = Output i of the DGP emulator
               RMSE = Root Mean Squared Error'
            }
            }
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    } else if ( method == "sampling" ){
      p_patch <- patchwork::wrap_plots(p_list) +
        patchwork::plot_annotation(
          title = 'Leave-One-Out (LOO) Cross Validation',
          caption = if ( style== 1) {
            if ( min_max ) {
              'Oi = Output i of the DGP emulator
               NRMSE = Normalized Root Median Squared Error
               CI = Credible Interval'
            } else {
              'Oi = Output i of the DGP emulator
               RMSE = Root Median Squared Error
               CI = Credible Interval'
            }
          } else {
            if ( min_max ) {
              'Oi = Output i of the DGP emulator
               NRMSE = Normalized Root Median Squared Error'
            } else {
              'Oi = Output i of the DGP emulator
               RMSE = Root Median Squared Error'
            }
          }
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")
    p_patch
    # For OOS
  } else if (!is.null(x_test) & !is.null(y_test)) {
    if ( isTRUE(verb) ) message("Post-processing OOS results ...", appendLF = FALSE)
    oos_res <- results$oos
    p_list <- list()
    if ( style==1 ){
      # If input is 1d
      if ( ncol(oos_res$x_test)==1 ){
        if ( "mean" %in% names(oos_res) ){
          method <- "mean_var"
        } else if ( "median" %in% names(oos_res) ){
          method <- "sampling"
        }
        # Extract training data points
        x_train <- results$constructor_obj$X[,1]
        y_train <- results$constructor_obj$Y
        rep <- results$constructor_obj$indices
        if ( !is.null(rep) ){
          rep <- rep + 1
          x_train <- x_train[rep]
        }
        # create a range for predictions
        x_min <- min(min(oos_res$x_test[,1]),min(x_train))
        x_max <- max(max(oos_res$x_test[,1]),max(x_train))
        x_range <- as.matrix(seq(x_min, x_max, length=500))
        if ( identical(cores,as.integer(1)) ){
          res <- results$emulator_obj$predict(x_range, method = method)
        } else {
          res <- results$emulator_obj$ppredict(x_range, method = method, core_num = cores)
        }
        if ( method=='sampling' ) {
          res_np <- pkg.env$np$array(res)
          quant <- pkg.env$np$transpose(pkg.env$np$quantile(res_np, c(0.025, 0.5, 0.975), axis=2L),c(0L,2L,1L))
        }

        for (l in 1:ncol(oos_res$y_test) ) {
          dat <- list()
          dat[["x_test"]] <- oos_res$x_test[,1]
          dat[["y_test"]] <- oos_res$y_test[,l]
          dat_train <- list('x_train' = x_train,'y_train' = y_train[,l])
          dat_range <- list()
          dat_range[["range"]] <- x_range[,1]
          if ( method == 'mean_var' ){
            dat_range[["mean"]] <- res[[1]][,l]
            std <- sqrt(res[[2]][,l])
            dat_range[["lower"]] <- dat_range$mean-2*std
            dat_range[["upper"]] <- dat_range$mean+2*std
          } else if ( method == 'sampling' ){
            dat_range[["median"]] <- quant[2,,l]
            dat_range[["lower"]] <- quant[1,,l]
            dat_range[["upper"]] <- quant[3,,l]
          }
          if ( min_max ) {
            p_list[[l]] <- plot_style_1_1d(as.data.frame(dat), as.data.frame(dat_range), as.data.frame(dat_train), method, type) +
              ggplot2::ggtitle(sprintf("O%i: NRMSE = %.2f%%", l, oos_res$nrmse[l]*100)) +
              ggplot2::theme(plot.title = ggplot2::element_text(size=10))
          } else {
            p_list[[l]] <- plot_style_1_1d(as.data.frame(dat), as.data.frame(dat_range), as.data.frame(dat_train), method, type) +
              ggplot2::ggtitle(sprintf("O%i: RMSE = %.6f", l, oos_res$rmse[l])) +
              ggplot2::theme(plot.title = ggplot2::element_text(size=10))
          }
        }
        # If input is at least 2d
      } else {
        if ( is.null(dim) ){
          rep1 <- pkg.env$np$unique(oos_res$x_test, return_index=TRUE, axis=0L)
          rep2 <- pkg.env$np$unique(oos_res$x_test, return_inverse=TRUE, axis=0L)
          idx <- seq(1,length(rep1[[2]]))[pkg.env$np$argsort(pkg.env$np$argsort(rep1[[2]]))+1][rep2[[2]]+1]
          isdup <- TRUE
        } else {
          idx <- oos_res$x_test[,dim]
          isdup <- FALSE
        }
        for (l in 1:ncol(oos_res$y_test) ) {
          dat <- list()
          dat[["idx"]] <- idx
          # extract mean or median
          if ( "mean" %in% names(oos_res) ){
            dat[["mean"]] <- oos_res$mean[,l]
            method <- "mean_var"
          } else if ( "median" %in% names(oos_res) ){
            dat[["median"]] <- oos_res$median[,l]
            method <- "sampling"
          }
          # extract other attributes
          dat[["lower"]] <- oos_res$lower[,l]
          dat[["upper"]] <- oos_res$upper[,l]
          dat[["y_validate"]] <- oos_res$y_test[,l]
          dat[["coverage"]] <- (dat[["y_validate"]]<=dat[["upper"]]) & (dat[["y_validate"]]>=dat[["lower"]])
          if ( min_max ) {
            p_list[[l]] <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
              ggplot2::ggtitle(sprintf("O%i: NRMSE = %.2f%%", l, oos_res$nrmse[l]*100)) +
              ggplot2::theme(plot.title = ggplot2::element_text(size=10))
          } else {
            p_list[[l]] <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
              ggplot2::ggtitle(sprintf("O%i: RMSE = %.6f", l, oos_res$rmse[l])) +
              ggplot2::theme(plot.title = ggplot2::element_text(size=10))
          }
        }
      }
    } else if ( style==2 ) {
      for (l in 1:ncol(oos_res$y_test) ) {
        dat <- list()
        # extract mean or median
        if ( "mean" %in% names(oos_res) ){
          dat[["mean"]] <- oos_res$mean[,l]
          method <- "mean_var"
        } else if ( "median" %in% names(oos_res) ){
          dat[["median"]] <- oos_res$median[,l]
          method <- "sampling"
        }
        dat[["y_validate"]] <- oos_res$y_test[,l]
        dat[["std"]] <- oos_res$std[,l]
        if ( min_max ) {
          p_list[[l]] <- plot_style_2(as.data.frame(dat), method, min_max, color) +
            ggplot2::ggtitle(sprintf("O%i: NRMSE = %.2f%%", l, oos_res$nrmse[l]*100)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        } else {
          p_list[[l]] <- plot_style_2(as.data.frame(dat), method, min_max, color) +
            ggplot2::ggtitle(sprintf("O%i: RMSE = %.6f", l, oos_res$rmse[l])) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        }
      }
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Plotting ...", appendLF = FALSE)
    if ( method == "mean_var" ){
      p_patch <- patchwork::wrap_plots(p_list) +
        patchwork::plot_annotation(
          title = 'Out-Of-Sample (OOS) Validation',
          caption = if ( style==1 ){
            if ( min_max ) {
              'Oi = Output i of the DGP emulator
               NRMSE = Normalized Root Mean Squared Error
               CI = Credible Interval'
            } else {
              'Oi = Output i of the DGP emulator
               RMSE = Root Mean Squared Error
               CI = Credible Interval'
            }
          } else {
            if ( min_max ) {
              'Oi = Output i of the DGP emulator
               NRMSE = Normalized Root Mean Squared Error'
            } else {
              'Oi = Output i of the DGP emulator
               RMSE = Root Mean Squared Error'
            }
          }
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    } else if ( method == "sampling" ){
      p_patch <- patchwork::wrap_plots(p_list) +
        patchwork::plot_annotation(
          title = 'Out-Of-Sample (OOS) Validation',
          caption = if ( style==1 ){
            if ( min_max ) {
              'Oi = Output i of the DGP emulator
               NRMSE = Normalized Root Median Squared Error
               CI = Credible Interval'
            } else {
              'Oi = Output i of the DGP emulator
               RMSE = Root Median Squared Error
               CI = Credible Interval'
            }
          } else {
            if ( min_max ) {
              'Oi = Output i of the DGP emulator
               NRMSE = Normalized Root Median Squared Error'
            } else {
              'Oi = Output i of the DGP emulator
               RMSE = Root Median Squared Error'
            }
          }
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")
    p_patch
  }
}

#' @rdname plot
#' @method plot lgp
#' @export
plot.lgp <- function(x, x_test = NULL, y_test = NULL, dim = NULL, method = 'mean_var', style = 1, min_max = TRUE, color = 'turbo', type = 'points', verb = TRUE, force = FALSE, cores = 1, threading = FALSE, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( style!=1&style!=2 ) stop("'style' must be either 1 or 2.", call. = FALSE)
  if ( type!='points'&type!='line' ) stop("'type' must be either 'points' or 'line'.", call. = FALSE)
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }
  if ( isTRUE(verb) ) message("Initializing ...", appendLF = FALSE)
  results <- validate(object = x, x_test = x_test, y_test = y_test, method = method, verb = FALSE, force = force, cores = cores, threading = threading)
  if ( isTRUE(verb) ) message(" done")

  if ( isTRUE(verb) ) message("Post-processing OOS results ...", appendLF = FALSE)
  oos_res <- results$oos
  if ( !is.list(oos_res$y_test) ) {
    y_test_list <- list(oos_res$y_test)
  } else {
    y_test_list <- oos_res$y_test
  }
  p_list <- list()
  if ( style==1 ){
    # check if the test input is 1d
    if ( !is.list(oos_res$x_test) ) {
      if ( ncol(oos_res$x_test)==1 ) {
        single_dim <- TRUE
      } else {
        single_dim <- FALSE
      }
    } else {
      for ( l in 1:length(oos_res$x_test) ){
        if ( l==1 ){
          if ( ncol(oos_res$x_test[[l]])==1 ) {
            single_dim <- TRUE
          } else {
            single_dim <- FALSE
          }
        } else {
          for ( k in 1:length(oos_res$x_test[[l]]) ){
            if ( is.null(oos_res$x_test[[l]][[k]]) ){
              single_dim <- single_dim * TRUE
            } else {
              single_dim <- single_dim * FALSE
            }
          }
        }
      }
    }
    # If input is 1d
    if ( isTRUE(single_dim) ){
      if ( "mean" %in% names(oos_res) ){
        method <- "mean_var"
      } else if ( "median" %in% names(oos_res) ){
        method <- "sampling"
      }

      # extract training input for range construction
      if ( results$emulator_obj$all_layer[[1]][[1]]$type == 'gp' ){
        x_train <- results$emulator_obj$all_layer[[1]][[1]]$structure$input[,1]
      } else if ( results$emulator_obj$all_layer[[1]][[1]]$type == 'dgp' ) {
        x_train <- results$emulator_obj$all_layer[[1]][[1]]$structure[[1]][[1]]$input[,1]
      }

      if ( is.list(oos_res$x_test) ) {
        x_min <- min(min(oos_res$x_test[[1]][,1]), min(x_train))
        x_max <- max(max(oos_res$x_test[[1]][,1]), max(x_train))
      } else {
        x_min <- min(min(oos_res$x_test[,1]), min(x_train))
        x_max <- max(max(oos_res$x_test[,1]), max(x_train))
      }

      x_range <- as.matrix(seq(x_min, x_max, length=500))
      if ( identical(cores,as.integer(1)) ){
        res <- results$emulator_obj$predict(x_range, method = method)
      } else {
        res <- results$emulator_obj$ppredict(x_range, method = method, core_num = cores)
      }

      counter <- 1
      for ( k in 1:length(oos_res$nrmse) ) {
        if ( method=='sampling' ) quant <- pkg.env$np$transpose(pkg.env$np$quantile(res[[k]], c(0.025, 0.5, 0.975), axis=2L),c(0L,2L,1L))
        for ( l in 1:length(oos_res$nrmse[[k]]) ) {
          dat <- list()
          if ( is.list(oos_res$x_test) ) {
            dat[["x_test"]] <- oos_res$x_test[[1]][,1]
          } else {
            dat[["x_test"]] <- oos_res$x_test[,1]
          }

          dat[["y_test"]] <- y_test_list[[k]][,l]
          dat_range <- list()
          dat_range[["range"]] <- x_range[,1]
          if ( method == 'mean_var' ){
            dat_range[["mean"]] <- res[[1]][[k]][,l]
            std <- sqrt(res[[2]][[k]][,l])
            dat_range[["lower"]] <- dat_range$mean-2*std
            dat_range[["upper"]] <- dat_range$mean+2*std
          } else if ( method == 'sampling' ){
            dat_range[["median"]] <- quant[2,,l]
            dat_range[["lower"]] <- quant[1,,l]
            dat_range[["upper"]] <- quant[3,,l]
          }
          if ( min_max ) {
            p_list[[counter]] <- plot_style_1_1d(as.data.frame(dat), as.data.frame(dat_range), NULL, method, type) +
              ggplot2::ggtitle(sprintf("E%iO%i: NRMSE = %.2f%%", k, l, oos_res$nrmse[[k]][l]*100)) +
              ggplot2::theme(plot.title = ggplot2::element_text(size=10))
          } else {
            p_list[[counter]] <- plot_style_1_1d(as.data.frame(dat), as.data.frame(dat_range), NULL, method, type) +
              ggplot2::ggtitle(sprintf("E%iO%i: RMSE = %.6f", k, l, oos_res$rmse[[k]][l])) +
              ggplot2::theme(plot.title = ggplot2::element_text(size=10))
          }
          counter <- counter + 1
        }
      }
      # If input is at least 2d
    } else {
      if ( is.null(dim) ) {
        idx <- seq(1,nrow(y_test_list[[1]]))
        isdup <- FALSE
      } else {
        if ( is.list(oos_res$x_test) ){
          if ( length(dim)==1 ){
            idx <- oos_res$x_test[[1]][,dim]
            dim <- c(1,dim)
          } else if ( length(dim)==3 ){
            if ( dim[1]<2 ) stop("The first element of 'dim' must be equal or greater than 2. See documentation for details.", call. = FALSE)
            target_emulator_input <- oos_res$x_test[[dim[1]]][[dim[2]]]
            if ( is.null(target_emulator_input) ) stop("The emulator specified by the first two elements of 'dim' has no global input.", call. = FALSE)
            idx <- target_emulator_input[,dim[3]]
          } else {
            stop("'dim' must be a vector of length 1 or 3. See documentation for details.", call. = FALSE)
          }
        } else {
          if ( length(dim)==1 ){
            idx <- oos_res$x_test[,dim]
            dim <- c(1,dim)
          } else {
            stop("'dim' must be a vector of length 1. See documentation for details.", call. = FALSE)
          }
        }
        isdup <- FALSE
      }

      counter <- 1
      for ( k in 1:length(oos_res$nrmse) ) {
        for ( l in 1:length(oos_res$nrmse[[k]]) ) {
          dat <- list()
          dat[["idx"]] <- idx
          # extract mean or median
          if ( "mean" %in% names(oos_res) ){
            dat[["mean"]] <- oos_res$mean[[k]][,l]
            method <- "mean_var"
          } else if ( "median" %in% names(oos_res) ){
            dat[["median"]] <- oos_res$median[[k]][,l]
            method <- "sampling"
          }
          # extract other attributes
          dat[["lower"]] <- oos_res$lower[[k]][,l]
          dat[["upper"]] <- oos_res$upper[[k]][,l]
          dat[["y_validate"]] <- y_test_list[[k]][,l]
          dat[["coverage"]] <- (dat[["y_validate"]]<=dat[["upper"]]) & (dat[["y_validate"]]>=dat[["lower"]])
          if ( min_max ) {
            p_list[[counter]] <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
              ggplot2::ggtitle(sprintf("E%iO%i: NRMSE = %.2f%%", k, l, oos_res$nrmse[[k]][l]*100)) +
              ggplot2::theme(plot.title = ggplot2::element_text(size=10))
          } else {
            p_list[[counter]] <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
              ggplot2::ggtitle(sprintf("E%iO%i: RMSE = %.6f", k, l, oos_res$rmse[[k]][l])) +
              ggplot2::theme(plot.title = ggplot2::element_text(size=10))
          }
          counter <- counter + 1
        }
      }
    }
  } else if ( style==2 ) {
    counter <- 1
    for ( k in 1:length(oos_res$nrmse) ) {
      for (l in 1:length(oos_res$nrmse[[k]]) ) {
        dat <- list()
        # extract mean or median
        if ( "mean" %in% names(oos_res) ){
          dat[["mean"]] <- oos_res$mean[[k]][,l]
          method <- "mean_var"
        } else if ( "median" %in% names(oos_res) ){
          dat[["median"]] <- oos_res$median[[k]][,l]
          method <- "sampling"
        }
        dat[["y_validate"]] <- y_test_list[[k]][,l]
        dat[["std"]] <- oos_res$std[[k]][,l]
        if ( min_max ) {
          p_list[[counter]] <- plot_style_2(as.data.frame(dat), method, min_max, color) +
            ggplot2::ggtitle(sprintf("E%iO%i: NRMSE = %.2f%%", k, l, oos_res$nrmse[[k]][l]*100)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        } else {
          p_list[[counter]] <- plot_style_2(as.data.frame(dat), method, min_max, color) +
            ggplot2::ggtitle(sprintf("E%iO%i: RMSE = %.6f", k, l, oos_res$rmse[[k]][l])) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        }
        counter <- counter + 1
      }
    }
  }
  if ( isTRUE(verb) ) Sys.sleep(0.5)
  if ( isTRUE(verb) ) message(" done")

  if ( isTRUE(verb) ) message("Plotting ...", appendLF = FALSE)
  if ( method == "mean_var" ){
    p_patch <- patchwork::wrap_plots(p_list) +
      patchwork::plot_annotation(
        title = 'Out-Of-Sample (OOS) Validation',
        caption = if ( style==1 ){
          if ( min_max ) {
            'EiOj = Output j of Emulator i in the final layer of the linked emulator
             NRMSE = Normalized Root Mean Squared Error
             CI = Credible Interval'
          } else {
            'EiOj = Output j of Emulator i in the final layer of the linked emulator
             RMSE = Root Mean Squared Error
             CI = Credible Interval'
          }
        } else {
          if ( min_max ) {
            'EiOj = Output j of Emulator i in the final layer of the linked emulator
             NRMSE = Normalized Root Mean Squared Error'
          } else {
            'EiOj = Output j of Emulator i in the final layer of the linked emulator
             RMSE = Root Mean Squared Error'
          }
        }
      ) +
      patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
  } else if ( method == "sampling" ){
    p_patch <- patchwork::wrap_plots(p_list) +
      patchwork::plot_annotation(
        title = 'Out-Of-Sample (OOS) Validation',
        caption = if ( style==1 ){
          if ( min_max ) {
            'EiOj = Output j of Emulator i in the final layer of the linked emulator
             NRMSE = Normalized Root Median Squared Error
             CI = Credible Interval'
          } else {
            'EiOj = Output j of Emulator i in the final layer of the linked emulator
             RMSE = Root Median Squared Error
             CI = Credible Interval'
          }
        } else {
          if ( min_max ) {
            'EiOj = Output j of Emulator i in the final layer of the linked emulator
             NRMSE = Normalized Root Median Squared Error'
          } else {
            'EiOj = Output j of Emulator i in the final layer of the linked emulator
             RMSE = Root Median Squared Error'
          }
        }
      ) +
      patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
  }
  if ( isTRUE(verb) ) Sys.sleep(0.5)
  if ( isTRUE(verb) ) message(" done")
  p_patch
}

#' @rdname plot
#' @method plot gp
#' @export
plot.gp <- function(x, x_test = NULL, y_test = NULL, dim = NULL, method = 'mean_var', style = 1, min_max = TRUE, color = 'turbo', type = 'points', verb = TRUE, force = FALSE, cores = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( style!=1&style!=2 ) stop("'style' must be either 1 or 2.", call. = FALSE)
  if ( type!='points'&type!='line' ) stop("'type' must be either 'points' or 'line'.", call. = FALSE)
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("The core number must be >= 1.", call. = FALSE)
  }

  if ( isTRUE(verb) ) message("Initializing ...", appendLF = FALSE)
  results <- validate(object = x, x_test = x_test, y_test = y_test, method = method, verb = FALSE, force = force, cores = cores)
  if ( isTRUE(verb) ) message(" done")

  dat <- list()
  # For LOO
  if ( is.null(x_test) & is.null(y_test) ){
    if ( isTRUE(verb) ) message("Post-processing LOO results ...", appendLF = FALSE)
    loo_res <- results$loo
    if ( style==1 ){
      # create indices
      if ( ncol(loo_res$x_train)==1 ){
        dat[["idx"]] <- loo_res$x_train[,1]
        isdup <- TRUE
      } else {
        if ( is.null(dim) ){
          rep1 <- pkg.env$np$unique(loo_res$x_train, return_index=TRUE, axis=0L)
          rep2 <- pkg.env$np$unique(loo_res$x_train, return_inverse=TRUE, axis=0L)
          idx <- seq(1,length(rep1[[2]]))[pkg.env$np$argsort(pkg.env$np$argsort(rep1[[2]]))+1][rep2[[2]]+1]
          dat[["idx"]] <- idx
          isdup <- TRUE
        } else {
          dat[["idx"]] <- loo_res$x_train[,dim]
          isdup <- FALSE
        }
      }
      # extract mean or median
      if ( "mean" %in% names(loo_res) ){
        dat[["mean"]] <- loo_res$mean[,1]
        method <- "mean_var"
      } else if ( "median" %in% names(loo_res) ){
        dat[["median"]] <- loo_res$median[,1]
        method <- "sampling"
      }
      # extract other attributes
      dat[["lower"]] <- loo_res$lower[,1]
      dat[["upper"]] <- loo_res$upper[,1]
      dat[["y_validate"]] <- loo_res$y_train[,1]
      dat[["coverage"]] <- (dat[["y_validate"]]<=dat[["upper"]]) & (dat[["y_validate"]]>=dat[["lower"]])
      if ( min_max ) {
        p <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
          ggplot2::ggtitle(sprintf('NRMSE = %.2f%%', loo_res$nrmse*100)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      } else {
        p <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
          ggplot2::ggtitle(sprintf('RMSE = %.6f', loo_res$rmse)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      }
    } else if ( style==2 ) {
      # extract mean or median
      if ( "mean" %in% names(loo_res) ){
        dat[["mean"]] <- loo_res$mean[,1]
        method <- "mean_var"
      } else if ( "median" %in% names(loo_res) ){
        dat[["median"]] <- loo_res$median[,1]
        method <- "sampling"
      }
      dat[["y_validate"]] <- loo_res$y_train[,1]
      dat[["std"]] <- loo_res$std[,1]
      if ( min_max ) {
        p <- plot_style_2(as.data.frame(dat), method, min_max, color) +
          ggplot2::ggtitle(sprintf('NRMSE = %.2f%%', loo_res$nrmse*100)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      } else {
        p <- plot_style_2(as.data.frame(dat), method, min_max, color) +
          ggplot2::ggtitle(sprintf('RMSE = %.6f', loo_res$rmse)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      }
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Plotting ...", appendLF = FALSE)
    if ( method == "mean_var" ){
      p_patch <- patchwork::wrap_plots(p) +
        patchwork::plot_annotation(
          title = 'Leave-One-Out (LOO) Cross Validation',
          caption = if ( style==1 ){
            if ( min_max ) {
              'NRMSE = Normalized Root Mean Squared Error
               CI = Credible Interval'
            } else {
              'RMSE = Root Mean Squared Error
               CI = Credible Interval'
            }
          } else {
            if ( min_max ) {
              'NRMSE = Normalized Root Mean Squared Error'
            } else {
              'RMSE = Root Mean Squared Error'
            }
          }
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    } else if ( method == "sampling" ){
      p_patch <- patchwork::wrap_plots(p) +
        patchwork::plot_annotation(
          title = 'Leave-One-Out (LOO) Cross Validation',
          caption = if ( style==1 ){
            if ( min_max ) {
              'NRMSE = Normalized Root Median Squared Error
               CI = Credible Interval'
            } else {
              'RMSE = Root Median Squared Error
               CI = Credible Interval'
            }
          } else {
            if ( min_max ) {
              'NRMSE = Normalized Root Median Squared Error'
            } else {
              'RMSE = Root Median Squared Error'
            }
          }
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")
    p_patch
    # For OOS
  } else if (!is.null(x_test) & !is.null(y_test)) {
    if ( isTRUE(verb) ) message("Post-processing OOS results ...", appendLF = FALSE)
    oos_res <- results$oos
    if ( style==1 ){
      # If input is 1d
      if ( ncol(oos_res$x_test)==1 ){
        dat[["x_test"]] <- oos_res$x_test[,1]
        dat[["y_test"]] <- oos_res$y_test[,1]
        # extract training data points
        dat_train <- list('x_train' = results$constructor_obj$X[,1], 'y_train' = results$constructor_obj$Y[,1])
        # construct a range for predictions
        dat_range <- list()
        if ( "mean" %in% names(oos_res) ){
          method <- "mean_var"
        } else if ( "median" %in% names(oos_res) ){
          method <- "sampling"
        }
        x_min <- min(min(oos_res$x_test[,1]), min(dat_train$x_train))
        x_max <- max(max(oos_res$x_test[,1]), max(dat_train$x_train))
        x_range <- as.matrix(seq(x_min, x_max, length=500))
        if ( identical(cores,as.integer(1)) ){
          res <- results$emulator_obj$predict(x_range, method = method, sample_size=500L)
        } else {
          res <- results$emulator_obj$ppredict(x_range, method = method, sample_size=500L, core_num = cores)
        }

        dat_range[["range"]] <- x_range[,1]
        if ( method == 'mean_var' ){
          dat_range[["mean"]] <- res[[1]][,1]
          std <- sqrt(res[[2]][,1])
          dat_range[["lower"]] <- dat_range$mean-2*std
          dat_range[["upper"]] <- dat_range$mean+2*std
        } else if ( method == 'sampling' ){
          quant <- t(pkg.env$np$quantile(res, c(0.025, 0.5, 0.975), axis=1L))
          dat_range[["median"]] <- quant[,2]
          dat_range[["lower"]] <- quant[,1]
          dat_range[["upper"]] <- quant[,3]
        }
        if ( min_max ) {
          p <- plot_style_1_1d(as.data.frame(dat), as.data.frame(dat_range), as.data.frame(dat_train), method, type) +
            ggplot2::ggtitle(sprintf('NRMSE = %.2f%%', oos_res$nrmse*100)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        } else {
          p <- plot_style_1_1d(as.data.frame(dat), as.data.frame(dat_range), as.data.frame(dat_train), method, type) +
            ggplot2::ggtitle(sprintf('RMSE = %.6f', oos_res$rmse)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        }
        # If input is at least 2d
      } else {
        if ( is.null(dim) ){
          rep1 <- pkg.env$np$unique(oos_res$x_test, return_index=TRUE, axis=0L)
          rep2 <- pkg.env$np$unique(oos_res$x_test, return_inverse=TRUE, axis=0L)
          idx <- seq(1,length(rep1[[2]]))[pkg.env$np$argsort(pkg.env$np$argsort(rep1[[2]]))+1][rep2[[2]]+1]
          dat[["idx"]] <- idx
          isdup <- TRUE
        } else {
          dat[["idx"]] <- oos_res$x_test[,dim]
          isdup <- FALSE
        }
        # extract mean or median
        if ( "mean" %in% names(oos_res) ){
          dat[["mean"]] <- oos_res$mean[,1]
          method <- "mean_var"
        } else if ( "median" %in% names(oos_res) ){
          dat[["median"]] <- oos_res$median[,1]
          method <- "sampling"
        }
        # extract other attributes
        dat[["lower"]] <- oos_res$lower[,1]
        dat[["upper"]] <- oos_res$upper[,1]
        dat[["y_validate"]] <- oos_res$y_test[,1]
        dat[["coverage"]] <- (dat[["y_validate"]]<=dat[["upper"]]) & (dat[["y_validate"]]>=dat[["lower"]])
        if ( min_max ) {
          p <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
            ggplot2::ggtitle(sprintf('NRMSE = %.2f%%', oos_res$nrmse*100)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        } else {
          p <- plot_style_1(as.data.frame(dat), method, dim, isdup) +
            ggplot2::ggtitle(sprintf('RMSE = %.6f', oos_res$rmse)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        }
      }
    } else if ( style==2 ) {
      # extract mean or median
      if ( "mean" %in% names(oos_res) ){
        dat[["mean"]] <- oos_res$mean[,1]
        method <- "mean_var"
      } else if ( "median" %in% names(oos_res) ){
        dat[["median"]] <- oos_res$median[,1]
        method <- "sampling"
      }
      dat[["y_validate"]] <- oos_res$y_test[,1]
      dat[["std"]] <- oos_res$std[,1]
      if ( min_max ) {
        p <- plot_style_2(as.data.frame(dat), method, min_max, color) +
          ggplot2::ggtitle(sprintf('NRMSE = %.2f%%', oos_res$nrmse*100)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      } else {
        p <- plot_style_2(as.data.frame(dat), method, min_max, color) +
          ggplot2::ggtitle(sprintf('RMSE = %.6f', oos_res$rmse)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      }
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Plotting ...", appendLF = FALSE)
    if ( method == "mean_var" ){
      p_patch <- patchwork::wrap_plots(p) +
        patchwork::plot_annotation(
          title = 'Out-Of-Sample (OOS) Validation',
          caption = if ( style==1 ){
            if ( min_max ) {
              'NRMSE = Normalized Root Mean Squared Error
               CI = Credible Interval'
            } else {
              'RMSE = Root Mean Squared Error
               CI = Credible Interval'
            }
          } else {
            if ( min_max ) {
              'NRMSE = Normalized Root Mean Squared Error'
            } else {
              'RMSE = Root Mean Squared Error'
            }
          }
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    } else if ( method == "sampling" ){
      p_patch <- patchwork::wrap_plots(p) +
        patchwork::plot_annotation(
          title = 'Out-Of-Sample (OOS) Validation',
          caption = if ( style==1 ){
            if ( min_max ) {
              'NRMSE = Normalized Root Median Squared Error
               CI = Credible Interval'
            } else {
              'RMSE = Root Median Squared Error
               CI = Credible Interval'
            }
          } else {
            if ( min_max ) {
              'NRMSE = Normalized Root Median Squared Error'
            } else {
              'RMSE = Root Median Squared Error'
            }
          }
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")
    p_patch
  }
}


plot_style_1 <- function(dat, method, dim, isdup) {
  if ( isTRUE(isdup) ){
    dup <- duplicated(dat$idx)
  } else {
    dup <- as.logical(rep(0, length(dat$idx)))
  }

  if ( is.null(dim) ){
    x_lab <- "Input position"
  } else {
    if ( length(dim)==1 ){
      x_lab <- sprintf("Input dimension %i", dim)
    } else if ( length(dim)==2 ){
      x_lab <- sprintf("Input dimension %i of emulators in layer %i", dim[2], dim[1])
    } else if ( length(dim)==3 ){
      x_lab <- sprintf("Global input dimension %i of emulator %i in layer %i", dim[3], dim[2], dim[1])
    }
  }

  coverage <- dat$coverage

  if ( all(coverage) ){
    p <- ggplot2::ggplot(dat, ggplot2::aes_(x=~idx, y=~y_validate, color = "Validation point inside CI"))
    if ( method=="sampling" ){
      p <- p +
        ggplot2::geom_pointrange(data=dat[!dup,], ggplot2::aes_(x=~idx, y=~median, ymin=~lower, ymax=~upper, color = "Median and 95% CI"), fatten = 1.5, size = 0.3) +
        ggplot2::scale_color_manual(name = "", values = c("Median and 95% CI"="#52854C", "Validation point outside CI"="#D55E00", "Validation point inside CI"="#E69F00"),
                                    limits = c("Median and 95% CI", "Validation point outside CI", "Validation point inside CI"))
    } else if ( method=="mean_var" ) {
      p <- p +
        ggplot2::geom_pointrange(data=dat[!dup,], ggplot2::aes_(x=~idx, y=~mean, ymin=~lower, ymax=~upper, color = "Mean and CI (+/-2SD)"), fatten = 1.5, size = 0.3) +
        ggplot2::scale_color_manual(name = "", values = c("Mean and CI (+/-2SD)"="#52854C", "Validation point outside CI"="#D55E00", "Validation point inside CI"="#E69F00"),
                                    limits = c("Mean and CI (+/-2SD)", "Validation point outside CI", "Validation point inside CI"))
    }

    p <- p +
      ggplot2::geom_point(size=0.8) +
      ggplot2::labs(x =x_lab, y = "Model output") +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 7),
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "blank", "blank"))))
  } else if ( all(!coverage) ){
    p <- ggplot2::ggplot(dat, ggplot2::aes_(x=~idx, y=~y_validate, color = "Validation point outside CI"))
    if ( method=="sampling" ){
      p <- p +
        ggplot2::geom_pointrange(data=dat[!dup,], ggplot2::aes_(x=~idx, y=~median, ymin=~lower, ymax=~upper, color = "Median and 95% CI"), fatten = 1.5, size = 0.3) +
        ggplot2::scale_color_manual(name = "", values = c("Median and 95% CI"="#52854C", "Validation point outside CI"="#D55E00", "Validation point inside CI"="#E69F00"),
                                    limits = c("Median and 95% CI", "Validation point outside CI", "Validation point inside CI"))
    } else if ( method=="mean_var" ) {
      p <- p +
        ggplot2::geom_pointrange(data=dat[!dup,], ggplot2::aes_(x=~idx, y=~mean, ymin=~lower, ymax=~upper, color = "Mean and CI (+/-2SD)"), fatten = 1.5, size = 0.3) +
        ggplot2::scale_color_manual(name = "", values = c("Mean and CI (+/-2SD)"="#52854C", "Validation point outside CI"="#D55E00", "Validation point inside CI"="#E69F00"),
                                    limits = c("Mean and CI (+/-2SD)", "Validation point outside CI", "Validation point inside CI"))
    }

    p <- p +
      ggplot2::geom_point(size=0.8) +
      ggplot2::labs(x =x_lab, y = "Model output") +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 7),
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "blank", "blank"))))
  } else {
    p <- ggplot2::ggplot(dat, ggplot2::aes_(x=~idx, y=~y_validate, color = ~coverage))
    if ( method=="sampling" ){
      p <- p +
        ggplot2::geom_pointrange(data=dat[!dup,], ggplot2::aes_(x=~idx, y=~median, ymin=~lower, ymax=~upper, color = "#52854C"), fatten = 1.5, size = 0.3) +
        ggplot2::scale_color_manual(name = "", labels = c("Median and 95% CI","Validation point outside CI", "Validation point inside CI"), values = c("#52854C", "#D55E00", "#E69F00"))
    } else if ( method=="mean_var" ) {
      p <- p +
        ggplot2::geom_pointrange(data=dat[!dup,], ggplot2::aes_(x=~idx, y=~mean, ymin=~lower, ymax=~upper, color = "#52854C"), fatten = 1.5, size = 0.3) +
        ggplot2::scale_color_manual(name = "", labels=c( "Mean and CI (+/-2SD)", "Validation point outside CI", "Validation point inside CI"), values = c("#52854C","#D55E00", "#E69F00"))
    }

    p <- p +
      ggplot2::geom_point(size=0.8) +
      ggplot2::labs(x =x_lab, y = "Model output") +
      ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 7),
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c("solid", "blank", "blank"))))
  }
  return(p)
}

plot_style_2 <- function(dat, method, min_max, color) {
  y_min <- min(dat$y_validate)
  y_max <- max(dat$y_validate)
  std_min <- min(dat$std)
  std_max <- max(dat$std)
  if (isTRUE(min_max)){
    if ( method=="sampling" ){
      p <- ggplot2::ggplot(dat, ggplot2::aes_(x=~(median-y_min)/(y_max-y_min), y=~(y_validate-y_min)/(y_max-y_min), color= ~(std-std_min)/(std_max-std_min))) +
        ggplot2::labs(x ="Normalized predictive median", y = "Normalized model output")
    } else if ( method=="mean_var" ) {
      p <- ggplot2::ggplot(dat, ggplot2::aes_(x=~(mean-y_min)/(y_max-y_min), y=~(y_validate-y_min)/(y_max-y_min), color= ~(std-std_min)/(std_max-std_min))) +
        ggplot2::labs(x ="Normalized predictive mean", y = "Normalized model output")
    }
  } else {
    if ( method=="sampling" ){
      p <- ggplot2::ggplot(dat, ggplot2::aes_(x=~median, y=~y_validate, color=~std)) +
        ggplot2::labs(x ="Predictive median", y = "Model output")
    } else if ( method=="mean_var" ) {
      p <- ggplot2::ggplot(dat, ggplot2::aes_(x=~mean, y=~y_validate, color=~std)) +
        ggplot2::labs(x ="Predictive mean", y = "Model output")
    }
  }

  p <- p +
    ggplot2::geom_abline(alpha=0.6,intercept = 0, slope = 1) +
    ggplot2::geom_point(alpha=0.8, size=1.5)

  if (isTRUE(min_max)){
    p <- p + ggplot2::scale_colour_viridis_c("Normalized Predictive SD", option = color, breaks=seq(0,1,0.2), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'))
  } else {
    p <- p + ggplot2::scale_colour_viridis_c("Predictive SD", option = color)
  }

  p <- p +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1, "cm"),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 7),
      legend.title.align=0.5
    ) +
    ggplot2::guides(colour = ggplot2::guide_colourbar(title.position="top"))

  return(p)
}

plot_style_1_1d <- function(dat1, dat2, dat3, method, type) {
  p <- ggplot2::ggplot(data=dat2, ggplot2::aes_(x=~range))

  if ( method=="sampling" ){
    p <- p +
      ggplot2::geom_ribbon(data=dat2, alpha=0.5, mapping=ggplot2::aes_(ymin=~lower, ymax=~upper, fill="95% CI"))

    if (type == 'line') {
      p <- p +
        ggplot2::geom_line(data=dat1, ggplot2::aes_(x=~x_test, y=~y_test, color = "Testing Function"), linetype="solid", size=0.5, alpha=0.9)
    }

    p <- p + ggplot2::geom_line(data=dat2, ggplot2::aes_(y=~median, color="Pred. Median"), linetype="dashed", size=0.5, alpha=0.9)

    if ( !is.null(dat3) ) p <- p + ggplot2::geom_point(data=dat3, ggplot2::aes_(x=~x_train, y=~y_train, color = "Training Point"), size=2, alpha=0.9, shape=19)

    if ( type == 'points' ){
      p <- p +
        ggplot2::geom_point(data=dat1, ggplot2::aes_(x=~x_test, y=~y_test, color = "Testing Point"), size=1.75, alpha=0.9, shape=17)
    }
    if (type == 'points') {
      if ( is.null(dat3) ){
        p <- p +
          ggplot2::scale_color_manual(name = "xx", values = c("Pred. Median"="#000000", "Testing Point"="#D55E00")) +
          ggplot2::scale_fill_manual(name = "xx", values=c("95% CI"="#999999"))
      } else {
        p <- p +
          ggplot2::scale_color_manual(name = "xx", values = c("Pred. Median"="#000000", "Training Point"="#0072B2", "Testing Point"="#D55E00")) +
          ggplot2::scale_fill_manual(name = "xx", values=c("95% CI"="#999999"))
      }
    } else if (type == 'line') {
      if ( is.null(dat3) ){
        p <- p +
          ggplot2::scale_color_manual(name = "xx", values = c("Pred. Median"="#000000", "Testing Function"="#D55E00")) +
          ggplot2::scale_fill_manual(name = "xx", values=c("95% CI"="#999999"))
      } else {
        p <- p +
          ggplot2::scale_color_manual(name = "xx", values = c("Pred. Median"="#000000", "Training Point"="#0072B2", "Testing Function"="#D55E00")) +
          ggplot2::scale_fill_manual(name = "xx", values=c("95% CI"="#999999"))
      }
    }
  } else if ( method=="mean_var" ) {
    p <- p +
      ggplot2::geom_ribbon(data=dat2, alpha=0.5, mapping=ggplot2::aes_(ymin=~lower, ymax=~upper, fill="CI (+/-2SD)"))

    if (type == 'line') {
      p <- p +
        ggplot2::geom_line(data=dat1, ggplot2::aes_(x=~x_test, y=~y_test, color = "Testing Function"), linetype="solid", size=0.5, alpha=0.9)
    }

    p <- p + ggplot2::geom_line(data=dat2, ggplot2::aes_(y=~mean, color="Pred. Mean"),linetype="dashed", size=0.5, alpha=0.9)

    if ( !is.null(dat3) ) p <- p + ggplot2::geom_point(data=dat3, ggplot2::aes_(x=~x_train, y=~y_train, color = "Training Point"), size=2, alpha=0.9, shape=19)

    if (type == 'points') {
      p <- p +
        ggplot2::geom_point(data=dat1, ggplot2::aes_(x=~x_test, y=~y_test, color = "Testing Point"), size=1.75, alpha=0.9, shape=17)
    }

    if (type == 'points') {
      if ( is.null(dat3) ){
        p <- p +
          ggplot2::scale_color_manual(name = "xx", values = c("Pred. Mean"="#000000", "Testing Point"="#D55E00")) +
          ggplot2::scale_fill_manual(name = "xx", values=c("CI (+/-2SD)"="#999999"))
      } else {
        p <- p +
          ggplot2::scale_color_manual(name = "xx", values = c("Pred. Mean"="#000000", "Training Point"="#0072B2", "Testing Point"="#D55E00")) +
          ggplot2::scale_fill_manual(name = "xx", values=c("CI (+/-2SD)"="#999999"))
      }
    } else if (type == 'line') {
      if ( is.null(dat3) ){
        p <- p +
          ggplot2::scale_color_manual(name = "xx", values = c("Pred. Mean"="#000000", "Testing Function"="#D55E00")) +
          ggplot2::scale_fill_manual(name = "xx", values=c("CI (+/-2SD)"="#999999"))
      } else {
        p <- p +
          ggplot2::scale_color_manual(name = "xx", values = c("Pred. Mean"="#000000", "Training Point"="#0072B2", "Testing Function"="#D55E00")) +
          ggplot2::scale_fill_manual(name = "xx", values=c("CI (+/-2SD)"="#999999"))
      }
    }
  }

  p <- p +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_blank(),
      legend.spacing.x = ggplot2::unit(4, "pt")
    ) +
    ggplot2::labs(x ="Input position", y = "Model output")

  if ( is.null(dat3) ){
    if (type == 'points') {
      p <- p +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c("dashed", "blank"), shape = c(NA, 17), size = c(0.5, 2), alpha = c(0.9, 0.9))))
    } else if (type == 'line') {
      p <- p +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c("dashed", "solid"), shape = c(NA, NA), size = c(0.5, 0.5), alpha = c(0.9, 0.9))))
    }
  } else {
    if (type == 'points') {
      p <- p +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c("dashed", "blank", "blank"), shape = c(NA, 17, 19), size = c(0.5, 2, 2), alpha = c(0.9, 0.9, 0.9))))
    } else if (type == 'line') {
      p <- p +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c("dashed", "solid", "blank"), shape = c(NA, NA, 19), size = c(0.5, 0.5, 2), alpha = c(0.9, 0.9, 0.9))))
    }
  }
  return(p)
}
