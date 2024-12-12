#' @title Validation plots of a constructed GP, DGP, or linked (D)GP emulator
#'
#' @description
#'
#' `r new_badge("updated")`
#'
#' This function draws validation plots of a GP, DGP, or linked (D)GP emulator.
#'
#' @param x can be one of the following emulator classes:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `lgp`.
#' @param x_test same as that of [validate()].
#' @param y_test same as that of [validate()].
#' @param dim `r new_badge("updated")` if `dim = NULL`, the index of an emulator's input within the design will be shown on the x-axis in validation plots. Otherwise, `dim` indicates
#'     which dimension of an emulator's input will be shown on the x-axis in validation plots:
#' * If `x` is an instance of the `gp` of `dgp` class, `dim` is an integer.
#' * `r lifecycle::badge("deprecated")` If `x` is an instance of the `lgp` class created by [lgp()] without specifying the `struc` argument in data frame form, `dim` can be:
#'   1. an integer referring to the dimension of the global input to emulators in the first layer of a linked emulator system; or
#'   2. a vector of three integers referring to the dimension (specified by the third integer) of the global input to an emulator
#'   (specified by the second integer) in a layer (specified by the first integer) that is not the first layer of a linked emulator
#'   system.
#'
#'   **This option for linked (D)GP emulators is deprecated and will be removed in the next release.**
#' * `r new_badge("new")` If `x` is an instance of the `lgp` class created by [lgp()] with argument `struc` in data frame form, `dim` is an integer referring
#'   to the dimension of the global input to the linked emulator system.
#'
#' This argument is only used when `style = 1`. Defaults to `NULL`.
#' @param method same as that of [validate()].
#' @param sample_size same as that of [validate()].
#' @param style either `1` or `2`, indicating two different plotting styles for validation.
#' @param min_max a bool indicating if min-max normalization will be used to scale the testing output, RMSE, predictive mean and std from the
#'     emulator. Defaults to `TRUE`. This argument is not applicable to DGP emulators with categorical likelihoods.
#' @param normalize `r new_badge("new")` a bool indicating if normalization will be used to scale the counts in validation plots of DGP emulators with categorical
#'    likelihoods when `style = 2`.  Defaults to `TRUE`.
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
#'     individual points when the input of the emulator is one-dimensional and `style = 1`. This argument is not applicable to DGP emulators with
#'     categorical likelihoods. Defaults to `'points'`
#' @param verb a bool indicating if trace information on plotting will be printed during execution.
#'     Defaults to `TRUE`.
#' @param M `r new_badge("new")` same as that of [validate()].
#' @param force same as that of [validate()].
#' @param cores same as that of [validate()].
#' @param ... N/A.
#'
#' @return A `patchwork` object.
#'
#' @note
#' * [plot()] calls [validate()] internally to obtain validation results for plotting. However, [plot()] will not export the
#'   emulator object with validation results. Instead, it only returns the plotting object. For small-scale validations (i.e., small
#'   training or testing data points), direct execution of [plot()] works well. However, for moderate- to large-scale validation,
#'   it is recommended to first run [validate()] to obtain and store validation results in the emulator object, and then supply the
#'   object to [plot()]. [plot()] checks the object's `loo` and `oos` slots prior to calling [validate()] and will not perform further calculation if the required information is already stored.
#' * [plot()] will only use stored OOS validation if `x_test` and `y_test` are identical to those used by [validate()] to produce the data contained in the object's `oos` slot, otherwise [plot()] will re-evaluate OOS validation before plotting.
#' * The returned [patchwork::patchwork] object contains the [ggplot2::ggplot2] objects. One can modify the included individual ggplots
#'   by accessing them with double-bracket indexing. See <https://patchwork.data-imaginist.com/> for further information.
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#' @md
#' @name plot
NULL
#' @importFrom rlang .data

#' @rdname plot
#' @method plot dgp
#' @export
plot.dgp <- function(x, x_test = NULL, y_test = NULL, dim = NULL, method = NULL, sample_size = 50, style = 1, min_max = TRUE, normalize = TRUE, color = 'turbo', type = 'points', verb = TRUE, M = 50, force = FALSE, cores = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( style!=1&style!=2 ) stop("'style' must be either 1 or 2.", call. = FALSE)
  if ( type!='points'&type!='line' ) stop("'type' must be either 'points' or 'line'.", call. = FALSE)
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }

  if ( isTRUE(verb) ) message("Validating and computing ...", appendLF = FALSE)
  results <- validate(object = x, x_test = x_test, y_test = y_test, method = method, sample_size = sample_size, verb = FALSE, M = M, force = force, cores = cores)
  if ( isTRUE(verb) ) message(" done")

  # For LOO
  if ( is.null(x_test) & is.null(y_test) ){
    if ( isTRUE(verb) ) message("Post-processing LOO results ...", appendLF = FALSE)
    loo_res <- results$loo
    if ( "label" %in% names(loo_res) ){
      is.categorical <- TRUE
    } else {
      is.categorical <- FALSE
    }
    p_list <- list()
    if ( style==1 ){
      # create indices
      if ( ncol(loo_res$x_train)==1 ){
        idx <- loo_res$x_train[,1]
        isdup <- TRUE
        if (!is.null(dim)) dim <- 1
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
      if (is.categorical){
        samples_df <- as.data.frame(loo_res$label)
        colnames(samples_df) <- paste0("sample_", 1:ncol(samples_df))
        dat <- cbind(data.frame(idx = idx, y_validate = loo_res$y_train[,1]), samples_df)
        p_list[[1]] <- plot_style_1_classify(dat, dim, isdup) +
          ggplot2::ggtitle(sprintf("Log Loss = %.4f", loo_res$log_loss)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      } else {
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
      }
    } else if ( style==2 ) {
      if (is.categorical){
        rep1 <- pkg.env$np$unique(loo_res$x_train, return_index=TRUE, axis=0L)
        rep2 <- pkg.env$np$unique(loo_res$x_train, return_inverse=TRUE, axis=0L)
        idx <- seq(1,length(rep1[[2]]))[pkg.env$np$argsort(pkg.env$np$argsort(rep1[[2]]))+1][rep2[[2]]+1]
        samples_df <- as.data.frame(loo_res$label)
        colnames(samples_df) <- paste0("sample_", 1:ncol(samples_df))
        dat <- cbind(data.frame(idx = idx, y_validate = loo_res$y_train[,1]), samples_df)
        p_list[[1]] <- plot_style_2_classify(dat, color, normalize) +
          ggplot2::ggtitle(sprintf("Log Loss = %.4f", loo_res$log_loss)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      } else {
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
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Plotting ...", appendLF = FALSE)
    if (is.categorical){
      p_patch <- patchwork::wrap_plots(p_list) +
        patchwork::plot_annotation(
          title = 'Leave-One-Out (LOO) Cross Validation'
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    } else {
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
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")
    p_patch
    # For OOS
  } else if (!is.null(x_test) & !is.null(y_test)) {
    if ( isTRUE(verb) ) message("Post-processing OOS results ...", appendLF = FALSE)
    oos_res <- results$oos
    if ( "label" %in% names(oos_res) ){
      is.categorical <- TRUE
    } else {
      is.categorical <- FALSE
    }
    p_list <- list()
    if ( style==1 ){
      if (is.categorical){
        if ( ncol(oos_res$x_test)==1 ){
          idx <- oos_res$x_test[,1]
          isdup <- TRUE
          if (!is.null(dim)) dim <- 1
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
        }
        samples_df <- as.data.frame(oos_res$label)
        colnames(samples_df) <- paste0("sample_", 1:ncol(samples_df))
        dat <- cbind(data.frame(idx = idx, y_validate = oos_res$y_test[,1]), samples_df)
        p_list[[1]] <- plot_style_1_classify(dat, dim, isdup) +
          ggplot2::ggtitle(sprintf("Log Loss = %.4f", oos_res$log_loss)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      } else {
        # If input is 1d
        if ( ncol(oos_res$x_test)==1 ){
          if ( "mean" %in% names(oos_res) ){
            method <- "mean_var"
          } else if ( "median" %in% names(oos_res) ){
            method <- "sampling"
          }
          if (!is.null(dim)) dim <- 1
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
            res <- results$emulator_obj$predict(x_range, method = method, m = M)
          } else {
            res <- results$emulator_obj$ppredict(x_range, method = method, m = M, core_num = cores)
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
      }
    } else if ( style==2 ) {
      if (is.categorical){
        rep1 <- pkg.env$np$unique(oos_res$x_test, return_index=TRUE, axis=0L)
        rep2 <- pkg.env$np$unique(oos_res$x_test, return_inverse=TRUE, axis=0L)
        idx <- seq(1,length(rep1[[2]]))[pkg.env$np$argsort(pkg.env$np$argsort(rep1[[2]]))+1][rep2[[2]]+1]
        samples_df <- as.data.frame(oos_res$label)
        colnames(samples_df) <- paste0("sample_", 1:ncol(samples_df))
        dat <- cbind(data.frame(idx = idx, y_validate = oos_res$y_test[,1]), samples_df)
        p_list[[1]] <- plot_style_2_classify(dat, color, normalize) +
          ggplot2::ggtitle(sprintf("Log Loss = %.4f", oos_res$log_loss)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      } else {
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
    }
    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")

    if ( isTRUE(verb) ) message("Plotting ...", appendLF = FALSE)
    if (is.categorical){
      p_patch <- patchwork::wrap_plots(p_list) +
        patchwork::plot_annotation(
          title = 'Out-Of-Sample (OOS) Cross Validation'
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    } else {
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
    }

    if ( isTRUE(verb) ) Sys.sleep(0.5)
    if ( isTRUE(verb) ) message(" done")
    p_patch
  }
}

#' @rdname plot
#' @method plot lgp
#' @export
plot.lgp <- function(x, x_test = NULL, y_test = NULL, dim = NULL, method = NULL, sample_size = 50, style = 1, min_max = TRUE, color = 'turbo', type = 'points', M = 50, verb = TRUE, force = FALSE, cores = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }

  if ( "metadata" %in% names(x$specs) ){
    if ( !("emulator_obj" %in% names(x)) ){
      stop("'object' is not activated for plotting. Please set `activate = TRUE` in `lgp()` to activate the emulator.", call. = FALSE)
    }
    is.df <- TRUE
  } else {
    is.df <- FALSE
  }

  if ( style!=1&style!=2 ) stop("'style' must be either 1 or 2.", call. = FALSE)
  if ( type!='points'&type!='line' ) stop("'type' must be either 'points' or 'line'.", call. = FALSE)
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }

  if ( isTRUE(verb) ) message("Validating and computing ...", appendLF = FALSE)
  results <- validate(object = x, x_test = x_test, y_test = y_test, method = method, sample_size = sample_size, verb = FALSE, M = M, force = force, cores = cores)
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
      if (!is.null(dim)) dim <- 1
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
        res <- results$emulator_obj$predict(x_range, method = method, m = M)
      } else {
        res <- results$emulator_obj$ppredict(x_range, method = method, m = M, core_num = cores)
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
            if ( dim[1]<2 ) stop("The first element of 'dim' must be greater than or equal to 2. See documentation for details.", call. = FALSE)
            target_emulator_input <- oos_res$x_test[[dim[1]]][[dim[2]]]
            if ( is.null(target_emulator_input) ) stop("The emulator specified by the first two elements of 'dim' has no global input.", call. = FALSE)
            idx <- target_emulator_input[,dim[3]]
          } else {
            stop("'dim' must be a vector of length 1 or 3. See documentation for details.", call. = FALSE)
          }
        } else {
          if ( length(dim)==1 ){
            idx <- oos_res$x_test[,dim]
            if (!is.df) dim <- c(1,dim)
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
plot.gp <- function(x, x_test = NULL, y_test = NULL, dim = NULL, method = NULL, sample_size = 50, style = 1, min_max = TRUE, color = 'turbo', type = 'points', verb = TRUE, M = 50, force = FALSE, cores = 1, ...) {
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  if ( style!=1&style!=2 ) stop("'style' must be either 1 or 2.", call. = FALSE)
  if ( type!='points'&type!='line' ) stop("'type' must be either 'points' or 'line'.", call. = FALSE)
  if( !is.null(cores) ) {
    cores <- as.integer(cores)
    if ( cores < 1 ) stop("'cores' must be >= 1.", call. = FALSE)
  }

  if ( isTRUE(verb) ) message("Validating and computing ...", appendLF = FALSE)
  results <- validate(object = x, x_test = x_test, y_test = y_test, method = method, sample_size = sample_size, verb = FALSE, M = M, force = force, cores = cores)
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
        if (!is.null(dim)) dim <- 1
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
        if (!is.null(dim)) dim <- 1
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
          res <- results$emulator_obj$predict(x_range, method = method, sample_size=500L, m = M)
        } else {
          res <- results$emulator_obj$ppredict(x_range, method = method, sample_size=500L, m = M, core_num = cores)
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

plot_style_1_classify <- function(dat, dim, isdup) {
  # Handle duplication logic
  if (isTRUE(isdup)) {
    dup <- duplicated(dat$idx)
  } else {
    dup <- as.logical(rep(0, length(dat$idx)))
  }

  # Set x-axis label based on dim input
  if (is.null(dim)) {
    x_lab <- "Input position"
  } else {
    x_lab <- sprintf("Input dimension %i", dim)
  }

  # Reshape the dataframe to long format (excluding duplicates)
  df_long <- reshape2::melt(dat[!dup, ], id.vars = c("idx", "y_validate"),
                            variable.name = "sample",
                            value.name = "predicted_value")

  if (is.character(dat$y_validate)) {
    dat$y_validate <- as.factor(dat$y_validate)
  }

  if (is.character(df_long$predicted_value)) {
    df_long$predicted_value <- as.factor(df_long$predicted_value)
  }

  prop_df <- dplyr::group_by(df_long, .data$idx, .data$predicted_value)
  prop_df <- dplyr::summarize(prop_df, count = dplyr::n(), .groups = 'drop')
  prop_df <- dplyr::group_by(prop_df, .data$idx)
  prop_df <- dplyr::mutate(prop_df, proportion = .data$count / sum(.data$count))
  prop_df <- dplyr::ungroup(prop_df)

  # Create the ggplot
  ggplot2::ggplot() +

    # Plot predicted value proportions from prop_df
    ggplot2::geom_tile(data = prop_df,
                       ggplot2::aes(x = .data$idx, y = .data$predicted_value, fill = .data$proportion),
                       height = 0.1,
                       alpha = 0.95) +   # Keep the height fixed

    #ggplot2::geom_point(data = prop_df,
    #                    ggplot2::aes(x = .data$idx, y = .data$predicted_value, color = .data$proportion), size = 2.5,
    #                    shape = 15, alpha = 0.85) +

    # Use a gradient color scale for proportions (continuous)
    ggplot2::scale_fill_viridis_c("Predicted Label Proportion", option = 'G',
                                   direction = -1,
                                   breaks=seq(0,1,0.2),
                                   labels=c('0.0','0.2','0.4','0.6','0.8','1.0'),
                                   limits = c(0, 1)) +

    # Plot true values (y_validate) from the original dataset 'dat'
    ggplot2::geom_point(data = dat,
                        ggplot2::aes(x = .data$idx, y = .data$y_validate, shape = "True labels"),
                        color = "#FFD700", size = 1.1) +

    # Use manual shape scale to differentiate true labels
    ggplot2::scale_shape_manual(values = c("True labels" = 16), name = "True Labels") +

    # Labels, axis, and theme settings
    ggplot2::labs(x = x_lab, y = "Model output") +

    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1, "cm"),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 7, hjust = 0.5)
    ) +

    # Use guides() for color and shape, hide the item text for True Labels
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(title.position = "top"),
      shape = ggplot2::guide_legend(override.aes = list(size = 3),
                                    title.position = "top",
                                    title = "True Labels",
                                    label = FALSE)  # Hide the item label
    )
}

plot_style_2_classify <- function(dat, color, normalize) {
  df_long <- reshape2::melt(dat, id.vars = c("idx", "y_validate"), variable.name = "sample", value.name = "predicted_value")

  # Group by idx to collect true labels and predicted values
  df_grouped <-
    dplyr::summarize(dplyr::group_by(df_long, .data$idx),
                     true_labels = list(unique(.data$y_validate)),          # Collect all true labels for each idx
                     predicted_values = list(.data$predicted_value) # Collect all predicted values for each idx
    )

  # Unnest the predicted values and compare each predicted value to the true labels
  df_expanded <-
    dplyr::mutate(dplyr::rowwise(tidyr::unnest(df_grouped, "predicted_values")), match = ifelse(.data$predicted_values %in% unlist(.data$true_labels), "Match", "No Match"))  # Compare each predicted value

  # Unnest the true labels so that we can attribute mismatches to all true labels of the idx
  df_confusion <-
    dplyr::mutate(tidyr::unnest(df_expanded, "true_labels"), final_label = ifelse(match == "No Match", .data$true_labels, .data$predicted_values))  # Mismatches contribute to all true labels

  # Filter out rows where match is marked as "Match" but true and predicted values are not identical
  df_confusion_filtered <-
    dplyr::filter( df_confusion, !(.data$match == "Match" & .data$predicted_values != .data$true_labels))

  # Create confusion matrix by counting how often each predicted value contributes to each true label
  conf_matrix <-
    dplyr::count(df_confusion_filtered, .data$true_labels, .data$predicted_values)

  # Get all unique true labels and predicted values to ensure zero-counts are included
  all_combinations <- expand.grid(
    true_labels = unique(df_confusion_filtered$true_labels),
    predicted_values = unique(df_confusion_filtered$predicted_values)
  )

  # Merge the actual confusion matrix with all possible combinations to include zero counts
  conf_matrix_complete <-
    tidyr::replace_na(dplyr::left_join(all_combinations, conf_matrix, by = c("true_labels", "predicted_values")), list(n = 0))  # Replace NA with 0 for missing counts

  if (normalize){
    conf_matrix_complete$n <- conf_matrix_complete$n / sum(conf_matrix_complete$n)
  }
  # Plot the confusion matrix using ggplot2
  p <- ggplot2::ggplot(conf_matrix_complete, ggplot2::aes_(x =~predicted_values, y =~true_labels, fill =~n)) +
    ggplot2::geom_tile(color = "black") +  # Add black border to each tile
    ggplot2::geom_text(ggplot2::aes_(
      label = ifelse(conf_matrix_complete$n %% 1 == 0,  # Check if the value is an integer
                     as.character(conf_matrix_complete$n),  # If integer, show as is
                     sprintf("%.2f", conf_matrix_complete$n)  # Otherwise, format to 2 decimal places
      )
    ), color = "white", show.legend = FALSE)+  # Add count labels
    ggplot2::labs(x = "Predicted label", y = "True label")

  if (normalize){
    p <- p + ggplot2::scale_fill_viridis_c("Normalized Count", option = color, breaks=seq(0,1,0.2), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'), limits = c(0, 1))
  } else {
    p <- p + ggplot2::scale_fill_viridis_c("Count", option = color)
  }

  p <- p +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5),  # Rotate x-axis labels for better readability
      axis.text.y = ggplot2::element_text(angle = 0, hjust = 0.5),  # Keep y-axis labels aligned
      axis.ticks = ggplot2::element_line(color = "black"),  # Add ticks to boundaries of cells
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1, "cm"),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 7),
      legend.title.align=0.5
    ) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(title.position="top"))

  return(p)
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
