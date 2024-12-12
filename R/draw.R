#' @title Validation and diagnostic plots for a sequential design
#'
#' @description
#'
#' `r new_badge("updated")`
#'
#' This function draws diagnostic and validation plots for a sequential design of a (D)GP emulator or a bundle of (D)GP emulators.
#'
#' @param object can be one of the following emulator classes:
#' * the S3 class `gp`.
#' * the S3 class `dgp`.
#' * the S3 class `bundle`.
#' @param type specifies the type of plot or visualization to generate:
#' - `"rmse"`: generates a trace plot of RMSEs, log-losses for DGP emulators with categorical likelihoods, or custom evaluation metrics specified via the `"eval"` argument in the `[design()]` function.
#' - `"design"`: shows visualizations of input designs created by the sequential design procedure.
#'
#' Defaults to `"rmse"`.
#' @param log a bool indicating whether to plot RMSEs, log-losses (for DGP emulators with categorical likelihoods), or custom evaluation metrics on a log scale when `type = "rmse"`.
#'     Defaults to `FALSE`.
#' @param emulator `r new_badge("updated")` an index or vector of indices of emulators packed in `object`. This argument is only used if `object` is an instance of the `bundle` class. When set to `NULL`, all
#'     emulators in the bundle are drawn. Defaults to `NULL`.
#' @param ... N/A.
#'
#' @return A `patchwork` object.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See design() for an example.
#' }
#' @md
#' @name draw
#' @export
draw <- function(object, ...){
  UseMethod("draw")
}

#' @rdname draw
#' @method draw gp
#' @export
draw.gp <- function(object, type = 'rmse', log = FALSE, ...){
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"gp") ){
    stop("'object' must be an instance of the 'gp' class.", call. = FALSE)
  }
  if ( !("design" %in% names(object)) ) stop("'object' must contain the 'design' slot created from design() to be used by draw().", call. = FALSE)

  if ( type == 'design' ) {
    design_data <- object$data$X
    total_N <- nrow(design_data)
    design_data <- as.data.frame(design_data)
    for (l in 1:length(design_data)) {
      names(design_data)[l] <- paste('X', l, sep="")
    }
    wave <- c()
    seq_N <- 0
    wave_N <- length(object$design)
    if ( "type" %in% names(object$design) ) wave_N <- wave_N - 1
    if ( "exclusion" %in% names(object$design) ) wave_N <- wave_N - 1
    if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) wave_N <- wave_N - 2
    for ( i in 1:wave_N ){
      Ni <- sum(object$design[[paste('wave',i,sep='')]]$enrichment)
      wave <- c(wave, rep(paste("wave",i,sep=""), Ni))
      seq_N <- seq_N + Ni
    }
    wave <- c(rep("Initial", total_N-seq_N), wave)
    design_data$Design <- wave

    p <- pair_design(design_data)

    p_patch <- patchwork::wrap_plots(p) +
      patchwork::plot_annotation(
      title = 'Sequential Design',
      caption =
        'Xi = Input dimension i of the GP emulator'
    )
  } else if ( type == 'rmse' ){
    total_N <- nrow(object$data$X)
    seq_N <- 0
    wave_N <- length(object$design)
    if ( "type" %in% names(object$design) ) wave_N <- wave_N - 1
    if ( "exclusion" %in% names(object$design) ) wave_N <- wave_N - 1
    cust <- object$design$type=='customized'
    if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) wave_N <- wave_N - 2
    for ( i in 1:wave_N ){
      Ni <- sum(object$design[[paste('wave',i,sep='')]]$enrichment)
      seq_N <- seq_N + Ni
    }
    init_N <- total_N - seq_N
    design_N <- c()
    rmse <- c()
    wave <- c()
    target <- c()
    target_wave <- c()
    for ( i in 1:wave_N ){
      Ni <- object$design[[paste('wave',i,sep='')]]$N
      if ( cust ){
        rmsei <- object$design[[paste('wave',i,sep='')]]$metric
      } else {
        rmsei <- object$design[[paste('wave',i,sep='')]]$rmse
      }
      Fi <- object$design[[paste('wave',i,sep='')]]$freq
      enrichi <- cumsum(object$design[[paste('wave',i,sep='')]]$enrichment)
      enrichi <- c(init_N, init_N + enrichi)
      step_Ni <- seq(0, Ni, Fi)
      if ( step_Ni[length(step_Ni)]!=Ni ) step_Ni <- c(step_Ni, Ni)
      design_Ni <- enrichi[step_Ni+1]
      design_N <- c(design_N, design_Ni)
      wave <- c(wave, rep(paste("wave",i,sep=" "), nrow(rmsei)))
      rmse <- rbind(rmse, rmsei)
      init_N <- design_Ni[length(design_Ni)]
      if ( 'target' %in% names(object$design[[paste('wave',i,sep='')]]) ){
        target <- c(target, object$design[[paste('wave',i,sep='')]]$target)
        target_wave <- c(target_wave, i)
      }
    }

    if ( !is.null(target) ){
      dat_target <- list()
      dat_target[["val"]] <- target
      dat_target[["Target"]] <- target_wave
      dat_target <- stats::aggregate(Target~val, dat_target, function(x) paste('wave',paste0(stats::na.omit(x), collapse = ","), sep=" "))
      dat_target <- dat_target[order(dat_target$val,decreasing = T),]
      dat_target <- as.data.frame(dat_target)
    } else {
      dat_target <- NULL
    }

    dat <- list()
    dat[["N"]] <- design_N
    dat[["rmse"]] <- rmse[,1]
    dat[["Design"]] <- wave
    p <- draw_seq_design(as.data.frame(dat),log = log, target = dat_target, cust = cust)
    p_patch <- patchwork::wrap_plots(p)
    if ( cust ){
      p_patch <- p_patch +
        patchwork::plot_annotation(
          title = 'Sequential Design Validation'
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    } else {
      p_patch <- p_patch +
        patchwork::plot_annotation(
          title = 'Sequential Design Validation',
          caption = 'RMSE = Root Mean Squared Error'
        ) +
        patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    }
  } else {
    stop("The provided 'type' is not supported.", call. = FALSE)
  }
  p_patch
}

#' @rdname draw
#' @method draw dgp
#' @export
draw.dgp <- function(object, type = 'rmse', log = FALSE, ...){
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"dgp") ){
    stop("'object' must be an instance of the 'dgp' class.", call. = FALSE)
  }
  if ( !("design" %in% names(object)) ) stop("'object' must contain the 'design' slot created from design() to be used by draw().", call. = FALSE)

  if ( type == 'design' ) {
    design_data <- object$data$X
    total_N <- nrow(design_data)
    design_data <- as.data.frame(design_data)
    for (l in 1:length(design_data)) {
      names(design_data)[l] <- paste('X', l, sep="")
    }
    wave <- c()
    seq_N <- 0
    wave_N <- length(object$design)
    if ( "type" %in% names(object$design) ) wave_N <- wave_N - 1
    if ( "exclusion" %in% names(object$design) ) wave_N <- wave_N - 1
    if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) wave_N <- wave_N - 2
    for ( i in 1:wave_N ){
      Ni <- sum(object$design[[paste('wave',i,sep='')]]$enrichment)
      wave <- c(wave, rep(paste("wave",i,sep=""), Ni))
      seq_N <- seq_N + Ni
    }
    wave <- c(rep("Initial", total_N-seq_N), wave)
    design_data$Design <- wave

    p <- pair_design(design_data)

    p_patch <- patchwork::wrap_plots(p) +
      patchwork::plot_annotation(
        title = 'Sequential Design',
        caption =
          'Xi = Input dimension i of the DGP emulator'
      )
  } else if ( type == 'rmse' ){
    if (object$constructor_obj$all_layer[[object$constructor_obj$n_layer]][[1]]$name == "Categorical") {
      is.categorical <- TRUE
    } else {
      is.categorical <- FALSE
    }
    total_N <- nrow(object$data$X)
    #output_D <- ncol(object$data$Y)
    seq_N <- 0
    wave_N <- length(object$design)
    if ( "type" %in% names(object$design) ) wave_N <- wave_N - 1
    if ( "exclusion" %in% names(object$design) ) wave_N <- wave_N - 1
    cust <- object$design$type=='customized'
    if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) wave_N <- wave_N - 2
    for ( i in 1:wave_N ){
      Ni <- sum(object$design[[paste('wave',i,sep='')]]$enrichment)
      seq_N <- seq_N + Ni
    }
    init_N <- total_N - seq_N
    design_N <- c()
    rmse <- c()
    wave <- c()
    target <- c()
    target_wave <- c()
    for ( i in 1:wave_N ){
      Ni <- object$design[[paste('wave',i,sep='')]]$N
      if ( cust ){
        rmsei <- object$design[[paste('wave',i,sep='')]]$metric
      } else {
        rmsei <- object$design[[paste('wave',i,sep='')]]$rmse
      }
      Fi <- object$design[[paste('wave',i,sep='')]]$freq
      enrichi <- cumsum(object$design[[paste('wave',i,sep='')]]$enrichment)
      enrichi <- c(init_N, init_N + enrichi)
      step_Ni <- seq(0, Ni, Fi)
      if ( step_Ni[length(step_Ni)]!=Ni ) step_Ni <- c(step_Ni, Ni)
      design_Ni <- enrichi[step_Ni+1]
      design_N <- c(design_N, design_Ni)
      wave <- c(wave, rep(paste("wave",i,sep=" "), nrow(rmsei)))
      rmse <- rbind(rmse, rmsei)
      init_N <- design_Ni[length(design_Ni)]
      if ( 'target' %in% names(object$design[[paste('wave',i,sep='')]]) ){
        target <- rbind(target, object$design[[paste('wave',i,sep='')]]$target)
        target_wave <- c(target_wave, i)
      }
    }
    output_D <- ncol(rmse)

    p_list <- list()
    for (l in 1:output_D ) {

      if ( !is.null(target) ){
        dat_target <- list()
        dat_target[["val"]] <- target[,l]
        dat_target[["Target"]] <- target_wave
        dat_target <- stats::aggregate(Target~val, dat_target, function(x) paste('wave',paste0(stats::na.omit(x), collapse = ","), sep=" "))
        dat_target <- dat_target[order(dat_target$val,decreasing = T),]
        dat_target <- as.data.frame(dat_target)
      } else {
        dat_target <- NULL
      }

      dat <- list()
      dat[["N"]] <- design_N
      dat[["rmse"]] <- rmse[,l]
      dat[["Design"]] <- wave
      p_list[[l]] <- draw_seq_design(as.data.frame(dat), log = log, target = dat_target, cust = cust, is.categorical = is.categorical)
      if ( output_D==ncol(object$data$Y) ){
        p_list[[l]] <- p_list[[l]] +
           ggplot2::ggtitle(sprintf("O%i", l)) +
           ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      }
    }
    p_patch <- patchwork::wrap_plots(p_list) +
      patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')
    if ( cust ){
      if ( output_D==ncol(object$data$Y) ){
        p_patch <- p_patch +
          patchwork::plot_annotation(
            title = 'Sequential Design Validation',
            caption = 'Oi = Output i of the DGP emulator'
          )
      } else {
        p_patch <- p_patch +
          patchwork::plot_annotation(
            title = 'Sequential Design Validation'
          )
      }
    } else {
      if ( is.categorical ){
        p_patch <- p_patch +
          patchwork::plot_annotation(
            title = 'Sequential Design Validation',
            caption = 'Oi = Output i of the DGP emulator'
          )
      } else {
        p_patch <- p_patch +
          patchwork::plot_annotation(
            title = 'Sequential Design Validation',
            caption = 'Oi = Output i of the DGP emulator
                   RMSE = Root Mean Squared Error'
          )
      }
    }
  } else {
    stop("The provided 'type' is not supported.", call. = FALSE)
  }
  p_patch
}


#' @rdname draw
#' @method draw bundle
#' @export
draw.bundle <- function(object, type = 'rmse', log = FALSE, emulator = NULL, ...){
  if ( is.null(pkg.env$dgpsi) ) {
    init_py(verb = F)
    if (pkg.env$restart) return(invisible(NULL))
  }
  #check class
  if ( !inherits(object,"bundle") ){
    stop("'object' must be an instance of the 'bundle' class.", call. = FALSE)
  }
  if ( !("design" %in% names(object)) ) stop("'object' must contain the 'design' slot created from design() to be used by draw().", call. = FALSE)

  if (is.null(emulator)){
    n_emulator <- length(grep("^emulator[0-9]+$", names(object), value = TRUE))
    emulator <- 1:n_emulator
  }

  if ( type == 'design' ) {
    p_list <- list()
    for (emu in 1:length(emulator)){
      design_data <- object$data$X[[paste('emulator', emulator[emu], sep="")]]
      total_N <- nrow(design_data)
      design_data <- as.data.frame(design_data)
      for (l in 1:length(design_data)) {
        names(design_data)[l] <- paste('X', l, sep="")
      }
      wave <- c()
      seq_N <- 0
      wave_N <- length(object$design)
      if ( "type" %in% names(object$design) ) wave_N <- wave_N - 1
      if ( "exclusion" %in% names(object$design) ) wave_N <- wave_N - 1
      if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) wave_N <- wave_N - 2
      for ( i in 1:wave_N ){
        Ni <- sum(object$design[[paste('wave',i,sep='')]]$enrichment[,emulator[emu]])
        wave <- c(wave, rep(paste("wave",i,sep=""), Ni))
        seq_N <- seq_N + Ni
      }
      wave <- c(rep("Initial", total_N-seq_N), wave)
      design_data$Design <- wave
      p_list[[emu]] <- pair_design(design_data) +
        ggplot2::ggtitle(sprintf("Emulator %i", emulator[emu])) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=10))

    }

    p_patch <- patchwork::wrap_plots(p_list) +
      patchwork::plot_layout() & ggplot2::theme(legend.position='bottom')

    p_patch <- p_patch +
      patchwork::plot_annotation(
        title = 'Sequential Design',
        caption = 'Xi = Input dimension i'
      )
  } else if ( type == 'rmse' ){
    p_list <- list()
    for (emu in 1:length(emulator)){
      is.categorical <- FALSE
      obj_k <- object[[paste('emulator',emulator[emu],sep='')]]
      if (inherits(obj_k,"dgp")){
        if (obj_k$constructor_obj$all_layer[[obj_k$constructor_obj$n_layer]][[1]]$name == "Categorical") {
          is.categorical <- TRUE
        }
      }

      total_N <- nrow(object$data$X[[paste('emulator', emulator[emu], sep="")]])
      seq_N <- 0
      wave_N <- length(object$design)
      if ( "type" %in% names(object$design) ) wave_N <- wave_N - 1
      if ( "exclusion" %in% names(object$design) ) wave_N <- wave_N - 1
      cust <- object$design$type=='customized'
      if ( "x_test" %in% names(object$design) & "y_test" %in% names(object$design) ) wave_N <- wave_N - 2
      for ( i in 1:wave_N ){
        Ni <- sum(object$design[[paste('wave',i,sep='')]]$enrichment[,emulator[emu]])
        seq_N <- seq_N + Ni
      }
      init_N <- total_N - seq_N
      design_N <- c()
      rmse <- c()
      wave <- c()
      target <- c()
      target_wave <- c()

      for ( i in 1:wave_N ){
        Ni <- object$design[[paste('wave',i,sep='')]]$N
        if ( cust ){
          if ( ncol(object$design[[paste('wave',i,sep='')]]$metric)==1 ) {
            emulator_tmp <- 1
          } else {
            emulator_tmp <- emulator[emu]
          }
        } else {
          emulator_tmp <- emulator[emu]
        }
        if ( cust ){
          rmsei <- object$design[[paste('wave',i,sep='')]]$metric[,emulator_tmp]
        } else {
          rmsei <- object$design[[paste('wave',i,sep='')]]$rmse[,emulator_tmp]
        }
        Fi <- object$design[[paste('wave',i,sep='')]]$freq
        enrichi <- cumsum(object$design[[paste('wave',i,sep='')]]$enrichment[,emulator[emu]])
        enrichi <- c(init_N, init_N + enrichi)
        step_Ni <- seq(0, Ni, Fi)
        if ( step_Ni[length(step_Ni)]!=Ni ) step_Ni <- c(step_Ni, Ni)
        design_Ni <- enrichi[step_Ni+1]
        design_N <- c(design_N, design_Ni)
        wave <- c(wave, rep(paste("wave",i,sep=" "), length(rmsei)))
        rmse <- c(rmse, rmsei)
        init_N <- design_Ni[length(design_Ni)]
        if ( 'target' %in% names(object$design[[paste('wave',i,sep='')]]) ){
          target <- c(target, object$design[[paste('wave',i,sep='')]]$target[emulator_tmp])
          target_wave <- c(target_wave, i)
        }
      }

      if ( !is.null(target) ){
        dat_target <- list()
        dat_target[["val"]] <- target
        dat_target[["Target"]] <- target_wave
        dat_target <- stats::aggregate(Target~val, dat_target, function(x) paste('wave',paste0(stats::na.omit(x), collapse = ","), sep=" "))
        dat_target <- dat_target[order(dat_target$val,decreasing = T),]
        dat_target <- as.data.frame(dat_target)
      } else {
        dat_target <- NULL
      }

      dat <- list()
      dat[["N"]] <- design_N
      dat[["rmse"]] <- rmse
      dat[["Design"]] <- wave
      p <- draw_seq_design(as.data.frame(dat), log = log, target = dat_target, cust = cust, is.categorical = is.categorical)

      if ( cust ){
        if (ncol(object$design[['wave1']]$metric)!=1) {
          p <- p +ggplot2::ggtitle( sprintf("Emulator %i", emulator[emu]) ) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=10))
        }
      } else {
        p <- p +ggplot2::ggtitle( sprintf("Emulator %i", emulator[emu]) ) +
          ggplot2::theme(plot.title = ggplot2::element_text(size=10))
      }
      p_list[[emu]] <- p
    }

    p_patch <- patchwork::wrap_plots(p_list) +
      patchwork::plot_layout(guides = 'collect') & ggplot2::theme(legend.position='bottom')

    if ( cust ){
      p_patch <- p_patch +
        patchwork::plot_annotation(
          title = 'Sequential Design Validation',
        )
    } else {
      if (is.categorical){
        p_patch <- p_patch +
          patchwork::plot_annotation(
            title = 'Sequential Design Validation'
          )
      } else {
      p_patch <- p_patch +
        patchwork::plot_annotation(
          title = 'Sequential Design Validation',
          caption = 'RMSE = Root Mean Squared Error'
        )
      }
    }
  } else {
    stop("The provided 'type' is not supported.", call. = FALSE)
  }
  p_patch
}

pair_design <- function(dat){
  c25 <- c("gray50", "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
           "#CAB2D6", "#FDBF6F", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
           "darkorange4", "brown")
  p <- ggplot2::ggplot(dat) +
    ggplot2::geom_point(ggplot2::aes_(x = ~.panel_x, y = ~.panel_y, shape=~Design, color =~Design)) +
    ggforce::facet_matrix(ggplot2::vars(names(dat)[1:(length(dat)-1)])) +
    ggplot2::scale_shape_manual(values = c(17,rep(16,length(unique(dat$Design))-1)), breaks=unique(dat$Design)) +
    ggplot2::scale_colour_manual(values = ggplot2::alpha(c25,0.7), breaks=unique(dat$Design)) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "right",
      legend.key.width = ggplot2::unit(0.5, "cm"),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 7, face = "bold"),
      legend.title.align=0.5
    )
  return(p)
}

draw_seq_design <- function(dat, log, target, cust, is.categorical = FALSE) {
  c24 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
           "#CAB2D6", "#FDBF6F", "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
           "darkorange4", "brown")
  p <- ggplot2::ggplot()

  if ( !is.null(target) ) {
    p <- p +
      ggplot2::geom_hline(data = target, mapping = ggplot2::aes_(yintercept=~val, group=~Target, linetype=~Target), alpha=0.8, color="gray20", size = 0.5) +
      ggplot2::scale_linetype_manual(values=c("dashed", "dotdash", "twodash", "dotted", "solid", "longdash"))
  }

  p <- p +
    ggplot2::geom_line(dat, mapping=ggplot2::aes_(x=~N, y=~rmse, group=~Design, color=~Design), alpha=0.8, stat = "unique") +
    ggplot2::geom_point(dat, mapping=ggplot2::aes_(x=~N, y=~rmse, group=~Design, color=~Design), size=1.5, alpha=0.8, stat = "unique") +
    ggplot2::scale_colour_manual(values = c24, breaks=unique(dat$Design))

  p <- p + ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.5, "cm"),
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 8.5)
    )

  if ( log ){
    p <- p + ggplot2::scale_y_continuous(trans='log10') +
      ggplot2::labs(x ="Number of design points", y = ifelse(cust, "Metric (log-scale)", ifelse(is.categorical, "Log-Loss (log-scale)", "RMSE (log-scale)" )))
  } else {
    p <- p +
      ggplot2::labs(x ="Number of design points", y = ifelse(cust, "Metric",  ifelse(is.categorical, "Log-Loss", "RMSE")))
  }
  return(p)
}





