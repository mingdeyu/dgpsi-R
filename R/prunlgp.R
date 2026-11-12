#' Select internal connections in a linked emulator
#'
#' Select input subsets for scalar component emulators using validation
#' predictions from the complete linked network. The initial network remains
#' eligible throughout the search.
#'
#' @param object An `lgp` object defining the initial network.
#' @param emulators The original list of fitted scalar `gp` or `dgp` components.
#'   Component IDs must match `object`. Their `data$X` columns use the original
#'   `To_Input` numbering. Multi-output components must first be split into
#'   scalar emulators.
#' @param x Numeric matrix of validation global inputs, in the original global
#'   column numbering. No observed internal outputs are supplied for prediction.
#' @param y Named list of validation response vectors, keyed by emulator ID.
#'   Only these outputs contribute to scoring and protection.
#' @param refit Function `function(model, keep, seed)` returning a newly fitted
#'   emulator with the same ID and response as `model`, trained on
#'   `model$data$X[, keep, drop = FALSE]`. `keep` refers to original columns.
#'   Set the fitting budget, kernel and other parameters explicitly in this
#'   function. It must not use validation responses or modify its input model.
#' @param weights Nonnegative scoring weights, named like `y`. Defaults to
#'   equal weights. At least one weight must be positive.
#' @param scale Positive output scales, named like `y`. Defaults to the sample
#'   standard deviations of component training responses. Supply a meaningful
#'   positive scale explicitly for constant training responses.
#' @param objective Optional function `function(rmse, scale, weights)` returning
#'   one finite numeric loss to minimise. Defaults to weighted mean NRMSE.
#' @param protect Optional nonnegative relative RMSE tolerance, either scalar
#'   or named like `y`. Each scored output must remain within this tolerance of
#'   its initial validation RMSE. A zero initial error permits only numerical
#'   tolerance. `NULL` applies no per-output protection.
#' @param groups Optional named character vector mapping emulator IDs to shared
#'   fitting groups. Copies in a group must have identical training data and
#'   original input dimensions. Deleting an input removes that dimension from
#'   every copy, including an origin copy supplied through `Global`.
#' @param search Either `"subsets"`, visiting each group and evaluating all
#'   original internal-input subsets, or `"backward"`, evaluating single-input
#'   deletions from the current network. Both use the system validation loss.
#' @param visit_order Group IDs in the desired visiting order for `"subsets"`.
#'   Defaults to first occurrence in topological order.
#' @param max_evals Maximum network evaluations, including the initial and
#'   all-internal-inputs-removed networks. Defaults to 200.
#' @param max_rounds Maximum full passes (subsets) or accepted deletions
#'   (backward). Defaults to one pass for subsets and unlimited for backward.
#' @param min_improvement Minimum absolute loss reduction for a non-tied search
#'   transition. Defaults to zero. The final choice compares all evaluated
#'   admissible networks, including smaller improvements.
#' @param tolerance Absolute numerical tolerance for score ties and protection.
#'   Ties prefer fewer internal connections, then the lexicographic retained-input key.
#' @param seed Positive integer controlling fitting, linking and prediction.
#' @param B Number of linked imputations; defaults to `object$specs$B`.
#' @param M Prediction conditioning-set size; defaults to 50.
#' @param workers Number of PSOCK processes for independent candidates.
#'   Functions supplied by the caller must be serialisable and self-contained.
#' @param threads Python numerical threads per process. Defaults to one.
#' @param worker_init Optional zero-argument function initialising Python in
#'   each worker, for example `function() init_py(dgpsi_ver = "dev")`.
#'   The default calls [init_py()] with its default environment.
#' @param verb Print accepted moves and completion information.
#'
#' @return A list of class `prunlgp`, containing `model` (an ordinary `lgp`),
#'   `initial_model` (the initial network activated with the same settings),
#'   `emulators`, `initial_structure`, `structure`, `kept_inputs`,
#'   `initial_score`, `score`, `initial_rmse`, `rmse`, `validation`, `history`,
#'   `settings` and `stop_reason`. `predict()` and `summary()` delegate to the
#'   selected linked model. The history includes failed and rejected candidates.
#'
#' @details
#' Every candidate is fitted from the original component training data. Other
#'   fitted parameters stay fixed, and predictions are recomputed through the
#'   entire network using native linked uncertainty propagation. Component
#'   copies can be tied through `groups`; candidate inputs are then the union
#'   of dimensions receiving internal predictions in any copy.
#'
#'   An additional empty-internal-input candidate is evaluated without changing
#'   the starting point of the search. A component must retain at least one
#'   input; a candidate removing its last input is recorded as failed. Global
#'   inputs outside grouped internal dimensions are always retained.
#'
#'   The returned network has the lowest admissible loss among evaluated
#'   networks, with ties resolved as above. This is a validation search, not a
#'   guarantee of a global optimum or lower test error. Use a separate test set.
#'   Common fitting seeds for a group and common linking/prediction seeds reduce
#'   random variation; different graph shapes need not consume identical random
#'   draws. R's random state is restored on exit. The Python random state and
#'   thread setting are set explicitly and are not restored.
#'
#'   Unlike [prune()], which operates inside a DGP, this function removes
#'   dependencies between scalar component emulators.
#' @examples
#' \dontrun{
#' # See vignette("prunlgp") for a complete executable example.
#' result <- prunlgp(network, components, x_validation, y_validation,
#'   refit = function(model, keep, seed) {
#'     set_seed(seed)
#'     gp(model$data$X[, keep, drop = FALSE], model$data$Y,
#'        name = "sexp", id = model$id, verb = FALSE)
#'   })
#' predict(result, x_test, full_layer = TRUE)
#' summary(result)
#' result$history
#' }
#' @md
#' @export
prunlgp <- function(object, emulators, x, y, refit, weights = NULL,
                    scale = NULL, objective = NULL, protect = NULL,
                    groups = NULL, search = c("subsets", "backward"),
                    visit_order = NULL, max_evals = 200L, max_rounds = NULL,
                    min_improvement = 0, tolerance = 1e-10, seed = 1L,
                    B = object$specs$B, M = 50L, workers = 1L,
                    threads = 1L, worker_init = NULL, verb = TRUE) {
  search <- match.arg(search)
  if (!inherits(object, "lgp")) stop("object must be an lgp emulator.")
  if (!is.function(refit)) stop("refit must be a function(model, keep, seed).")
  if (!is.null(worker_init) && !is.function(worker_init)) stop("worker_init must be a function.")
  if (is.null(objective)) objective <- function(rmse, scale, weights)
    sum(weights * rmse / scale) / sum(weights)
  if (!is.function(objective)) stop("objective must be a function.")
  integer_arg <- function(z, name, lower = 1) {
    if (length(z) != 1L || !is.finite(z) || z < lower ||
        z != floor(z) || z > .Machine$integer.max)
      stop(name, " must be an integer >= ", lower, ".")
    as.integer(z)
  }
  seed <- integer_arg(seed, "seed"); workers <- integer_arg(workers, "workers")
  threads <- integer_arg(threads, "threads"); B <- integer_arg(B, "B")
  M <- integer_arg(M, "M"); max_evals <- integer_arg(max_evals, "max_evals")
  if (is.null(max_rounds)) max_rounds <- if (search == "subsets") 1 else Inf
  if (!identical(max_rounds, Inf)) max_rounds <- integer_arg(max_rounds, "max_rounds")
  for (z in list(min_improvement, tolerance))
    if (length(z) != 1L || !is.finite(z) || z < 0) stop("Invalid loss tolerance.")
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) old_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_rng) assign(".Random.seed", old_rng, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  ids <- vapply(emulators, function(m) as.character(m$id), character(1))
  if (!length(ids) || anyDuplicated(ids) || any(ids == "Global"))
    stop("Components need unique IDs other than Global.")
  names(emulators) <- ids
  for (m in emulators) {
    if (!inherits(m, c("gp", "dgp")) || ncol(as.matrix(m$data$Y)) != 1L)
      stop("Use scalar gp/dgp components; split multi-output emulators first.")
    if (!is.numeric(m$data$X) || !is.numeric(m$data$Y) ||
        any(!is.finite(m$data$X)) || any(!is.finite(m$data$Y)))
      stop("Component training data must be finite numeric matrices.")
  }
  struc <- object$specs$struc
  if (!setequal(ids, object$specs$metadata$Emulator))
    stop("emulators must match the initial object's component IDs.")
  # Reuse the package's dimension, coverage and acyclicity checks.
  struc <- validate_emulator_data(struc, emulators)
  struc$From_Emulator <- as.character(struc$From_Emulator)
  struc$To_Emulator <- as.character(struc$To_Emulator)
  struc <- struc[order(match(struc$To_Emulator, ids), struc$To_Input), ]
  rownames(struc) <- NULL
  x <- as.matrix(x)
  if (!is.numeric(x) || !nrow(x) || any(!is.finite(x)) ||
      ncol(x) < max(struc$From_Output[struc$From_Emulator == "Global"]))
    stop("x must contain all original global input columns and finite values.")
  if (!is.list(y) || !length(y) || is.null(names(y)) || anyDuplicated(names(y)) ||
      any(!names(y) %in% ids)) stop("y must be a named list keyed by component ID.")
  y <- lapply(y, function(z) {
    if (!is.numeric(z) || length(z) != nrow(x) || any(!is.finite(z)))
      stop("Each y vector must match the validation rows and be finite.")
    as.numeric(z)
  })
  scored <- names(y)
  named_values <- function(z, default, name, positive = FALSE) {
    if (is.null(z)) z <- default
    if (length(z) == 1L && is.null(names(z))) z <- stats::setNames(rep(z, length(scored)), scored)
    if (is.null(names(z)) || anyDuplicated(names(z)) || !setequal(names(z), scored))
      stop(name, " must be named like y.")
    z <- z[scored]
    if (!is.numeric(z) || any(!is.finite(z)) || any(if (positive) z <= 0 else z < 0))
      stop(name, " contains invalid values; constant training responses need an explicit scale.")
    z
  }
  weights <- named_values(weights, stats::setNames(rep(1, length(scored)), scored), "weights")
  if (sum(weights) <= 0) stop("At least one weight must be positive.")
  scale <- named_values(scale, vapply(emulators[scored], function(m)
    stats::sd(as.numeric(m$data$Y)), numeric(1)), "scale", TRUE)
  if (!is.null(protect)) protect <- named_values(protect, NULL, "protect")
  if (is.null(groups)) groups <- stats::setNames(ids, ids)
  if (!is.character(groups) || anyNA(groups) || any(!nzchar(groups)) ||
      is.null(names(groups)) || anyDuplicated(names(groups)) || !setequal(names(groups), ids))
    stop("groups must map every component ID to a nonempty group name.")
  groups <- groups[ids]
  members <- split(ids, groups)
  for (g in members) for (id in g) {
    a <- emulators[[g[1]]]$data; b <- emulators[[id]]$data
    if (!identical(a, b)) stop("Copies in a fitting group must have identical training data.")
  }
  inputs <- lapply(emulators, function(m) seq_len(ncol(m$data$X)))
  slots <- lapply(members, function(g) sort(unique(struc$To_Input[
    struc$To_Emulator %in% g & struc$From_Emulator != "Global"])))
  topo <- object$specs$metadata
  default_order <- unique(unname(groups[as.character(topo$Emulator[order(topo$Layer, topo$Pos_in_Layer)])]))
  if (is.null(visit_order)) visit_order <- default_order
  if (anyDuplicated(visit_order) || !setequal(visit_order, names(members)))
    stop("visit_order must contain each fitting group exactly once.")
  set_seed(seed); set_thread_num(threads)
  state <- list(models = lapply(emulators, serialize, light = FALSE), struc = struc,
    ids = ids, inputs = inputs, members = members, x = x, y = y, weights = weights,
    scale = scale, objective = objective, refit = refit, seed = seed, B = B, M = M,
    threads = threads)
  pool <- NULL
  if (workers > 1L) {
    pool <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(pool), add = TRUE)
    libs <- .libPaths()
    initialise_worker <- function(paths, initialise) {
      .libPaths(paths); loadNamespace("dgpsi")
      if (is.null(initialise)) dgpsi::init_py(verb = FALSE) else initialise()
      NULL
    }
    # Avoid loading a namespace from the worker's default library before its
    # library paths have been set (the closure environment is serialised too).
    environment(initialise_worker) <- baseenv()
    parallel::clusterCall(pool, initialise_worker, libs, worker_init)
  }
  history <- list(); results <- list(); cache <- new.env(parent = emptyenv())
  evaluate <- function(candidates, stage) {
    keys <- vapply(candidates, .prunlgp_key, character(1))
    fresh <- !duplicated(keys) & !vapply(keys, exists, logical(1), envir = cache, inherits = FALSE)
    todo <- which(fresh)
    todo <- utils::head(todo, max(0L, max_evals - length(history)))
    if (length(todo)) {
      jobs <- candidates[todo]
      ans <- if (is.null(pool)) lapply(jobs, .prunlgp_evaluate, state = state)
        else parallel::parLapply(pool, jobs, function(job, state)
          utils::getFromNamespace(".prunlgp_evaluate", "dgpsi")(job, state), state = state)
      for (j in seq_along(todo)) {
        k <- todo[j]; r <- ans[[j]]; n <- length(history) + 1L
        r$key <- keys[k]; r$keep <- candidates[[k]]; r$index <- n
        results[[n]] <<- r; assign(keys[k], n, envir = cache)
        history[[n]] <<- data.frame(evaluation = n, stage = stage, key = keys[k],
          loss = if (r$ok) r$loss else NA_real_, internal_links = r$n_edges,
          admissible = FALSE, accepted = FALSE, selected = FALSE,
          error = if (r$ok) "" else r$error, stringsAsFactors = FALSE)
        for (id in names(y)) history[[n]][[paste0("RMSE_", id)]] <<-
          if (r$ok) unname(r$rmse[id]) else NA_real_
      }
    }
    unique(vapply(keys[vapply(keys, exists, logical(1), envir = cache, inherits = FALSE)],
      get, integer(1), envir = cache, inherits = FALSE))
  }
  initial <- evaluate(list(inputs), "initial")[1]
  if (!results[[initial]]$ok) stop("Initial network failed: ", results[[initial]]$error)
  base <- results[[initial]]$rmse
  admissible <- function(i) results[[i]]$ok && (is.null(protect) ||
    all(results[[i]]$rmse <= base * (1 + protect) + tolerance))
  better <- function(i, j) {
    if (!admissible(i)) return(FALSE)
    a <- results[[i]]; b <- results[[j]]
    if (a$loss < b$loss - max(tolerance, min_improvement)) return(TRUE)
    abs(a$loss - b$loss) <= tolerance &&
      (a$n_edges < b$n_edges || (a$n_edges == b$n_edges && a$key < b$key))
  }
  choose <- function(indices, current) {
    for (i in indices) if (better(i, current)) current <- i
    current
  }
  empty <- inputs
  for (g in names(members)) for (id in members[[g]]) empty[[id]] <- setdiff(inputs[[id]], slots[[g]])
  evaluate(list(empty), "empty")
  current <- initial; round <- 0L; reason <- "no_improvement"
  repeat {
    if (length(history) >= max_evals) {reason <- "evaluation_budget"; break}
    if (round >= max_rounds) {reason <- "round_budget"; break}
    round <- round + 1L; before <- current
    if (search == "subsets") {
      for (g in visit_order) {
        n <- length(slots[[g]])
        # Enumerate lazily in batches so the budget also bounds memory use.
        if (!n) next
        if (n > 30L) stop("Subset search supports at most 30 candidate inputs per group; use backward.")
        start <- 0; incumbent <- current; context <- results[[current]]$keep
        while (start < 2^n && length(history) < max_evals) {
          masks <- seq.int(start, min(2^n - 1, start + max(1L, workers) - 1L))
          candidates <- lapply(masks, function(mask) {
            k <- context; retained <- slots[[g]][as.logical(intToBits(as.integer(mask))[seq_len(n)])]
            for (id in members[[g]]) k[[id]] <- sort(c(setdiff(inputs[[id]], slots[[g]]), retained))
            k
          })
          incumbent <- choose(evaluate(candidates, paste0("round", round, ":", g)), incumbent)
          start <- start + length(masks)
        }
        current <- incumbent
        history[[current]]$accepted <- TRUE
        if (verb) message("Round ", round, ", group ", g, ": loss ", signif(results[[current]]$loss, 6))
        if (length(history) >= max_evals) break
      }
    } else {
      candidates <- list(); context <- results[[current]]$keep
      for (g in visit_order) for (slot in intersect(context[[members[[g]][1]]], slots[[g]])) {
        k <- context
        for (id in members[[g]]) k[[id]] <- setdiff(k[[id]], slot)
        candidates[[length(candidates) + 1L]] <- k
      }
      if (!length(candidates)) {reason <- "no_internal_inputs"; break}
      current <- choose(evaluate(candidates, paste0("round", round)), current)
      history[[current]]$accepted <- TRUE
    }
    if (current == before) break
  }
  eligible <- which(vapply(seq_along(results), admissible, logical(1)))
  losses <- vapply(results[eligible], `[[`, numeric(1), "loss")
  tied <- eligible[losses <= min(losses) + tolerance]
  best <- tied[order(vapply(results[tied], `[[`, numeric(1), "n_edges"),
                     vapply(results[tied], `[[`, character(1), "key"))][1]
  for (i in seq_along(history)) history[[i]]$admissible <- admissible(i)
  history[[best]]$selected <- TRUE
  chosen <- results[[best]]
  if (length(history) >= max_evals) reason <- "evaluation_budget"
  if (verb) message("Selected ", chosen$n_edges, " internal links after ", length(history), " evaluations.")
  structure(list(model = deserialize(chosen$model),
    initial_model = deserialize(results[[initial]]$model),
    emulators = lapply(chosen$models, deserialize), initial_structure = struc,
    structure = chosen$struc, kept_inputs = chosen$keep,
    initial_score = results[[initial]]$loss, score = chosen$loss,
    initial_rmse = base, rmse = chosen$rmse, validation = chosen$prediction,
    history = do.call(rbind, history), stop_reason = reason,
    settings = list(seed = seed, B = B, M = M, workers = workers, threads = threads,
      search = search, visit_order = visit_order, groups = groups, weights = weights,
      scale = scale, protect = protect, max_evals = max_evals, max_rounds = max_rounds,
      min_improvement = min_improvement, tolerance = tolerance)), class = "prunlgp")
}

.prunlgp_key <- function(keep) paste(vapply(keep, function(k)
  paste(sprintf("%05d", k), collapse = ","), character(1)), collapse = "|")

.prunlgp_structure <- function(struc, keep) {
  use <- vapply(seq_len(nrow(struc)), function(i)
    struc$To_Input[i] %in% keep[[struc$To_Emulator[i]]], logical(1))
  out <- struc[use, , drop = FALSE]
  out$To_Input <- vapply(seq_len(nrow(out)), function(i)
    match(out$To_Input[i], keep[[out$To_Emulator[i]]]), integer(1))
  rownames(out) <- NULL
  out
}

.prunlgp_evaluate <- function(keep, state) {
  struc <- .prunlgp_structure(state$struc, keep)
  n_edges <- sum(struc$From_Emulator != "Global")
  tryCatch({
    if (any(lengths(keep) == 0L)) stop("A component would have zero inputs.")
    set_thread_num(state$threads)
    models <- lapply(state$models, deserialize)
    for (g in seq_along(state$members)) {
      members <- state$members[[g]]; id <- members[1]
      if (identical(keep[[id]], state$inputs[[id]])) next
      s <- as.integer((as.double(state$seed) + 104729 * g) %% 2147483646 + 1)
      set_seed(s)
      fitted <- state$refit(models[[id]], keep[[id]], s)
      if (!identical(class(fitted), class(models[[id]])) || as.character(fitted$id) != id ||
          !isTRUE(all.equal(unname(fitted$data$X), unname(models[[id]]$data$X[, keep[[id]], drop = FALSE]))) ||
          !isTRUE(all.equal(unname(fitted$data$Y), unname(models[[id]]$data$Y))))
        stop("refit must preserve the ID/response and model class, and use exactly the retained training columns.")
      for (copy in members) models[[copy]] <- set_id(fitted, copy)
    }
    set_seed(state$seed)
    linked <- lgp(struc, models, B = state$B, activate = TRUE, verb = FALSE)
    set_seed(state$seed)
    pred <- stats::predict(linked, state$x, method = "mean_var", full_layer = TRUE, M = state$M)$results
    means <- do.call(c, unname(pred$mean)); variances <- do.call(c, unname(pred$var))
    rmse <- vapply(names(state$y), function(id) {
      mu <- as.numeric(means[[id]]); va <- as.numeric(variances[[id]])
      if (length(mu) != length(state$y[[id]]) || length(va) != length(mu) ||
          any(!is.finite(c(mu, va))) || any(va < -1e-10)) stop("Invalid predictive moments.")
      sqrt(mean((mu - state$y[[id]])^2))
    }, numeric(1))
    loss <- state$objective(rmse, state$scale, state$weights)
    if (!is.numeric(loss) || length(loss) != 1L || !is.finite(loss)) stop("objective must return one finite loss.")
    list(ok = TRUE, loss = loss, rmse = rmse, n_edges = n_edges, struc = struc,
      prediction = list(mean = means, var = variances), model = serialize(linked, light = FALSE),
      models = lapply(models, serialize, light = FALSE))
  }, error = function(e) list(ok = FALSE, n_edges = n_edges, error = conditionMessage(e)))
}

#' @rdname prunlgp
#' @param ... Arguments forwarded to the selected model's method.
#' @export
predict.prunlgp <- function(object, ...) stats::predict(object$model, ...)

#' @rdname prunlgp
#' @export
summary.prunlgp <- function(object, ...) summary(object$model, ...)

#' @rdname prunlgp
#' @export
print.prunlgp <- function(x, ...) {
  cat("Pruned linked emulator\n",
      "  Validation loss: ", signif(x$initial_score, 6), " -> ", signif(x$score, 6), "\n",
      "  Internal links: ", sum(x$initial_structure$From_Emulator != "Global"), " -> ",
      sum(x$structure$From_Emulator != "Global"), "\n",
      "  Evaluations: ", nrow(x$history), " (", x$stop_reason, ")\n", sep = "")
  invisible(x)
}
