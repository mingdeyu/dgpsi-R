# Run from a clone of this repository with the fork installed.
if (!requireNamespace("dgpsi", quietly=TRUE) ||
    !"prunlgp" %in% getNamespaceExports("dgpsi")) {
  stop("Install this fork first; see PRUNLGP.md.")
}
if (!requireNamespace("rmarkdown", quietly=TRUE)) {
  stop("Install rmarkdown before running the tutorial.")
}
args <- commandArgs(trailingOnly=FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
if (length(file_arg)) {
  root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else {
  source_file <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  root <- if (length(source_file)) dirname(dirname(normalizePath(tail(source_file, 1)[[1]]))) else getwd()
}
input <- file.path(root, "vignettes", "prunlgp.Rmd")
if (!file.exists(input)) stop("Run this script from the repository containing vignettes/prunlgp.Rmd.")
output <- file.path(root, "tutorial-output")
dir.create(output, showWarnings=FALSE)
Sys.setenv(DGPSI_RUN_PRUNLGP_VIGNETTE="true")
if (!nzchar(Sys.getenv("DGPSI_TUTORIAL_BACKEND"))) {
  Sys.setenv(DGPSI_TUTORIAL_BACKEND="dev")
}
run <- new.env(parent=globalenv())
elapsed <- system.time(rmarkdown::render(input, output_dir=output,
  knit_root_dir=root, envir=run, quiet=FALSE))[["elapsed"]]
write.csv(run$scores, file.path(output, "errors.csv"), row.names=FALSE)
write.csv(run$connection_summary, file.path(output, "connections.csv"), row.names=FALSE)
write.csv(run$selected$history, file.path(output, "search.csv"), row.names=FALSE)
write.csv(run$selected$structure, file.path(output, "structure.csv"), row.names=FALSE)
write.csv(data.frame(partition=rep(c("validation", "test"), each=2),
  model=rep(c("Unpruned", "Pruned"), 2),
  NRMSE=c(run$selected$initial_score, run$selected$score, unname(run$system_test))),
  file.path(output, "system-errors.csv"), row.names=FALSE)
write.csv(rbind(run$interval_summary(run$before, "Unpruned"),
                run$interval_summary(run$after, "Pruned")),
          file.path(output, "intervals.csv"), row.names=FALSE)
stopifnot(nrow(run$scores)==3L, all(is.finite(as.matrix(run$scores[,-1]))),
  sum(run$connection_summary$initial)==6L,
  sum(run$connection_summary$retained)==sum(run$selected$structure$From_Emulator!="Global"),
  run$selected$score <= run$selected$initial_score + run$selected$settings$tolerance)
capture.output(sessionInfo(), file=file.path(output, "session-info.txt"))
saveRDS(list(scores=run$scores, summary=run$connection_summary,
  history=run$selected$history, initial_score=run$selected$initial_score,
  score=run$selected$score, settings=run$selected$settings,
  initial_test=run$before, pruned_test=run$after, elapsed_seconds=elapsed),
  file.path(output,"results.rds"))
cat(sprintf("\nTutorial completed in %.1f seconds.\n%s\n", elapsed,
            file.path(output,"prunlgp.html")))
