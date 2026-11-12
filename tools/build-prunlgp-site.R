# Build the fork website without changing the upstream pkgdown configuration.
args <- commandArgs(trailingOnly=FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
root <- if(length(file_arg)) dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1])))) else getwd()
repo <- Sys.getenv("GITHUB_REPOSITORY", "zhuangchaoqun/dgpsi-R")
owner <- strsplit(repo, "/", fixed=TRUE)[[1]][1]
override <- list(url=paste0("https://", owner, ".github.io/dgpsi-R"),
  development=list(mode="release"),
  repo=list(url=list(
    home=paste0("https://github.com/", repo, "/"),
    source=paste0("https://github.com/", repo, "/blob/prunlgp/"),
    issue=paste0("https://github.com/", repo, "/issues/"),
    user="https://github.com/")))
Sys.setenv(DGPSI_RUN_PRUNLGP_VIGNETTE="false")
pkgdown::build_site(pkg=root, examples=FALSE, preview=FALSE,
  new_process=FALSE, install=FALSE, override=override)
Sys.setenv(DGPSI_RUN_PRUNLGP_VIGNETTE="true", DGPSI_TUTORIAL_BACKEND="dev")
pkgdown::build_article("prunlgp", pkg=root, new_process=FALSE,
  override=override, quiet=FALSE)
file.create(file.path(root,"docs",".nojekyll"))
cat("Site built in ",file.path(root,"docs"),"\n",sep="")
