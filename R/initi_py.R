pkg.env <- new.env(parent = emptyenv())
pkg.env$dgpsi <- NULL
pkg.env$py_buildin <- NULL
pkg.env$np <- NULL
pkg.env$copy <- NULL
pkg.env$py_gc <- NULL
pkg.env$restart <- FALSE
pkg.env$thread_num <- NULL
pkg.env$base64 <- NULL
pkg.env$dill <- NULL

#' @title 'python' environment initialization
#'
#' @description This function initializes the 'python' environment for the package.
#'
#' @param py_ver a string that gives the 'python' version to be installed. If `py_ver = NULL`, the default 'python'
#'     version '3.9.13' will be installed.
#' @param dgpsi_ver a string that gives the 'python' version of 'dgpsi' to be used. If `dgpsi_ver = NULL`,
#' * the latest 'python' version of 'dgpsi' will be used, if the package is installed from CRAN;
#' * the development 'python' version of 'dgpsi' will be used, if the package is installed from GitHub.
#' @param reinstall a bool that indicates whether to reinstall the 'python' version of 'dgpsi' specified
#'    in `dgpsi_ver` if it has already been installed. This argument is useful when the development version
#'    of the R package is installed and one may want to regularly update the development 'python' version
#'    of 'dgpsi'. Defaults to `FALSE`.
#' @param uninstall a bool that indicates whether to uninstall the 'python' version of 'dgpsi' specified
#'    in `dgpsi_ver` if it has already been installed. This argument is useful when the 'python' environment
#'    is corrupted and one wants to completely uninstall and reinstall it. Defaults to `FALSE`.
#' @param verb a bool indicating if trace information will be printed during function execution.
#'     Defaults to `TRUE`.
#'
#' @return No return value, called to install required 'python' environment.
#'
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#'
#' @md
#' @export
init_py <- function(py_ver = NULL, dgpsi_ver = NULL, reinstall = FALSE, uninstall = FALSE, verb = TRUE) {
  if ( is.null(py_ver) ) py_ver <- '3.9.13'
  if ( is.null(dgpsi_ver) ) {
    ##For devel version
    #dgpsi_ver <- c('cython>=0.29.30', 'dill>=0.3.2, <=0.3.5.1', 'jupyter>=1.0.0', 'matplotlib-base>=3.2.1', 'numba >=0.51.2',
    #               'numpy >=1.18.2', 'pathos ==0.2.9', 'multiprocess ==0.70.13', 'psutil >=5.8.0', 'pybind11 >=2.10.0', 'pythran >=0.11.0',
    #               'scikit-build >=0.15.0', 'scikit-learn >=0.22.0', 'scipy >=1.4.1', 'tqdm >=4.50.2', 'tabulate >=0.8.7', 'faiss-cpu >=1.7.4')
    #env_name <- 'dgp_si_R_2_4_0_9000'
    ##For release version
    dgpsi_ver <- 'dgpsi==2.5.0'
    env_name <- 'dgp_si_R_2_5_0'
  } else {
    env_name <- paste('dgp_si_R_', gsub(".", "_", dgpsi_ver,fixed=TRUE), sep = "")
    dgpsi_ver <- paste('dgpsi==', dgpsi_ver, sep = "")
  }
  #Check if there is any conda binary installed, if not, request to install it.
  #restart <- FALSE
  if (is.null(tryCatch(reticulate::conda_binary(), error = function(e) NULL))){
    ans <- readline(prompt="I am unable to find a conda binary. Do you want me to install it for you? (Y/N) ")
    #If the user would like to have the conda binary to be installed
    if ( tolower(ans)=='y'|tolower(ans)=='yes' ){
      message("Installing the Conda binary...")
      reticulate::install_miniconda()
      conda_path <- reticulate::conda_binary()
      install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
      pkg.env$restart <- TRUE
    } else{
      #stop("Please first install Miniforge, Miniconda, or Anaconda, and then re-initialize the Python environment.", call. = FALSE)
      stop("Please first install Miniforge, Miniconda, or Anaconda, and then re-load the package.", call. = FALSE)
    }
  } else {
    conda_path <- reticulate::conda_binary()
    no_dgpsi <- inherits(tryCatch(reticulate::conda_python(envname = env_name, conda = conda_path), error = identity), "error")
    if (no_dgpsi){
      if (any(grepl('^dgp_si_R', reticulate::conda_list(conda = conda_path)$name))){
        #conda_list <- reticulate::conda_list(conda = conda_path)$name
        #dgpsi_list <- conda_list[grepl('^dgp_si_R', conda_list)]
        #cat("I found Python environment(s) for other versions of the package.")
        #ans <- readline(prompt="Do you want me to remove them? (Y/N) ")
        #if ( tolower(ans)=='y'|tolower(ans)=='yes' ){
        #  message(sprintf("Removing Python environment(s): %s.\n", paste(dgpsi_list, collapse = ', ')))
        #  for (item in dgpsi_list) {
        #    reticulate::conda_remove(envname = item, conda = conda_path)
        #  }
        #  message("Done.")
        #}
        install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
        pkg.env$restart <- TRUE
      } else {
        ans <- readline(prompt="Is this your first time using the package? (Y/N) ")
        if ( tolower(ans)=='n'|tolower(ans)=='no' ){
              message("I am unable to find the required Python environment. It may be because your conda binary has changed.")
              cat("I am re-setting it for you now ...")
              install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
              pkg.env$restart <- TRUE
              } else {
                install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver)
                pkg.env$restart <- TRUE
              }
      }
    } else {
      if ( uninstall ){
        reticulate::conda_remove(envname = env_name, conda = conda_path)
        pkg.env$restart <- TRUE
        #message("Uninstallation finished. Please restart R and run 'init_py()' to reinstall the Python environment.")
        message("Uninstallation finished. Please restart R.")
      } else {
        if (isTRUE(reinstall)) {
          install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver, reinsatll = TRUE)
          #if (grepl('9000',env_name)) {
          #  reticulate::conda_install(envname = env_name, packages = c("git+https://github.com/mingdeyu/DGP.git") , conda = conda_path, pip = TRUE, pip_options = c('--no-deps', '--force-reinstall'))
          #} else {
          #  reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver) , conda = conda_path)
          #}
          #if (Sys.info()[["sysname"]] == 'Linux'){
          #  libstdc_path <- paste(gsub("bin.*$", "", conda_path), 'envs/', env_name, '/lib/libstdc++.so.6.0.30', sep='')
          #  libstdc_sys_path <- "/usr/lib/x86_64-linux-gnu/libstdc++.so.6"
          #  system(paste("sudo rm",libstdc_sys_path))
          #  system(paste("sudo ln -s", libstdc_path, libstdc_sys_path))
          #}
          #message("Installation finished. Please restart R.")
          pkg.env$restart <- TRUE
        }
      }
    }
  }

  if ( isFALSE(pkg.env$restart) ){
    if ( verb ) message("Connecting to Python ...", appendLF = FALSE)
    warning_error_handler(with_warning_handler(reticulate::use_condaenv(condaenv = env_name, conda = conda_path, required = TRUE)))
    if ( verb ) message(" done")

    if ( verb ) message("Importing required Python modules ...", appendLF = FALSE)
    assign('dgpsi', reticulate::import("dgpsi"), pkg.env)
    assign('py_buildin', reticulate::import_builtins(), pkg.env)
    assign('np', reticulate::import("numpy"), pkg.env)
    assign('copy', reticulate::import("copy"), pkg.env)
    assign('py_gc', reticulate::import("gc"), pkg.env)
    assign('dill', reticulate::import("dill"), pkg.env)
    assign('base64', reticulate::import("base64"), pkg.env)
    pkg.env$thread_num <- pkg.env$dgpsi$get_thread()
    if ( verb ) message(" done")
    Sys.sleep(0.5)
    if ( verb ) message("The Python environment for 'dgpsi' is successfully loaded.")
  }
}

install_dgpsi <- function(env_name, py_ver, conda_path, dgpsi_ver, reinsatll = FALSE) {
  if (!reinsatll) message(sprintf("Setting up the Python environment for %s ...\n", dgpsi_ver))
  if (!reinsatll) reticulate::conda_create(envname = env_name, python_version = py_ver, conda = conda_path)
  if (reinsatll) {
    message("Re-installing the required Python packages ...")
  } else {
    message("Installing the required Python packages ...")
  }
  if (Sys.info()[["sysname"]] == "Darwin" & Sys.info()[["machine"]] == "arm64"){
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*accelerate"') , conda = conda_path)
  } else if (grepl("Intel",benchmarkme::get_cpu()$model_name)){
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*mkl"') , conda = conda_path)
    reticulate::conda_install(envname = env_name, packages = c('icc_rt') , channel = c("numba"), conda = conda_path)
  } else {
    reticulate::conda_install(envname = env_name, packages = c(dgpsi_ver) , conda = conda_path)
  }
  if (grepl('9000',env_name)) {
    if (reinsatll) {
      reticulate::conda_install(envname = env_name, packages = c("git+https://github.com/mingdeyu/DGP.git") , conda = conda_path, pip = TRUE, pip_options = c('--no-deps', '--force-reinstall'))
    } else {
      reticulate::conda_install(envname = env_name, packages = c("git+https://github.com/mingdeyu/DGP.git") , conda = conda_path, pip = TRUE, pip_options = c('--no-deps'))
    }
  }
  #if (Sys.info()[["sysname"]] == 'Linux' & !any(grepl("libstdc++.so.6.0.3",list.files("/usr/lib/x86_64-linux-gnu/"), fixed = TRUE))){
  #  cat("The required file 'libstdc++.so.6.0.30' or above is missing from /usr/lib/x86_64-linux-gnu/.")
  #  ans <- readline(prompt="To proceed, would you like to grant sudo permissions to resolve the issue? (Y/N) ")
  #  if ( tolower(ans)=='y'|tolower(ans)=='yes' ){
  #    libstdc_path <- paste(gsub("bin.*$", "", conda_path), 'envs/', env_name, '/lib/libstdc++.so.6.0.3*', sep='')
  #    system(paste("sudo cp", libstdc_path, "/usr/lib/x86_64-linux-gnu/"))
  #    libstdc_sys_path <- "/usr/lib/x86_64-linux-gnu/libstdc++.so.6"
  #    system(paste("sudo rm",libstdc_sys_path))
  #    system(paste("sudo ln -s", "/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.3*", libstdc_sys_path))
  #  } else {
  #    stop("Please link /usr/lib/x86_64-linux-gnu/libstdc++.so.6 to 'libstdc++.so.6.0.30' or above and run init_py(reinstall = T).", call. = FALSE)
  #  }
  #}

  if (Sys.info()[["sysname"]] == 'Linux') {
    # Retrieve conda environment information
    conda_env_path <- reticulate::conda_list(conda = conda_path)
    conda_dgpsi_path <- conda_env_path$python[conda_env_path$name == env_name]
    libstdc_path <- paste(gsub("bin.*$", "", conda_dgpsi_path), 'lib', sep = '')

    # Construct the export command
    export_command <- paste0("export LD_LIBRARY_PATH=", libstdc_path, ":$LD_LIBRARY_PATH")

    # Detect the user's shell
    shell <- Sys.getenv("SHELL")
    if (grepl("bash", shell)) {
      rc_file <- "~/.bashrc"
    } else if (grepl("zsh", shell)) {
      rc_file <- "~/.zshrc"
    } else {
      # Default to ~/.bashrc if the shell is unrecognized
      rc_file <- "~/.bashrc"
    }

    # Inform the user
    message("To use the package properly, we need to update your LD_LIBRARY_PATH.")
    permission <- readline(prompt = "Can we automatically update this in your shell configuration file? (Y/N) ")

    if (tolower(permission) == 'y' || tolower(permission) == 'yes') {
      rc_file_path <- path.expand(rc_file)

      # Check if the rc file exists, if not create it
      if (!file.exists(rc_file_path)) {
        #message(paste0(rc_file, " does not exist. Creating a new ", rc_file, " file."))
        is.create <- file.create(rc_file_path)
      } else {
        # Read the content of the rc file
        rc_content <- readLines(rc_file_path)

        # Identify lines that match the existing LD_LIBRARY_PATH export command
        matching_lines <- grepl("LD_LIBRARY_PATH=.*dgp_si_R", rc_content)

        # Remove any existing export commands with the pattern "dgp_si_R"
        rc_content_cleaned <- rc_content[!matching_lines]

        # Write back the cleaned content (without the old export commands)
        writeLines(rc_content_cleaned, rc_file_path)
      }

      # Append the new export command
      cat(export_command, file = rc_file_path, append = TRUE, sep = "\n")

      message(paste0("The path has been updated in your ", rc_file, "."))
      message(paste0("You may need to restart your terminal or run 'source ", rc_file, "' to apply the updates."))
    } else {
      # Print the command for manual addition
      message("Please manually add the following line to your ", rc_file, ":")
      cat(export_command, "\n")
    }
  }

  message("Installation finished. Please restart R.")
}

with_warning_handler <- function(...)
{
  Sys.unsetenv("RETICULATE_PYTHON")
  withCallingHandlers(..., warning = function(w)
  { condition <- conditionMessage(w)
  reg1 <- "Previous request to"
  reg2 <- "will be ignored. It is superseded by request to"
  if(grepl(reg1, condition) & grepl(reg2, condition)) invokeRestart("muffleWarning")
  })
}

warning_error_handler <- function(...){
  tryCatch(...,
           error = function(r)
           { condition <- conditionMessage(r)
           reg <- "failed to initialize requested version of Python"
           if(grepl(reg, condition)) {
             cat("NOTE: please clear your R workspace and delete the workspace image file '.RData' before restarting the R session.")
           } else {
             message(paste("ERROR in", condition))
           }
           })
}
