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
#' @param py_ver a string that gives the 'python' version to be installed. Supported versions are 3.10, 3.11, and 3.12.
#'    If `py_ver = NULL`, the default 'python' version '3.10' will be installed.
#' @param dgpsi_ver a string that gives the 'python' version of 'dgpsi' to be used. If `dgpsi_ver = NULL`,
#' * the latest 'python' version of 'dgpsi' will be used, if the package is installed from CRAN;
#' * the development 'python' version of 'dgpsi' will be used, if the package is installed from GitHub.
#' @param conda optional path to a conda binary. If `NULL`, [init_py()]
#'   attempts to locate conda automatically; see [reticulate::conda_binary()] for the search order.
#'   If provided, the path is used directly and must point to an existing conda executable.
#' @param reinstall a bool that indicates whether to reinstall the 'python' version of 'dgpsi' specified
#'    in `dgpsi_ver` if it has already been installed. This argument is useful when the development version
#'    of the R package is installed and one may want to regularly update the development 'python' version
#'    of 'dgpsi'. Defaults to `FALSE`.
#' @param uninstall a bool that indicates whether to uninstall the 'python' version of 'dgpsi' specified
#'    in `dgpsi_ver` if it has already been installed. This argument is useful when the 'python' environment
#'    is corrupted and one wants to completely uninstall and reinstall it. Defaults to `FALSE`.
#' @param verb a bool indicating if trace information will be printed during function execution.
#'     Defaults to `TRUE`.
#' @param show_config a bool indicating whether to print the Python backend
#'   configuration when `verb = TRUE`. Defaults to `FALSE`.
#'
#' @return No return value, called to install required 'python' environment.
#' @note
#' On Linux, `init_py()` may ask for permission during installation. To
#' auto-accept these prompts in non-interactive settings,
#' set the environment variable `DGPSI_AUTO_YES` to
#' `"TRUE"` before calling [init_py()] for the first time. Accepted values
#' include `"TRUE"`, `"true"`, `"yes"`, `"y"`, and
#' `"1"`. If unset, [init_py()] will continue to prompt interactively.
#' @details See further examples and tutorials at <`r get_docs_url()`>.
#' @examples
#' \dontrun{
#'
#' # See gp(), dgp(), or lgp() for an example.
#' }
#'
#' @md
#' @export
init_py <- function(py_ver = NULL, dgpsi_ver = 'dev', conda = NULL, reinstall = FALSE, uninstall = FALSE, verb = TRUE, show_config = FALSE) {
  if ( is.null(py_ver) ) py_ver <- '3.10'
  if ( is.null(dgpsi_ver) ) {
    ##For devel version
    dgpsi_ver <- c('dill>=0.3.2', 'matplotlib-base>=3.2.1', 'numba >=0.51.2',
                   'numpy >=1.18.2', 'pathos >=0.2.9', 'multiprocess >=0.70.13', 'psutil >=5.8.0',
                   'scikit-learn >=0.22.0', 'scipy >=1.4.1', 'tqdm >=4.50.2', 'tabulate >=0.8.7', 'faiss-cpu >=1.7.4', 'tbb', 'pip')
    env_name <- 'dgp_si_R_2_6_0_9000'

    ##For release version
    #dgpsi_ver <- 'dgpsi==2.6.0'
    #env_name <- 'dgp_si_R_2_6_0'
  } else {
    ##For devel version
    if (dgpsi_ver=='dev'){
      dgpsi_ver <- c('dill>=0.3.2', 'matplotlib-base>=3.2.1', 'numba >=0.51.2',
                     'numpy >=1.18.2', 'pathos >=0.2.9', 'multiprocess >=0.70.13', 'psutil >=5.8.0',
                     'scikit-learn >=0.22.0', 'scipy >=1.4.1', 'tqdm >=4.50.2', 'tabulate >=0.8.7', 'faiss-cpu >=1.7.4', 'tbb', 'pip')
      env_name <- 'dgp_si_R_2_6_0_dev'
    } else {
      dgpsi_ver <- paste('dgpsi==', dgpsi_ver, sep = "")
      env_name <- paste('dgp_si_R_', gsub(".", "_", dgpsi_ver,fixed=TRUE), sep = "")
    }

    ##For release version
    # dgpsi_ver <- paste('dgpsi==', dgpsi_ver, sep = "")
    # env_name <- paste('dgp_si_R_', gsub(".", "_", dgpsi_ver,fixed=TRUE), sep = "")
  }
  #Check if there is any conda binary installed, if not, request to install it.
  #restart <- FALSE
  Sys.setenv(CONDA_PLUGINS_AUTO_ACCEPT_TOS = "yes")
  auto_yes <- identical(Sys.getenv("GITHUB_ACTIONS"), "true") ||
    tolower(Sys.getenv("DGPSI_AUTO_YES", "")) %in% c("1", "true", "yes", "y")
  conda_missing <- if (is.null(conda)) {
    is.null(tryCatch(reticulate::conda_binary(), error = function(e) NULL))
  } else {
    tryCatch(
      { reticulate::conda_version(conda = conda); FALSE },
      error = function(e) TRUE
    )
  }
  if (conda_missing){
    if (is.null(conda)) {
      ans <- if (auto_yes) {
        "y"
      } else {
        readline("I am unable to find a conda binary. Do you want me to install it for you? (Y/N) ")
      }
    } else {
      stop("I am unable to find the specified conda binary. Please check the path supplied to `conda`.", call. = FALSE)
    }

    #If the user would like to have the conda binary to be installed
    if ( tolower(trimws(ans))=='y'|tolower(trimws(ans))=='yes' ){
      message("Installing the Conda binary...")
      reticulate::install_miniconda()
      conda_path <- reticulate::conda_binary()
      install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver, auto_yes)
      pkg.env$restart <- TRUE
    } else{
      #stop("Please first install Miniforge, Miniconda, or Anaconda, and then re-initialize the Python environment.", call. = FALSE)
      stop("Please first install Miniforge, Miniconda, or Anaconda, and then re-load the package.", call. = FALSE)
    }
  } else {
    conda_path <- if (is.null(conda)) {
      reticulate::conda_binary()
    } else {
      conda
    }
    no_dgpsi <- inherits(tryCatch(reticulate::conda_python(envname = env_name, conda = conda_path), error = identity), "error")
    if (no_dgpsi){
      install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver, auto_yes)
      pkg.env$restart <- TRUE
    } else {
      if ( uninstall ){
        reticulate::conda_remove(envname = env_name, conda = conda_path)

        if (Sys.info()[["sysname"]] == 'Linux') {
          shell <- Sys.getenv("SHELL")
          if (grepl("bash", shell)) {
            rc_file <- "~/.bashrc"
          } else if (grepl("zsh", shell)) {
            rc_file <- "~/.zshrc"
          } else {
            rc_file <- "~/.bashrc"
          }

          rc_file_path <- path.expand(rc_file)

          if (file.exists(rc_file_path)) {
            rc_content <- readLines(rc_file_path)

            matching_lines <- grepl(
              paste0("^[[:space:]]*export[[:space:]]+(R_)?LD_LIBRARY_PATH[[:space:]]*=.*", env_name),
              rc_content
            )

            rc_content_cleaned <- rc_content[!matching_lines]
            writeLines(rc_content_cleaned, rc_file_path)
          }
        }

        pkg.env$restart <- TRUE
        #message("Uninstallation finished. Please restart R and run 'init_py()' to reinstall the Python environment.")
        message("Uninstallation finished. Please restart R.")
      } else {
        if (isTRUE(reinstall)) {
          install_dgpsi(env_name, py_ver, conda_path, dgpsi_ver, auto_yes, reinsatll = TRUE)
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
    if ( verb ) {
      message("The Python environment for 'dgpsi' is successfully loaded.")
      if (show_config) {
        cfg <- reticulate::py_config()
        message(
          "\n[dgpsi Python configuration]\n",
          "  Conda binary      : ", conda_path, "\n",
          "  Conda environment : ", env_name, "\n",
          "  Python path       : ", cfg$python, "\n",
          "  Python version    : ", cfg$version_string
        )
      }
    }
  }
}

install_dgpsi <- function(env_name, py_ver, conda_path, dgpsi_ver, auto_yes, reinsatll = FALSE) {
  if (!reinsatll) message(sprintf("Setting up the Python environment for %s ...\n", dgpsi_ver))
  if (reinsatll) {
    message("Re-installing the required Python packages ...")
    if (!grepl("9000|dev", env_name)){
      reticulate::conda_remove(envname = env_name, conda = conda_path)

      if (Sys.info()[["sysname"]] == 'Linux') {
        shell <- Sys.getenv("SHELL")
        if (grepl("bash", shell)) {
          rc_file <- "~/.bashrc"
        } else if (grepl("zsh", shell)) {
          rc_file <- "~/.zshrc"
        } else {
          rc_file <- "~/.bashrc"
        }

        rc_file_path <- path.expand(rc_file)

        if (file.exists(rc_file_path)) {
          rc_content <- readLines(rc_file_path)

          matching_lines <- grepl(
            paste0("^[[:space:]]*export[[:space:]]+(R_)?LD_LIBRARY_PATH[[:space:]]*=.*", env_name),
            rc_content
          )

          rc_content_cleaned <- rc_content[!matching_lines]
          writeLines(rc_content_cleaned, rc_file_path)
        }
      }
      reinsatll <- FALSE
    }
  } else {
    message("Installing the required Python packages ...")
  }
  if (Sys.info()[["sysname"]] == "Darwin" & Sys.info()[["machine"]] == "arm64"){
    macos_version <- system("sw_vers -productVersion", intern = TRUE)
    version_nums <- as.numeric(strsplit(macos_version, "\\.")[[1]])
    current_version <- version_nums[1] + version_nums[2] / 10
    if (current_version>=13.3){
      if (!reinsatll) reticulate::conda_create(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*newaccelerate"'),
                                               python_version = py_ver, conda = conda_path, forge = TRUE, additional_create_args = c('--strict-channel-priority'))
    } else {
      if (!reinsatll) reticulate::conda_create(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*accelerate"'),
                                               python_version = py_ver, conda = conda_path, forge = TRUE, additional_create_args = c('--strict-channel-priority'))
    }
  } else if ( isTRUE(grepl("Intel",benchmarkme::get_cpu()$model_name)) ){
    if (Sys.info()[["sysname"]] %in% c("Windows", "Linux")) dgpsi_ver <- c(dgpsi_ver, 'intel-cmplr-lib-rt')
    if (!reinsatll) reticulate::conda_create(envname = env_name, packages = c(dgpsi_ver, '"libblas=*=*mkl"', 'mkl>=2022'),
                                             python_version = py_ver, conda = conda_path, forge = TRUE, additional_create_args = c('--strict-channel-priority'))
  } else {
    if (!reinsatll) reticulate::conda_create(envname = env_name, packages = c(dgpsi_ver),
                                             python_version = py_ver, conda = conda_path, forge = TRUE, additional_create_args = c('--strict-channel-priority'))
  }
  if (grepl("9000|dev", env_name)) {
    git_path = "git+https://github.com/mingdeyu/DGP.git"
    if (grepl("dev", env_name)){
      git_path = paste0(git_path, '@dev')
    }
    if (reinsatll) {
      reticulate::conda_install(
        envname  = env_name,
        packages = "tbb",
        conda    = conda_path,
        forge = TRUE
      ) # only this dev version - remove after the next release
      reticulate::conda_install(envname = env_name, packages = c(git_path) , conda = conda_path, pip = TRUE, pip_options = c('--no-deps', '--force-reinstall'))
    } else {
      reticulate::conda_install(envname = env_name, packages = c(git_path) , conda = conda_path, pip = TRUE, pip_options = c('--no-deps'))
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
    libstdc_path <- normalizePath(
      paste0(gsub("bin.*$", "", conda_dgpsi_path), "lib"),
      winslash = "/",
      mustWork = FALSE
    )

    # Construct the export command
    export_command <- sprintf(
      paste0(
        "export R_LD_LIBRARY_PATH=\"",
        "%s",
        "$(R RHOME 2>/dev/null | sed 's#^#:#; s#$#/lib#')",
        "${R_LD_LIBRARY_PATH:+:${R_LD_LIBRARY_PATH}}",
        "\""
      ),
      libstdc_path
    )

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
    message("To use the package properly, we need to update your R_LD_LIBRARY_PATH.")
    permission <- if (auto_yes) { "y" } else { readline(prompt = "Can we automatically update this in your shell configuration file? (Y/N) ") }

    if (tolower(trimws(permission)) == 'y' || tolower(trimws(permission)) == 'yes') {
      rc_file_path <- path.expand(rc_file)

      # Check if the rc file exists, if not create it
      if (!file.exists(rc_file_path)) {
        #message(paste0(rc_file, " does not exist. Creating a new ", rc_file, " file."))
        is.create <- file.create(rc_file_path)
      } else {
        # Read the content of the rc file
        rc_content <- readLines(rc_file_path)

        # Identify lines that match the existing R_LD_LIBRARY_PATH export command
        matching_lines <- grepl(
          "^[[:space:]]*export[[:space:]]+(R_)?LD_LIBRARY_PATH[[:space:]]*=.*dgp_si_R",
          rc_content
        )

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
