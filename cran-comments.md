# Resubmission
This is a resubmission. Following the comments (see below) from CRAN in the first submission, various changes (see responses below) have been made in this version (`2.1.5`).

## Responses to comments from CRAN

- Comment 1:

  **CRAN**: Please omit the redundant "R" at the start of your title and the description.

  **Maintainer's Response**: Done.

- Comment 2:

  **CRAN**: Please always write package names, software names and API (application programming interface) names in single quotes in title and description. e.g: --> 'python' Please note that package names are case sensitive.
  
  **Maintainer's Response**: All package names, software names, and API mentioned in title and description are now in single quotes. 

- Comment 3:

  **CRAN**: Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar). 
  
  Missing Rd-tags:
  * init_py.Rd: \arguments,  \value
  * trace_plot.Rd: \value
  * write.Rd: \value
  
  **Maintainer's Response**: All exported methods now have \arguments and \value in their .Rd files. The structure of the output (class) from exported methods are now all carefully explained.

- Comment 4:

  **CRAN**: Please add small executable examples in your Rd-files to illustrate the use of the exported function but also enable automatic testing.
  
  **Maintainer's Response**: Small examples are added to key functions `gp()`, `dgp()`, and `lgp()` of the package. Since the examples involve the implementations of other exported functions, we did not repeat the same examples in different Rd-files. Instead, we highlight in `Examples` sections of other exported functions that their usages can be found in `Examples` sections of `gp()`, `dgp()`, and `lgp()`. Since the package depends on the underlying 'python' implementation, the added small examples are wrapped in `\dontrun{}`, following other similar packages, such as [`{reticulate}`](https://github.com/rstudio/reticulate), [`{keras}`](https://github.com/rstudio/keras), etc. All small examples were carefully tested before being included in Rd-files. 

## Test environments

- Local Test: macOS-arm64 (release)
- GitHub Actions: macOS-latest (release)
- GitHub Actions: windows-latest (release)
- GitHub Actions: ubuntu-latest (devel)
- GitHub Actions: ubuntu-latest (release)
- GitHub Actions: ubuntu-latest (oldrel-1)
- R-Hub: windows-x86_64-devel (r-devel)
- R-Hub: ubuntu-gcc-release (r-release)
- R-Hub: fedora-clang-devel (r-devel)
- win-builder: windows-x86_64-w64-mingw32 (r-devel)

## R CMD check results

- There were no ERRORs, WARNINGs and NOTEs on Local Test and GitHub Actions.

- There were no ERRORs, WARNINGs and 5 NOTEs on R-Hub:

  * On windows-x86_64-devel (r-devel)
        
    ```
    checking CRAN incoming feasibility ... [21s] NOTE
    Maintainer: 'Deyu Ming <deyu.ming.16@ucl.ac.uk>'

    New submission
    ```

    **Maintainer's Comment**: New submission. This note can be safely ignored.

  * On windows-x86_64-devel (r-devel)
    
    ```
    checking for detritus in the temp directory ... NOTE
    Found the following files/directories:
    'lastMiKTeXException'
    ```

    **Maintainer's Comment**: This note only appeared in this build. It did not appear in any other builds. In addition, As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX. Thus, this note can be likely ignored.

  * On ubuntu-gcc-release (r-release)
        
    ```
    checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Deyu Ming <deyu.ming.16@ucl.ac.uk>’

    New submission
    ```

    **Maintainer's Comment**: New submission. This note can be safely ignored.

  * On fedora-clang-devel (r-devel)
        
    ```
    checking CRAN incoming feasibility ... [7s/23s] NOTE
    Maintainer: ‘Deyu Ming <deyu.ming.16@ucl.ac.uk>’

    New submission
    ```

    **Maintainer's Comment**: New submission. This note can be safely ignored.

  * On fedora-clang-devel (r-devel)
             
    ```
    checking HTML version of manual ... NOTE
    Skipping checking HTML validation: no command 'tidy' found
    ```
            
    **Maintainer's Comment**: Not sure about the cause of this note. Some relevant discussions are give                   [here](https://groups.google.com/g/r-sig-mac/c/7u_ivEj4zhM) and it seems that the test Server does not have `tidy` installed. This note only appeared in this build, and did not appear in any other builds. Thus, this note can be likely ignored.

- There were no ERRORs, WARNINGs and 1 NOTE on win-builder:
        
  ```
  checking CRAN incoming feasibility ... [10s] NOTE
  Maintainer: 'Deyu Ming <deyu.ming.16@ucl.ac.uk>'

  New submission

  Possibly misspelled words in DESCRIPTION:
    Guillas (15:56, 16:34)
  ```

  **Maintainer's Comment**: The note on New submission can be safely ignored. The possible misspelled word "Guillas" is the surname of a co-author of articles cited in DESCRIPTION. This note did not appear in any other R CMD checks because `.aspell/` directory has been created under the package root to mitigate this note as per instructions given [here](http://dirk.eddelbuettel.com/blog/2017/08/10/#008_aspell_cran_incoming) and examples given in [`{edmdata}`](https://github.com/tmsalab/edmdata). Thus, this note can be safely ignored.
