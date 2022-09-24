This is the new submission of 'dgpsi' 2.1.4

## Test environments

-   Local Test: macOS-arm64 (release)
-   GitHub Actions: macOS-latest (release)
-   GitHub Actions: windows-latest (release)
-   GitHub Actions: ubuntu-latest (devel)
-   GitHub Actions: ubuntu-latest (release)
-   GitHub Actions: ubuntu-latest (oldrel-1)
-   R-Hub: windows-x86_64-devel (r-devel)
-   R-Hub: ubuntu-gcc-release (r-release)
-   R-Hub: fedora-clang-devel (r-devel)
-   win-builder: windows-x86_64-w64-mingw32 (r-devel)

## R CMD check results

-   There were no ERRORs, WARNINGs and NOTEs on Local Test and GitHub Actions.

-   There were no ERRORs, WARNINGs and 5 NOTEs on R-Hub:

    -   On windows-x86_64-devel (r-devel)
        
        ```
        checking CRAN incoming feasibility ... [24s] NOTE
        Maintainer: 'Deyu Ming <deyu.ming.16@ucl.ac.uk>'

        New submission
        ```

        **Maintainer's Comment: New submission. This note can be safely ignored.**

    -   On windows-x86_64-devel (r-devel)
    
        ```
        checking for detritus in the temp directory ... NOTE
        Found the following files/directories:
          'lastMiKTeXException'
        ```

        **Maintainer's Comment: This note only appeared in this build. It did not appear in any other builds. In addition, As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX. Thus, this note can be likely ignored.**

    -   On ubuntu-gcc-release (r-release)
        
        ```
        checking CRAN incoming feasibility ... NOTE
        Maintainer: ‘Deyu Ming <deyu.ming.16@ucl.ac.uk>’

        New submission
        ```

        **Maintainer's Comment: New submission. This note can be safely ignored.**

    -   On fedora-clang-devel (r-devel)
        
        ```
        checking CRAN incoming feasibility ... [7s/22s] NOTE
        Maintainer: ‘Deyu Ming <deyu.ming.16@ucl.ac.uk>’

        New submission
        ```

        **Maintainer's Comment: New submission. This note can be safely ignored.**

    -   On fedora-clang-devel (r-devel)
             
        ```
        checking HTML version of manual ... NOTE
        Skipping checking HTML validation: no command 'tidy' found
        ```
            
        **Maintainer's Comment: Not sure about the cause of this note. Some relevant discussions are given [here](https://groups.google.com/g/r-sig-mac/c/7u_ivEj4zhM) and it seems that the test Server does not have `tidy` installed. This note only appeared in this build, and did not appear in any other builds. Thus, this note can be likely ignored.**

-   There were no ERRORs, WARNINGs and 1 NOTE on win-builder:
        
    ```
    checking CRAN incoming feasibility ... [10s] NOTE
    Maintainer: 'Deyu Ming <deyu.ming.16@ucl.ac.uk>'

    New submission

    Possibly misspelled words in DESCRIPTION:
      Guillas (15:56, 16:34)
    ```

    **Maintainer's Comment: The note on New submission can be safely ignored. The possible misspelled word "Guillas" is the surname of a co-author of articles cited in DESCRIPTION. This note did not appear in any other R CMD checks because `.aspell/` directory has been created under the package root to mitigate this note as per instructions given [here](http://dirk.eddelbuettel.com/blog/2017/08/10/#008_aspell_cran_incoming) and examples given in [`{edmdata}`](https://github.com/tmsalab/edmdata). Thus, this note can be safely ignored.**
