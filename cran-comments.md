This is the submission of 'dgpsi' 2.4.0.

Various updates have been made in this version. A detailed list of changes that have been made since the last version 2.3.0 is contained in NEWS.md.  

## Test environments

- Local Test: macOS-arm64 (release)
- GitHub Actions: macOS-latest (release)
- GitHub Actions: windows-latest (release)
- GitHub Actions: ubuntu-latest (devel)
- GitHub Actions: ubuntu-latest (release)
- GitHub Actions: ubuntu-latest (oldrel-1)
- win-builder: windows-x86_64-w64-mingw32 (r-devel)
- win-builder: windows-x86_64-w64-mingw32 (r-release)
- R-Hub: windows-x86_64-devel (r-devel)
- R-Hub: fedora-clang-devel (r-devel)

## R CMD check results

- There were no ERRORs, WARNINGs and NOTEs on Local Test, GitHub Actions, and win-builder.

- There were no ERRORs, WARNINGs and 3 NOTEs on R-Hub:

  * On windows-x86_64-devel (r-devel)
    ```
    checking for non-standard things in the check directory ... NOTE
    Found the following files/directories:
    ''NULL''
    ```
    
    **Maintainer's Comment**: This note only appeared in this build. It did not appear in any other builds. In addition, As noted in [R-hub issue #560](https://github.com/r-hub/rhub/issues/560), this could be a R-hub issue. Thus, this note can be ignored.
    
    ```
    checking for detritus in the temp directory ... NOTE
    Found the following files/directories:
    'lastMiKTeXException'
    ```

    **Maintainer's Comment**: This note only appeared in this build. It did not appear in any other builds. In addition, As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX. Thus, this note can be likely ignored.

  * On fedora-clang-devel (r-devel)
             
    ```
    checking HTML version of manual ... NOTE
    Skipping checking HTML validation: no command 'tidy' found
    ```
            
    **Maintainer's Comment**: Not sure about the cause of this note. Some relevant discussions are give                   [here](https://groups.google.com/g/r-sig-mac/c/7u_ivEj4zhM) and it seems that the test Server does not have `tidy` installed. This note only appeared in this build, and did not appear in any other builds. Thus, this note can be likely ignored.
