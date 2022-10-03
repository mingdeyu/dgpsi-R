# dgpsi 2.1.5-9000 (development version)

- A bug is found in multi-core predictions in `predict()` when `object` is an instance of `lgp` class and `x` is a list. This bug needs a modification in the underlying Python implementation and will be fixed in the next `dgpsi` release version `2.1.6`. The single-core predictions in the current release version `2.1.5` are not affected by the bug.  

# dgpsi 2.1.5

- Initial release of the R interface to the Python package `dgpsi v2.1.5`.
