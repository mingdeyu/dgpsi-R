# Validation

The contribution was applied to upstream R dev commit
`12319fb6fc55df7136e23b8d36f26e02d6117047`. Tests and the tutorial were run on
Windows 11 with R 4.4.2, Python 3.10.20 and Python dgpsi dev commit
`480bb30260d3a1c22af1ddffe1a11df3e4de0820`.

## Checks

- Native integration suite: 11 tests and 41 expectations passed, with no failures
  or skips. Serial/parallel consistency and both GP and DGP prediction are covered.
- The installed package's complete tutorial script ran successfully.
- Source build and `R CMD check --no-manual --no-build-vignettes` completed
  with Status OK. The native tests and executed tutorial were checked separately.
- The example uses five scalar DGPs, 96 training rows, 48 validation rows and
  96 test rows. Training uses 40 iterations, `ess_burn=10` and five imputations.
- No seed, variable or test-period search was performed for this contribution.

## Example results

| Quantity | Initial | Selected |
|---|---:|---:|
| Internal connections | 6 | 4 |
| Validation mean NRMSE | 0.464404 | 0.452438 |
| Test mean NRMSE | 0.456788 | 0.462622 |

The selected graph removes both internal inputs of b while a and c retain theirs.
Validation mean NRMSE decreases by 2.58%; test mean NRMSE increases by 1.28%.
The validation RMSE of b decreases by 8.79%, while its test RMSE increases by
4.93%. This is a validation-selection example, not a claim of uniformly better
test predictions.

The source function, reference documentation and tutorial all use the same
selection rule. The initial network remains eligible throughout the search.

## Compatibility changes

The upstream dev linked constructor reorders component inputs differently from
the earlier master-based prototype. That implementation is preserved. Only the
one-layer linked prediction path and restoration of serialization metadata are
patched. The training burn-in is explicit because dev changed its default
from ten to five.

Earlier prototype results are not reused for this dev-based example. The small
benchmark evaluates filtered, standardised responses and is not a prospective
field evaluation.

## Repeating the checks

Follow [PRUNLGP.md](PRUNLGP.md) for installation, the complete tutorial and native
tests. Default package checks skip Python integration tests unless
`DGPSI_RUN_PRUNLGP_TESTS=true`. The executed tutorial is a separate check.
