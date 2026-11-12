# Using prunlgp

This branch adds validation-based internal connection selection to
[Deyu Ming's dgpsi-R package](https://github.com/mingdeyu/dgpsi-R/tree/dev).
Each component has one output. Removing an input refits that receiving model,
then recomputes predictions through the complete linked network.

- [Executed tutorial](https://zhuangchaoqun.github.io/dgpsi-R/articles/prunlgp.html)
- [Function reference](https://zhuangchaoqun.github.io/dgpsi-R/reference/prunlgp.html)
- [Source](https://github.com/zhuangchaoqun/dgpsi-R/blob/prunlgp/R/prunlgp.R)

## 1. Install the R package

Use R 4.0 or newer. The checks for this contribution used R 4.4.2 and Python 3.10.
Install Rtools on Windows if an R dependency needs compilation. RStudio includes
Pandoc; otherwise install Pandoc to render the tutorial.

In a fresh R session:

```r
install.packages(c("remotes", "rmarkdown", "testthat"))
remotes::install_github(
  "zhuangchaoqun/dgpsi-R", ref="prunlgp",
  dependencies=NA, upgrade="never", build_vignettes=FALSE
)
library(dgpsi)
stopifnot("prunlgp" %in% getNamespaceExports("dgpsi"))
```

If you previously installed `mingdeyu/dgpsi-R@dev`, run the same installation
command to replace it with this branch, then restart R. The package name remains
`dgpsi`.

Sourcing only `R/prunlgp.R` is insufficient: the function also needs registered
methods, prediction support for a network with no internal connections, and
linked-model serialization fixes. Installing the full branch includes these
changes together.

## 2. Initialise Python

```r
library(dgpsi)
init_py(dgpsi_ver="dev")
set_thread_num(1)
```

On first use, `init_py()` installs or locates its Conda environment and Python
dependencies. Complete any installation prompt and restart R if requested.
Run the two commands again in the fresh session. Initial environment setup can
take considerably longer than the example itself.

If Python has already been loaded from another environment, restart R before
calling `init_py()`. For Linux shared-library errors, follow the diagnostics
printed by `init_py()` and the upstream installation guide.

### Reproduce the published example

The published run uses the upstream R dev commit
`12319fb6fc55df7136e23b8d36f26e02d6117047` plus this contribution, and Python
dgpsi dev commit `480bb30260d3a1c22af1ddffe1a11df3e4de0820`.
The R example explicitly sets `ess_burn=10` instead of relying on a
development-version default. Record package versions when comparing results.

The example uses the same data and fitting budget throughout selection.
Different upstream revisions can change fitted predictions and selected graphs;
see [VALIDATION.md](VALIDATION.md) for the measured results and software versions.

## 3. Run the tutorial from start to finish

Clone this branch:

```sh
git clone --branch prunlgp https://github.com/zhuangchaoqun/dgpsi-R.git
cd dgpsi-R
Rscript tools/run-prunlgp-tutorial.R
```

Alternatively, open the clone in RStudio and run:

```r
source("tools/run-prunlgp-tutorial.R")
```

The script trains five scalar DGPs, constructs six candidate internal connections,
selects inputs using validation data, predicts the held-out test period and renders
the English tutorial. It writes the following files to `tutorial-output/`:

| File | Contents |
|---|---|
| `prunlgp.html` | Tutorial, interactive networks and error plots |
| `errors.csv` | Per-output validation and test RMSE |
| `connections.csv` | Retained and removed inputs for each output |
| `search.csv` | All evaluated candidates, including failures |
| `structure.csv` | Selected network's actual connection table |
| `system-errors.csv` | Mean normalised RMSE |
| `intervals.csv` | Coverage, interval width and interval score |
| `session-info.txt` | R package versions |
| `results.rds` | Numeric results and selection settings |

No previous model cache is needed. Allow several minutes for the example after
Python setup; the first run also compiles numerical routines. Runtime depends
on the machine and backend version.

A successful run ends with `Tutorial completed`. Its internal checks confirm
three scored outputs, six initial dependencies, finite errors and a selected
validation objective no higher than the initial objective. The full candidate
history is retained even when no deletion improves validation performance.

## 4. Use the function with another network

```r
fit <- prunlgp(
  network, components,
  x=x_validation,
  y=list(a=y_a_validation, b=y_b_validation),
  refit=function(model, keep, seed) {
    dgpsi::set_seed(seed)
    dgpsi::dgp(
      model$data$X[, keep, drop=FALSE], model$data$Y,
      id=model$id, name="sexp", depth=2, node=2,
      N=40, B=5, ess_burn=10, share=FALSE,
      bounds=c(0.1, 10), prior="ga", nugget_est=TRUE, verb=FALSE
    )
  },
  search="subsets", max_rounds=2, max_evals=32, seed=57
)
prediction <- predict(fit, x_test, full_layer=TRUE)
summary(fit)
fit$history
```

Supply one fitted model per output and the original training matrices.
`keep` contains indices into each model's original input columns. External
inputs remain available; internal validation inputs come from upstream predictions.

The default objective is mean RMSE normalised by training-response standard
deviations. Set `protect=0` to disallow any scored output's validation RMSE
from exceeding its initial value, or `protect=0.01` to allow a 1% tolerance.
Those constraints define different searches. They do not guarantee lower test error.

For larger input sets, `search="backward"` evaluates individual deletions.
For parallel candidate evaluation, use `workers=2, threads=1` and provide
`worker_init=function() dgpsi::init_py(dgpsi_ver="dev", verb=FALSE)`.
Use self-contained callbacks that can be sent to another R process.

## 5. Run tests

From the repository root, after installing this branch:

```r
Sys.setenv(DGPSI_RUN_PRUNLGP_TESTS="true", DGPSI_TEST_BACKEND="dev")
library(dgpsi)
testthat::test_file("tests/testthat/test-prunlgp.R")
```

The native tests cover output-specific deletion, grouped copies, upstream
information propagation, the empty-link case, serialization, reproducibility,
parallel candidates and invalid refits. They are opt-in because they require
a working Python environment.

For package checks:

```sh
R CMD build --no-manual .
R CMD check --no-manual --no-build-vignettes dgpsi_2.6.0-9000.tar.gz
```

Default package checks do not train the tutorial. Run the tutorial script
separately to verify its computations. See [VALIDATION.md](VALIDATION.md)
for the recorded versions and results.

## 6. Update after upstream review

Keep the `prunlgp` branch so shared source links remain valid. After the
contribution is merged upstream:

```sh
git remote add upstream https://github.com/mingdeyu/dgpsi-R.git
git fetch upstream
git switch prunlgp
git merge upstream/dev
git push origin prunlgp
```

Add the remote only if it does not already exist. Resolve any reported conflicts,
then repeat the tests and tutorial before pushing. The fork's Pages site uses
the same `articles/prunlgp.html` and `reference/prunlgp.html` paths.
An upstream website has its own addresses; merging does not redirect the fork's links.

Rebuild the website after updating the code:

```r
install.packages("pkgdown")
```

```sh
Rscript tools/build-prunlgp-site.R
git fetch origin gh-pages
git worktree add ../dgpsi-R-pages origin/gh-pages
cp -R docs/. ../dgpsi-R-pages/
git -C ../dgpsi-R-pages add .
git -C ../dgpsi-R-pages commit -m "Update prunlgp documentation"
git -C ../dgpsi-R-pages push origin HEAD:gh-pages
```

The shell commands above use Bash (or Git Bash on Windows). Create the Pages
worktree once and reuse it for later updates. Keep GitHub Pages configured to
publish the root of `gh-pages`. The site build executes the tutorial and requires
the same Python environment as the example.

The sample is de-identified operational data with neutral variable names,
causal filtering and training-only scaling. Raw measurements, dates, physical
variable mappings and inverse transforms are not distributed.
