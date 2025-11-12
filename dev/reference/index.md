# Package index

## Initalization

The function to initiate the underlying Python environment.

- [`init_py()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/init_py.md)
  : 'python' environment initialization

## GP, DGP, and Linked (D)GP Emulations

Functions to train, validate, and make predictions from GP, DGP, and
Linked (D)GP emulators.

- [`gp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/gp.md)
  **\[updated\]** : Gaussian process emulator construction
- [`dgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/dgp.md)
  **\[updated\]** : Deep Gaussian process emulator construction
- [`lgp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/lgp.md) :
  Linked (D)GP emulator construction
- [`continue()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/continue.md)
  : Continue training a DGP emulator
- [`predict(`*`<dgp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  [`predict(`*`<lgp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  [`predict(`*`<gp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/predict.md)
  : Prediction from GP, DGP, or linked (D)GP emulators
- [`validate()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/validate.md)
  : Validate a constructed GP, DGP, or linked (D)GP emulator
- [`plot(`*`<dgp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  [`plot(`*`<lgp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  [`plot(`*`<gp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/plot.md)
  : Validation plots of a constructed GP, DGP, or linked (D)GP emulator
- [`prune()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/prune.md)
  : Static pruning of a DGP emulator

## Sequential Design

Functions to implement sequential designs for (D)GP and bundles of (D)GP
emulators.

- [`design()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/design.md)
  : Sequential design of a (D)GP emulator or a bundle of (D)GP emulators
- [`alm()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/alm.md) :
  Locate the next design point(s) for a (D)GP emulator or a bundle of
  (D)GP emulators using Active Learning MacKay (ALM)
- [`mice()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/mice.md) :
  Locate the next design point for a (D)GP emulator or a bundle of (D)GP
  emulators using MICE
- [`vigf()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/vigf.md) :
  Locate the next design point for a (D)GP emulator or a bundle of (D)GP
  emulators using VIGF
- [`update()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/update.md)
  : Update a GP or DGP emulator
- [`draw()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/draw.md) :
  Validation and diagnostic plots for a sequential design
- [`pack()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/pack.md) :
  Pack GP and DGP emulators into a bundle
- [`unpack()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/unpack.md)
  : Unpack a bundle of (D)GP emulators

## Helpers

- [`summary(`*`<gp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/summary.md)
  [`summary(`*`<dgp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/summary.md)
  [`summary(`*`<lgp>`*`)`](http://mingdeyu.github.io/dgpsi-R/dev/reference/summary.md)
  : Summary of a constructed GP, DGP, or linked (D)GP emulator
- [`write()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/write.md)
  : Save the constructed emulator
- [`read()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/read.md) :
  Load the stored emulator
- [`set_seed()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_seed.md)
  : Random seed generator
- [`set_imp()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_imp.md)
  : Reset number of imputations for a DGP emulator
- [`window()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/window.md)
  : Trim the sequence of hyperparameter estimates within a DGP emulator
- [`trace_plot()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/trace_plot.md)
  : Trace plot for DGP hyperparameters
- [`nllik()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/nllik.md)
  : Calculate the predictive negative log-likelihood
- [`set_thread_num()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_thread_num.md)
  : Set the number of threads
- [`get_thread_num()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/get_thread_num.md)
  : Get the number of threads
- [`set_vecchia()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_vecchia.md)
  : Add or remove the Vecchia approximation
- [`set_id()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/set_id.md)
  : Set Emulator ID
- [`serialize()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/serialize.md)
  : Serialize the constructed emulator
- [`deserialize()`](http://mingdeyu.github.io/dgpsi-R/dev/reference/deserialize.md)
  : Restore the serialized emulator
