# dioscRi 1.0.0

Manuscript-accompanying release. This is the version of the package
used to produce the analyses in *dioscRi enables transferable
prediction of clinical outcomes in multi-parameter cytometry data*.
For the manuscript code, analyses, and reproduction workflow see the
companion repository,
[ecool50/dioscRi_manuscript](https://github.com/ecool50/dioscRi_manuscript).

## Core capabilities

* MMD-VAE normalisation for high-parameter cytometry data
  (`trainVAEModel`, `decodeSamples`), applicable to unseen samples
  within a compatible marker panel without retraining.
* Frobenius-norm reference-sample selection
  (`computeReferenceSample`): chooses training samples whose
  covariance lies near the cohort centroid so the MMD-VAE learns
  broadly applicable normalisation rather than fitting outlier or
  single-batch structure.
* LDA / KNN / multinomial cell-type classifier
  (`trainCellTypeClassifier`) for transferring cluster or
  manually-gated labels from training to test samples.
* Sample-level features (`computeFeatures`): cell-type proportions
  (with optional logit transform) and per-cell-type marker means.
* Hierarchical overlapping group LASSO (`fitModel`,
  `generateTree`, `generateGroups`) for interpretable clinical
  prediction, with both data-driven and hand-built hierarchies
  supported.
* Visualisation and interpretation: `plotAUC`, `plotSigFeatures`,
  `visualiseModelTree`, `getSigFeatures`.

## Behaviour and API

* `trainVAEModel()` now accepts `earlyStop`, `patience`, and
  `verbose` arguments. Defaults: `earlyStop = TRUE, patience = 15,
  verbose = 1`. The manuscript vignettes use `earlyStop = FALSE,
  verbose = 0L` to match the published 100-epoch training run.
* `trainCellTypeClassifier()` supports three back-ends via the
  `model` argument: `"lda"` (default), `"knn"`, `"multinom"`. The
  caret dependency was removed in favour of direct calls to
  `MASS::lda` with auto-subsampling for large cohorts.
* `fitModel()` defaults: alpha grid is the 10 values
  {0.1, 0.2, ..., 1.0}, lambda is a fixed 10,000-point uniform grid
  over [0.05, 1]; selection is per-alpha elbow on BIC, then global
  minimum over alpha.

## Bug fixes

* `plotAUC()` now imports `ggtitle` from ggplot2 (previously
  relied on the caller having ggplot2 attached).
* `computeFeatures()` now uses `tidyr::replace_na` explicitly
  (previously relied on the caller having tidyr attached).
* `inst/CITATION` syntax error (stray `e` at end of line) fixed so
  `citation("dioscRi")` works.

## Documentation

* New pkgdown site at
  [ecool50.github.io/dioscRi](https://ecool50.github.io/dioscRi).
* Vignettes:
  * Quickstart on the bundled toy dataset.
  * BioHEART-CT walkthrough, unsupervised-clustering pipeline.
  * BioHEART-CT walkthrough, manually gated cell types.
  * Adding clinical covariates (age, sex) to the model.

## Datasets evaluated

The package was evaluated on four cytometry datasets in the
manuscript (raw and processed inputs archived on
[Zenodo](https://zenodo.org/records/15694581)):

* BioHEART-CT (CyTOF; CAD prediction).
* Wagner et al. breast cancer (CyTOF; tumour vs. non-tumour).
* CMV multi-study (CyTOF; SDY519 held out for evaluation).
* Mathew et al. COVID-19 (flow cytometry; CD8+ non-naive T cells).

## Pinned dependencies

| Component             | Version          |
|-----------------------|------------------|
| R                     | 4.5.0            |
| `keras3` (R)          | 1.2.0            |
| `tensorflow` (Python) | 2.16.2           |
| Python                | 3.10.15 (`r-reticulate` virtualenv) |

A scripted helper that installs these matched versions is available
at `revision/scripts/helpers/setup_environment.R` in the companion
manuscript repository.
