# dioscRi v1.0.0

**Manuscript-accompanying release.** This is the version of the `dioscRi` R package used to produce all analyses in:

> *dioscRi enables transferable prediction of clinical outcomes in multi-parameter cytometry data*

For the manuscript code, analyses, and reproduction workflow, see the companion repository: [ecool50/dioscRi_manuscript](https://github.com/ecool50/dioscRi_manuscript).

## What the package provides

- **MMD-VAE normalization** — a Maximum Mean Discrepancy variational autoencoder that learns a normalization for high-parameter cytometry data and can be applied to unseen samples (within a compatible marker panel) without retraining.
- **Frobenius-norm reference-sample selection** — chooses training samples whose covariance lies near the cohort centroid, biasing the latent space toward broadly applicable representations.
- **Hierarchical overlapping group LASSO** — a structured-regularization model that incorporates biologically or empirically derived cell-type hierarchies for interpretable clinical prediction.
- **Helpers** — model visualization and post-hoc interpretation utilities.

## Pinned dependencies

| Component | Version |
|---|---|
| R | 4.5.0 |
| `keras3` (R) | 1.2.0 |
| `tensorflow` (Python) | 2.16.2 |
| Python | 3.10.15 (virtualenv `r-reticulate`) |

A scripted helper that installs these matched versions is available at `revision/scripts/helpers/setup_environment.R` in the companion repository.

## Datasets evaluated

The package was evaluated on four cytometry datasets (raw + processed inputs archived on [Zenodo](https://zenodo.org/records/15694581)):

- BioHEART-CT (CyTOF; CAD prediction)
- Wagner et al. breast cancer (CyTOF; tumor vs. non-tumor)
- CMV multi-study (CyTOF; SDY519 held out for evaluation)
- Mathew et al. COVID-19 (flow cytometry; CD8+ non-naive T cells)

## License

GPL-3.0
