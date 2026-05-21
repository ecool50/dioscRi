dioscRi
======================================================
A deep learning framework that combines an `MMD-VAE` with hierarchical `Group-Lasso` for clinical prediction
in `high parameter cytometry assays`.

<img src=https://raw.githubusercontent.com/ecool50/dioscRi/main/inst/dioscRi_overview.jpg align="middle" height="500" width="1000">

Overview
--------

**dioscRi** provides predictive modelling of `high parameter cytometry assays`.
This pipeline uses the `MMMD-VAE` architecture coupled with `Group-Lasso` that incorporates cell type hierarchies for predicting clinical outcomes.
We also provide functions for model visualisation and interpretation.

Installation
--------
Before installing this package, [tensorflow - 2.16.2](https://tensorflow.rstudio.com) and [keras3 - 1.2.0](https://keras.rstudio.com) must be installed in Python and connected to R.

If you would like the most up-to-date features, install the development version from GitHub.
```
# install.packages("devtools")
devtools::install_github("https://github.com/ecool50/dioscRi/")
library(dioscRi)
```

Pinned versions used in the manuscript:

| Component | Version |
|---|---|
| R | 4.5.0 |
| `keras3` (R) | 1.2.0 |
| `tensorflow` (Python) | 2.16.2 |
| Python | 3.10.15 (virtualenv `r-reticulate`) |

A scripted helper that installs these matched versions is available at
`revision/scripts/helpers/setup_environment.R` in the
[companion manuscript repository](https://github.com/ecool50/dioscRi_manuscript).

Reproducing the manuscript
--------
The manuscript analysis code, response materials, and step-by-step reproduction workflow
live in the companion repository:
[ecool50/dioscRi_manuscript](https://github.com/ecool50/dioscRi_manuscript). The release tag
`v1.0.0` (both repositories) is the immutable reference for this manuscript.

### Submitting an issue or feature request

`dioscRi` is still under active development. We would greatly appreciate any and 
all feedback related to the package.

* R package related issues should be raised [here](https://github.com/ecool50/dioscRi/issues).
* For general questions and feedback, please contact us directly via [ewil3501@uni.sydney.edu.au](mailto:ewil3501@uni.sydney.edu.au).


## Author

* **Elijah Willie**
* **Ellis Patrick**  - [@TheEllisPatrick](https://twitter.com/TheEllisPatrick)
