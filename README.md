dioscRi
======================================================
<img src=https://raw.githubusercontent.com/ecool50/dioscRi/main/inst/dioscRi_overview.png align="middle" height="200" width="200">


A deep learning framework that combines an `MMD-VAE` with hierarchical `Group-Lasso` for clinical prediction
in `high parameter cytometry assays`.

Overview
--------

**dioscRi** provides predictive modelling of `high parameter cytometry assays`.
This pipeline uses the `MMMD-VAE` architecture coupled with `Group-Lasso` that incorporates cell type hierarchies for predicting clinical outcomes.
We also provide functions for model visualisation and interpretation.

Installation
--------
Before installing this package, [tensorflow - 2.16.2](https://tensorflow.rstudio.com) and [keras3 - 1.2.0](https://keras.rstudio.com) must be installed in Python and connected to R .

If you would like the most up-to-date features, install the development version from GitHub.
```
# install.packages("devtools")
devtools::install_github("https://github.com/ecool50/dioscRi/")
library(dioscRi)
```
### Submitting an issue or feature request

`dioscRi` is still under active development. We would greatly appreciate any and 
all feedback related to the package.

* R package related issues should be raised [here](https://github.com/ecool50/dioscRi/issues).
* For general questions and feedback, please contact us directly via [ewil3501@uni.sydney.edu.au](mailto:ewil3501@uni.sydney.edu.au).


## Author

* **Elijah Willie**
* **Ellis Patrick**  - [@TheEllisPatrick](https://twitter.com/TheEllisPatrick)
