# Bayesian nonparametric modeling of heterogeneous populations of networks

This repository provides the companion code for the published article:

> Barile, F., Lunagómez, S., & Nipoti, B. (2025). *Bayesian Nonparametric Modeling of Heterogeneous Populations of Networks*. **Bayesian Analysis**. [DOI: 10.1214/26-BA1588](https://doi.org/10.1214/26-BA1588)


## Installation

The `DPmCER` R package can be installed directly from GitHub using the `remotes` package:

```r
# install.packages("remotes")
remotes::install_github("FrancescoBarile/DPmCER")
```

Once installed, the package can be loaded with:

```r
library(DPmCER)
```

## Reproducing the simulations
All experiments on synthetic data described in **Section 4 of the main article** and **Section 4 of the Supplementary Material** are implemented in this repository.

After installing the `DPmCER` package, the simulation studies can be reproduced by running the R scripts available in:

```text
reproduction/simulations/
```



## Reference

If you use the code or the `DPmCER` package, please cite:

```bibtex
@article{10.1214/26-BA1588,
  author = {Francesco Barile and Sim{\'o}n Lunag{\'o}mez and Bernardo Nipoti},
  title = {{Bayesian Nonparametric Modeling of Heterogeneous Populations of Networks}},
  journal = {Bayesian Analysis},
  publisher = {International Society for Bayesian Analysis},
  pages = {1--27},
  year = {2025},
  doi = {10.1214/26-BA1588},
  url = {https://doi.org/10.1214/26-BA1588}
}
```
