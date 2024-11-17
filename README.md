

<!-- README.md is generated from README.Rmd. Please edit that file -->



# env.MRmega (Environment-adjusted MR-MEGA)

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11047160.svg)](https://doi.org/10.5281/zenodo.11047160)
<!-- badges: end -->

This is an R package for performing environment-adjusted MR-MEGA analysis aggregating all involved GWAS data. 

Website available: <https://SiruRooney.github.io/env.MRmega/>

For more details, please see:

S Wang, OO Ojewunmi,  A Kamiza, M Ramsay, AP Morris, T Chikowore, S Fatumo, JL Asimit. Accounting for heterogeneity due to environmental sources in meta-analysis of genome-wide association studies. *Communications Biology*. 2024 Nov 14;7(1):1512. https://doi.org/10.1038/s42003-024-07236-9 [PDF](https://rdcu.be/d0hIh)


## System requirements
env.MRmega could be installed with ease on versions of R > 4.3.0.


## Installation

Before installing R package of env.MRmega, the following packages from CRAN are required:

```r
install.packages("data.table")
install.packages("doMC")
install.packages("doParallel")
install.packages("dplyr")
install.packages("iterators")
install.packages("parallel")
install.packages("R.utils")
```

You can install the GitHub version of env.MRmega:

``` r
# install.packages("devtools")
devtools::install_github("SiruRooney/env.MRmega")
```

