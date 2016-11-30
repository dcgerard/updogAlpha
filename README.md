
<!-- README.md is generated from README.Rmd. Please edit that file -->
updog: Using Parental Data for Offspring Genotyping
===================================================

[![Linux Build Status](https://travis-ci.org/dcgerard/updog.svg?branch=master)](https://travis-ci.org/dcgerard/updog) [![Windows Build Status](https://ci.appveyor.com/api/projects/status/github/dcgerard/updog?branch=master&svg=true)](https://ci.appveyor.com/project/dcgerard/updog) [![Coverage Status](https://img.shields.io/codecov/c/github/dcgerard/updog/master.svg)](https://codecov.io/github/dcgerard/updog?branch=master)

Suppose you have sequencing data from possibly polyploid individuals that you want to genotype. If you know that these individuals are siblings, then you can borrow strength between the individuals to estimate each individual's genotype. This is the main idea behind `updog`. This package will allow you to

-   Estimate the parental genotypes using just the offspring genotypes.
-   Incorporate this information in estimating the offspring genotypes.
-   Include parental sequence data, if available, to improve the estimates of the offspring genotypes.
-   Incorporate sequencing error in genotype estimation for polyploid individuals.
-   Estimate the sequencing error rate for a particular SNP.
-   Estimate and incorporate overdispersion in the sequencing process.
-   Automatically detect and remove outliers.

Please report any bugs/issues [here](https://github.com/dcgerard/updog/issues).

Installation
============

To install, run the following code in R:

``` r
install.packages("devtools")
devtools::install_github("dcgerard/updog")
```

Citation
========

Please cite this package as

> Gerard, David. 2016. *updog: Using Parental Data for Offspring Genotyping*. <https://github.com/dcgerard/updog>.

Or, using BibTex:

``` tex
@Manual{gerard2016updog,
    title = {{updog}: Using Parental Data for Offspring Genotyping},
    author = {David Gerard},
    year = {2016},
    note = {R package version 0.1.0},
    url = {https://github.com/dcgerard/updog},
}
```
