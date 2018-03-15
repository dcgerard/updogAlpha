
<!-- README.md is generated from README.Rmd. Please edit that file -->
updog: Using Parental Data for Offspring Genotyping
===================================================

[![Linux Build Status](https://travis-ci.org/dcgerard/updog.svg?branch=master)](https://travis-ci.org/dcgerard/updog) [![Windows Build Status](https://ci.appveyor.com/api/projects/status/github/dcgerard/updog?branch=master&svg=true)](https://ci.appveyor.com/project/dcgerard/updog) <!-- [![Coverage Status](https://img.shields.io/codecov/c/github/dcgerard/updog/master.svg)](https://codecov.io/github/dcgerard/updog?branch=master) --> <!-- [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0) -->

This package will fit an empirical Bayesian procedure to genotype autopolyploid individuals using reduced-representation next-generation sequencing (NGS) data, such as genotyping by sequencing (GBS) (Elshire 2011) or restriction site-associated DNA sequencing (RAD-Seq) (Baird 2008). For such NGS data there exist other methods for genotyping --- see for example [ebg](https://github.com/pblischak/polyploid-genotyping) (Blischak, Kubatko, and Wolfe 2017) and [TET](http://www.g3journal.org/lookup/suppl/doi:10.1534/g3.117.039008/-/DC1) (Maruki and Lynch 2017). `updog` adds to this field by:

-   Incorporating hierarchical structure into the genotyping procedure by either assuming the population consists of siblings (either an F1 or S1 population), or the population is in Hardy-Weinberg equilibrium. Though see Serang (2012) and Li (2011) for other discussion on this.
-   Estimating sequencing error for each SNP.
-   Accounting for read-mapping bias.
-   Incorporating overdispersion, a key indicator for the quality of a SNP.
-   Automatically accounting for outlying observations.
-   Beautifully plotting the NGS data (see `plot.updog`).
-   Running a simple goodness of fit test (see `summary.updog`).
-   Simulating under a fitted model (see `rupdog`).

We've included a few SNP's from the data of Shirasawa et al. (2017) to show off the features of `updog`. See `snpdata`.

A vignette is available [here](https://dcgerard.github.io/updog/articles/smells_like_updog.html).

Please report any bugs/issues [here](https://github.com/dcgerard/updog/issues).

Installation
============

To install, run the following code in R:

``` r
# install.packages("devtools")
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

References
==========

Baird, Paul D. AND Atwood, Nathan A. AND Etter. 2008. “Rapid SNP Discovery and Genetic Mapping Using Sequenced RAD Markers.” *PLOS ONE* 3 (10). Public Library of Science: 1–7. doi:[10.1371/journal.pone.0003376](https://doi.org/10.1371/journal.pone.0003376).

Blischak, Paul, Laura Kubatko, and Andrea Wolfe. 2017. “SNP Genotyping and Parameter Estimation in Polyploids Using Low-Coverage Sequencing Data.” *bioRxiv*. Cold Spring Harbor Labs Journals. doi:[10.1101/120261](https://doi.org/10.1101/120261).

Elshire, Jeffrey C. AND Sun, Robert J. AND Glaubitz. 2011. “A Robust, Simple Genotyping-by-Sequencing (GBS) Approach for High Diversity Species.” *PLOS ONE* 6 (5). Public Library of Science: 1–10. doi:[10.1371/journal.pone.0019379](https://doi.org/10.1371/journal.pone.0019379).

Li, Heng. 2011. “A Statistical Framework for SNP Calling, Mutation Discovery, Association Mapping and Population Genetical Parameter Estimation from Sequencing Data.” *Bioinformatics* 27 (21): 2987. doi:[10.1093/bioinformatics/btr509](https://doi.org/10.1093/bioinformatics/btr509).

Maruki, Takahiro, and Michael Lynch. 2017. “Genotype Calling from Population-Genomic Sequencing Data.” *G3: Genes, Genomes, Genetics* 7 (5). G3: Genes, Genomes, Genetics: 1393–1404. doi:[10.1534/g3.117.039008](https://doi.org/10.1534/g3.117.039008).

Serang, Marcelo AND Garcia, Oliver AND Mollinari. 2012. “Efficient Exact Maximum a Posteriori Computation for Bayesian SNP Genotyping in Polyploids.” *PLOS ONE* 7 (2). Public Library of Science: 1–13. doi:[10.1371/journal.pone.0030906](https://doi.org/10.1371/journal.pone.0030906).

Shirasawa, Kenta, Masaru Tanaka, Yasuhiro Takahata, Daifu Ma, Qinghe Cao, Qingchang Liu, Hong Zhai, et al. 2017. “A High-Density SNP Genetic Map Consisting of a Complete Set of Homologous Groups in Autohexaploid Sweetpotato (Ipomoea Batatas).” *Scientific Reports* 7. Nature Publishing Group. doi:[10.1038/srep44207](https://doi.org/10.1038/srep44207).
