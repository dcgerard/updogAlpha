% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{obj_offspring_weights_reparam}
\alias{obj_offspring_weights_reparam}
\title{Reparameterization of \code{\link{obj_offspring_weights}}.}
\usage{
obj_offspring_weights_reparam(ocounts, osize, weight_vec, ploidy, prob_geno, s,
  ell, r)
}
\arguments{
\item{ocounts}{The observed counts of the refernce
allele for each individual.}

\item{osize}{The observed number of reads for each
individuals.}

\item{weight_vec}{A vector of numerics between 0 and 1. They don't have to sum to 1.}

\item{ploidy}{An integer. The ploidy of the species. This is assumed
to be the same for each individual.}

\item{prob_geno}{The allele frequencies of the genotypes. See \code{\link{get_prob_geno}}.}

\item{s}{Same as \code{log(bias_val)} in \code{\link{obj_offspring}}.}

\item{ell}{We have \code{seq_error = expit(ell)} from \code{\link{obj_offspring}}.}

\item{r}{Same as \code{log((1.0 - od_param) / od_param)} from \code{\link{obj_offspring}}.}
}
\description{
It doesn't make any sense to have outlier = true since the weights are just for the EM already.
}
\author{
David Gerard
}
