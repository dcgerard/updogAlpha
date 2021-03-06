% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_all.R
\name{grad_wrapp_all}
\alias{grad_wrapp_all}
\title{wrapper for \code{\link{grad_offspring_weights}}}
\usage{
grad_wrapp_all(parvec, ocounts, osize, weight_vec, ploidy, p1geno, p2geno,
  allele_freq = 0.5, p1counts = NULL, p1size = NULL, p1weight = NULL,
  p2counts = NULL, p2size = NULL, p2weight = NULL, bound_bias = FALSE,
  update_bias_val = TRUE, update_seq_error = TRUE, update_od_param = TRUE,
  seq_error_mean = -4.7, seq_error_sd = 1, bias_val_mean = 0,
  bias_val_sd = 0.7, model = c("f1", "s1", "hw", "uniform"),
  is_seq_ninf = FALSE)
}
\arguments{
\item{parvec}{A vector of three elements, s, ell, and r. We have s = log(bias_val) = log(d),
ell = logit(seq_error) = logit(eps), and r = - logit(od_param) = - logit(tau).}

\item{ocounts}{The observed counts of the refernce
allele for each individual.}

\item{osize}{The observed number of reads for each
individuals.}

\item{weight_vec}{A vector of weights obtained via the E-step.}

\item{ploidy}{An integer. The ploidy of the species. This is assumed
to be the same for each individual.}

\item{p1geno}{The genotype of parent 1 if \code{model = "f1"} or \code{model = "s1"}.}

\item{p2geno}{The genotype of parent 2 if \code{model = "f1"}. This needs to be null if \code{model = "s1"}}

\item{allele_freq}{The allele-frequency if \code{model = "hw"}}

\item{p1counts}{The number of reference alleles observed from parent 1.}

\item{p1size}{The number of reads observed from parent 1.}

\item{p1weight}{The posterior probability that parent 1 is not an outlier.}

\item{p2counts}{The number of reference alleles observed from parent 2.}

\item{p2size}{The number of reads observed from parent 2.}

\item{p2weight}{The posterior probability that parent 2 is not an outlier.}

\item{bound_bias}{A logical. Should we bound the bias parameter \code{parvec[1]} by a somewhat arbitrary value
(\code{TRUE}) or not (\code{FALSE})?}

\item{update_bias_val}{A logical. If \code{FALSE}, then the first position
of the returned gradient will be zero.}

\item{update_seq_error}{A logical. If \code{FALSE}, then the second position
of the returned gradient will be zero.}

\item{update_od_param}{A logical. If \code{FALSE}, then the third position
of the returned gradient will be zero.}

\item{seq_error_mean}{The mean of the logit-normal prior on the sequencing error rate, which corresponds
to \code{parvec[2]}. Set \code{seq_error_sd = Inf} to have no penalty on the sequencing error rate.
The default is -4.7, which roughly corresponds to a mean sequencing error value of 0.009.}

\item{seq_error_sd}{The standard deviation of the logit-normal prior on the sequencing error rate, which
corresponds to \code{parvec[2]}. The default is 1, which at three standard deviations is about a sequencing
error rate of 0.15.}

\item{bias_val_mean}{The prior mean on the log of \code{bias_val} (corresponding to \code{parvec[1]}).
Set \code{bias_val_sd = Inf} to have no penalty on the bias parameter.}

\item{bias_val_sd}{The prior standard deviation on the log of \code{bias_val}
(corresponding to \code{parvec[1]}).
Set \code{bias_val_sd = Inf} to have no penalty on the bias parameter.}

\item{model}{The model for the genotype distribution. Should we assume an F1 population (\code{"f1"}),
an S1 population (\code{"s1"}),
Hardy-Weinberg equilibrium (\code{"hw"}), or a uniform distribution (\code{"uniform"})?}

\item{is_seq_ninf}{Hacky way to get around fact that optim won't allow \code{-Inf} in
\code{par}, even if \code{gr} will always be zero for that parameter.
How dare they not anticipate my unique situation!}
}
\description{
wrapper for \code{\link{grad_offspring_weights}}
}
\author{
David Gerard
}
