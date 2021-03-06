% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_all.R
\name{fn_cov_mle}
\alias{fn_cov_mle}
\title{Wrapper for \code{\link{obj_offspring_reparam}}.}
\usage{
fn_cov_mle(par, ocounts, osize, ploidy, prob_geno, p1counts = NULL,
  p1size = NULL, p2counts = NULL, p2size = NULL, p1geno = NULL,
  p2geno = NULL, bias_val_mean = 0, bias_val_sd = 0.7,
  seq_error_mean = -4.7, seq_error_sd = 1, model = c("f1", "s1", "hw",
  "uniform"))
}
\arguments{
\item{par}{A numeric vector of length 4. The elements are \code{s},
\code{ell}, \code{r}, and the logit of \code{out_prop} from
\code{\link{obj_offspring_reparam}}.}

\item{ocounts}{The observed counts of the refernce
allele for each individual.}

\item{osize}{The observed number of reads for each
individuals.}

\item{ploidy}{An integer. The ploidy of the species. This is assumed
to be the same for each individual.}

\item{prob_geno}{The allele frequencies of the genotypes. See \code{\link{get_prob_geno}}.}

\item{p1counts}{A vector of non-negative integers. The ith element
is the number of reads of the reference allele in the ith
sample of parent 1.  If \code{NULL} then the prior
probabilities on parent 1's genotype will default to uniform.}

\item{p1size}{A vector of positive integers. The ith element is the
total number of reads in the ith sample of parent 1.  If
\code{NULL} then the prior probabilities on parent 1's genotype
will default to uniform.}

\item{p2counts}{A vector of non-negative integers. The ith element
is the number of reads of the reference allele in the ith
sample of parent 2.  If \code{NULL} then the prior
probabilities on parent 2's genotype will default to uniform.}

\item{p2size}{A vector of positive integers. The ith element is the
total number of reads in the ith sample of parent 2.  If
\code{NULL} then the prior probabilities on parent 2's genotype
will default to uniform.}

\item{p1geno}{The initial value of the first parental genotype.}

\item{p2geno}{The initial value of the second parental genotype.}

\item{bias_val_mean}{The prior mean on the log of \code{bias_val} (corresponding to \code{parvec[1]}).
Set \code{bias_val_sd = Inf} to have no penalty on the bias parameter.}

\item{bias_val_sd}{The prior standard deviation on the log of \code{bias_val}
(corresponding to \code{parvec[1]}).
Set \code{bias_val_sd = Inf} to have no penalty on the bias parameter.}

\item{seq_error_mean}{The mean of the logit-normal prior on the sequencing error rate, which corresponds
to \code{parvec[2]}. Set \code{seq_error_sd = Inf} to have no penalty on the sequencing error rate.
The default is -4.7, which roughly corresponds to a mean sequencing error value of 0.009. If you want
to constain \code{seq_error} to be zero, you need to set \code{update_seq_error = FALSE},
\code{seq_error = 0}, and\code{seq_error_mean = -Inf}.}

\item{seq_error_sd}{The standard deviation of the logit-normal prior on the sequencing error rate, which
corresponds to \code{parvec[2]}. The default is 1, which at three standard deviations is about a sequencing
error rate of 0.15. Set \code{seq_error_sd = Inf} to have no penalty on the sequencing error rate.}

\item{model}{The model for the genotype distribution. Do we assume an
F1 population (\code{"f1"}), an S1 population (\code{"s1"}), Hardy-Weinberg equilibrium (\code{"hw"}),
or a uniform distribution (\code{"uniform"}).}
}
\description{
This exists to calculate the Hessian of the parameters at the end of
\code{\link{updog_update_all}}.
}
\author{
David Gerard
}
