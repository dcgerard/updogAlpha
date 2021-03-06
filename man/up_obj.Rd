% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updog.R
\name{up_obj}
\alias{up_obj}
\title{The objective function}
\usage{
up_obj(pival, p1geno, p2geno, alpha, beta, ocounts, osize, dbinommat, qarray,
  r1vec, r2vec)
}
\arguments{
\item{pival}{A proportion. The probability of not being a mistake.}

\item{p1geno}{The genotype of parent 1}

\item{p2geno}{The genotype of parent 2}

\item{alpha}{The alpha parameter in \code{\link{dbetabinom}}.}

\item{beta}{The beta parameter in \code{\link{dbetabinom}}.}

\item{ocounts}{A vector of non-negative integers. The ith element
is the number of reads of the reference allele in the ith child.}

\item{osize}{A vector of positive integers. The ith element is the
total number of reads for the ith child.}

\item{dbinommat}{A matrix. The rows index the genotype and the
columns index the offspring.  These are the binomial
probabilities of the count given the genotype.}

\item{qarray}{An three-way array of proportions. The (i, j, k)th
element is the probability of an offspring having k - 1
reference alleles given that parent 1 has i - 1 refrerence
alleles and parent 2 has j - 1 reference alleles. Each
dimension of the array is \code{ploidy + 1}.}

\item{r1vec}{A vector of prior probabilities whose ith element is
the prior probability that parent 1 has i -1 copies of allele
A.}

\item{r2vec}{A vector of prior probabilities whose ith element is
the prior probability that parent 2 has i -1 copies of allele
A.}
}
\description{
The objective function
}
\author{
David Gerard
}
