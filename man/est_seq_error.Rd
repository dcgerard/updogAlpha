% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_seq.R
\name{est_seq_error}
\alias{est_seq_error}
\title{Estimate sequencing error from counts.}
\usage{
est_seq_error(ncounts, ssize, ploidy = NULL, eps = NULL)
}
\arguments{
\item{ncounts}{A vector of non-negative integers. The ith element
is the number of counts of the ith sample.}

\item{ssize}{A vector of positive integers. The ith element is the
total number of counts of the ith sample.}

\item{ploidy}{The ploidy of the species. Cannot be \code{NULL} if \code{eps} is \code{NULL}.}

\item{eps}{To estimate the sequencing erorr, we will only use points such that \code{ncounts / ssize >= eps}.
If \code{ploidy} is specified and \code{eps = NULL}, then \code{eps} will be chosen by a heuristic.}
}
\description{
Assuming that the reference allele is in great enough quantities so that there are enough points
}
\author{
David Gerard
}
