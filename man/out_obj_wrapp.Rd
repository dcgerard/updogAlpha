% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_all.R
\name{out_obj_wrapp}
\alias{out_obj_wrapp}
\title{Wrapper for \code{\link{outlier_obj}}.}
\usage{
out_obj_wrapp(parvec, ocounts, osize, weight_vec, p1counts = NULL,
  p1size = NULL, p1weight = NULL, p2counts = NULL, p2size = NULL,
  p2weight = NULL, min_disp = 0)
}
\arguments{
\item{parvec}{A numeric vector of length two. The first element is the outlier mean.
The second element is the outlier overdispersion parameter.}

\item{ocounts}{The offspring counts of the reference allele.}

\item{osize}{The offspring counts of reads.}

\item{weight_vec}{The probability a point is an outlier.}

\item{p1counts}{The parent 1 counts of the reference allele.}

\item{p1size}{The parent 1 counts of reads.}

\item{p1weight}{The probability parent 1 is an outlier.}

\item{p2counts}{The parent 2 counts of the reference allele.}

\item{p2size}{The parent 2 counts of reads.}

\item{p2weight}{The probability parent 2 is an outlier.}

\item{min_disp}{The bound on the overdispersion parameter (can't be too small).}
}
\description{
Wrapper for \code{\link{outlier_obj}}.
}
\author{
David Gerard
}
