% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dbeta_dtau_withxi}
\alias{dbeta_dtau_withxi}
\title{Derivative of overdispersion parameter with mean already calculated}
\usage{
dbeta_dtau_withxi(x, n, xi, tau)
}
\arguments{
\item{x}{The observed counts of reference allele.}

\item{n}{The total number of observed counts.}

\item{xi}{The mean proportion.}

\item{tau}{The overdispersion parameter.}
}
\description{
Derivative of overdispersion parameter with mean already calculated
}
\seealso{
dbeta_dtau
}
\author{
David Gerard
}
