% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{plot_beta_dist_gg}
\alias{plot_beta_dist_gg}
\title{Plot the beta distribution using ggplot2.}
\usage{
plot_beta_dist_gg(alpha = NULL, beta = NULL, mu = NULL, rho = NULL)
}
\arguments{
\item{alpha}{The shape1 parameter.}

\item{beta}{The shape2 parameter.}

\item{mu}{The mean parameter.}

\item{rho}{The overdispersion parameter.}
}
\description{
This is the same as \code{\link{plot_beta_dist}}, except we use \code{ggplot2} to do the plotting.
See \code{\link{plot_beta_dist}} for details.
}
\seealso{
\code{\link{plot_beta_dist}} for the R Base graphics version of this function.
}
\author{
David Gerard
}
