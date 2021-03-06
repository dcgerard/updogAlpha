% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dxi_df}
\alias{dxi_df}
\title{Returns derivative to f / (d * (1 - f) + f)}
\usage{
dxi_df(d, f)
}
\arguments{
\item{d}{The bias term. Must be greater than 0.}

\item{f}{The sequencing-error adjusted probability of the refernence allele.}
}
\description{
Returns derivative to f / (d * (1 - f) + f)
}
\author{
David Gerard
}
