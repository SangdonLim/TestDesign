% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/theta_functions.R
\name{computeEAPFromPosterior}
\alias{computeEAPFromPosterior}
\title{(Internal) Convert posterior densities into an EAP estimate}
\usage{
computeEAPFromPosterior(posterior, theta_grid)
}
\arguments{
\item{posterior}{a named list posterior densities and likelihoods.}

\item{theta_grid}{a vector containing quadrature points corresponding to the above.}
}
\value{
\code{\link{computeEAPFromPosterior}} returns a named list containing an EAP theta estimate.
}
\description{
\code{\link{computeEAPFromPosterior}} is an internal function for converting posterior densities into an EAP estimate.
}
\keyword{internal}
