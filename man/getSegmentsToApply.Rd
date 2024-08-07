% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exposure_control_functions.R
\name{getSegmentsToApply}
\alias{getSegmentsToApply}
\title{(Internal) Make a vector for segment-dimensioned matrix update}
\usage{
getSegmentsToApply(n_segment, segments)
}
\arguments{
\item{n_segment}{the total number of segments.}

\item{segments}{a vector of target segments.}
}
\value{
\code{\link{getSegmentsToApply}} returns a vector.
}
\description{
\code{\link{getSegmentsToApply}} is an internal function for
creating a vector for the purpose of updating various matrices that use theta segments as a dimension.
}
\keyword{internal}
