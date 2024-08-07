% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exposure_control_functions.R
\name{updateDiagnosticExposureRecord}
\alias{updateDiagnosticExposureRecord}
\title{(Internal) Update diagnostic exposure record}
\usage{
updateDiagnosticExposureRecord(
  diagnostic_exposure_record,
  j,
  exposure_record,
  config,
  simulation_constants
)
}
\arguments{
\item{diagnostic_exposure_record}{a named list containing diagnostic exposure record.}

\item{j}{the numeric index of the examinee.}

\item{exposure_record}{a named list containing exposure records.}

\item{config}{a \code{\linkS4class{config_Shadow}} object. Use \code{\link{createShadowTestConfig}} for this.}

\item{simulation_constants}{a named list containing simulation constants.}
}
\value{
\code{\link{updateDiagnosticExposureRecord}} returns an updated diagnostic exposure record.
}
\description{
\code{\link{updateDiagnosticExposureRecord}} is an internal function for updating diagnostic exposure record.
This function does not do actual exposure control;
this function is used for updating diagnostic exposure record after exposure control is performed.
}
\keyword{internal}
