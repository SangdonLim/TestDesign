% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loading_functions.R
\docType{class}
\name{item_attrib-class}
\alias{item_attrib-class}
\alias{loadItemAttrib}
\title{Load item attributes}
\usage{
loadItemAttrib(object, pool)
}
\arguments{
\item{object}{item attributes. Can be a \code{\link{data.frame}} or the file path of a .csv file.
The content should at least include an 'ID' column that matches with item IDs (the 'ID' column) of the \code{\linkS4class{item_pool}} object.}

\item{pool}{an \code{\linkS4class{item_pool}} object. Use \code{\link{loadItemPool}} for this.}
}
\value{
\code{\link{loadItemAttrib}} returns an \code{\linkS4class{item_attrib}} object.

\itemize{
  \item{\code{data}} a \code{\link{data.frame}} containing item attributes.
}
}
\description{
\code{\link{loadItemAttrib}} is a data loading function for creating an \code{\linkS4class{item_attrib}} object.
\code{\link{loadItemAttrib}} can read item attributes from a \code{\link{data.frame}} or a .csv file.
}
\examples{
## Read from data.frame:
itempool_science   <- loadItemPool(itempool_science_data)
itemattrib_science <- loadItemAttrib(itemattrib_science_data, itempool_science)

## Read from file: write to tempdir() for illustration and clean afterwards
f <- file.path(tempdir(), "itemattrib_science.csv")
write.csv(itemattrib_science_data, f, row.names = FALSE)
itemattrib_science <- loadItemAttrib(f, itempool_science)
file.remove(f)

}
\seealso{
\code{\link{dataset_science}}, \code{\link{dataset_reading}}, \code{\link{dataset_fatigue}}, \code{\link{dataset_bayes}} for examples.
}
