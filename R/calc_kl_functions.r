#' @include calc_prob_functions.r
NULL

#' Calculate Kullback-Leibler information
#'
#' \code{\link{calcKL}} is a function to calculate Kullback-Leibler information.
#'
#' @param object an \code{\linkS4class{item_pool}} object.
#' @param theta the theta value to use. This must be a single value.
#' @param theta_width the theta width to use to define the range of integration. The range is \code{c(theta - theta_width, theta + theta_width)}.
#' @param quadrature_unit the quadrature spacing to use for integration. Smaller numbers yield more accurate values.
#'
#' @return
#' \describe{
#'   \item{\code{\linkS4class{item_pool}} object:}{\code{\link{calcKL}} returns a length-\emph{ni} vector of information values.}
#' }
#' \describe{
#'   \item{\emph{notations}}{\itemize{
#'     \item{\emph{ni} denotes the number of items in the \code{\linkS4class{item_pool}} object.}
#'   }}
#' }
#'
#' @examples
#' info_pool <- calcKL(itempool_science, 0, 1, 10)
#'
#' @template 1pl-ref
#' @template 2pl-ref
#' @template 3pl-ref
#' @template pc-ref
#' @template gpc-ref
#' @template gr-ref
#'
#' @export
#' @docType methods
#' @rdname calcKL-methods
setGeneric(
  name = "calcKL",
  def = function(object, theta, theta_width, quadrature_unit) {
    standardGeneric("calcKL")
  }
)

#' @rdname calcKL-methods
#' @aliases calcKL,item_pool,matrix-method
setMethod(
  f = "calcKL",
  signature = c("item_pool", "matrix"),
  definition = function(object, theta, theta_width, quadrature_unit) {

    ni <- object@ni
    NCAT <- object@NCAT
    info <- numeric(ni)

    theta_q <- seq(theta - theta_width, theta + theta_width, by = quadrature_unit)

    for (i in 1:ni) {
      ncat_thisitem <- NCAT[i]
      prob_origin <- calcProb(object@parms[[i]], theta)
      prob_area   <- calcProb(object@parms[[i]], theta_q)
      info_thisitem <- 0
      for (k in 1:ncat_thisitem) {
        info_thisitem <- info_thisitem +
          sum(prob_origin[k] * log(prob_origin[k] / prob_area[, k]))
      }
      info[i] <- info[i] + info_thisitem
    }

    return(info)

  }

)
