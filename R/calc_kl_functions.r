#' @include calc_prob_functions.r
NULL

#' Calculate Kullback-Leibler information
#'
#' \code{\link{calcKL}} is a function to calculate Kullback-Leibler information.
#'
#' @param object an \code{\linkS4class{item_pool}} object.
#' @param theta theta values to use.
#' @param theta_width the theta width to use to define the range of integration. The range is \code{c(theta - theta_width, theta + theta_width)}.
#' @param quadrature_unit the quadrature spacing to use for integration. Smaller numbers yield more accurate values.
#' @param use_ellipse if \code{TRUE}, use an ellipse domain for integration. if \code{FALSE}, use a rectangular domain for integration.
#'
#' @return
#' \describe{
#'   \item{\code{\linkS4class{item_pool}} object:}{\code{\link{calcKL}} returns a (\emph{nx}, \emph{ni}) matrix of information values.}
#' }
#' \describe{
#'   \item{\emph{notations}}{\itemize{
#'     \item{\emph{nx} denotes the number of theta values in the 'theta' argument.}
#'     \item{\emph{ni} denotes the number of items in the \code{\linkS4class{item_pool}} object.}
#'   }}
#' }
#'
#' @examples
#' info_pool <- calcKL(itempool_science, seq(-3, 3, 1), 1, 10, FALSE)
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
  def = function(object, theta, theta_width, quadrature_unit, use_ellipse) {
    standardGeneric("calcKL")
  }
)

#' @rdname calcKL-methods
#' @aliases calcKL,item_pool,matrix-method
setMethod(
  f = "calcKL",
  signature = c("item_pool", "numeric"),
  definition = function(object, theta, theta_width, quadrature_unit, use_ellipse) {
    theta <- matrix(theta, , 1)
    return(calcKL(object, theta, theta_width, quadrature_unit, use_ellipse))
  }
)

#' @rdname calcKL-methods
#' @aliases calcKL,item_pool,matrix-method
setMethod(
  f = "calcKL",
  signature = c("item_pool", "matrix"),
  definition = function(object, theta, theta_width, quadrature_unit, use_ellipse) {

    ni <- object@ni
    NCAT <- object@NCAT
    info <- numeric(ni)

    n_theta <- nrow(theta)
    info <- matrix(0, n_theta, ni)

    for (theta_idx in 1:n_theta) {

      this_theta <- theta[theta_idx, , drop = FALSE]

      theta_q <- getThetaQuadrature(this_theta, theta_width, quadrature_unit, use_ellipse)

      for (i in 1:ni) {
        ncat_thisitem <- NCAT[i]
        prob_origin <- calcProb(object@parms[[i]], this_theta)
        prob_area   <- calcProb(object@parms[[i]], theta_q)
        info_thisitem <- 0
        for (k in 1:ncat_thisitem) {
          info_thisitem <- info_thisitem +
          mean(prob_origin[k] * log(prob_origin[k] / prob_area[, k]))
        }
        info[theta_idx, i] <- info[theta_idx, i] + info_thisitem
      }

    }

    return(info)

  }
)
