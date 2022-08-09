#' @include calc_prob_functions.r
NULL

#' Calculate Kullback-Leibler information
#'
#' \code{\link{calcKL}} is a function to calculate Kullback-Leibler information.
#'
#' @param object an \code{\link{item}} or an \code{\linkS4class{item_pool}} object.
#' @param theta theta values to use.
#' @param theta_width the theta width to use to define the range of integration. The range is \code{c(theta - theta_width, theta + theta_width)}.
#' @param quadrature_unit the quadrature spacing to use for integration. Smaller numbers yield more accurate values.
#'
#' @return
#' \describe{
#'   \item{\code{\link{item}} object:}{\code{\link{calcKL}} returns a (\emph{nq}, \emph{1}) matrix of information values.}
#'   \item{\code{\linkS4class{item_pool}} object:}{\code{\link{calcKL}} returns a (\emph{nq}, \emph{ni}) matrix of information values.}
#' }
#' \describe{
#'   \item{\emph{notations}}{\itemize{
#'     \item{\emph{nq} denotes the number of theta values.}
#'     \item{\emph{ni} denotes the number of items in the \code{\linkS4class{item_pool}} object.}
#'   }}
#' }
#'
#' @examples
#' item_1      <- new("item_1PL", difficulty = 0.5)
#' item_2      <- new("item_2PL", slope = 1.0, difficulty = 0.5)
#' item_3      <- new("item_3PL", slope = 1.0, difficulty = 0.5, guessing = 0.2)
#' item_4      <- new("item_PC", threshold = c(-1, 0, 1), ncat = 4)
#' item_5      <- new("item_GPC", slope = 1.2, threshold = c(-0.8, -1.0, 0.5), ncat = 4)
#' item_6      <- new("item_GR", slope = 0.9, category = c(-1, 0, 1), ncat = 4)
#'
#' info_item_1 <- calcKL(item_1, seq(-3, 3, 1), 1, 10)
#' info_item_2 <- calcKL(item_2, seq(-3, 3, 1), 1, 10)
#' info_item_3 <- calcKL(item_3, seq(-3, 3, 1), 1, 10)
#' info_item_4 <- calcKL(item_4, seq(-3, 3, 1), 1, 10)
#' info_item_5 <- calcKL(item_5, seq(-3, 3, 1), 1, 10)
#' info_item_6 <- calcKL(item_6, seq(-3, 3, 1), 1, 10)
#' info_pool   <- calcKL(itempool_science, seq(-3, 3, 1), 1, 10)
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
calcKL <- function(current.theta, d) {
  KL <- numeric(ni)
  interval <- seq(current.theta - d, current.theta + d, length.out = 10)
  available <- items.available
  if (content.balancing) {
    available <- items.available & (content.cat == .GetNextContent())
  }
  for (i in 1:ni) {
    if (available[i]) {
      ncat <- NCAT[i]
      p <- as.vector(TestDesign::calcProb(item.pool@parms[[i]], current.theta))
      p.interval <- TestDesign::calcProb(item.pool@parms[[i]], interval)
      for (k in 1:ncat) {
        KL[i] <- KL[i] + sum(p[k] * log(p[k] / p.interval[, k]))
      }
    }
  }
  return(KL)
}
