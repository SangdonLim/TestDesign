#' @include calc_prob_functions.r
NULL

#' Calculate mutual information
#'
#' \code{\link{calcMI}} is a function for calculating mutual information.
#'
#' @param object an \code{\link{item}} or an \code{\linkS4class{item_pool}} object.
#' @param posterior the posterior density of estimated theta.
#'
#' @return
#' \describe{
#'   \item{\code{\link{item}} object:}{\code{\link{calcMI}} returns a (\emph{nq}, \emph{1}) matrix of information values.}
#'   \item{\code{\linkS4class{item_pool}} object:}{\code{\link{calcMI}} returns a (\emph{nq}, \emph{ni}) matrix of information values.}
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
#' info_item_1 <- calcMI(item_1, seq(-3, 3, 1))
#' info_item_2 <- calcMI(item_2, seq(-3, 3, 1))
#' info_item_3 <- calcMI(item_3, seq(-3, 3, 1))
#' info_item_4 <- calcMI(item_4, seq(-3, 3, 1))
#' info_item_5 <- calcMI(item_5, seq(-3, 3, 1))
#' info_item_6 <- calcMI(item_6, seq(-3, 3, 1))
#' info_pool   <- calcMI(itempool_science, seq(-3, 3, 1))
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
#' @rdname calcMI-methods
calcMI <- function(pos) {
  MI <- numeric(ni)
  available <- items.available
  if (content.balancing) {
    available <- items.available & (content.cat == .GetNextContent())
  }
  for (i in 1:ni) {
    if (available[i]) {
      ncat <- NCAT[i]
      p <- numeric(ncat)
      posterior.k <- matrix(NA, ncat, length(pos))
      for (k in 1:ncat) {
        posterior.k[k,] <- pos * pp[, i, k]
        p[k] <- sum(posterior.k[k, ])
      }
      p <- p / sum(p)
      for (k in 1:ncat) {
        MI[i] <- MI[i] + sum(posterior.k[k, ] * log(posterior.k[k, ] / (pos * p[k])))
      }
    }
  }
  return(MI)
}
