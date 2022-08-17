#' @include calc_prob_functions.r
NULL

#' Calculate mutual information
#'
#' \code{\link{calcMI}} is a function for calculating mutual information.
#'
#' @param prob_grid the response probability of a \code{\linkS4class{item_pool}} object over quadrature points.
#' @param posterior the posterior density of estimated theta. Must be from the same quadrature points used for \code{prob_grid}.
#'
#' @return
#' \describe{
#'   \item{\code{\linkS4class{item_pool}} object:}{\code{\link{calcMI}} returns a (\emph{nq}, \emph{ni}) matrix of information values.}
#' }
#' \describe{
#'   \item{\emph{notations}}{\itemize{
#'     \item{\emph{ni} denotes the number of items in the \code{\linkS4class{item_pool}} object.}
#'   }}
#' }
#'
#' @examples
#' info_pool <- calcMI(itempool_science, seq(-3, 3, 1))
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
calcMI <- function(prob_grid, posterior) {

  ni <- length(prob_grid)
  info <- numeric(ni)

  for (i in 1:ni) {
    ncat_thisitem <- ncol(prob_grid[[i]])
    p <- numeric(ncat_thisitem)
    posterior_k <- matrix(NA, ncat_thisitem, length(posterior))
    for (k in 1:ncat_thisitem) {
      posterior_k[k,] <- posterior * prob_grid[[i]][, k]
      p[k] <- sum(posterior_k[k, ])
    }
    p <- p / sum(p)
    for (k in 1:ncat_thisitem) {
      info[i] <- info[i] + sum(posterior_k[k, ] * log(posterior_k[k, ] / (posterior * p[k])))
    }
  }

  return(info)

}
