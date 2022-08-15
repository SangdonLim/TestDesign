#' @include calc_prob_functions.r
NULL

#' Calculate mutual information
#'
#' \code{\link{calcMI}} is a function for calculating mutual information.
#'
#' @param object an \code{\linkS4class{item_pool}} object.
#' @param posterior the posterior density of estimated theta.
#' @param theta_q the theta quadrature points corresponding to the posterior density vector.
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
calcMI <- function(object, posterior, theta_q) {

  ni <- object@ni
  NCAT <- object@NCAT
  info <- numeric(ni)

  pp <- calcProb(object, theta_q)

  for (i in 1:ni) {
    ncat_thisitem <- NCAT[i]
    p <- numeric(ncat_thisitem)
    posterior_k <- matrix(NA, ncat_thisitem, length(posterior))
    for (k in 1:ncat_thisitem) {
      posterior_k[k,] <- posterior * pp[[i]][, k]
      p[k] <- sum(posterior_k[k, ])
    }
    p <- p / sum(p)
    for (k in 1:ncat_thisitem) {
      info[i] <- info[i] + sum(posterior_k[k, ] * log(posterior_k[k, ] / (posterior * p[k])))
    }
  }

  return(info)

}
