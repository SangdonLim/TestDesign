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
    p <- posterior %*% prob_grid[[i]]
    p <- p / sum(p)
    p <- p[1, ]
    posterior_k <- t(posterior * prob_grid[[i]])
    info_k <- diag(posterior_k %*% log(t(t(prob_grid[[i]]) / p)))
    info[i] <- sum(info_k)
  }

  return(info)

}
