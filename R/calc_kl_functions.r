#' @include calc_prob_functions.r
NULL

.CalcKL <- function(current.theta, d) {
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
