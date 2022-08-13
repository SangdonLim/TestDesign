.CalcMI <- function(pos) {
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
