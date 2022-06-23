# Multidimensional MAP estimation demo

library(TestDesign)
library(mvnfast)

# M2PL -------------------------------------------------------------------------

plot(
  0, 0, type = "n", xlim = c(-3, 3), ylim = c(-3, 3),
  main = "M2PL",
  xlab = "Dimension 1",
  ylab = "Dimension 2"
)

for (t1 in seq(-3, 3, .5)) {
for (t2 in seq(-3, 3, .5)) {

  set.seed(1)
  ni <- 10000
  ipar <- matrix(NA, ni, 4)
  ipar[, 1:2] <- abs(rmvn(ni, rep(1, 2), diag(1, 2)))
  ipar[, 3]   <- rmvn(ni, 0, 1)
  ipar[, 4]   <- 999

  true_theta <- matrix(c(t1, t2), 1, 2)

  response <- rep(NA, ni)

  set.seed(1)
  for (i in 1:ni) {
    p <- p_m_2pl(true_theta, ipar[i, 1:2], ipar[i, 3])
    r <- runif(1)
    if (p > r) response[i] <- 1 else response[i] <- 0
  }

  theta_est_fisher <- estimate_theta_map(
    ipar,
    item_model = matrix(102, 1, ni),
    nd = 2,
    n_pars = matrix(3, 1, ni),
    response = matrix(response, 1),
    start_theta = matrix(c(0, 0), 1, 2),
    sigma = diag(1, 2),
    use_Fisher = TRUE
  )
  theta_est_hessian <- estimate_theta_map(
    ipar,
    item_model = matrix(102, 1, ni),
    nd = 2,
    n_pars = matrix(3, 1, ni),
    response = matrix(response, 1),
    start_theta = matrix(c(0, 0), 1, 2),
    sigma = diag(1, 2),
    use_Fisher = FALSE
  )

  points(t1, t2, col = "blue")
  points(
    theta_est_fisher$theta[, 1], theta_est_fisher$theta[, 2],
    col = "red", cex = 1
  )
  points(
    theta_est_hessian$theta[, 1], theta_est_hessian$theta[, 2],
    pch = 21, col = "magenta", bg = "magenta", cex = 0.25
  )

}}

# M3PL -------------------------------------------------------------------------

plot(
  0, 0, type = "n", xlim = c(-3, 3), ylim = c(-3, 3),
  main = "M3PL",
  xlab = "Dimension 1",
  ylab = "Dimension 2"
)

for (t1 in seq(-3, 3, .5)) {
for (t2 in seq(-3, 3, .5)) {

  set.seed(1)
  ni <- 10000
  ipar <- matrix(NA, ni, 4)
  ipar[, 1:2] <- abs(rmvn(ni, rep(1, 2), diag(1, 2)))
  ipar[, 3]   <- rmvn(ni, 0, 1)
  ipar[, 4]   <- runif(ni, 0, 0.2)

  true_theta <- matrix(c(t1, t2), 1, 2)

  response <- rep(NA, ni)

  set.seed(1)
  for (i in 1:ni) {
    p <- p_m_3pl(true_theta, ipar[i, 1:2], ipar[i, 3], ipar[i, 4])
    r <- runif(1)
    if (p > r) response[i] <- 1 else response[i] <- 0
  }

  theta_est_fisher <- estimate_theta_map(
    ipar,
    item_model = matrix(103, 1, ni),
    nd = 2,
    n_pars = matrix(4, 1, ni),
    response = matrix(response, 1),
    start_theta = matrix(c(0, 0), 1, 2),
    sigma = diag(1, 2),
    use_Fisher = TRUE
  )
  theta_est_hessian <- estimate_theta_map(
    ipar,
    item_model = matrix(103, 1, ni),
    nd = 2,
    n_pars = matrix(4, 1, ni),
    response = matrix(response, 1),
    start_theta = matrix(c(0, 0), 1, 2),
    sigma = diag(1, 2),
    use_Fisher = FALSE
  )

  points(t1, t2, col = "blue")
  points(
    theta_est_fisher$theta[, 1], theta_est_fisher$theta[, 2],
    col = "red", cex = 1
  )
  points(
    theta_est_hessian$theta[, 1], theta_est_hessian$theta[, 2],
    pch = 21, col = "magenta", bg = "magenta", cex = 0.25
  )

}}

# MGPC -------------------------------------------------------------------------

plot(
  0, 0, type = "n", xlim = c(-3, 3), ylim = c(-3, 3),
  main = "MGPC",
  xlab = "Dimension 1",
  ylab = "Dimension 2"
)

for (t1 in seq(-3, 3, .5)) {
for (t2 in seq(-3, 3, .5)) {

  set.seed(1)
  ni <- 10000
  ipar <- matrix(NA, ni, 4)
  ipar[, 1:2] <- abs(rmvn(ni, rep(1, 2), diag(1, 2)))
  ipar[, 3]   <- rmvn(ni, 0, 1)
  ipar[, 4]   <- rmvn(ni, 0, 1)

  true_theta <- matrix(c(t1, t2), 1, 2)

  response <- rep(NA, ni)

  set.seed(1)
  for (i in 1:ni) {
    p <- p_m_gpc(true_theta, ipar[i, 1:2], ipar[i, 3:4])
    p <- cumsum(p)
    r <- runif(1)
    response[i] <- sum(r > p)
  }

  theta_est_fisher <- estimate_theta_map(
    ipar,
    item_model = matrix(105, 1, ni),
    nd = 2,
    n_pars = matrix(4, 1, ni),
    response = matrix(response, 1),
    start_theta = matrix(c(0, 0), 1, 2),
    sigma = diag(1, 2),
    use_Fisher = TRUE
  )
  theta_est_hessian <- estimate_theta_map(
    ipar,
    item_model = matrix(105, 1, ni),
    nd = 2,
    n_pars = matrix(4, 1, ni),
    response = matrix(response, 1),
    start_theta = matrix(c(0, 0), 1, 2),
    sigma = diag(1, 2),
    use_Fisher = FALSE
  )

  points(t1, t2, col = "blue")
  points(
    theta_est_fisher$theta[, 1], theta_est_fisher$theta[, 2],
    col = "red", cex = 1
  )
  points(
    theta_est_hessian$theta[, 1], theta_est_hessian$theta[, 2],
    pch = 21, col = "magenta", bg = "magenta", cex = 0.25
  )

}}

# MGR --------------------------------------------------------------------------

plot(
  0, 0, type = "n", xlim = c(-3, 3), ylim = c(-3, 3),
  main = "MGR",
  xlab = "Dimension 1",
  ylab = "Dimension 2"
)

for (t1 in seq(-3, 3, .5)) {
for (t2 in seq(-3, 3, .5)) {

  set.seed(1)
  ni <- 10000
  ipar <- matrix(NA, ni, 4)
  ipar[, 1:2] <- abs(rmvn(ni, rep(1, 2), diag(1, 2)))
  ipar[, 3]   <- rmvn(ni, -1, 1)
  ipar[, 4]   <- rmvn(ni, 1, 1)
  dd <- apply(ipar[, 3:4], 1, sort)
  ipar[, 4:3] <- t(dd)

  true_theta <- matrix(c(t1, t2), 1, 2)

  response <- rep(NA, ni)

  set.seed(1)
  for (i in 1:ni) {
    p <- p_m_gr(true_theta, ipar[i, 1:2], ipar[i, 3:4])
    p <- cumsum(p)
    r <- runif(1)
    response[i] <- sum(r > p)
  }

  theta_est_fisher <- estimate_theta_map(
    ipar,
    item_model = matrix(106, 1, ni),
    nd = 2,
    n_pars = matrix(4, 1, ni),
    response = matrix(response, 1),
    start_theta = matrix(c(0, 0), 1, 2),
    sigma = diag(1, 2),
    use_Fisher = TRUE
  )
  theta_est_hessian <- estimate_theta_map(
    ipar,
    item_model = matrix(106, 1, ni),
    nd = 2,
    n_pars = matrix(4, 1, ni),
    response = matrix(response, 1),
    start_theta = matrix(c(0, 0), 1, 2),
    sigma = diag(1, 2),
    use_Fisher = FALSE
  )

  points(t1, t2, col = "blue")
  points(
    theta_est_fisher$theta[, 1], theta_est_fisher$theta[, 2],
    col = "red", cex = 1
  )
  points(
    theta_est_hessian$theta[, 1], theta_est_hessian$theta[, 2],
    pch = 21, col = "magenta", bg = "magenta", cex = 0.25
  )

}}
