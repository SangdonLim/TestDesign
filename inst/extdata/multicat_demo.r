library(TestDesign)

prior_mean  <- rep(0, 3)
prior_sigma <- matrix(.3, 3, 3)
diag(prior_sigma) <- 1

cfg_1 <- createShadowTestConfig(
  MIP = list(solver = "GUROBI"),
  item_selection = list(method = "DIRINFO-45"),
  interim_theta = list(
    method = "MAP",
    shrinkage_correction = TRUE,
    prior_par = list(
      prior_mean = prior_mean,
      prior_sigma = prior_sigma
    )
  ),
  final_theta = list(
    method = "MAP",
    shrinkage_correction = TRUE,
    prior_par = list(
      prior_mean = prior_mean,
      prior_sigma = prior_sigma
    )
  ),
  exposure_control = list(
    method = "NONE"
  ),
  theta_grid = as.matrix(expand.grid(
    seq(-4, 4, .25),
    seq(-4, 4, .25),
    seq(-4, 4, .25)
  ))
)

cfg_2 <- createShadowTestConfig(
  MIP = list(solver = "GUROBI"),
  content_balancing = list(method = "NONE"),
  item_selection = list(method = "D-OPTIMAL"),
  interim_theta = list(
    method = "MAP",
    shrinkage_correction = TRUE,
    prior_par = list(
      prior_mean = prior_mean,
      prior_sigma = prior_sigma
    )
  ),
  final_theta = list(
    method = "MAP",
    shrinkage_correction = TRUE,
    prior_par = list(
      prior_mean = prior_mean,
      prior_sigma = prior_sigma
    )
  ),
  exposure_control = list(
    method = "NONE"
  ),
  theta_grid = as.matrix(expand.grid(
    seq(-4, 4, .25),
    seq(-4, 4, .25),
    seq(-4, 4, .25)
  ))
)

cfg_3 <- createShadowTestConfig(
  MIP = list(solver = "GUROBI"),
  item_selection = list(
    method = "MKL",
    KL_width = 1,
    KL_quadrature_unit = 0.25
  ),
  interim_theta = list(
    method = "MAP",
    shrinkage_correction = TRUE,
    prior_par = list(
      prior_mean = prior_mean,
      prior_sigma = prior_sigma
    )
  ),
  final_theta = list(
    method = "MAP",
    shrinkage_correction = TRUE,
    prior_par = list(
      prior_mean = prior_mean,
      prior_sigma = prior_sigma
    )
  ),
  exposure_control = list(
    method = "NONE"
  ),
  theta_grid = as.matrix(expand.grid(
    seq(-4, 4, .25),
    seq(-4, 4, .25),
    seq(-4, 4, .25)
  ))
)

cfg_4 <- createShadowTestConfig(
  MIP = list(solver = "GUROBI"),
  item_selection = list(
    method = "MMI"
  ),
  interim_theta = list(
    method = "MAP",
    shrinkage_correction = TRUE,
    prior_par = list(
      prior_mean = prior_mean,
      prior_sigma = prior_sigma
    )
  ),
  final_theta = list(
    method = "MAP",
    shrinkage_correction = TRUE,
    prior_par = list(
      prior_mean = prior_mean,
      prior_sigma = prior_sigma
    )
  ),
  exposure_control = list(
    method = "NONE"
  ),
  theta_grid = as.matrix(expand.grid(
    seq(-4, 4, .25),
    seq(-4, 4, .25),
    seq(-4, 4, .25)
  ))
)

itempool    <- loadItemPool(itempool_mat_data)
itemattrib  <- loadItemAttrib(itempool_mat_data[, 1, drop = FALSE], itempool)
constraints_data <- constraints_science_data[1, ]
constraints_data$LB <- 50
constraints_data$UB <- 50
constraints <- loadConstraints(constraints_data, itempool, itemattrib)

n_examinees <- 1000
set.seed(1)
true_theta <- matrix(
  mvnfast::rmvn(n_examinees, prior_mean, prior_sigma),
  n_examinees, 3
)
set.seed(1)
resp_data <- simResp(itempool, true_theta)

# using shadowtests (directional information)
set.seed(1)
o1 <- Shadow(cfg_1, constraints, true_theta, resp_data)
# D-optimality
set.seed(1)
o2 <- Shadow(cfg_2, constraints, true_theta, resp_data)
# Kullback-Liebler
set.seed(1)
o3 <- Shadow(cfg_3, constraints, true_theta, resp_data)
# Mutual Information
set.seed(1)
o4 <- Shadow(cfg_4, constraints, true_theta, resp_data)

final_theta_est_1 <- lapply(o1@output, function(x) x@final_theta_est)
final_theta_est_1 <- do.call(rbind, final_theta_est_1)
final_theta_est_2 <- lapply(o2@output, function(x) x@final_theta_est)
final_theta_est_2 <- do.call(rbind, final_theta_est_2)
final_theta_est_3 <- lapply(o3@output, function(x) x@final_theta_est)
final_theta_est_3 <- do.call(rbind, final_theta_est_3)
final_theta_est_4 <- lapply(o4@output, function(x) x@final_theta_est)
final_theta_est_4 <- do.call(rbind, final_theta_est_4)

write.csv(true_theta, "inst/extdata/true_theta.csv", row.names = FALSE)
write.csv(final_theta_est_1, "inst/extdata/final_theta_est_dirinfo.csv", row.names = FALSE)
write.csv(final_theta_est_2, "inst/extdata/final_theta_est_doptim.csv", row.names = FALSE)
write.csv(final_theta_est_3, "inst/extdata/final_theta_est_kl.csv", row.names = FALSE)
write.csv(final_theta_est_4, "inst/extdata/final_theta_est_mi.csv", row.names = FALSE)

true_theta <- read.csv("inst/extdata/true_theta.csv")
final_theta_est_1 <- read.csv("inst/extdata/final_theta_est_dirinfo.csv")
final_theta_est_2 <- read.csv("inst/extdata/final_theta_est_doptim.csv")
final_theta_est_3 <- read.csv("inst/extdata/final_theta_est_kl.csv")
final_theta_est_4 <- read.csv("inst/extdata/final_theta_est_mi.csv")

rmse_1 <- sqrt(apply((final_theta_est_1 - true_theta) ** 2, 2, mean))
rmse_2 <- sqrt(apply((final_theta_est_2 - true_theta) ** 2, 2, mean))
rmse_3 <- sqrt(apply((final_theta_est_3 - true_theta) ** 2, 2, mean))
rmse_4 <- sqrt(apply((final_theta_est_4 - true_theta) ** 2, 2, mean))

rmse_1
rmse_2
rmse_3
rmse_4

theta_range <- c(-4, 4)

for (d in 1:3) {

  png(sprintf("inst/extdata/multicat_dim%s.png", d), width = 400, height = 400)

  plot(
    0, 0,
    type = "n",
    xlim = theta_range, ylim = theta_range,
    xlab = "True theta",
    ylab = "(Estimate - True)",
    main = sprintf("Dimension %s", d)
  )
  abline(0, 0, lty = 3)

  points(true_theta[, d], final_theta_est_1[, d] - true_theta[, d], col = rgb(0.0, 0.0, 1.0, 0.5), pch = 21, cex = 0.5)
  points(true_theta[, d], final_theta_est_2[, d] - true_theta[, d], col = rgb(1.0, 0.0, 0.0, 0.5), pch = 21, cex = 0.5)
  points(true_theta[, d], final_theta_est_3[, d] - true_theta[, d], col = rgb(1.0, 0.5, 0.0, 0.5), pch = 21, cex = 0.5)
  points(true_theta[, d], final_theta_est_4[, d] - true_theta[, d], col = rgb(0.0, 0.5, 1.0, 0.5), pch = 21, cex = 0.5)

  legend(
    "topleft",
    c(
      sprintf("shadow-dirinfo RMSE = %1.3f", rmse_1[d]),
      sprintf("D-optimality RMSE = %1.3f", rmse_2[d]),
      sprintf("shadow-KL RMSE = %1.3f", rmse_3[d]),
      sprintf("shadow-MI RMSE = %1.3f", rmse_4[d])
    ),
    col = c(
      rgb(0.0, 0.0, 1.0),
      rgb(1.0, 0.0, 0.0),
      rgb(1.0, 0.5, 0.0),
      rgb(0.5, 0.5, 1.0)
    ),
    pch = 1,
    pt.cex = 0.5
  )

  dev.off()

}
