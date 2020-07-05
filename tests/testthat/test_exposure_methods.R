test_that("exposure control works", {

  skip_on_cran()
  skip_on_travis()

  set.seed(1)
  true_theta <- runif(100, -3.5, 3.5)
  resp_bayes <- makeTest(itempool_bayes, info_type = "FISHER", true_theta = true_theta)@data

  cfg <- createShadowTestConfig(
    MIP = list(solver = "LPSOLVE"),
    exposure_control = list(method = "ELIGIBILITY")
  )
  set.seed(1)
  solution <- Shadow(cfg, constraints_bayes, true_theta, data = resp_bayes)
  exposure_rate <- solution@exposure_rate[, 2]

  expect_lte(
    max(exposure_rate), 0.35
  )

  cfg <- createShadowTestConfig(
    MIP = list(solver = "LPSOLVE"),
    exposure_control = list(
      method = "BIGM",
      M = 100
    )
  )
  set.seed(1)
  solution <- Shadow(cfg, constraints_bayes, true_theta, data = resp_bayes)
  exposure_rate <- solution@exposure_rate[, 2]

  expect_lte(
    max(exposure_rate), 0.35
  )

  cfg <- createShadowTestConfig(
    MIP = list(solver = "LPSOLVE"),
    exposure_control = list(method = "BIGM-BAYESIAN"),
    interim_theta = list(method = "EB"))
  set.seed(1)
  solution <- Shadow(cfg, constraints_bayes, true_theta, data = resp_bayes)
  exposure_rate <- solution@exposure_rate[, 2]

  expect_lte(
    max(exposure_rate), 0.35
  )

  cfg <- createShadowTestConfig(
    MIP = list(solver = "LPSOLVE"),
    exposure_control = list(method = "BIGM-BAYESIAN"),
    interim_theta = list(method = "FB"))
  set.seed(1)
  solution <- Shadow(cfg, constraints_bayes, true_theta, data = resp_bayes)
  exposure_rate <- solution@exposure_rate[, 2]

  expect_lt(
    max(exposure_rate), 0.35
  )

})