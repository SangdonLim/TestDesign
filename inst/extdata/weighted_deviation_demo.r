library(TestDesign)

cfg <- createShadowTestConfig(
  content_balancing = list(
    method = "WEIGHTED-DEVIATION"
  )
)

constraints <- constraints_science
constraints@constraints$WEIGHT <- c(
  rep(20, 7), rep(1, 13), rep(20, 3), rep(1, 13)
)

set.seed(1)
true_theta <- rnorm(100)
o <- Shadow(cfg, constraints, true_theta, seed = 1)

summary(o)
