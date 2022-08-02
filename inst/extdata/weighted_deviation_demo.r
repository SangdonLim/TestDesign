library(TestDesign)

cfg <- createShadowTestConfig(
  content_balancing = list(
    method = "WEIGHTED-DEVIATION"
  )
)

set.seed(1)
true_theta <- rnorm(100)
o <- Shadow(cfg, constraints_science, true_theta, seed = 1)

summary(o)
