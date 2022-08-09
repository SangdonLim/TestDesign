
theta <- seq(-4, 4, .1)
theta_width <- 1
n_theta_grid <- 10
item <- itempool_science[1]
info <- calcKL(item, theta, theta_width, n_theta_grid)

png("KL_demo.png", width = 600, height = 600)

plot(
  theta, info, type = "n", ylim = c(0, max(info)),
  xlab = "Theta",
  ylab = "Information",
  main = "Kullback-Leibler Information"
)
grid()
lines(theta, info, col = "blue")
text(-4, 0,
  sprintf(
    "Integration interval: ±%s\nQuadrature points: %s",
    theta_width, n_theta_grid
  ),
  adj = c(0, 0)
)

text(4, 0,
  sprintf(
    "Model: %s\nParameters: %s",
    item@model,
    paste0(round(item@ipar, 3), collapse = " ")
  ),
  adj = c(1, 0)
)
box(lwd = 1)

dev.off()
