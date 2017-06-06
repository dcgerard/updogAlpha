context("diagnostics")

test_that("rupdog works", {
  uobj <- readRDS(file = "updog_obj.RDS")
  ## plot(uobj, plot_beta = FALSE, use_colorblind = FALSE)
  samp_obj <- rupdog(obj = uobj)
  ## plot(samp_obj, plot_beta = FALSE, use_colorblind = FALSE)
}
)


test_that("rbetabinom_mu_rho works", {
  set.seed(10)
  n <- 100000
  mu   <- c(0, 1, 0, 1, 0, 1, 0.5, 0.5, 0.5)
  rho  <- c(0, 0, 1, 1, 0.5, 0.5, 0, 1, 0.5)
  size <- rep(n, length = length(mu))

  sout <- rbetabinom_mu_rho(mu = mu, rho = rho, size = size)

  expect_equal(sout[1:6], c(0, n, 0, n, 0, n))
  expect_equal(sout[7] / n, 0.5, tol = 2 * sqrt(0.5 ^ 2 / n))

}
)
