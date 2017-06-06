context("diagnostics")

test_that("rupdog works", {
  uobj <- readRDS(file = "updog_obj.RDS")
  llike_me <- dupdog(uobj)
  expect_equal(uobj$llike, llike_me)


  plot(uobj, plot_beta = FALSE, use_colorblind = FALSE)
  samp_obj <- rupdog(obj = uobj)
  samp_obj$llike
  plot(samp_obj, plot_beta = FALSE, use_colorblind = FALSE)
}
)

test_that("no NA in warn_dat.RDS", {
  uobj <- readRDS(file = "warn_dat.RDS")
  pl <- plot_geno(ocounts = uobj$input$ocounts, osize = uobj$input$osize,
                  p1counts = uobj$input$p1counts, p1size = uobj$input$p1size,
                  p2counts = uobj$input$p2counts, p2size = uobj$input$p2size,
                  ploidy = uobj$input$ploidy, ogeno = uobj$ogeno,
                  seq_error = uobj$seq_error,
                  bias_val = uobj$bias_val,
                  prob_ok = uobj$prob_ok, maxpostprob = uobj$maxpostprob,
                  p1geno = uobj$p1geno, p2geno = uobj$p2geno,
                  use_colorblind = FALSE)
  print(pl)

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
