context("alpha error")


test_that("alpha_err_dat.RDS data will fit", {
  obj <- readRDS("alpha_err_dat.RDS")


  expect_equal(dbetabinom(x = 10, size = 10, alpha = 1, beta = 0), 1)
  expect_equal(dbetabinom(x = 9, size = 10, alpha = 1, beta = 0), 0)
  expect_equal(dbetabinom(x = 0, size = 10, alpha = 0, beta = 1), 1)
  expect_equal(dbetabinom(x = 1, size = 10, alpha = 0, beta = 1), 0)
  expect_error(dbetabinom(x = 1, size = 10, alpha = 0, beta = 0))

  uout <- updog(ocounts = obj$ocounts, osize = obj$osize, ploidy = obj$ploidy, seq_error = 0,
                update_rho = TRUE, update_pi = FALSE, overdispersion = TRUE, update_geno = FALSE)
}
)
