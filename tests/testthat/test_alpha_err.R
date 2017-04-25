context("alpha error")


test_that("alpha_err_dat.RDS data will fit", {
  obj <- readRDS("alpha_err_dat.RDS")

  expect_error(uout <- updog(ocounts = obj$ocounts, osize = obj$osize, ploidy = obj$ploidy))
}
)
