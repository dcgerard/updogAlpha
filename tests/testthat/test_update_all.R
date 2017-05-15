context("test updog_update_all")

test_that("updog_update_all works ok", {
  dat          <- readRDS(file = "example_snp.RDS")
  ocounts      <- dat$ocounts
  osize        <- dat$osize
  ploidy       <- 6
  tol          <- 10 ^ -4
  maxiter      <- 500
  print_update <- TRUE

}
)


test_that("bad_dat.RDS will not return NaN", {
  obj <- readRDS("bad_dat.RDS")
  val <- obj_offspring_weights_reparam(ocounts = obj$ocounts, osize = obj$osize,
                                       weight_vec = obj$weight_vec,
                                       ploidy = obj$ploidy, p1geno = obj$p1geno, p2geno = obj$p2geno,
                                       s = obj$parvec[1], ell = obj$parvec[2], r = obj$parvec[3])
  expect_equal(val, -Inf)
}
)

test_that("expit works fro very large numbers", {
  obj <- expit(1500)
  expect_equal(obj, 1)

  obj <- expit(-1500)
  expect_equal(obj, 0)
}
)
