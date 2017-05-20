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


test_that("non_llike_increase_dat can be fit", {
  skip("need to fix")
  obj <- readRDS("non_llike_increase_dat.RDS")

  llike_old <- obj_offspring(ocounts = obj$ocounts, osize = obj$osize, ploidy = obj$ploidy,
                             p1geno = obj$p1geno, p2geno = obj$p2geno, bias_val = obj$bias_val,
                             seq_error = obj$seq_error, od_param = obj$od_param, outlier = TRUE,
                             out_prop = obj$out_prop, out_mean = obj$out_mean, out_disp = obj$out_disp)

  obj$weight_vec <- get_out_prop(ocounts = obj$ocounts, osize = obj$osize,
                                 ploidy = obj$ploidy, p1geno = obj$p1geno,
                                 p2geno = obj$p2geno, d = obj$bias_val,
                                 eps = obj$seq_error, tau = obj$od_param,
                                 out_prop = obj$out_prop, out_mean = obj$out_mean,
                                 out_disp = obj$out_disp)

  ## out weights
  obj$out_prop <- mean(obj$weight_vec)

  llike_new <- obj_offspring(ocounts = obj$ocounts, osize = obj$osize, ploidy = obj$ploidy,
                             p1geno = obj$p1geno, p2geno = obj$p2geno, bias_val = obj$bias_val,
                             seq_error = obj$seq_error, od_param = obj$od_param, outlier = TRUE,
                             out_prop = obj$out_prop, out_mean = obj$out_mean, out_disp = obj$out_disp)

  expect_true(llike_old <= llike_new)
  llike_old <- llike_new


  ## good
  obj$s   <- log(obj$bias_val)
  obj$ell <- log(obj$seq_error / (1 - obj$seq_error))
  obj$r   <- log((1 - obj$od_param) / obj$od_param)
  obj$parvec <- c(obj$s, obj$ell, obj$r)
  gout <- update_good(parvec = obj$parvec, ocounts = obj$ocounts, osize = obj$osize,
                      weight_vec = 1 - obj$weight_vec, ploidy = obj$ploidy)
  obj$bias_val  <- gout$bias
  obj$seq_error <- gout$seq_error
  obj$od_param  <- gout$od

  llike_new <- obj_offspring(ocounts = obj$ocounts, osize = obj$osize, ploidy = obj$ploidy,
                             p1geno = obj$p1geno, p2geno = obj$p2geno, bias_val = obj$bias_val,
                             seq_error = obj$seq_error, od_param = obj$od_param, outlier = TRUE,
                             out_prop = obj$out_prop, out_mean = obj$out_mean, out_disp = obj$out_disp)

  expect_true(llike_old <= llike_new)
  llike_old <- llike_new

  ## out dist
  oout <- stats::optim(par = c(obj$out_mean, obj$out_disp),
                       fn = out_obj_wrapp, gr = out_grad_wrapp,
                       ocounts = obj$ocounts, osize = obj$osize,
                       weight_vec = obj$weight_vec, min_disp = 0,
                       method = "BFGS",
                       control = list(fnscale = -1))
  obj$out_mean <- oout$par[1]
  obj$out_disp <- oout$par[2]

  llike_new <- obj_offspring(ocounts = obj$ocounts, osize = obj$osize, ploidy = obj$ploidy,
                             p1geno = obj$p1geno, p2geno = obj$p2geno, bias_val = obj$bias_val,
                             seq_error = obj$seq_error, od_param = obj$od_param, outlier = TRUE,
                             out_prop = obj$out_prop, out_mean = obj$out_mean, out_disp = obj$out_disp)

  expect_true(llike_old <= llike_new)
  llike_old
  llike_new
}
)


test_that("update_* doesn't update parameters set to FALSE", {
  osize <- rbinom(n = 20, size = 70, prob = 0.5)
  ocounts <- rbinom(n = 20, size = osize, prob = 0.5)
  uout <- updog_vanilla(ocounts = ocounts, osize = osize,
                        ploidy = 4, update_bias_val = FALSE, bias_val = 0.9)
  expect_equal(uout$bias_val, 0.9)

  uout <- updog_vanilla(ocounts = ocounts, osize = osize,
                        ploidy = 4, update_seq_error = FALSE, seq_error = 0.3)
  expect_equal(uout$seq_error, 0.3)

  uout <- updog_vanilla(ocounts = ocounts, osize = osize,
                        ploidy = 4, update_od_param = FALSE, od_param = 0.1)
  expect_equal(uout$od_param, 0.1)
}
)
