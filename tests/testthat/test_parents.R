context("test parents")

test_that("obj_parent works ok", {

  pcounts <- 10
  psize <- 12
  pgeno <- 2
  ploidy <- 4
  bias_val <- 0.9
  seq_error <- 0.01
  od_param <- 0.01
  outlier = TRUE
  out_prop = 0.01
  out_mean = 0.5
  out_disp = 1/3
  cppout <- obj_parent(pcounts = pcounts, psize = psize, pgeno = pgeno,
                       ploidy = ploidy, bias_val = bias_val,
                       seq_error = seq_error, od_param = od_param,
                       outlier = outlier, out_prop = out_prop,
                       out_mean = out_mean, out_disp = out_disp)

  mu <- pbias(prob = pgeno / ploidy, bias = bias_val, seq_error = seq_error)
  good_obj <- (1 - out_prop) * dbetabinom_mu_rho_cpp_double(x = pcounts, size = psize, mu = mu, rho = od_param, return_log = FALSE)
  out_obj <- out_prop * dbetabinom_mu_rho_cpp_double(x = pcounts, size = psize, mu = out_mean, rho = out_disp, return_log = FALSE)
  obj <- good_obj + out_obj
  expect_equal(log(obj), cppout)

  s <- log(bias_val)
  ell <- log(seq_error / (1 - seq_error))
  r <- log((1 - od_param) / od_param)

  orig_out <- obj_parent(pcounts = pcounts, psize = psize, pgeno = pgeno,
                         ploidy = ploidy, bias_val = bias_val,
                         seq_error = seq_error, od_param = od_param,
                         outlier = FALSE, out_prop = out_prop,
                         out_mean = out_mean, out_disp = out_disp)

  reparam_out <- obj_parent_reparam(pcounts = pcounts, psize = psize,
                                    pgeno = pgeno, ploidy = ploidy,
                                    s = s, ell = ell, r = r)
  expect_equal(orig_out, reparam_out)
}
)


