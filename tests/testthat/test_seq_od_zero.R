context("od_param or seq_error set to 0")

test_that("seq_error set to 0 works ok", {
  skip("not finished yet")
  dat <- readRDS("example_snp.RDS")
  ploidy <- 6
  uout <- updog_vanilla(ocounts = dat$ocounts, osize = dat$osize, ploidy = ploidy,
                        update_outprop = FALSE, update_bias_val = FALSE, update_od_param = FALSE,
                        update_outdisp = FALSE, update_outmean = FALSE, update_seq_error = FALSE,
                        update_pgeno = FALSE, p1geno = 4, p2geno = 3, seq_error = 0,
                        seq_error_mean = -Inf)

  input_list <- readRDS(file = "obj_wrapp_all_input.RDS")
  input_list$parvec <- c(0.000, -Inf, 4.595)
  input_list$seq_error_mean <- -Inf
  do.call(what = obj_wrapp_all, args = input_list) ## finite objective




  do.call(what = grad_wrapp_all, args = input_list)

  temp <- grad_offspring_mat(ocounts = 310, osize = 653, ploidy = 6,
                     prob_geno = c(0.00, 0.01, 0.12, 0.37, 0.37, 0.12, 0.01),
                     s = 0, ell = -Inf, r = 4.5)

  dbeta_ds(310, 653, 0, -Inf, 0, exp(-4.5))
  dbeta_ds(310, 653, 0, -Inf, 0.01, exp(-4.5))
  dbeta_ds(310, 653, 0, -Inf, 1, exp(-4.5))
  dbeta_dl(310, 653, 1, -Inf, 0, exp(-4.5))
  dbeta_dl(310, 653, 1, -Inf, 0.01, exp(-4.5))
  dbeta_dl(310, 653, 1, -Inf, 1, exp(-4.5))
  dbeta_dr_ell(310, 653, 1, -Inf, 0, exp(-4.5))
  dbeta_dr_ell(310, 653, 1, -Inf, 0.1, exp(-4.5))
  dbeta_dr_ell(310, 653, 1, -Inf, 1, exp(-4.5))

  x <- 310
  n <- 653
  s <- 0
  ell <- -Inf
  r <- -4.5
  p <- 0
  h <- exp(r)
  d <- exp(s)
  tau <- 1 / (h + 1)
  eps <- expit(ell)
  xi <- pbias_double(p, d, eps)
  f <- (1 - p) * eps + p * (1 - eps)
  dbdxi <- dbeta_dprop(x, n, xi, tau)
  dbdxi
  dxidd <- dxi_dd(d, f)
  dxidd
  ddds <- dd_ds(s)
  ddds
  dxidf <- dxi_df(d, f)
  dxidf
  dfdeps <- df_deps(eps, p)
  dfdeps
  depsdl <- deps_dell(ell)
  depsdl

  dbeta_dprop(x, n, xi, tau)
  dbeta_dr(x, n, xi, r) ## bad because of dbdh

  dbdh <- dbeta_dh(x, n, xi, h) ## bad
  dbdh

  dhdr <- dh_dr(r)
  dhdr

  ### important ones --------------------------------------------------
  dbeta_dh(2, 4, 0, 0.2)
  dbeta_dh(2, 4, 0.01, 0.01)
  dbeta_dh(2, 4, 0.01, 0)

  dbeta_dprop(2, 4, 0, 0.1)
  dbeta_dprop(2, 4, 0.01, 0.1)
  dbeta_dprop(2, 4, 0.01, 0)

}
)
