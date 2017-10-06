context("od_param or seq_error set to 0")

test_that("od_param set to 0 works ok", {
  skip("not finished yet")
  dat <- readRDS("example_snp.RDS")
  ploidy <- 6
  uout <- updog_vanilla(ocounts = dat$ocounts, osize = dat$osize, ploidy = ploidy,
                        update_outprop = FALSE, update_bias_val = TRUE, update_od_param = FALSE,
                        update_outdisp = FALSE, update_outmean = FALSE, update_seq_error = TRUE,
                        update_pgeno = TRUE, p1geno = 4, p2geno = 3, od_param = 0, seq_error = 0.01, out_prop = 0)

  err_list <- readRDS("err_list.RDS")

  for (index in 1:10000) {
    grad <- grad_wrapp_all(parvec = err_list$par, ocounts = err_list$ocounts, osize = err_list$osize,
                           weight_vec = err_list$weight_vec, ploidy = err_list$ploidy, p1geno = err_list$p1geno,
                           p2geno = err_list$p2geno, allele_freq = 2, p1counts = err_list$p1counts,
                           p1size = err_list$p1size, p1weight = err_list$p1weight, p2counts = err_list$p2counts,
                           p2size = err_list$p2size, p2weight = err_list$p2weight, bound_bias = err_list$bound_bias,
                           update_bias_val = TRUE, update_seq_error = TRUE,
                           update_od_param = err_list$update_od_param, seq_error_mean = err_list$seq_error_mean,
                           seq_error_sd = err_list$seq_error_sd, bias_val_mean = err_list$bias_val_mean,
                           bias_val_sd = err_list$bias_val_sd, model = err_list$model)
    grad
    err_list$par <- err_list$par + grad / 10000

    obj <- obj_wrapp_all(parvec = err_list$par, ocounts = err_list$ocounts, osize = err_list$osize,
                  weight_vec = err_list$weight_vec, ploidy = err_list$ploidy, p1geno = err_list$p1geno,
                  p2geno = err_list$p2geno, allele_freq = 2, p1counts = err_list$p1counts,
                  p1size = err_list$p1size, p1weight = err_list$p1weight, p2counts = err_list$p2counts,
                  p2size = err_list$p2size, p2weight = err_list$p2weight, bound_bias = err_list$bound_bias,
                  update_bias_val = TRUE, update_seq_error = TRUE,
                  update_od_param = err_list$update_od_param, seq_error_mean = err_list$seq_error_mean,
                  seq_error_sd = err_list$seq_error_sd, bias_val_mean = err_list$bias_val_mean,
                  bias_val_sd = err_list$bias_val_sd, model = err_list$model, bound_od = FALSE)
    cat(obj, "\n")

  }




  input_list <- readRDS(file = "obj_wrapp_all_input.RDS")
  input_list$parvec <- c(0.000, -Inf, 4.595)
  input_list$seq_error_mean <- -Inf
  do.call(what = obj_wrapp_all, args = input_list) ## finite objective




  do.call(what = grad_wrapp_all, args = input_list)

  temp <- grad_offspring_mat(ocounts = 310, osize = 653, ploidy = 6,
                     prob_geno = c(0.00, 0.01, 0.12, 0.37, 0.37, 0.12, 0.01),
                     s = 0, ell = -Inf, r = 4.5)

  dbeta_ds(310, 653, 0, -Inf, 0, exp(-4.5))
  dbeta_ds(310, 653, 0, -Inf, 0.1, exp(-4.5))
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
  dbeta_dh(2, 4, 1, 0.2)
  dbeta_dh(2, 4, 0.01, 0.01)
  dbeta_dh(2, 4, 0.1, 0)

  dbeta_dprop(0, 4, 0.1, 0.1)
  dbeta_dprop(2, 4, 0.01, 0.1)
  dbeta_dprop(4, 4, 0.1, 0)
  dbeta_dprop(4, 4, 0, 0.1)

  museq <- seq(0, 1, length = 100)
  dbbseq <- rep(NA, length = length(museq))
  for (index in 1:length(museq)) {
    dbbseq[index] <- dbetabinom_mu_rho_cpp(4, 4, museq[index], 0)
  }
  plot(museq, dbbseq, type = "l")

}
)
