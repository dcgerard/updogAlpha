context("test llike increases")

test_that("problem_params.RDS no longer problem", {
  skip("for now, because problem_params.RDS not loading correctly.")
  params <- readRDS(file = "problem_params.RDS")

  llike_old <- obj_wrapp_all(parvec = params$parvec, ocounts = params$ocounts, osize = params$osize,
                weight_vec = params$weight_vec, ploidy = params$ploidy, p1geno = params$p1geno,
                p2geno = params$p2geno, allele_freq = params$allele_freq, p1counts = params$p1counts,
                p1size = params$p1size, p1weight = params$p1weight, p2counts = params$p2counts,
                p2size = params$p2size, p2weight = params$p2weight, bound_bias = params$bound_bias,
                bound_od = FALSE, update_bias_val = params$update_bias_val,
                update_seq_error = params$update_seq_error, update_od_param = params$update_od_param,
                seq_error_mean = params$seq_error_mean, seq_error_sd = params$seq_error_sd,
                bias_val_mean = params$bias_val_mean, bias_val_sd = params$bias_val_sd, model = params$model)

  gout <- update_good(parvec = params$parvec, ocounts = params$ocounts, osize = params$osize,
                      weight_vec = params$weight_vec, ploidy = params$ploidy, p1geno = params$p1geno,
                      p2geno = params$p2geno, allele_freq = params$allele_freq, p1counts = params$p1counts,
                      p1size = params$p1size, p1weight = params$p1weight, p2counts = params$p2counts,
                      p2size = params$p2size, p2weight = params$p2weight, bound_bias = params$bound_bias,
                      update_bias_val = params$update_bias_val,
                      update_seq_error = params$update_seq_error, update_od_param = params$update_od_param,
                      seq_error_mean = params$seq_error_mean, seq_error_sd = params$seq_error_sd,
                      bias_val_mean = params$bias_val_mean, bias_val_sd = params$bias_val_sd, model = params$model)

  params_new <- params
  params_new$p1geno <- gout$p1geno
  params_new$p2geno <- gout$p2geno

  s   <- log(gout$bias)
  ell <- log(gout$seq_error / (1 - gout$seq_error))
  r   <- log((1 - gout$od) / gout$od)

  params_new$parvec <- c(s, ell, r)
  params_new$allele_freq <- gout$allele_freq

  llike_new <- obj_wrapp_all(parvec = params_new$parvec, ocounts = params_new$ocounts, osize = params_new$osize,
                      weight_vec = params_new$weight_vec, ploidy = params_new$ploidy, p1geno = params_new$p1geno,
                      p2geno = params_new$p2geno, allele_freq = params_new$allele_freq, p1counts = params_new$p1counts,
                      p1size = params_new$p1size, p1weight = params_new$p1weight, p2counts = params_new$p2counts,
                      p2size = params_new$p2size, p2weight = params_new$p2weight, bound_bias = params_new$bound_bias,
                      update_bias_val = params_new$update_bias_val,
                      update_seq_error = params_new$update_seq_error, update_od_param = params_new$update_od_param,
                      seq_error_mean = params_new$seq_error_mean, seq_error_sd = params_new$seq_error_sd,
                      bias_val_mean = params_new$bias_val_mean, bias_val_sd = params_new$bias_val_sd, model = params_new$model)

  expect_true(llike_new > llike_old)
}
)

test_that("bad_vals works", {
  bad_vals <- readRDS("bad_vals.RDS")
  l1 <- obj_offspring(ocounts = bad_vals$ocounts, osize = bad_vals$osize, ploidy = bad_vals$ploidy,
                prob_geno = bad_vals$prob_geno, bias_val = bad_vals$bias_val,
                seq_error = bad_vals$seq_error, od_param = bad_vals$od_param, outlier = TRUE,
                out_prop = 0, out_mean = bad_vals$out_mean, out_disp = bad_vals$out_disp)
  l2 <- obj_offspring(ocounts = bad_vals$ocounts, osize = bad_vals$osize, ploidy = bad_vals$ploidy,
                      prob_geno = bad_vals$prob_geno, bias_val = bad_vals$bias_val,
                      seq_error = bad_vals$seq_error, od_param = bad_vals$od_param, outlier = FALSE,
                      out_prop = 0.1, out_mean = bad_vals$out_mean, out_disp = bad_vals$out_disp)
  expect_equal(l1, l2)
}
)
