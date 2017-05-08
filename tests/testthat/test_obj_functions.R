context("objective functions")


test_that("obj_offspring works", {
  set.seed(8)

  ocounts <- 1:10
  osize <- ocounts + stats::rbinom(10, 6, 0.5)
  ploidy <- 4
  bias_val <- 1
  seq_error <- 0
  od_param <- 0
  outlier <- FALSE
  p1geno <- 1
  p2geno <- 2

  obj_offspring(ocounts = ocounts, osize = osize,
                ploidy = ploidy, p1geno = p1geno,
                p2geno = p2geno, outlier = TRUE)


}
)

test_that("logsumexp works", {
  set.seed(9)
  xx <- matrix(log(abs(rnorm(110))), nrow = 10)
  expect_equal(log(rowSums(exp(xx))), logsumexp(xx))
}
)


test_that("up_bb_obj and up_obj and obj_offspring can be reconciled", {


  ocounts <- rbinom(n = 2, size = 20, prob = 1/2)
  osize    <- rep(20, 2)
  rho <- 1/3
  ploidy <- 4
  r1vec <- rep(1/(ploidy + 1), ploidy + 1)
  r2vec <- r1vec
  pk <- (0:ploidy) / ploidy
  p1geno <- 1
  p2geno <- 2
  qarray <- get_q_array_cpp(ploidy = ploidy)

  ## When no outliers -----------------------------------------------------------
  ## need to add back prior to up_bb_obj
  rout <- up_bb_obj(pival = 1, p1geno = p1geno, p2geno = p2geno, rho = rho, out_mu = 1/2, out_rho = 1/3,
                    ocounts = ocounts, osize = osize, qarray = qarray, r1vec = r1vec, r2vec = r2vec,
                    pk = pk) - log(r1vec[p1geno + 1]) - log(r2vec[p2geno + 1])
  cppout <- obj_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy, p1geno = p1geno,
                          p2geno = p2geno, bias_val = 1, seq_error = 0, od_param = rho,
                          outlier = FALSE, out_prop = 0, out_mean = 1/2, out_disp = 1/3)
  expect_equal(rout, cppout)

  ## When outliers ---------------------------------------------------------------
  rout <- up_bb_obj(pival = 1/2, p1geno = p1geno, p2geno = p2geno, rho = rho, out_mu = 1/2, out_rho = 1/3,
                    ocounts = ocounts, osize = osize, qarray = qarray, r1vec = r1vec, r2vec = r2vec,
                    pk = pk) - log(r1vec[p1geno + 1]) - log(r2vec[p2geno + 1])
  cppout <- obj_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy, p1geno = p1geno,
                          p2geno = p2geno, bias_val = 1, seq_error = 0, od_param = rho,
                          outlier = TRUE, out_prop = 1/2, out_mean = 1/2, out_disp = 1/3)
  expect_equal(rout, cppout)

}
)

test_that("obj_offspring_weights works ok.", {
  ocounts <- rbinom(n = 2, size = 20, prob = 1/2)
  osize    <- rep(20, 2)
  rho <- 1/3
  ploidy <- 4
  r1vec <- rep(1/(ploidy + 1), ploidy + 1)
  r2vec <- r1vec
  pk <- (0:ploidy) / ploidy
  p1geno <- 1
  p2geno <- 2
  bias_val <- 0.9
  seq_error <- 0.01
  weight_vec <- rep(1, 2)

  cppout <- obj_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy, p1geno = p1geno,
                          p2geno = p2geno, bias_val = bias_val, seq_error = seq_error, od_param = rho,
                          outlier = FALSE)
  cppout2 <- obj_offspring_weights(ocounts = ocounts, osize = osize, ploidy = ploidy, p1geno = p1geno,
                           p2geno = p2geno, bias_val = bias_val, seq_error = seq_error, od_param = rho,
                           outlier = FALSE, weight_vec = weight_vec)

  expect_equal(cppout, cppout2)
}
)
