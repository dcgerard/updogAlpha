context("alpha error")


test_that("alpha_err_dat.RDS data will fit", {
  obj <- readRDS("alpha_err_dat.RDS")


  expect_equal(dbetabinom(x = 10, size = 10, alpha = 1, beta = 0), 1)
  expect_equal(dbetabinom(x = 9, size = 10, alpha = 1, beta = 0), 0)
  expect_equal(dbetabinom(x = 0, size = 10, alpha = 0, beta = 1), 1)
  expect_equal(dbetabinom(x = 1, size = 10, alpha = 0, beta = 1), 0)
  expect_error(dbetabinom(x = 1, size = 10, alpha = 0, beta = 0))

}
)


test_that("input throws error", {
  expect_error(pbias(-1, 1, 1))
  expect_error(pbias( 1, -1, 1))
  expect_error(pbias( 1, 1, -1))

  expect_equal(pbias(c(0.75, 0.5), 0.5, 0.01),
               c(pbias_double(0.75, 0.5, 0.01),
                 pbias_double(0.5, 0.5, 0.01)))
}
)


test_that("dbetabinom_cpp same as dbetabinom", {
  set.seed(1)
  a <- rchisq(1000, 1)
  b <- rchisq(1000, 1)
  expect_equal(dbetabinom(round(a), round(a) + round(b), 1, 1, TRUE),
               dbetabinom_cpp(round(a), round(a) + round(b), 1, 1, TRUE))
  expect_equal(dbetabinom_mu_rho(round(a), round(a) + round(b), 0.1, 0.3, TRUE),
               dbetabinom_mu_rho_cpp(round(a), round(a) + round(b), 0.1, 0.3, TRUE))
}
)

test_that("dbetabinom_mu_rho_cpp same as dbinom when rho = 0", {
  expect_equal(dbetabinom_mu_rho_cpp(x = 4, size = 10,
                                     mu = 0.3, rho = 0),
               dbinom(x = 4, size = 10, prob = 0.3))
}
)

test_that("dhyper_rcpp works", {
  expect_equal(dhyper_cpp(0:4, 4, 5, 4),
               dhyper(0:4, 4, 5, 4))
}
)

test_that("get_q_array_cpp works", {
  expect_true(all(get_q_array(2) - get_q_array_cpp(2) < 10 ^ -6))
  expect_true(all(get_q_array(4) - get_q_array_cpp(4) < 10 ^ -6))
  expect_true(all(get_q_array(6) - get_q_array_cpp(6) < 10 ^ -6))
  expect_true(all(get_q_array(8) - get_q_array_cpp(8) < 10 ^ -6))
}
)


test_that("dbetabinom_mu_rho_cpp_double returns indicators", {
  expect_equal(dbetabinom_mu_rho_cpp_double(0, 2, 0, 1/2, TRUE), 0)
  expect_equal(dbetabinom_mu_rho_cpp_double(1, 2, 0, 1/2, TRUE), -Inf)
  expect_equal(dbetabinom_mu_rho_cpp_double(0, 2, 0, 1/2, FALSE), 1)
  expect_equal(dbetabinom_mu_rho_cpp_double(1, 2, 0, 1/2, FALSE), 0)

  expect_equal(dbetabinom_mu_rho_cpp_double(2, 2, 1, 1/2, TRUE), 0)
  expect_equal(dbetabinom_mu_rho_cpp_double(1, 2, 1, 1/2, TRUE), -Inf)
  expect_equal(dbetabinom_mu_rho_cpp_double(2, 2, 1, 1/2, FALSE), 1)
  expect_equal(dbetabinom_mu_rho_cpp_double(1, 2, 1, 1/2, FALSE), 0)
}
)

test_that("dbetabinom_mu_rho_cpp_double for small od is about the same as binomial", {
  expect_equal(dbetabinom_mu_rho_cpp(x = 2, size = 5, mu = 0.5, rho = 4^-16),
               dbetabinom_mu_rho_cpp(x = 2, size = 5, mu = 0.5, rho = 0), tol = 10^-4)
}
)

test_that("get_pvec works", {
  ploidy <- 4
  d <- 1.2
  seq_error <- 0.01

  p <- 0:ploidy / ploidy
  q <- p * (1 - seq_error) + (1 - p) * seq_error
  xi <- q / (d * (1 - q) + q)
  expect_equal(xi, get_pvec(ploidy = ploidy, bias_val = d, seq_error = seq_error))
}
)

test_that("colSums_cpp works", {
  set.seed(61)
  X <- matrix(stats::rnorm(110), ncol = 11)

  expect_equal(colSums(X), colSums_cpp(X))
}
)
