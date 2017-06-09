context("test hw")


test_that("get_prob_geno works", {
  expect_equal(get_prob_geno(4, "f1", 1, 1, 0.5), c(1/4, 1/2, 1/4, 0, 0))
  expect_equal(get_prob_geno(4, "hw", 1, 1, 0.5), stats::dbinom(0:4, 4, 1/2))
  expect_equal(get_prob_geno(4, "uniform", 1, 1, 0.5), rep(1/5, 5))
  expect_error(get_prob_geno(4, "hello", 1, 1, 0.5))
}
)
