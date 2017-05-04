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
