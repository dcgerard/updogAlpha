context("Test Hypergeometric")

test_that("get_q_array is same as segreg_poly", {
  ploidy <- 6

  qout <- get_q_array(ploidy = ploidy)

  alt_q_array <- array(NA, dim = rep(ploidy + 1, 3))
  for (i in 0:ploidy) {
    for(j in 0:ploidy){
      alt_q_array[i + 1, j + 1, ] <- segreg_poly(m = ploidy, dP = i, dQ = j)
    }
  }

  q53 <- segreg_poly(m = ploidy, dP = 4, dQ = 2)
  q33 <- segreg_poly(m = ploidy, dP = 2, dQ = 2)

  expect_true(all(abs(q53 - qout[5, 3, ]) < 10 ^ -14))
  expect_true(all(abs(q33 - qout[3, 3, ]) < 10 ^ -14))

  expect_true(all(abs(alt_q_array - qout) < 10 ^ -14))
}
)


test_that("dbetabinom is correct", {
  alpha <- 1
  beta <- 3

  mu <- alpha / (alpha + beta)
  rho <- 1 / (1 + alpha + beta)
  a <- dbetabinom(x = 11, size = 23, alpha = alpha, beta = beta)
  if (requireNamespace("VGAM", quietly = TRUE)) {
    a <- dbetabinom(x = 11, size = 23, alpha = alpha, beta = beta)
    b <- VGAM::dbetabinom(x = 11, size = 23, prob = mu, rho = rho)
    expect_equal(a, b)
  }
  c <- dbetabinom_mu_rho(x = 11, size = 23, mu = mu, rho = rho)

  expect_equal(a, c)

  a <- dbetabinom(x = 11, size = 23, alpha = alpha, beta = beta, log = TRUE)
  c <- dbetabinom_mu_rho(x = 11, size = 23, mu = mu, rho = rho, log = TRUE)
  expect_equal(a, c)
}
)
