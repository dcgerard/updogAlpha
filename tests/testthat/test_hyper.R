context("Test Hypergeometric")

test_that("get_q_array is same as segreg_poly", {
  ploidy <- 6

  qout <- get_q_array(ploidy = ploidy)

  q53 <- segreg_poly(m = ploidy, dP = 4, dQ = 2)
  q33 <- segreg_poly(m = ploidy, dP = 2, dQ = 2)

  expect_true(all(abs(q53 - qout[5, 3, ]) < 10 ^ -14))
  expect_true(all(abs(q33 - qout[3, 3, ]) < 10 ^ -14))
}
)
