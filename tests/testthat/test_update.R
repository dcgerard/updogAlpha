context("updog")


test_that("bin_post works", {
  set.seed(18)
  ncounts <- c(11, 13)
  ssize   <- c(21, 37)
  prior   <- 6

  postprob <- bin_post(ncounts = ncounts, ssize = ssize, prior = prior)
  expect_equal(sum(postprob), 1)

  ## test zeros
  postprob1 <- bin_post(ncounts = 0, ssize = 2, prior = 2)
  postprob2 <- bin_post(ncounts = 2, ssize = 2, prior = 2)
  temp <- ((1 - seq(0, 2) / 2) ^ 2) / 3
  expect_true(all(abs(temp / sum(temp) - postprob1) < 10 ^ -14))
  expect_true(all(abs(temp / sum(temp) - postprob2[3:1]) < 10 ^ -14))
}
)

test_that("updog works", {
  set.seed(734)
  ocounts  <- c(11, 13, 0, 11)
  osize    <- c(21, 37, 2, 11)
  p1counts <- c(11, 9, 7, 7)
  p1size   <- c(23, 29, 17, 15)
  p2counts <- c(7, 9, 11)
  p2size   <- c(11, 9, 11)
  ploidy   <- 6

  uout <- updog(ocounts = ocounts, osize = osize, p1counts = p1counts,
                p1size = p1size, p2counts = p2counts, p2size = p2size,
                ploidy = ploidy)

  expect_true(all(abs(rowSums(uout) - 1) < 10 ^ -14))

}
)
