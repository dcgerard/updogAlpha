context("Sequencing Error")

test_that("est_seq works", {
  set.seed(1225)
  seq_error <- 0.01
  allA <- stats::rbinom(n = 100, size = 11:110, prob = 1 - seq_error)
  ncounts <- c(allA, stats::rbinom(n = 100, size = 11:110, prob = 5 / 6))
  ssize <- rep(11:110, times = 2)
  ploidy <- 6
  eps <- NULL

  est_seq_error(ncounts = ncounts, ssize = ssize, ploidy = ploidy)


}
)

test_that("Estimating Sequencing Error Using Posterior", {
  set.seed(701)
  seq_error <- 0.01
  allA <- stats::rbinom(n = 100, size = 11:110, prob = min(jitter(rep(1 - seq_error, 100)), 1))
  ncounts <- c(allA, stats::rbinom(n = 100, size = 11:110, prob = abs(jitter(rep(5 / 6, 100)))))
  ssize <- rep(11:110, times = 2)
  ploidy <- 6
  eps <- NULL

}
)




