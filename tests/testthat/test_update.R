context("updog")


test_that("bin_post works", {
  set.seed(18)
  ncounts <- c(11, 13)
  ssize   <- c(21, 37)
  prior   <- 6

  postprob <- bin_post(ncounts = ncounts, ssize = ssize, prior = prior)
  expect_equal(sum(postprob), 1)

  ## test zeros
  postprob1 <- bin_post(ncounts = 0, ssize = 2, prior = 2, seq_error = 0)
  postprob2 <- bin_post(ncounts = 2, ssize = 2, prior = 2, seq_error = 0)
  temp <- ((1 - seq(0, 2) / 2) ^ 2) / 3
  expect_true(all(abs(temp / sum(temp) - postprob1) < 10 ^ -14))
  expect_true(all(abs(temp / sum(temp) - postprob2[3:1]) < 10 ^ -14))


}
)

test_that("updog works", {
  set.seed(734)
  ocounts  <- c(11, 13, 0, 100)
  osize    <- c(21, 37, 2, 101)
  p1counts <- c(11, 9, 7, 7)
  p1size   <- c(23, 29, 17, 15)
  p2counts <- c(7, 9, 11)
  p2size   <- c(11, 9, 11)
  ploidy   <- 6
  seq_error <- 0.01

  uout <- updog(ocounts = ocounts, osize = osize, p1counts = p1counts,
                p1size = p1size, p2counts = p2counts, p2size = p2size,
                ploidy = ploidy, do_eb = FALSE, seq_error = 0.01, overdispersion = FALSE)

  uout2 <- updog(ocounts = ocounts, osize = osize,
                ploidy = ploidy, do_eb = FALSE, seq_error = 0.01, overdispersion = FALSE)

  expect_true(all(abs(colSums(uout$opostprob) - 1) < 10 ^ -14))

  uout$m_opostprob
  uout2$m_opostprob

  # uout <- updog(ocounts = ocounts, osize = osize, p1counts = p1counts,
  #               p1size = p1size, p2counts = p2counts, p2size = p2size,
  #               ploidy = ploidy, do_mcmc = TRUE, iterate = TRUE)
  #
  # plot(uout$opostprob, uout$m_opostprob)
  # abline(0, 1)
  #
  # plot(uout$p1postprob, uout$m_p1postprob)
  # abline(0, 1)
  #
  # plot(uout$p2postprob, uout$m_p2postprob)

}
)


test_that("up_fix increases the up_obj at every iteration", {
    dat <- readRDS("dat.RDS")
    seq_error <- 0.01
    ocounts <- dat[1:(nrow(dat) - 12), 1]
    osize <- rowSums(dat[1:(nrow(dat) - 12), ])
    p1counts <- dat[(nrow(dat) - 12 + 1):(nrow(dat) - 12 + 6), 1]
    p1size <- rowSums(dat[(nrow(dat) - 12 + 1):(nrow(dat) - 12 + 6), ])
    p2counts <- dat[(nrow(dat) - 12 + 7):nrow(dat), 1]
    p2size <- rowSums(dat[(nrow(dat) - 12 + 7):nrow(dat), ])

    ploidy <- 6
    prior_vec <- rep(1/(ploidy + 1), length = ploidy + 1)
    r1vec <- bin_post(ncounts = p1counts, ssize = p1size, prior = prior_vec)
    r2vec <- bin_post(ncounts = p2counts, ssize = p2size, prior = prior_vec)


    pk <- seq(0, ploidy) / ploidy ## the possible probabilities
    pk <- (1 - seq_error) * pk + seq_error * (1 - pk)
    pival <- 0.9
    p1geno <- which.max(r1vec)
    p2geno <- which.max(r2vec)
    alpha <- 0.01
    beta <- 0.01
    dbinommat <-  mapply(FUN = stats::dbinom, ocounts, osize, MoreArgs = list(prob = pk))
    ldbinommat <-  mapply(FUN = stats::dbinom, ocounts, osize, MoreArgs = list(prob = pk, log = TRUE))
    qarray <- get_q_array(ploidy)

    objout <- up_obj(pival = pival, p1geno = p1geno, p2geno = p2geno, alpha = alpha,
                     beta = beta, ocounts = ocounts, osize = osize,
                     dbinommat = dbinommat, qarray = qarray, r1vec = r1vec,
                     r2vec = r2vec)

    itermax <- 100
    objvec <- rep(NA, length = itermax)
    for(iter_index in 1:itermax) {
        objout_old <- objout
        fout <- up_fix(pival = pival, p1geno = p1geno, p2geno = p2geno, alpha = alpha,
                       beta = beta, ocounts = ocounts, osize = osize,
                       ldbinommat = ldbinommat, qarray = qarray, r1vec = r1vec,
                       r2vec = r2vec, update_beta = TRUE, update_pi = TRUE, update_geno = TRUE)
        p1geno <- fout$p1geno
        p2geno <- fout$p2geno
        pival <- fout$pival
        beta <- fout$beta
        alpha <- fout$alpha

        objout <- up_obj(pival = pival, p1geno = p1geno, p2geno = p2geno, alpha = alpha,
                         beta = beta, ocounts = ocounts, osize = osize,
                         dbinommat = dbinommat, qarray = qarray, r1vec = r1vec,
                         r2vec = r2vec)

        objvec[iter_index] <- objout

        expect_true(objout - objout_old > - 10 ^ -12)
    }

    r1vec <- rep(1 / (ploidy + 1), length = ploidy + 1)
    r2vec <- rep(1 / (ploidy + 1), length = ploidy + 1)
    maxout <- updog_maximize(ocounts = ocounts, osize = osize, qarray = qarray,
                   r1vec = r1vec, r2vec = r2vec, pk = pk, pival = 0.99,
                   alpha = 0.1, beta = 0.1,
                   est_fudge = TRUE,
                   tol = 10 ^ -4, itermax = 1000,
                   update_geno = TRUE, update_pi = TRUE, update_beta = TRUE)

    plot_geno(ocounts = ocounts, osize = osize, ploidy = ploidy, p1counts = p1counts,
              p2counts = p2counts, p1size = p1size, p2size = p2size, col = maxout$ogeno, prob_ok = maxout$theta)

    plot(sort(apply(maxout$opostprob, 2, max)))

}
)
