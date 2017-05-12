context("beta-binomial")

test_that("up_bb_fix increases likelihood", {
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
    p1geno <- which.max(r1vec) - 1
    p2geno <- which.max(r2vec) - 1
    out_mu <- 0.5
    out_rho <- 0.8
    rho <- 0.001
    qarray <- get_q_array(ploidy)

    objnew <- up_bb_obj(pival = pival, p1geno = p1geno, p2geno = p2geno,
                        rho = rho, out_mu = out_mu, out_rho = out_rho,
                        ocounts = ocounts, osize = osize, qarray = qarray,
                        r1vec = r1vec, r2vec = r2vec, pk = pk)

    itermax <- 1
    objvec <- rep(NA, length = itermax)
    for (iterindex in 1:itermax) {
        objold <- objnew
        fout <- up_bb_fix(pival = pival, p1geno = p1geno, p2geno = p2geno,
                          rho = rho, out_mu = out_mu, out_rho = out_rho,
                          ocounts = ocounts, osize = osize, qarray = qarray,
                          r1vec = r1vec, r2vec = r2vec, pk = pk, update_geno = FALSE,
                          update_pival = TRUE, update_rho = TRUE, update_outlier = TRUE)
        pival   <- fout$pival
        rho     <- fout$rho
        out_mu  <- fout$out_mu
        out_rho <- fout$out_rho

        objnew <- up_bb_obj(pival = pival, p1geno = p1geno, p2geno = p2geno,
                            rho = rho, out_mu = out_mu, out_rho = out_rho,
                            ocounts = ocounts, osize = osize, qarray = qarray,
                            r1vec = r1vec, r2vec = r2vec, pk = pk)

        objvec[iterindex] <- objnew

        expect_true(objnew - objold >= -10^-8)
    }


    plot(objvec)

}
)


test_that("bb_post works", {

  set.seed(1238)
  ploidy  <- 6
  ncounts <- rbinom(n = 11, size = 10, prob = 0.5)
  ssize   <- rbinom(n = 11, size = 100, prob = 0.5)
  prior   <- rep(1 / (ploidy + 1), length = ploidy + 1)
  seq_error <- 0.01
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)
  rho <- 0.01

  bbout <- bb_post(ncounts = ncounts, ssize = ssize, prior = prior, seq_error = 0.01)


}
)


test_that("dbetabinom_mu_rho_cpp_double works", {
  x <- 4
  n <- 6
  mu <- 0.5
  rho <- 1/3

  expect_equal(dbetabinom_mu_rho_cpp(x, n, mu, rho, FALSE),
               dbetabinom_mu_rho_cpp_double(x, n, mu, rho, FALSE))
}
)


test_that("bbpost_double and bb_simple_post give same results", {
  x <- 4
  n <- 6
  ploidy <- 6
  p1geno <- 2
  p2geno <- 1
  seq_error <- 0.01
  bias_val <- 0.9
  od_param <- 0.01

  bbout <- bbpost_double(x = x, n = n, ploidy = ploidy, p1geno = p1geno, p2geno = p2geno,
                         seq_error = seq_error, bias_val = bias_val, od_param = od_param,
                         outlier = FALSE)
  bbout_r <- bb_simple_post(ncounts = x, ssize = n, ploidy = ploidy, p1geno = p1geno, p2geno = p2geno,
                            seq_error = seq_error, bias_val = bias_val, od_param = od_param)

  bbvec <- bbpost_tot(c(x, x), c(n, n), ploidy = ploidy, p1geno = p1geno, p2geno = p2geno,
                      seq_error = seq_error, bias_val = bias_val, od_param = od_param,
                      outlier = FALSE)

  expect_equal(bbout, bbout_r)
  expect_equal(bbout, bbvec[1, , drop = TRUE])
  expect_equal(bbout, bbvec[2, , drop = TRUE])

  ## Test outlier prob
  temp <- get_out_prop(ocounts = c(x, x), osize = c(n, n), ploidy = ploidy, p1geno = p1geno,
                       p2geno = p2geno, d = bias_val, eps = seq_error, tau = od_param,
                       out_prop = 0.1, out_mean = 0.5,
                       out_disp = 1/3)
  expect_equal(temp[1], temp[2])
}
)
