context("gradients")

test_that("dbeta_dprop works", {
  x <- 4
  n <- 6
  tau <- 1/3
  xi <- 1/4
  h <- (1 - tau) / tau

  ## compare against R version
  comp1 <- choose(n, x) * gamma(n - x + (1 - xi) * h) * gamma(h) /
    (gamma(n + h) * gamma(xi * h) * gamma((1 - xi) * h)) * h * digamma(x + xi * h) * gamma(x + xi * h)

  comp2 <- choose(n, x) * gamma(x + xi * h) * gamma(h) /
    (gamma(n + h) * gamma(xi * h) * gamma((1 - xi) * h)) * h * digamma(n - x + (1 - xi) * h) * gamma(n - x + (1 - xi) * h)

  comp3 <- choose(n, x) * gamma(x + xi * h) * gamma(n - x + (1 - xi) * h) * gamma(h) /
    (gamma(n + h) * gamma((1 - xi) * h)) * digamma(xi * h) * h / gamma(xi * h)

  comp4 <- choose(n, x) * gamma(x + xi * h) * gamma(n - x + (1 - xi) * h) * gamma(h) /
    (gamma(n + h) * gamma(xi * h)) * digamma((1 - xi) * h) * h / gamma((1 - xi) * h)

  rderiv <- comp1 - comp2 - comp3 + comp4

  cderiv <- dbeta_dprop(x = x, n = n, xi = xi, tau = tau)
  expect_equal(rderiv, cderiv)


  ## Compare against numeric derivative
  myenv <- new.env()
  assign("tau", tau, envir = myenv)
  assign("x", x, envir = myenv)
  assign("n", n, envir = myenv)
  assign("xi", xi, envir = myenv)
  nout <- stats::numericDeriv(quote(dbetabinom_mu_rho_cpp(x, n, xi, tau, FALSE)), "xi", myenv)
  expect_equal(attr(nout, "gradient")[1, 1], cderiv)

}
)





test_that("dbeta_dh works", {
  x  <- 4
  n  <- 6
  h  <- 2
  p <- 0.5
  d <- 0.9
  eps <- 0.2
  xi <- pbias_double(p, d, eps)
  ell <- log(eps / (1 - eps))

  ## R version -----------------------------------------------
  dense <- dbetabinom_mu_rho_cpp(x, n, xi, 1 / (h + 1),
                                         return_log = FALSE)

  comp1 <-  xi * digamma(x + xi * h)
  comp2 <- (1 - xi) * digamma(n - x + (1 - xi) * h)
  comp3 <- digamma(h)
  comp4 <- digamma(n + h)
  comp5 <- xi * digamma(xi * h)
  comp6 <- (1 - xi) * digamma((1 - xi) * h)

  rderiv <- (comp1 + comp2 + comp3 - comp4 - comp5 - comp6) * dense
  cderiv <- dbeta_dh(x, n, xi, h)

  expect_equal(cderiv, rderiv)

  ## Numerical version ---------------------------------------

  db_wrapper <- function(x, n, xi, h) {
    tau <- 1 / (h + 1)
    dbetabinom_mu_rho_cpp(x, n, xi, tau, return_log = FALSE)
  }

  myenv <- new.env()
  assign("h", h, envir = myenv)
  assign("x", x, envir = myenv)
  assign("n", n, envir = myenv)
  assign("xi", xi, envir = myenv)
  nout <- stats::numericDeriv(quote(db_wrapper(x, n, xi, h)), "h", myenv)
  expect_equal(cderiv, attr(nout, "gradient")[1, 1])
}
)


test_that("dbeta_dr works", {
  x   <- 4
  n   <- 6
  h   <- 2
  r   <- log(h)
  p   <- 0.5
  d   <- 0.9
  eps <- 0.2
  xi  <- pbias_double(p, d, eps)
  ell <- log(eps / (1 - eps))

  ## Numerical version ---------------------------------------

  db_wrapper <- function(x, n, xi, r) {
    tau <- 1 / (exp(r) + 1)
    dbetabinom_mu_rho_cpp(x, n, xi, tau, return_log = FALSE)
  }

  myenv <- new.env()
  assign("r", r, envir = myenv)
  assign("x", x, envir = myenv)
  assign("n", n, envir = myenv)
  assign("xi", xi, envir = myenv)
  nout <- stats::numericDeriv(quote(db_wrapper(x, n, xi, r)), "r", myenv)
  cderiv <- dbeta_dr(x, n, xi, r)
  expect_equal(cderiv, attr(nout, "gradient")[1, 1])

  cderiv2 <- dbeta_dr_ell(x = x, n = n, d = d, ell = ell, p = p, r = r)
  expect_equal(cderiv2, cderiv)

}
)

test_that("dxi_df works ok", {

  fwrapper <- function(d, f) {
    return(f / (d * (1 - f) + f));
  }

  d <- 1/3
  f <- 1/6

  myenv <- new.env()
  assign("d", d, envir = myenv)
  assign("f", f, envir = myenv)
  nout <- stats::numericDeriv(quote(fwrapper(d, f)), "f", myenv)
  cderiv <- dxi_df(d, f)

  expect_equal(attr(nout, "gradient")[1, 1], cderiv)
}
)

test_that("dxi_dd works ok", {

  fwrapper <- function(d, f) {
    return(f / (d * (1 - f) + f));
  }

  d <- 1/3
  f <- 1/6

  myenv <- new.env()
  assign("d", d, envir = myenv)
  assign("f", f, envir = myenv)
  nout <- stats::numericDeriv(quote(fwrapper(d, f)), "d", myenv)
  cderiv <- dxi_dd(d, f)

  expect_equal(attr(nout, "gradient")[1, 1], cderiv)
}
)


test_that("df_deps works ok", {

  fwrapper <- function(eps, p) {
    return(p * (1 - eps) + (1 - p) * eps);
  }

  p <- 1/5
  eps <- 1/11

  myenv <- new.env()
  assign("p", p, envir = myenv)
  assign("eps", eps, envir = myenv)
  nout <- stats::numericDeriv(quote(fwrapper(eps, p)), "eps", myenv)
  cderiv <- df_deps(eps = eps, p = p)

  expect_equal(attr(nout, "gradient")[1, 1], cderiv)
}
)

test_that("deps_dell works ok", {

  fwrapper <- function(ell) {
    return(exp(ell) / (1 + exp(ell)))
  }

  ell <- 6

  myenv <- new.env()
  assign("ell", ell, envir = myenv)
  nout <- stats::numericDeriv(quote(fwrapper(ell)), "ell", myenv)
  cderiv <- deps_dell(ell)

  expect_equal(attr(nout, "gradient")[1, 1], cderiv)
}
)


test_that("dbeta_dl works ok", {
  set.seed(114)
  x <- 4
  n <- 6
  d <- 3/2
  ell <- 1
  p <- 1/3
  tau <- 1 / 3
  h <- (1 - tau) / tau

  beta_wrap <- function(x, n, d, ell, p, tau) {
    eps <- exp(ell) / (1 + exp(ell))
    xi <- pbias_double(prob = p, bias = d, seq_error = eps)
    dbetabinom_mu_rho_cpp(x = x, size = n, mu = xi, rho = tau, return_log = FALSE)
  }

  myenv <- new.env()
  assign("x", x, envir = myenv)
  assign("n", n, envir = myenv)
  assign("d", d, envir = myenv)
  assign("ell", ell, envir = myenv)
  assign("p", p, envir = myenv)
  assign("tau", tau, envir = myenv)
  nout <- stats::numericDeriv(quote(beta_wrap(x, n, d, ell, p, tau)), "ell", myenv)
  cderiv <- dbeta_dl(x, n, d, ell, p, h)
  expect_equal(attr(nout, "gradient")[1, 1], cderiv)
}
)

test_that("dbeta_ds works ok", {
  set.seed(534)
  x <- 4
  n <- 6
  d <- 3/2
  s <- log(d)
  ell <- 1
  p <- 1/3
  tau <- 1 / 3
  h <- (1 - tau) / tau

  beta_wrap <- function(x, n, s, ell, p, tau) {
    eps <- exp(ell) / (1 + exp(ell))
    d = exp(s)
    xi <- pbias_double(prob = p, bias = d, seq_error = eps)
    dbetabinom_mu_rho_cpp(x = x, size = n, mu = xi, rho = tau, return_log = FALSE)
  }

  myenv <- new.env()
  assign("x", x, envir = myenv)
  assign("n", n, envir = myenv)
  assign("s", s, envir = myenv)
  assign("ell", ell, envir = myenv)
  assign("p", p, envir = myenv)
  assign("tau", tau, envir = myenv)
  nout <- stats::numericDeriv(quote(beta_wrap(x, n, s, ell, p, tau)), "s", myenv)
  cderiv <- dbeta_ds(x, n, s, ell, p, h)
  expect_equal(attr(nout, "gradient")[1, 1], cderiv, tol = 10 ^ -5)
}
)

test_that("the grad_offspring_mat works", {
  set.seed(9)
  osize   <- stats::rbinom(n = 2, size = 70, prob = 0.4)
  ocounts <- stats::rbinom(n = 2, size = osize, prob = 0.4)
  ploidy <- 6
  p1geno <- 2
  p2geno <- 2
  d <- 2
  s <- log(d)
  eps <- 0.01
  ell <- log(eps / (1 - eps))
  tau <- 0.1
  h <- (1 - tau) / (tau)
  r = log(h)

  prob_geno <- get_prob_geno(ploidy = ploidy, model = "f1", p1geno = p1geno, p2geno = p2geno, allele_freq = 0.5)

  gdmat <- grad_offspring_mat(ocounts = ocounts, osize = osize, ploidy = ploidy, prob_geno,
                              s = s, ell = ell, r = r)

  pvec <- 0:ploidy / ploidy
  qout <- get_q_array(ploidy)
  tempsum <- 0
  for (index in 1:length(pvec)) {
    tempsum <- tempsum + dbeta_ds(x = ocounts[1], n = osize[1], s = s, ell = ell, p = pvec[index], h = h) *
      qout[p1geno + 1, p2geno + 1,index]
  }

  ldenom_vec = obj_offspring_vec(ocounts, osize,
                                 ploidy, prob_geno = prob_geno,
                                 d, eps, tau, FALSE, 0, 1.0 / 2.0, 1.0 / 3.0)

  expect_equal(tempsum * exp(-1 * ldenom_vec[1]), gdmat[1,1])
}
)


test_that("the grad_offspring works", {
  set.seed(9)
  osize   <- stats::rbinom(n = 1, size = 70, prob = 0.4)
  ocounts <- stats::rbinom(n = 1, size = osize, prob = 0.4)
  ploidy <- 4
  p1geno <- 2
  p2geno <- 1
  d <- 2
  s <- log(d)
  eps <- 0.01
  ell <- log(eps / (1 - eps))
  tau <- 0.1
  h <- (1 - tau) / (tau)
  r <- log(h)

  prob_geno <- get_prob_geno(ploidy = ploidy, model = "f1", p1geno = p1geno, p2geno = p2geno, allele_freq = 0.5)

  ## Numerical implementation ---------------------------------------
  tempfunc <- function(s, ell, r) {
    obj_offspring_reparam(ocounts = ocounts, osize = osize, ploidy = ploidy, prob_geno = prob_geno,
                          s = s, ell = ell, r = r)
  }

  myenv <- new.env()
  assign("s", s, envir = myenv)
  assign("ell", ell, envir = myenv)
  assign("r", r, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(s, ell, r)), c("s", "ell", "r"), myenv)

  ## Rcpp version ---------------------------------------------------
  cderiv <- grad_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy,
                           prob_geno = prob_geno, s = s, ell = ell, r = r)

  expect_equal(c(attr(nout, "gradient")), c(cderiv), tol = 10 ^ -5)
}
)

test_that("the grad_offspring_weights works", {
  set.seed(9)
  nsamp <- 4
  osize   <- stats::rbinom(n = nsamp, size = 70, prob = 0.4)
  ocounts <- stats::rbinom(n = nsamp, size = osize, prob = 0.4)
  ploidy <- 4
  p1geno <- 2
  p2geno <- 1
  d <- 2
  s <- log(d)
  eps <- 0.01
  ell <- log(eps / (1 - eps))
  tau <- 0.1
  h <- (1 - tau) / (tau)
  r <- log(h)
  weight_vec <- stats::runif(nsamp)

  prob_geno <- get_prob_geno(ploidy = ploidy, model = "f1", p1geno = p1geno, p2geno = p2geno, allele_freq = 0.5)

  ## Numerical implementation ---------------------------------------
  tempfunc <- function(s, ell, r) {
    obj_offspring_weights_reparam(ocounts = ocounts, osize = osize,
                                  weight_vec = weight_vec,
                                  ploidy = ploidy,
                                  prob_geno = prob_geno,
                                  s = s, ell = ell, r = r)
  }

  myenv <- new.env()
  assign("s", s, envir = myenv)
  assign("ell", ell, envir = myenv)
  assign("r", r, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(s, ell, r)), c("s", "ell", "r"), myenv)

  ## Rcpp version ---------------------------------------------------
  cderiv <- grad_offspring_weights(ocounts = ocounts, osize = osize,
                                   weight_vec = weight_vec,
                                   ploidy = ploidy,
                                   prob_geno = prob_geno,
                                   s = s, ell = ell, r = r)

  expect_equal(c(attr(nout, "gradient")), c(cderiv), tol = 10 ^ -6)
}
)


######
## Using original parameterization
######

test_that("dbeta_deps works ok", {
  set.seed(9)
  x <- 10
  n <- 12
  d <- 0.1
  eps <- 0.1
  p <- 0.25
  tau <- 0.1

  tempfunc <- function(eps) {
    xi <- pbias(prob = p, bias = d, seq_error = eps)
    dbetabinom_mu_rho_cpp(x = x, size = n, mu = xi, rho = tau, return_log = FALSE)
  }

  myenv <- new.env()
  assign("eps", eps, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(eps)), c("eps"), myenv)
  cderiv <- dbeta_deps(x = x, n = n, d = d, eps = eps, p = p, tau = tau)

  expect_equal(attr(nout, "gradient")[1, 1], cderiv, tol = 10 ^ -6)
}
)


test_that("dbeta_dd works ok", {
  set.seed(9)
  x <- 10
  n <- 12
  d <- 0.1
  eps <- 0.1
  p <- 0.25
  tau <- 0.1

  tempfunc <- function(d) {
    xi <- pbias(prob = p, bias = d, seq_error = eps)
    dbetabinom_mu_rho_cpp(x = x, size = n, mu = xi, rho = tau, return_log = FALSE)
  }

  myenv <- new.env()
  assign("d", d, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(d)), c("d"), myenv)
  cderiv <- dbeta_dd(x = x, n = n, d = d, eps = eps, p = p, tau = tau)

  expect_equal(attr(nout, "gradient")[1, 1], cderiv, tol = 10 ^ -6)
}
)


test_that("dbeta_dtau works ok", {
  set.seed(9)
  x <- 10
  n <- 12
  d <- 0.1
  eps <- 0.1
  p <- 0.25
  tau <- 0.1

  tempfunc <- function(tau) {
    xi <- pbias(prob = p, bias = d, seq_error = eps)
    dbetabinom_mu_rho_cpp(x = x, size = n, mu = xi, rho = tau, return_log = FALSE)
  }

  myenv <- new.env()
  assign("tau", tau, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(tau)), c("tau"), myenv)
  cderiv <- dbeta_dtau(x = x, n = n, d = d, eps = eps, p = p, tau = tau)

  expect_equal(attr(nout, "gradient")[1, 1], cderiv, tol = 10 ^ -6)
}
)



test_that("the grad_offspring_original works", {
  set.seed(9)
  osize   <- stats::rbinom(n = 2, size = 70, prob = 0.4)
  ocounts <- stats::rbinom(n = 2, size = osize, prob = 0.4)
  ploidy <- 4
  p1geno <- 2
  p2geno <- 1
  d <- 2
  eps <- 0.01
  tau <- 0.1

  prob_geno <- get_prob_geno(ploidy = ploidy, model = "f1", p1geno = p1geno, p2geno = p2geno, allele_freq = 0.5)

  ## Numerical implementation ---------------------------------------
  tempfunc <- function(d, eps, tau) {
    obj_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy,
                  prob_geno = prob_geno, bias_val = d, seq_error = eps, od_param = tau,
                  outlier = FALSE)
  }

  myenv <- new.env()
  assign("d", d, envir = myenv)
  assign("eps", eps, envir = myenv)
  assign("tau", tau, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(d, eps, tau)), c("d", "eps", "tau"), myenv)

  ## Rcpp version ---------------------------------------------------
  cderiv <- grad_offspring_original(ocounts = ocounts, osize = osize, ploidy = ploidy,
                                    prob_geno = prob_geno,
                                    d = d, eps = eps, tau = tau)

  expect_equal(c(attr(nout, "gradient")), c(cderiv), tol = 10 ^ -5)
}
)


test_that("the grad_offspring_weights_original works", {
  set.seed(9)
  osize   <- stats::rbinom(n = 2, size = 70, prob = 0.4)
  ocounts <- stats::rbinom(n = 2, size = osize, prob = 0.4)
  ploidy <- 4
  p1geno <- 2
  p2geno <- 1
  d <- 2
  eps <- 0.01
  tau <- 0.1
  weight_vec <- stats::runif(2)

  prob_geno <- get_prob_geno(ploidy = ploidy, model = "f1", p1geno = p1geno, p2geno = p2geno, allele_freq = 0.5)

  ## Numerical implementation ---------------------------------------
  tempfunc <- function(d, eps, tau) {
    obj_offspring_weights(ocounts = ocounts, weight_vec = weight_vec,
                          osize = osize, ploidy = ploidy,
                          prob_geno = prob_geno, bias_val = d, seq_error = eps, od_param = tau,
                          outlier = FALSE)
  }

  myenv <- new.env()
  assign("d", d, envir = myenv)
  assign("eps", eps, envir = myenv)
  assign("tau", tau, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(d, eps, tau)), c("d", "eps", "tau"), myenv)

  ## Rcpp version ---------------------------------------------------
  cderiv <- grad_offspring_weights_original(ocounts = ocounts, osize = osize,
                                            weight_vec = weight_vec, ploidy = ploidy,
                                            prob_geno = prob_geno,
                                            d = d, eps = eps, tau = tau)

  expect_equal(c(attr(nout, "gradient")), c(cderiv), tol = 10 ^ -4)
}
)


test_that("outlier_grad matches outlier_obj", {

  ocounts <- c(12, 14, 16)
  osize   <- c(20, 30, 40)
  weight_vec <- c(0.1, 0.2, 0.3)
  out_mean <- 0.4
  out_disp <- 0.1

  tempfunc <- function(out_mean, out_disp) {
    outlier_obj(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                out_mean = out_mean, out_disp = out_disp)
  }

  myenv <- new.env()
  assign("out_mean", out_mean, envir = myenv)
  assign("out_disp", out_disp, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(out_mean, out_disp)), c("out_mean", "out_disp"), myenv)
  cderiv <- outlier_grad(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                         out_mean = out_mean, out_disp = out_disp)
  expect_equal(c(attr(nout, "gradient")), cderiv, tol = 10 ^ -6)
}
)


test_that("grad_parent matchs obj_parent", {
  set.seed(10)
  pcounts   <- 11
  psize     <- 20
  ploidy    <- 4
  bias_val  <- 0.9
  seq_error <- 0.1
  od_param  <- 0.9
  pgeno     <- 1
  weight    <- 0.8

  tempfunc <- function(bias_val, seq_error, od_param) {
    obj_parent(pcounts = pcounts, psize = psize, pgeno = pgeno,
               ploidy = ploidy, bias_val = bias_val,
               seq_error = seq_error, od_param = od_param,
               outlier = FALSE, weight = weight)
  }

  myenv <- new.env()
  assign("bias_val", bias_val, envir = myenv)
  assign("seq_error", seq_error, envir = myenv)
  assign("od_param", od_param, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(bias_val, seq_error, od_param)), c("bias_val", "seq_error", "od_param"), myenv)
  cderiv <- grad_parent(pcounts = pcounts, psize = psize, pgeno = pgeno,
                        ploidy = ploidy, bias_val = bias_val, seq_error = seq_error,
                        od_param = od_param, weight = weight)
  expect_equal(c(attr(nout, "gradient")), cderiv, tol = 10 ^ -6)

}
)


test_that("grad_parent_reparam matches obj_parent_reparam", {
  set.seed(10)
  pcounts   <- 11
  psize     <- 20
  ploidy    <- 4
  bias_val  <- 0.9
  seq_error <- 0.1
  od_param  <- 0.9
  pgeno     <- 1
  weight    <- 0.8

  s <- log(bias_val)
  ell <- log(seq_error / (1 - seq_error))
  r <- log((1 - od_param) / od_param)

  ## test equality of obj functions while I have the reparameterizations ---
  o1 <- obj_parent_reparam(pcounts = pcounts, psize = psize, ploidy = ploidy, pgeno = pgeno, s = s, ell = ell, r = r, weight = weight)
  o2 <- obj_parent(pcounts = pcounts, psize = psize, ploidy = ploidy, pgeno = pgeno, bias_val = bias_val, seq_error = seq_error, od_param = od_param, weight = weight)
  expect_equal(o1, o2)

  ## now test gradient ---

  tempfunc <- function(s, ell, r) {
    obj_parent_reparam(pcounts = pcounts, psize = psize, ploidy = ploidy,
                       pgeno = pgeno, s = s, ell = ell, r = r, weight = weight)
  }

  myenv <- new.env()
  assign("s", s, envir = myenv)
  assign("ell", ell, envir = myenv)
  assign("r", r, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(s, ell, r)), c("s", "ell", "r"), myenv)
  cderiv <- grad_parent_reparam(pcounts = pcounts, psize = psize, pgeno = pgeno,
                                ploidy = ploidy, s = s, ell = ell, r = r,
                                weight = weight)
  expect_equal(c(attr(nout, "gradient")), cderiv, tol = 10 ^ -5)
}
)


test_that("grad_wrapp_all matches obj_wrap_all", {
  parvec <- c(1, 0, 1)
  ocounts <- c(10, 11)
  osize <- c(20, 21)
  weight_vec <- c(0.5, 0.7)
  p1counts <- 12
  p1size <- 13
  p1weight <- 0.9
  p2counts <- 20
  p2size <- 41
  p2weight <- 0.4
  ploidy <- 4
  p1geno <- 4
  p2geno <- 2

  tempfunc <- function(parvec) {
    obj_wrapp_all(parvec = parvec, ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                  ploidy = ploidy, p1geno = p1geno, p2geno = p2geno, p1counts = p1counts,
                  p1size = p1size, p1weight = p1weight, p2counts = p2counts, p2size = p2size,
                  p2weight = p2weight, bound_bias = FALSE)
  }

  myenv <- new.env()
  assign("parvec", parvec, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(parvec)), "parvec", myenv)
  cderiv <- grad_wrapp_all(parvec = parvec, ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                           ploidy = ploidy, p1geno = p1geno, p2geno = p2geno, p1counts = p1counts,
                           p1size = p1size, p1weight = p1weight, p2counts = p2counts, p2size = p2size,
                           p2weight = p2weight, bound_bias = FALSE)
  expect_equal(cderiv, c(attr(nout, "gradient")), tol = 10 ^ -6)

}
)

test_that("out_grad_parent matches dbetabinom_mu_rho_cpp_double", {
  pcounts  <- 13
  psize    <- 21
  out_mean <- 0.8
  out_disp <- 0.4
  weight   <- 0.8

  tempfunc <- function(out_mean, out_disp) {
    weight * dbetabinom_mu_rho_cpp_double(x = pcounts, size = psize, mu = out_mean, rho = out_disp, return_log = TRUE)
  }

  myenv <- new.env()
  assign("out_mean", out_mean, envir = myenv)
  assign("out_disp", out_disp, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(out_mean, out_disp)), c("out_mean", "out_disp"), myenv)
  cderiv <- out_grad_parent(pcounts = pcounts, psize = psize,
                            out_mean = out_mean, out_disp = out_disp,
                            weight = weight)
  expect_equal(cderiv, c(attr(nout, "gradient")), tol = 10 ^ -5)
}
)

test_that("out_grad_wrapp matches out_obj_wrapp", {
  parvec     <- c(0.5, 0.3)
  ocounts    <- c(11, 12)
  osize      <- c(20, 30)
  weight_vec <- c(0.1, 0.2)
  p1counts   <- 12
  p1size     <- 20
  p1weight   <- 0.3
  p2counts   <- 18
  p2size     <- 44
  p2weight   <- 0.2
  p1geno     <- 2
  p2geno     <- 2
  ploidy     <- 4

  tempfunc <- function(parvec) {
    out_obj_wrapp(parvec = parvec, ocounts = ocounts, osize = osize,
                  weight_vec = weight_vec, p1counts = p1counts,
                  p1size = p1size, p1weight = p1weight,
                  p2counts = p2counts, p2size = p2size, p2weight = p2weight)
  }

  myenv <- new.env()
  assign("parvec", parvec, envir = myenv)
  nout <- stats::numericDeriv(quote(tempfunc(parvec)), "parvec", myenv)
  cderiv <- out_grad_wrapp(parvec = parvec, ocounts = ocounts, osize = osize,
                           weight_vec = weight_vec, p1counts = p1counts,
                           p1size = p1size, p1weight = p1weight,
                           p2counts = p2counts, p2size = p2size, p2weight = p2weight)
  expect_equal(cderiv, c(attr(nout, "gradient")), tol = 10 ^ -6)
}
)


test_that("get_cov_mle works", {
  skip("need to load snpdat")
  data("snpdat")
  subdat    <- snpdat[snpdat$snp == "SNP1", ]
  ocounts   <- subdat$counts
  osize     <- subdat$size

  p1geno <- 5
  p2geno <- 5
  model  <- "f1"
  ploidy <- 6

  old_hess <- structure(c(-602.766513627832, -5.23302863060546, -26.1415245533606,
                          -5.23302863060546, -5.13155719179381, -1.95466333012144, -26.1415245533606,
                          -1.95466333012144, -17.585650065388), .Dim = c(3L, 3L))
  prob_geno <- get_prob_geno(ploidy = ploidy, model = model, p1geno = p1geno,
                             p2geno = p2geno, allele_freq = 0)
  bias_val  <- 1.022
  seq_error <- 0.001102
  od_param  <- 0.01264
  out_prop  <- 0.02
  bias_val_mean  <- 0
  bias_val_sd    <- 1
  seq_error_mean <- -4.7
  seq_error_sd   <- 1

  p1counts <- NULL
  p1size <- NULL
  p2counts <- NULL
  p2size <- NULL

  covout <- get_cov_mle(ocounts = ocounts, osize = osize, ploidy = ploidy, prob_geno = prob_geno, bias_val = bias_val, seq_error = seq_error, od_param = od_param, out_prop = out_prop, p1counts = p1counts, p1size = p1size, p1geno = p1geno, p2counts = p2counts, p2size = p2size, p2geno = p2geno, bias_val_mean = bias_val_mean, bias_val_sd = bias_val_sd, seq_error_mean = seq_error_mean, seq_error_sd = seq_error_sd, model = model)
  covout$cov
  old_cov <- -1 * solve(old_hess)
  expect_true(all(abs(old_cov - covout$cov[1:3, 1:3]) < 10 ^ -1))
}
)
