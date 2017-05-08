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
  nout <- stats::numericDeriv(quote(updog:::dbetabinom_mu_rho_cpp(x, n, xi, tau, FALSE)), "xi", myenv)
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
  dense <- updog:::dbetabinom_mu_rho_cpp(x, n, xi, 1 / (h + 1),
                                         return_log = FALSE)

  comp1 <-  xi * digamma(x + xi * h)
  comp2 <- (1 - xi) * digamma(n - x + (1 - xi) * h)
  comp3 <- digamma(h)
  comp4 <- digamma(n + h)
  comp5 <- xi * digamma(xi * h)
  comp6 <- (1 - xi) * digamma((1 - xi) * h)

  rderiv <- (comp1 + comp2 + comp3 - comp4 - comp5 - comp6) * dense
  cderiv <- dbeta_dh(x, n, xi, h)

  cderiv2 <- dbeta_dh_ell(x = x, n = n, d = d, ell = ell, p = p, h = h)
  expect_equal(cderiv, rderiv)
  expect_equal(cderiv, cderiv2)

  ## Numerical version ---------------------------------------

  db_wrapper <- function(x, n, xi, h) {
    tau <- 1 / (h + 1)
    updog:::dbetabinom_mu_rho_cpp(x, n, xi, tau, return_log = FALSE)
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

test_that("dbeta_dd works ok", {
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
  nout <- stats::numericDeriv(quote(beta_wrap(x, n, d, ell, p, tau)), "d", myenv)
  cderiv <- dbeta_dd(x, n, d, ell, p, h)
  expect_equal(attr(nout, "gradient")[1, 1], cderiv)
}
)


