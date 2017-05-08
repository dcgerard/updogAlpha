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
  xi <- 1/4

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
  expect_equal(cderiv, rderiv)

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




