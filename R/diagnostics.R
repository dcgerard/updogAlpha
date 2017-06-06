## Diagnostic functions


#' Sample counts using the parameter values in an \code{updog} object.
#'
#' @param obj An object of class \code{updog} returned from \code{\link{updog}}.
#'
#' @author David Gerard
#'
#' @export
#'
rupdog <- function(obj) {
  assertthat::are_equal(class(obj), "updog")

  ## Possible mean probabilities ---
  pvec <- get_pvec(ploidy = obj$input$ploidy,
                   bias_val = obj$bias_val,
                   seq_error = obj$seq_error)

  ## segregation probabilities ---
  seg_probs <- get_q_array_cpp(ploidy = obj$input$ploidy)[obj$p1geno + 1, obj$p2geno + 1, ]

  ## sample offspring genotypes --
  rogeno <- sample(x = 0:obj$input$ploidy, size = length(obj$input$ocounts),
                   prob = seg_probs, replace = TRUE)



}


#' Draw from beta-binomial distribution parameterized
#' by mean and overdispersion parameter.
#'
#' @param mu The mean of the underlying beta distribution. This is
#'     equal to \eqn{a / (a + b)} where \eqn{a} and \eqn{b} are the two
#'     shape parameters of the usual parameterization of the beta. This can
#'     be a vector of values.
#' @param rho The over-dispersion parameter of the underlying beta
#'     distribution. This is equal to \eqn{1 / (a + b + 1)} where \eqn{a}
#'     and \eqn{b} are the two shape parameters of the usual parameterization
#'     of the beta. This can be a vector of values.
#' @param size The sample size. This can equal a vector of observations.
#'
#' @author David Gerard
#'
#' @export
#'
rbetabinom_mu_rho <- function(mu, rho, size) {
  tol <- 10 ^ -12
  ## check input --------------------------------------------
  assertthat::assert_that(all(mu >= 0), all(mu <= 1))
  assertthat::assert_that(all(rho >= 0), all(rho <= 1))
  assertthat::assert_that(all(size >= 0))
  assertthat::are_equal(length(mu), length(rho), length(size))

  ## sample probabilities ---

  which_normal <- (! mu < tol) & !(mu > 1 - 10 ^ -12) & !(rho < tol) & !(rho > 1 - tol)
  a <- mu * (1 - rho) / rho
  b <- (1 - mu) * (1 - rho) / rho
  psamp <- rep(NA, length = length(mu))
  psamp[which_normal] <- stats::rbeta(n = sum(which_normal), shape1 = a[which_normal], shape2 = b[which_normal])
  psamp[rho < tol]     <- mu[rho < tol]
  which_extreme <- rho > 1 - tol
  psamp[which_extreme] <- stats::rbinom(n = sum(which_extreme), size = 1, prob = mu[which_extreme])
  psamp[mu < tol]      <- 0
  psamp[mu > 1 - tol]  <- 1

  ## sample counts ---
  return(stats::rbinom(n = length(psamp), size = size, prob = psamp))
}





