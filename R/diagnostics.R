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
  ogeno <- sample(x = 0:obj$input$ploidy, size = length(obj$input$ocounts),
                   prob = seg_probs, replace = TRUE)


  ## Sample outliers ---
  prob_out <- stats::rbinom(n = length(ogeno), size = 1, prob = obj$out_prop)
  prob_ok  <- 1 - prob_out

  ## Collate outlier and good infor ---
  ogeno[prob_out == 1] <- NA
  mean_vec <- pvec[ogeno + 1]
  mean_vec[is.na(mean_vec)] <- obj$out_mean
  od_vec   <- rep(obj$od_param, length(mean_vec))
  od_vec[prob_out == 1] <- obj$out_disp

  ocounts <- rbetabinom_mu_rho(mu = mean_vec, rho = od_vec, size = obj$input$osize)


  new_obj <- obj
  new_obj$input$ocounts <- ocounts
  new_obj$ogeno         <- ogeno
  new_obj$prob_ok       <- prob_ok
  new_obj$prob_out      <- prob_out
  new_obj$num_iter      <- NULL
  new_obj$convergence   <- NULL

  return(new_obj)
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





