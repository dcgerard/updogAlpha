## Diagnostic functions


#' Samples counts using the parameter values in an \code{updog} object.
#'
#' @param obj An object of class \code{updog} returned from \code{\link{updog}}.
#'
#' @author David Gerard
#'
#' @export
#'
rupdog <- function(obj) {
  assertthat::are_equal(class(obj), "updog")
  if (obj$input$model == "s1") {
    if (obj$p1geno != obj$p2geno & obj$p2geno >= 0) {
      warning("p2geno != p1geno and model = s1. Setting p2geno <- p1geno.")
    }
    obj$p2geno <- obj$p1geno
  }

  ## Possible mean probabilities ---
  pvec <- get_pvec(ploidy = obj$input$ploidy,
                   bias_val = obj$bias_val,
                   seq_error = obj$seq_error)

  ## segregation probabilities ---
  seg_probs <- get_prob_geno(ploidy = obj$input$ploidy, model = obj$input$model, p1geno = obj$p1geno,
                p2geno = obj$p2geno, allele_freq = obj$allele_freq)

  ## sample offspring genotypes --
  ogeno <- sample(x = 0:obj$input$ploidy, size = length(obj$input$osize),
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

  ## Now sample parents if have them ---
  if (!is.null(obj$input$p1counts) & !is.null(obj$input$p1size) & (obj$input$model == "f1" | obj$input$model == "s1")) {
    p1_isout <- sample(x = c(TRUE, FALSE), prob = c(obj$out_prop, 1 - obj$out_prop), size = 1)
    if (p1_isout) {
      p1counts <- rbetabinom_mu_rho(mu = obj$out_mean, rho = obj$out_disp, size = obj$input$p1size)
    } else {
      p1counts <- rbetabinom_mu_rho(mu = pvec[obj$p1geno + 1], rho = obj$od_param, size = obj$input$p1size)
    }

    new_obj$input$p1counts <- p1counts
    new_obj$p1_prob_out    <- p1_isout * 1
  }
  if (!is.null(obj$input$p2counts) & !is.null(obj$input$p2size) & obj$input$model == "f1") {
    p2_isout <- sample(x = c(TRUE, FALSE), prob = c(obj$out_prop, 1 - obj$out_prop), size = 1)
    if (p2_isout) {
      p2counts <- rbetabinom_mu_rho(mu = obj$out_mean, rho = obj$out_disp, size = obj$input$p2size)
    } else {
      p2counts <- rbetabinom_mu_rho(mu = pvec[obj$p2geno + 1], rho = obj$od_param, size = obj$input$p2size)
    }
    new_obj$input$p2counts <- p2counts
    new_obj$p2_prob_out    <- p2_isout * 1
  }

  new_obj$llike <- dupdog(new_obj)

  return(new_obj)
}

#' Returns log-likelihood of \code{updog} object.
#'
#' @inheritParams rupdog
#'
#' @author David Gerard
#'
#' @export
#'
dupdog <- function(obj) {
  assertthat::are_equal(class(obj), "updog")

  if (obj$input$model == "s1") {
    if (obj$p1geno != obj$p2geno & obj$p2geno >= 0) {
      warning("p2geno != p1geno and model = s1. Setting p2geno <- p1geno.")
    }
    obj$p2geno <- obj$p1geno
  }
  ## get genotype frequencies --------------------------------------------------------
  prob_geno <- get_prob_geno(ploidy = obj$input$ploidy, model = obj$input$model, p1geno = obj$p1geno, p2geno = obj$p2geno, allele_freq = obj$allele_freq)

  llike_new <- obj_offspring(ocounts   = obj$input$ocounts,
                             osize     = obj$input$osize,
                             ploidy    = obj$input$ploidy,
                             prob_geno = prob_geno,
                             bias_val  = obj$bias_val,
                             seq_error = obj$seq_error,
                             od_param  = obj$od_param,
                             outlier   = TRUE,
                             out_prop  = obj$out_prop,
                             out_mean  = obj$out_mean,
                             out_disp  = obj$out_disp)
  if (!is.null(obj$input$p1counts) & !is.null(obj$input$p1size) & (obj$input$model == "f1" | obj$input$model == "s1")) { ## add parent 1 if available
    llike_new <- llike_new + obj_parent(pcounts   = obj$input$p1counts,
                                        psize     = obj$input$p1size,
                                        ploidy    = obj$input$ploidy,
                                        pgeno     = obj$p1geno,
                                        bias_val  = obj$bias_val,
                                        seq_error = obj$seq_error,
                                        od_param  = obj$od_param,
                                        outlier   = TRUE,
                                        out_prop  = obj$out_prop,
                                        out_mean  = obj$out_mean,
                                        out_disp  = obj$out_disp)
  }
  if (!is.null(obj$input$p2counts) & !is.null(obj$input$p2size) & (obj$input$model == "f1")) { ## add parent 2 if available
    llike_new <- llike_new + obj_parent(pcounts   = obj$input$p2counts,
                                        psize     = obj$input$p2size,
                                        ploidy    = obj$input$ploidy,
                                        pgeno     = obj$p2geno,
                                        bias_val  = obj$bias_val,
                                        seq_error = obj$seq_error,
                                        od_param  = obj$od_param,
                                        outlier   = TRUE,
                                        out_prop  = obj$out_prop,
                                        out_mean  = obj$out_mean,
                                        out_disp  = obj$out_disp)
  }
  return(llike_new)
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


#' Returns the distribution of proportion of individuals genotyped correctly using an the oracle estimator
#' when the parameter values equal the estimated values.
#'
#' @param obj An object of class \code{updog}.
#' @param itermax The number of repetitions to draw.
#'
#' @return A vector of length \code{itermax} which contains samples from the distribution of proportion of
#'     individuals genotyped correctly by an oracle estimator.
#'
#' @export
#'
#' @author David Gerard
screen_oracle <- function(obj, itermax = 1000) {
  if (obj$input$model != "hw") {
    obj$allele_freq <- 0.5
  }
  if (obj$input$model == "s1") {
    obj$p2geno <- obj$p1geno
  } else if (obj$model == "uniform") {
    obj$p1geno <- 0
    obj$p2geno <- 0
  }

  prop_correct_vec <- rep(NA, length = itermax)
  for (index in 1:itermax) {
    robj <- rupdog(obj)
    prob_geno <- get_prob_geno(ploidy = robj$input$ploidy, model = robj$input$model,
                               p1geno = robj$p1geno, p2geno = robj$p2geno,
                               allele_freq = robj$allele_freq)
    postmat <- bbpost_tot(ocounts = robj$input$ocounts, osize = robj$input$osize,
                          ploidy = robj$input$ploidy, prob_geno = prob_geno,
                          seq_error = robj$seq_error,
                          bias_val = robj$bias_val,
                          od_param = robj$od_param,
                          outlier = TRUE, out_prop = robj$out_prop,
                          out_mean = robj$out_mean,
                          out_disp = robj$out_disp)
    ogeno <- apply(postmat, 1, which.max) - 1
    prop_correct_vec[index] <- mean(ogeno == robj$ogeno)
  }
  return(prop_correct_vec)
}





