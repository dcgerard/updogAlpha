## Method for overdispersion using beta-binomial

#' EM algorithm for over-dispersed mixture model.
#'
#' @inheritParams updog_maximize
#' @param rho The overdispersion parameter. Needs to be between 0 and 1.
#'     If it is exactly 0, you should be using \code{\link{updog_maximize}}.
#' @param update_rho A logical. Should we we update the overdispersion parameter (\code{TRUE}) or not (\code{FALSE})?
#' @param out_mu The outlier mean.
#' @param out_rho The outlier overdispersion parameter.
#' @param update_pival A logical. Should we update the mixing proportion (\code{TRUE}) or not (\code{FALSE})?
#' @param update_outlier A logical. Should we update the outlier distribution (\code{TRUE}) or not (\code{FALSE})?
#'
#' @author David Gerard
up_max_bb <- function(ocounts, osize, qarray, r1vec, r2vec, pk, pival = 0.99,
                      out_mu = 0.5, out_rho = 0.8,
                      rho = 0.001,
                      est_fudge = TRUE,
                      tol = 10 ^ -4, itermax = 1000,
                      update_geno = FALSE, update_pival = FALSE, update_outlier = FALSE,
                      update_rho = FALSE) {

  ## check input ---------------------------------------------------------------
  ploidy <- length(r1vec) - 1
  assertthat::are_equal(ploidy + 1, length(r2vec))
  assertthat::are_equal(length(ocounts), length(osize))
  assertthat::assert_that(all(ocounts >= 0))
  assertthat::assert_that(all(osize >= ocounts))
  assertthat::assert_that(all(abs(dim(qarray) - ploidy - 1) < 10 ^ -14))
  assertthat::are_equal(sum(r1vec), 1)
  assertthat::are_equal(sum(r2vec), 1)
  assertthat::assert_that(all(pk >= 0))
  assertthat::assert_that(all(pk <= 1))
  assertthat::assert_that(is.logical(est_fudge))
  assertthat::assert_that(pival >= 0, pival <= 1)
  assertthat::assert_that(tol > 0)
  assertthat::assert_that(itermax > 0)
  assertthat::assert_that(is.logical(update_rho))
  assertthat::assert_that(rho > 0)
  assertthat::assert_that(rho < 1)

  if (sum(max(r1vec) == r1vec) > 1) {
    p1geno <- round(ploidy / 2)
  } else {
    p1geno <- which.max(r1vec) - 1
  }

  if (sum(max(r2vec) == r2vec) > 1) {
    p2geno <- round(ploidy / 2)
  } else {
    p2geno <- which.max(r2vec) - 1
  }

  objnew <- up_bb_obj(pival = pival, p1geno = p1geno, p2geno = p2geno,
                      rho = rho, out_mu = out_mu, out_rho = out_rho,
                      ocounts = ocounts, osize = osize, qarray = qarray,
                      r1vec = r1vec, r2vec = r2vec, pk = pk)

  err <- tol + 1
  iterindex <- 1
  while (iterindex < itermax & err > tol) {
    objold <- objnew
    fout <- up_bb_fix(pival = pival, p1geno = p1geno, p2geno = p2geno,
                      rho = rho, out_mu = out_mu, out_rho = out_rho,
                      ocounts = ocounts, osize = osize, qarray = qarray,
                      r1vec = r1vec, r2vec = r2vec, pk = pk, update_geno = update_geno,
                      update_pival = update_pival, update_rho = update_rho,
                      update_outlier = update_outlier)

    pival   <- fout$pival
    rho     <- fout$rho
    out_mu  <- fout$out_mu
    out_rho <- fout$out_rho
    p1geno  <- fout$p1geno
    p2geno  <- fout$p2geno

    objnew <- up_bb_obj(pival = pival, p1geno = p1geno, p2geno = p2geno,
                        rho = rho, out_mu = out_mu, out_rho = out_rho,
                        ocounts = ocounts, osize = osize, qarray = qarray,
                        r1vec = r1vec, r2vec = r2vec, pk = pk)

    err <- abs(objnew / objold - 1)
    assertthat::assert_that(objnew - objold >= -10^-8)
    iterindex <- iterindex + 1
  }

  ## calculate posterior summaries ------------------------------------------
  dbbmat <- matrix(NA, nrow = ploidy + 1, ncol = length(ocounts))
  for (kindex in 1:(ploidy + 1)) {
    dbbmat[kindex, ] <- dbetabinom_mu_rho(x = ocounts, size = osize,
                                          mu = pk[kindex], rho = rho, log = FALSE)
  }

  dbb_out <- dbetabinom_mu_rho(x = ocounts, size = osize, mu = out_mu, rho = out_rho,
                               log = FALSE)

  avec <- qarray[p1geno + 1, p2geno + 1, ]

  unprobmat <- sweep(x = pival * dbbmat, MARGIN = 2, STATS = (1 - pival) * dbb_out, FUN = `+`) * avec
  postprob <- sweep(x = unprobmat, MARGIN = 2, STATS = colSums(unprobmat), FUN = `/`)

  ogeno <- apply(postprob, 2, which.max) - 1

  ## return values -----------------------------------------------------------
  return_list           <- list()
  return_list$pival     <- pival
  return_list$rho       <- rho
  return_list$out_mu    <- out_mu
  return_list$out_rho   <- out_rho
  return_list$p1geno    <- p1geno
  return_list$p2geno    <- p2geno
  return_list$prob_ok   <- fout$theta
  return_list$opostprob <- postprob
  return_list$ogeno     <- ogeno
  return(return_list)
}

#' Objective function for overdispersed binomial mixture model
#'
#' @inheritParams updog_maximize
#' @param ocounts A vector of non-negative integers. The observed number of A's.
#' @param osize A vector of positive integers. The observed total counts.
#' @param rho The overdispersion parameter.
#' @param out_mu The outlier mean.
#' @param out_rho The outlier overdispersion parameter.
#' @param p1geno The genotype of parent 1. From 0 to ploidy.
#' @param p2geno The genotype of parent 2. From 0 to ploidy.
#'
#' @author David Gerard
#'
up_bb_obj <- function(pival, p1geno, p2geno, rho, out_mu, out_rho,
                      ocounts, osize, qarray, r1vec, r2vec, pk) {

    ## check input -----------------------------------------------------------
    ploidy <- length(pk) - 1
    assertthat::assert_that(all(ocounts >= 0))
    assertthat::assert_that(all(osize >= ocounts))
    assertthat::assert_that(all(pk >= 0))
    assertthat::assert_that(all(pk <= 1))
    assertthat::assert_that(rho >= 0)
    assertthat::assert_that(rho <= 1)
    assertthat::assert_that(pival >= 0)
    assertthat::assert_that(pival <= 1)
    assertthat::assert_that(out_mu >= 0)
    assertthat::assert_that(out_mu <= 1)
    assertthat::assert_that(out_rho >= 0)
    assertthat::assert_that(out_rho <= 1)
    assertthat::assert_that(p1geno %in% 0:ploidy)
    assertthat::assert_that(p2geno %in% 0:ploidy)
    assertthat::assert_that(all(abs(dim(qarray) - ploidy - 1) < 10 ^ -14))
    assertthat::assert_that(all(r1vec >= 0))
    assertthat::assert_that(all(r1vec <= 1))
    assertthat::assert_that(all(r2vec >= 0))
    assertthat::assert_that(all(r2vec <= 1))
    assertthat::are_equal(length(r1vec), ploidy + 1)
    assertthat::are_equal(length(r2vec), ploidy + 1)

    ## bbmat ------------------------------------------------------------------
    dbbmat <- matrix(NA, nrow = ploidy + 1, ncol = length(ocounts))
    for (kindex in 1:(ploidy + 1)) {
        dbbmat[kindex, ] <- dbetabinom_mu_rho(x = ocounts, size = osize,
                                              mu = pk[kindex], rho = rho, log = FALSE)
    }

    dbb_out <- dbetabinom_mu_rho(x = ocounts, size = osize, mu = out_mu, rho = out_rho,
                                 log = FALSE)

    avec <- qarray[p1geno + 1, p2geno + 1, ]

    lval <- sum(log(pival * colSums(avec * dbbmat) + (1 - pival) * dbb_out)) +
        log(r1vec[p1geno + 1]) + log(r2vec[p2geno + 1])

    return(lval)

}

#' EM iteration for overdispersed binomial mixture model.
#'
#' @inheritParams updog_maximize
#' @inheritParams up_bb_obj
#' @param update_pival A logical. Should we update the mixture proportion
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_rho A logical. Should we update the overdispersion
#'     parameters (\code{TRUE}) or not (\code{FALSE})?
#' @param update_geno A logical. Should we update the parental
#'     genotypes (\code{TRUE}) or not (\code{FALSE})?
#' @param update_outlier A logical. Should we update the outlier
#'     distribution (\code{TRUE}) or not (\code{FALSE})?
#' @param p1geno The current genotype of parent 1. Must be in \code{0:ploidy}.
#' @param p2geno The current genotype of parent 1. Must be in \code{0:ploidy}.
#' @param upper_rho The upperbound on the overdispersion parameter. Defaults to 1 / (ploidy + 1) which corrsponds
#'     to bounding smallest probabilities to have finite density at 0 and 1.
#'
#' @author David Gerard
#'
up_bb_fix <- function(pival, p1geno, p2geno, rho, out_mu, out_rho,
                      ocounts, osize, qarray, r1vec, r2vec, pk,
                      update_pival = TRUE, update_rho = TRUE,
                      update_geno = FALSE, update_outlier = TRUE, upper_rho = NULL) {

    ## check input -----------------------------------------------------------
    ploidy <- length(pk) - 1
    assertthat::assert_that(all(ocounts >= 0))
    assertthat::assert_that(all(osize >= ocounts))
    assertthat::assert_that(all(pk >= 0))
    assertthat::assert_that(all(pk <= 1))
    assertthat::assert_that(rho >= 0)
    assertthat::assert_that(rho <= 1)
    assertthat::assert_that(pival >= 0)
    assertthat::assert_that(pival <= 1)
    assertthat::assert_that(out_mu >= 0)
    assertthat::assert_that(out_mu <= 1)
    assertthat::assert_that(out_rho >= 0)
    assertthat::assert_that(out_rho <= 1)
    assertthat::assert_that(p1geno %in% 0:ploidy)
    assertthat::assert_that(p2geno %in% 0:ploidy)
    assertthat::assert_that(all(abs(dim(qarray) - ploidy - 1) < 10 ^ -14))
    assertthat::assert_that(all(r1vec >= 0))
    assertthat::assert_that(all(r1vec <= 1))
    assertthat::assert_that(all(r2vec >= 0))
    assertthat::assert_that(all(r2vec <= 1))
    assertthat::are_equal(length(r1vec), ploidy + 1)
    assertthat::are_equal(length(r2vec), ploidy + 1)

    if (is.null(upper_rho)) {
      upper_rho <- 1 / (ploidy + 1)
    }

    assertthat::assert_that(upper_rho >= 10^-6)
    assertthat::assert_that(upper_rho <= 1)

    ## bbmat ------------------------------------------------------------------
    dbbmat <- matrix(NA, nrow = ploidy + 1, ncol = length(ocounts))
    for (kindex in 1:(ploidy + 1)) {
        dbbmat[kindex, ] <- dbetabinom_mu_rho(x = ocounts, size = osize,
                                              mu = pk[kindex], rho = rho, log = FALSE)
    }

    dbb_out <- dbetabinom_mu_rho(x = ocounts, size = osize, mu = out_mu, rho = out_rho,
                                 log = FALSE)

    avec <- qarray[p1geno + 1, p2geno + 1, ]

    ## local probability of being an outlier
    num <- pival * colSums(dbbmat * avec)
    denom <- num + (1 - pival) * dbb_out
    theta <- num / denom

    ## update pival ----------------------------------------------------------
    if (update_pival) {
        pival <- mean(theta)
    }

    ## update outlier model---------------------------------------------------
    if (update_outlier) {
        if (all(1 - theta < 10 ^ -10)) {
            break
        }

        eps_out <- 10^-8
        oout <- stats::optim(par = c(out_mu, out_rho), fn = outlier_bb_obj,
                             method = "L-BFGS-B", lower = c(eps_out, eps_out),
                             upper = c(1 - eps_out, 1 - eps_out),
                             control = list(fnscale = -1, maxit = 20),
                             ocounts = ocounts, osize = osize, theta = theta)

        out_mu <- oout$par[1]
        out_rho <- oout$par[2]
    }

    ## jointly update parental genotypes and variance inflation parameter
    if (!update_geno & update_rho) {
        eps_out <- 10^-6
        oout <- stats::optim(par = rho, fn = good_bb_obj,
                             method = "Brent", lower = eps_out, upper = upper_rho,
                             control = list(fnscale = -1), ocounts = ocounts,
                             osize = osize, theta = theta, pk = pk, avec = avec)
        rho <- oout$par
    } else if (update_geno & !update_rho) {
      marg_lik_mat <- matrix(NA, nrow = ploidy + 1, ncol = ploidy + 1)
      for (ell1 in 0:ploidy) {
        for(ell2 in 0:ploidy) {
          if (r1vec[ell1 + 1] <= 0 | r2vec[ell2 + 1] <= 0) {
            marg_lik_mat[ell1 + 1, ell2 + 1] <- -Inf
          } else {
            current_avec <- qarray[ell1 + 1, ell2 + 1, ]
            marg_lik_mat[ell1 + 1, ell2 + 1] <- good_bb_obj(rho = rho, ocounts = ocounts, osize = osize,
                                                            theta = theta, avec = current_avec, pk = pk) +
              log(r1vec[ell1 + 1]) + log(r2vec[ell2 + 1])
          }
        }
      }
      best_indices <- which(marg_lik_mat == max(marg_lik_mat), arr.ind = TRUE)
      p1geno <- best_indices[1] - 1
      p2geno <- best_indices[2] - 1
    } else if (update_geno & update_rho) {
      eps_out <- 10 ^ -6
      marg_lik_mat <- matrix(NA, nrow = ploidy + 1, ncol = ploidy + 1)
      out_rho_mat  <- matrix(NA, nrow = ploidy + 1, ncol = ploidy + 1)
      for (ell1 in 0:ploidy) {
        for (ell2 in 0:ploidy) {
          if (r1vec[ell1 + 1] <= 0 | r2vec[ell2 + 1] <= 0) {
            marg_lik_mat[ell1 + 1, ell2 + 1] <- -Inf
          } else {
            current_avec <- qarray[ell1 + 1, ell2 + 1, ]
            oout <- stats::optim(par = rho, fn = good_bb_obj,
                                 method = "Brent", lower = eps_out, upper = upper_rho,
                                 control = list(fnscale = -1), ocounts = ocounts,
                                 osize = osize, theta = theta, pk = pk, avec = current_avec)
            marg_lik_mat[ell1 + 1, ell2 + 1] <- oout$value + log(r1vec[ell1 + 1]) + log(r2vec[ell2 + 1])
            out_rho_mat[ell1 + 1, ell2 + 1] <- oout$par
          }
        }
      }
      best_indices <- which(marg_lik_mat == max(marg_lik_mat), arr.ind = TRUE)
      p1geno <- best_indices[1] - 1
      p2geno <- best_indices[2] - 1
      rho <- out_rho_mat[best_indices[1], best_indices[2]]
    }

    return_list         <- list()
    return_list$pival   <- pival
    return_list$out_mu  <- out_mu
    return_list$out_rho <- out_rho
    return_list$rho     <- rho
    return_list$theta   <- theta
    return_list$p1geno  <- p1geno
    return_list$p2geno  <- p2geno
    return(return_list)
}

#' Objective function for outlier model.
#'
#' @inheritParams updog
#' @param par A 2-vector. The first element is the mean of the
#'      beta-binomial. The second is the overdispersion parameter.
#' @param theta A vector of proportions. The posterior probabilities of each point being ok.
#'
#' @author David Gerard
#'
outlier_bb_obj <- function(par, ocounts, osize, theta) {
    out_mu  <- par[1]
    out_rho <- par[2]

    if (out_mu < 0 | out_mu > 1 | out_rho < 0 | out_rho > 1) {
        return(-Inf)
    }
    assertthat::assert_that(out_mu >= 0)
    assertthat::assert_that(out_mu <= 1)
    assertthat::assert_that(out_rho >= 0)
    assertthat::assert_that(out_rho <= 1)
    assertthat::assert_that(all(ocounts >= 0))
    assertthat::assert_that(all(osize >= ocounts))
    assertthat::assert_that(all(theta >= 0))
    assertthat::assert_that(all(theta <= 1))

    lvec <- dbetabinom_mu_rho(x = ocounts, size = osize, mu = out_mu,
                              rho = out_rho, log = TRUE)

    return(sum(lvec * (1 - theta)))
}

#' Objective function for non-outlier model
#'
#' @inheritParams updog
#' @inheritParams up_bb_obj
#' @param avec The current slice of \code{q_array}.
#' @param pk The parameterspace.
#' @param theta A vector of proportions. The posterior probabilities of each point being ok.
#'
#' @author David Gerard
#'
good_bb_obj <- function(rho, ocounts, osize, theta, avec, pk) {
    ploidy <- length(avec) - 1
    assertthat::are_equal(length(pk), ploidy + 1)
    assertthat::are_equal(length(ocounts), length(osize))
    assertthat::are_equal(length(ocounts), length(theta))
    assertthat::assert_that(all(ocounts >= 0))
    assertthat::assert_that(all(osize >= ocounts))
    assertthat::assert_that(rho >= 0)
    assertthat::assert_that(rho <= 1)

    dbbmat <- matrix(NA, nrow = ploidy + 1, ncol = length(ocounts))
    for (kindex in 1:(ploidy + 1)) {
        dbbmat[kindex, ] <- dbetabinom_mu_rho(x = ocounts, size = osize,
                                              mu = pk[kindex], rho = rho, log = FALSE)
    }

    return(sum(theta * log(colSums(dbbmat * avec))))
}

#' Beta-binomial density function.
#'
#' Note that it is possible here for \code{alpha = 0} or for \code{beta = 0} (but not both).
#' In which case, this function becomes the indicator function for \code{x == 0}
#' (when \code{alpha = 0}) or \code{x == size} (when \code{beta = 0}).
#'
#' @param x The observed counts.
#' @param size The number of trials.
#' @param alpha The alpha parameter of the underlying beta.
#' @param beta The beta parameter of the underlying beta
#' @param log A logical. Should we return the log of the density
#'     (\code{TRUE}) or not (\code{FALSE})?
#'
#' @export
#'
#' @author David Gerard
#'
dbetabinom <- function(x, size, alpha, beta, log = FALSE) {

  ## check input
  assertthat::assert_that(all(x >= 0))
  assertthat::assert_that(all(size >= x))
  assertthat::assert_that(alpha >= 0)
  assertthat::assert_that(beta >= 0)
  assertthat::assert_that(is.logical(log))

  if (alpha > 0 & beta > 0) {
    ldense <- lchoose(size, x) + lbeta(x + alpha, size - x + beta) - lbeta(alpha, beta)
  }
  else if (alpha == 0 & beta > 0) {
    ldense <- log(x == 0)
  } else if (beta == 0 & alpha > 0) {
    ldense <- log(x == size)
  } else {
    stop("alpha and beta cannot both be zero.")
  }

  if (log) {
    return(ldense)
  } else {
    return(exp(ldense))
  }
}

#' beta-binomial density function in terms of mean and correlation.
#'
#' @inheritParams dbetabinom
#'
#' @param mu The mean proportion. Equivalent to \code{alpha / (alpha + beta)} in \code{\link{dbetabinom}}.
#' @param rho The correlation/overdispersion parameter.
#'     Equivalent to \code{1 / (1 + alpha + beta)} in \code{\link{dbetabinom}}.
dbetabinom_mu_rho <- function(x, size, mu, rho, log = FALSE) {
  assertthat::assert_that(all(mu >= 0))
  assertthat::assert_that(all(mu <= 1))
  assertthat::assert_that(all(rho > 0))
  assertthat::assert_that(all(rho < 1))
  assertthat::assert_that(length(mu) == length(rho) | length(mu) == 1 | length(rho) == 1)

  alpha <- mu * (1 - rho) / rho
  beta  <- (1 - mu) * (1 - rho) / rho
  return(dbetabinom(x = x, size = size, alpha = alpha, beta = beta, log = log))
}

#' Posterior inference in beta-binomial model.
#'
#' The main difference between this function and \code{bin_post} is that here we also estimate an
#' overdispersion parameter by maximum marginal likelihood.
#'
#' @inheritParams bin_post
#' @param log A logical. Should we return the log probabilties (\code{TRUE}) or not (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @export
#'
#'
bb_post <- function(ncounts, ssize, prior, seq_error = NULL, log = FALSE) {

  if (abs(sum(prior) - 1) < 10 ^ -14) {
    ploidy <- length(prior) - 1
  } else {
    ploidy <- prior
    prior  <- rep(1 / (ploidy + 1), length = ploidy + 1)
  }

  if (is.null(seq_error)) {
    seq_error <- est_seq_error(ncounts = ncounts, ssize = ssize, ploidy = ploidy)
  }

  assertthat::are_equal(length(ncounts), length(ssize))
  assertthat::assert_that(seq_error >= 0)
  assertthat::assert_that(seq_error <= 1)
  assertthat::assert_that(all(ncounts >= 0))
  assertthat::assert_that(all(ssize >= ncounts))

  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)

  rho <- 0.01
  oout <- stats::optim(par = rho, fn = dbb_post_obj, ncounts = ncounts,
                       ssize = ssize, prior = prior, pk = pk,
                       control = list(fnscale = -1), method = "Brent",
                       lower = 0, upper = 1)
  rho <- oout$par

  dbbmat <- matrix(NA, nrow = ploidy + 1, ncol = length(ncounts))
  for(index in 1:(ploidy + 1)) {
    dbbmat[index, ] <- dbetabinom_mu_rho(x = ncounts, size = ssize, mu = pk[index], rho = rho, log = TRUE)
  }
  tmat <- rowSums(dbbmat) + prior
  tmat_max <- max(tmat)
  ldenom <- log(sum(exp(tmat - tmat_max))) + tmat_max

  lprob <- tmat - ldenom

  if (log) {
    return(list(lprob = lprob, rho = rho))
  }  else {
    probvec <- exp(lprob)
    assertthat::are_equal(sum(probvec), 1)
    return(list(prob = probvec, rho = rho))
  }
}


dbb_post_obj <- function(rho, ncounts, ssize, prior, pk) {

  ploidy <- length(prior) - 1
  assertthat::are_equal(ploidy, length(pk) - 1)
  assertthat::are_equal(length(ncounts), length(ssize))
  assertthat::assert_that(all(ncounts >= 0))
  assertthat::assert_that(all(ssize >= ncounts))
  assertthat::assert_that(rho >= 0)
  assertthat::assert_that(rho <= 1)

  dbbmat <- matrix(NA, nrow = ploidy + 1, ncol = length(ncounts))
  for(index in 1:(ploidy + 1)) {
    dbbmat[index, ] <- dbetabinom_mu_rho(x = ncounts, size = ssize, mu = pk[index], rho = rho, log = TRUE)
  }

  tmat <- rowSums(dbbmat) + prior

  tmat_max <- max(tmat)

  return(log(sum(exp(tmat - tmat_max))) + tmat_max)
}




