
#' This function just updates everything. No options allowed!
#'
#' This is the same as assuming a uniform prior on the parental genotypes,
#' then estimating these genotypes by maximum marginal likelihood.
#' Though this implementation does not allow for an outlier model.
#'
#' I mostly created this to check my gradient and objective function
#' implementations.
#'
#' @inheritParams updog
#' @param print_val A logical. Should we print the updates?
#'
#'
#' @author David Gerard
#'
#' @export
#'
updog_update_all <- function(ocounts, osize, ploidy, print_val = TRUE) {

  obj_wrapp_all <- function(parvec, ocounts, osize, ploidy, p1geno, p2geno) {
    obj_offspring_reparam(ocounts = ocounts, osize = osize,
                          ploidy = ploidy, p1geno = p1geno,
                          p2geno = p2geno,
                          s = parvec[1], ell = parvec[2], r = parvec[3])
  }

  grad_wrapp_all <- function(parvec, ocounts, osize, ploidy, p1geno, p2geno) {
    gout <- grad_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy,
                           p1geno = p1geno, p2geno = p2geno, s = parvec[1],
                           ell = parvec[2], r = parvec[3])
    return(gout)
  }

  ## starting values ------------------------------------------
  seq_error <- 0.01
  od_param  <- 0.01
  bias_val  <- 0.9

  s   <- log(bias_val)
  ell <- log(seq_error / (1 - seq_error))
  r   <- log((1 - od_param) / od_param)

  parvec <- c(s, ell, r)
  best_par <- c(0, 0, 0)
  best_p1 <- 0
  best_p2 <- 0
  best_llike <- -Inf

  for (p1geno in 0:ploidy) {
    for (p2geno in 0:p1geno) {
      oout <- stats::optim(par = parvec, fn = obj_wrapp_all, gr = grad_wrapp_all,
                           ocounts = ocounts, osize = osize, ploidy = ploidy,
                           p1geno = p1geno, p2geno = p2geno, method = "BFGS",
                           control = list(fnscale = -1, maxit = 1000))
      if (print_val) {
        cat("P1 P2 LL:", p1geno, p2geno, oout$value, "\n")
        cat("   B S O:", oout$par, "\n\n")
      }
      if (best_llike < oout$value) {
        best_par    <- oout$par
        best_llike <- oout$value
        best_p1    <- p1geno
        best_p2    <- p2geno
      }
    }
  }

  return_list           <- list()
  return_list$p1geno    <- best_p1
  return_list$p2geno    <- best_p2
  return_list$bias      <- exp(best_par[1])
  return_list$seq_error <- expit(best_par[2])
  return_list$od        <- 1 / (1 + exp(best_par[3]))
  return_list$llike     <- best_llike
  return(return_list)
}


#####################################################################################################
## Code to get posterior summaries after having done optimization -----------------------------------
#####################################################################################################

#' Simple posterior inference for beta-binomial when ncounts and ssize are only length 1.
#'
#' @inheritParams bb_post
#' @param p1geno The first parental genotype.
#' @param p2geno The second parental genotype.
#' @param bias_val The bias parameter. A value of 1 means no bias.
#' @param od_param The over-dispersion parameter. A value of 0 means no overdispersion.
#'     A value of 1 means total overdispersion.
#' @param ploidy The ploidy of the species.
#'
#' @author David Gerard
#'
#' @export
#'
bb_simple_post <- function(ncounts, ssize, ploidy, p1geno, p2geno, seq_error = 0, bias_val = 1, od_param = 0) {
  ## Check input --------------------------------------------------------------
  assertthat::assert_that(od_param >= 0, od_param < 1)
  assertthat::assert_that(bias_val > 0)
  assertthat::assert_that(seq_error >= 0, seq_error <= 1)
  assertthat::are_equal(length(ncounts), length(ssize), 1)
  assertthat::assert_that(ncounts >= 0, ncounts <= ssize)

  pvec <- get_pvec(ploidy = ploidy, bias_val = bias_val, seq_error = seq_error)
  qarray <- get_q_array(ploidy = ploidy)

  prob_vec <- rep(NA, length = ploidy + 1)
  for(index in 1:(ploidy + 1)) {
    prob_vec[index] <- dbetabinom_mu_rho_cpp(x = ncounts, size = ssize, mu = pvec[index], rho = od_param, return_log = FALSE) *
      qarray[p1geno + 1, p2geno + 1, index]
  }
  prob_vec <- prob_vec / sum(prob_vec)
  return(prob_vec)
}


#' This is a vanilla version of updog that will optimize everything, assuming a uniform prior on the parental
#' genotypes and not assuming that we have parental sequence data, and then return the posterior summaries.
#'
#' @inheritParams updog
#' @param print_val Should we print the updates?
#'
#' @author David Gerard
#'
#' @export
#'
updog_vanilla <- function(ocounts, osize, ploidy, print_val) {

  ## Get the best parameters
  parout <- updog_update_all(ocounts, osize, ploidy, print_val = print_val)
  postmat <- matrix(NA, ncol = ploidy + 1, nrow = length(ocounts))
  for (index in 1:length(ocounts)) {
    postmat[index, ] <- bb_simple_post(ncounts = ocounts[index], ssize = osize[index],
                                       ploidy = ploidy, p1geno = parout$p1geno,
                                       p2geno = parout$p2geno, seq_error = parout$seq_error, bias_val = parout$bias,
                                       od_param = parout$od)
  }
  return(list(par = parout, post = postmat))
}


