
#' wrapper for \code{\link{obj_offspring_weights_reparam}}
#' @inheritParams obj_offspring_vec
#' @param parvec A vector of three elements, s, ell, and r. We have s = log(bias_val) = log(d),
#' ell = logit(seq_error) = logit(eps), and r = 1 / logit(od_param) = 1 / logit(tau).
#' @param weight_vec A vector of weights obtained via the E-step.
#'
#' @author David Gerard
#'
obj_wrapp_all <- function(parvec, ocounts, osize, weight_vec,
                          ploidy, p1geno, p2geno) {
  ## Check if second to last pk is too large -----------------------------------------
  eps <- expit(parvec[2])
  mind <- eps / (1 - eps) + 0.05 ## ad hoc bound
  if (exp(parvec[1]) < mind) {
    return(-Inf)
  }

  ## Normal value ------------------------------------------------------------------------
  val <- obj_offspring_weights_reparam(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                                        ploidy = ploidy, p1geno = p1geno,
                                        p2geno = p2geno,
                                        s = parvec[1], ell = parvec[2], r = parvec[3])
  return(val)
}

#' wrapper for \code{\link{grad_offspring_weights}}
#'
#' @inheritParams obj_offspring_vec
#' @inheritParams obj_wrapp_all
#'
#' @author David Gerard
#'
grad_wrapp_all <- function(parvec, ocounts, osize, weight_vec, ploidy, p1geno, p2geno) {
  gout <- grad_offspring_weights(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                                 ploidy = ploidy, p1geno = p1geno, p2geno = p2geno, s = parvec[1],
                                 ell = parvec[2], r = parvec[3])
  return(gout)
}


#' Update the OK points
#'
#' @inheritParams obj_offspring_vec
#' @inheritParams obj_wrapp_all
#'
#' @author David Gerard
#'
update_good <- function(parvec, ocounts, osize, weight_vec, ploidy,
                        p1geno = NULL, p2geno = NULL) {
  best_par <- parvec
  if (is.null(p1geno) | is.null(p2geno)) { ## try all p1geno/p2geno combos --
    best_p1 <- 0
    best_p2 <- 0
    best_llike <- -Inf

    for (p1geno in 0:ploidy) {
      for (p2geno in 0:p1geno) {
        ## If start position has -Inf, start somewhere else
        val <- obj_offspring_weights_reparam(ocounts = ocounts, osize = osize,
                                             weight_vec = weight_vec,
                                             ploidy = ploidy, p1geno = p1geno,
                                             p2geno = p2geno,
                                             s = parvec[1], ell = parvec[2],
                                             r = parvec[3])
        if (val == -Inf) {
          start_vec <- c(0, -4.5, 4.5) ## about (1, 0.01, 0.01) for (bias_val, seq_error, od_param)
        } else {
          start_vec <- parvec
        }
        oout <- stats::optim(par = start_vec, fn = obj_wrapp_all, gr = grad_wrapp_all,
                             ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                             ploidy = ploidy, p1geno = p1geno, p2geno = p2geno, method = "BFGS",
                             control = list(fnscale = -1, maxit = 1000))
        if (oout$convergence != 0) {
          warning(oout$message)
        }
        ## cat(p1geno, p2geno, ":", oout$value, "\n")
        if (best_llike < oout$value) {
          best_par   <- oout$par
          best_llike <- oout$value
          best_p1    <- p1geno
          best_p2    <- p2geno
        }
      }
    }
  } else { ## run only user-provided parental genotypes ----
    ## If start position has -Inf, start somewhere else
    val <- obj_offspring_weights_reparam(ocounts = ocounts, osize = osize,
                                         weight_vec = weight_vec,
                                         ploidy = ploidy, p1geno = p1geno,
                                         p2geno = p2geno,
                                         s = parvec[1], ell = parvec[2],
                                         r = parvec[3])
    if (val == -Inf) {
      start_vec <- c(0, -4.5, 4.5) ## about (1, 0.01, 0.01) for (bias_val, seq_error, od_param)
    } else {
      start_vec <- parvec
    }
    oout <- stats::optim(par = start_vec, fn = obj_wrapp_all, gr = grad_wrapp_all,
                         ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                         ploidy = ploidy, p1geno = p1geno, p2geno = p2geno, method = "BFGS",
                         control = list(fnscale = -1, maxit = 1000))
    best_p1 <- p1geno
    best_p2 <- p2geno
    best_llike <- oout$value
    best_par   <- oout$par
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

out_obj_wrapp <- function(obj, ocounts, osize, weight_vec) {
  outlier_obj(ocounts = ocounts, osize = osize, weight_vec = weight_vec, out_mean = obj[1], out_disp = obj[2])
}

out_grad_wrapp <- function(obj, ocounts, osize, weight_vec) {
  outlier_grad(ocounts = ocounts, osize = osize, weight_vec = weight_vec, out_mean = obj[1], out_disp = obj[2])
}



#' This function just updates everything. No options allowed!
#'
#' This is the same as assuming a uniform prior on the parental genotypes,
#' then estimating these genotypes by maximum marginal likelihood.
#' Though this implementation does not allow for an outlier model.
#'
#'
#' @inheritParams updog
#' @param print_val A logical. Should we print the updates?
#' @param tol The stopping criterion
#' @param maxiter The maximum number of iterations
#' @param commit_num The number of consecutive iterations where the parental
#'     genotypes do not change before we commit to those parental genotypes.
#' @param update_outmean A logical. Should we update \code{out_mean}?
#'
#' @author David Gerard
#'
#' @export
#'
updog_update_all <- function(ocounts, osize, ploidy, print_val = TRUE,
                             tol = 10 ^ -4, maxiter = 500,
                             commit_num = 4, update_outmean = FALSE) {

  ## starting values ------------------------------------------
  seq_error  <- 0.01
  od_param   <- 0.01
  bias_val   <- 1
  out_prop   <- 0.01
  out_mean   <- 0.5
  out_disp   <- 1/3
  weight_vec <- rep(out_prop, length = length(ocounts))
  p1geno     <- 3
  p2geno     <- 3
  non_mono   <- 0
  non_mono_max <- 2 # non-monotonicity leniancy.

  llike_new <- -Inf

  index <- 1
  err <- tol + 1
  parental_count <- 0 ## counts number of consecutive times parental genotypes do not change.
  while ((index <= maxiter) & err > tol) {
    llike_old <- llike_new

    ## E-step ------------
    if (index > 1) {
      weight_vec <- get_out_prop(ocounts = ocounts, osize = osize,
                                 ploidy = ploidy, p1geno = p1geno,
                                 p2geno = p2geno, d = bias_val,
                                 eps = seq_error, tau = od_param,
                                 out_prop = out_prop, out_mean = out_mean,
                                 out_disp = out_disp)
    }


    ## Update out_prop ---
    out_prop <- mean(weight_vec)

    ## reparameterization --
    s   <- log(bias_val)
    ell <- log(seq_error / (1 - seq_error))
    r   <- log((1 - od_param) / od_param)
    parvec <- c(s, ell, r)

    ## update good --------
    if (parental_count >= 4) {
      gout <- update_good(parvec = parvec, ocounts = ocounts, osize = osize,
                          weight_vec = 1 - weight_vec, ploidy = ploidy,
                          p1geno = p1geno, p2geno = p2geno)
    } else {
      gout <- update_good(parvec = parvec, ocounts = ocounts, osize = osize,
                          weight_vec = 1 - weight_vec, ploidy = ploidy)
    }

    bias_val  <- gout$bias
    seq_error <- gout$seq_error
    od_param  <- gout$od


    ## update parental_count if no change
    if (p1geno == gout$p1geno & p2geno == gout$p2geno) {
      parental_count <- parental_count + 1
    } else {
      parental_count <- 0
    }
    p1geno <- gout$p1geno
    p2geno <- gout$p2geno

    ## update bad --------
    if (update_outmean) {
      oout <- stats::optim(par = c(out_mean, out_disp),
                           fn = out_obj_wrapp, gr = out_grad_wrapp,
                           ocounts = ocounts, osize = osize,
                           weight_vec = weight_vec, method = "BFGS",
                           control = list(fnscale = -1))
      out_mean <- oout$par[1]
      out_disp <- oout$par[2]
    } else {
      oout <- stats::optim(par = out_disp, fn = outlier_obj,
                           method = "Brent", control = list(fnscale = -1),
                           lower = 0, upper = 1 - 10 ^ -3,
                           ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                           out_mean = out_mean)
      out_disp <- oout$par
    }


    ## Calculate log-likelihood and update err and index -------
    llike_new <- obj_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy,
                               p1geno = p1geno, p2geno = p2geno, bias_val = bias_val,
                               seq_error = seq_error, od_param = od_param, outlier = TRUE,
                               out_prop = out_prop, out_mean = out_mean, out_disp = out_disp)
    err <- abs(llike_new - llike_old)

    if (print_val) {
      cat("    Log-Likelihood:", llike_new, "\n")
      cat("Parental Genotypes:", gout$p1geno, gout$p2geno, "\n")
      cat("              Bias:", bias_val, "\n")
      cat("  Sequencing Error:", seq_error, "\n")
      cat("   Over-dispersion:", od_param, "\n")
      cat("Outlier Proportion:", out_prop, "\n")
      cat("      Outlier Mean:", out_mean, "\n")
      cat("        Outlier OD:", out_disp, "\n\n")
    }
    if (index > 1) {
      if (llike_new - llike_old < -10 ^ -6) {
        non_mono <- non_mono + 1
      }
      if (non_mono >= non_mono_max) {
        stop("Your likelihood is not increasing.")
      }
    }

    index <- index + 1
  }

  return_list <- list()
  return_list$bias_val    <- bias_val
  return_list$seq_error   <- seq_error
  return_list$od_param    <- od_param
  return_list$p1geno      <- p1geno
  return_list$p2geno      <- p2geno
  return_list$out_prop    <- out_prop
  return_list$out_mean    <- out_mean
  return_list$out_disp    <- out_disp
  return_list$prob_out    <- weight_vec
  return_list$num_iter    <- index
  return_list$convergence <- (index >= maxiter) * 1
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
#' @param update_outmean A logical. Should we update \code{out_mean}?
#'
#' @author David Gerard
#'
#' @export
#'
updog_vanilla <- function(ocounts, osize, ploidy, print_val = FALSE, update_outmean = FALSE) {
  assertthat::assert_that(is.logical(print_val))
  assertthat::assert_that(is.logical(update_outmean))

  ## Get the best parameters
  parout <- updog_update_all(ocounts, osize, ploidy, print_val = print_val,
                             update_outmean = update_outmean)
  parout$postmat <- bbpost_tot(ocounts = ocounts, osize = osize,
                               ploidy = ploidy, p1geno = parout$p1geno,
                               p2geno = parout$p2geno,
                               seq_error = parout$seq_error,
                               bias_val = parout$bias_val,
                               od_param = parout$od_param,
                               outlier = TRUE, out_prop = parout$out_prop,
                               out_mean = parout$out_mean,
                               out_disp = parout$out_disp)
  parout$ogeno <- apply(parout$postmat, 1, which.max) - 1
  parout$ogeno[abs(parout$prob_out - 1) < 10 ^ -3] <- NA   ## Put NA for ogeno when prob_out is almost 1 ----
  parout$maxpostprob <- parout$postmat[cbind(1:nrow(parout$postmat), parout$ogeno + 1)]
  parout$prob_ok <- 1 - parout$prob_out
  parout$input <- list()
  parout$input$ocounts <- ocounts
  parout$input$osize <- osize
  parout$input$ploidy <- ploidy
  class(parout) <- "updog"
  return(parout)
}


