
#' wrapper for \code{\link{obj_offspring_weights_reparam}}
#' @inheritParams obj_offspring_vec
#' @param parvec A vector of three elements, s, ell, and r. We have s = log(bias_val) = log(d),
#' ell = logit(seq_error) = logit(eps), and r = 1 / logit(od_param) = 1 / logit(tau).
#' @param weight_vec A vector of weights obtained via the E-step.
#' @param bound_bias A logical. Should we bound the bias parameter \code{parvec[1]} by a somewhat arbitrary value
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param bound_od A logical. Should we bound the overdisperion parameter
#'     \code{parvec[3]} by a very small number (\code{TRUE}) or not
#'     (\code{FALSE}). This is mostly because my code for super super small
#'     values of \code{parvec[3]} can be unstable.
#' @param update_bias_val A logical. Not used. Here for compatability
#'     with \code{\link{grad_wrapp_all}}.
#' @param update_seq_error A logical. Not used. Here for compatability
#'     with \code{\link{grad_wrapp_all}}.
#' @param update_od_param A logical. Not used. Here for compatability
#'     with \code{\link{grad_wrapp_all}}.
#' @param p1counts The number of reference alleles observed from parent 1.
#' @param p1size The number of reads observed from parent 1.
#' @param p1weight The posterior probability that parent 1 is not an outlier.
#' @param p2counts The number of reference alleles observed from parent 2.
#' @param p2size The number of reads observed from parent 2.
#' @param p2weight The posterior probability that parent 2 is not an outlier.
#'
#' @author David Gerard
#'
obj_wrapp_all <- function(parvec, ocounts, osize, weight_vec,
                          ploidy, p1geno, p2geno,
                          p1counts = NULL, p1size = NULL, p1weight = NULL,
                          p2counts = NULL, p2size = NULL, p2weight = NULL,
                          bound_bias = TRUE,
                          bound_od = TRUE,
                          update_bias_val = TRUE,
                          update_seq_error = TRUE,
                          update_od_param = TRUE) {
  ## Check if second to last pk is too large -----------------------------------------
  if (bound_bias) {
    eps <- expit(parvec[2])
    mind <- eps / (1 - eps) + 0.05 ## ad hoc bound
    if (exp(parvec[1]) < mind) {
      return(-Inf)
    }
  }
  if (bound_od) {
    if (parvec[3] > 18) {
      return(-Inf)
    }
  }


  ## Normal value ------------------------------------------------------------------------
  val <- obj_offspring_weights_reparam(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                                        ploidy = ploidy, p1geno = p1geno,
                                        p2geno = p2geno,
                                        s = parvec[1], ell = parvec[2], r = parvec[3])

  if (!is.null(p1counts) & !is.null(p1size) & !is.null(p1weight)) { ## add contribution from parent 1
    val <- val + obj_parent_reparam(pcounts = p1counts, psize = p1size, ploidy = ploidy,
                                    pgeno = p1geno, s = parvec[1], ell = parvec[2],
                                    r = parvec[3], weight = p1weight)
  }
  if (!is.null(p2counts) & !is.null(p2size) & !is.null(p2weight)) { ## add contribution from parent 2
    val <- val + obj_parent_reparam(pcounts = p2counts, psize = p2size, ploidy = ploidy,
                                    pgeno = p2geno, s = parvec[1], ell = parvec[2],
                                    r = parvec[3], weight = p2weight)
  }

  return(val)
}

#' wrapper for \code{\link{grad_offspring_weights}}
#'
#' @inheritParams obj_offspring_vec
#' @inheritParams obj_wrapp_all
#' @param update_bias_val A logical. If \code{FALSE}, then the first position
#'     of the returned gradient will be zero.
#' @param update_seq_error A logical. If \code{FALSE}, then the second position
#'     of the returned gradient will be zero.
#' @param update_od_param A logical. If \code{FALSE}, then the third position
#'     of the returned gradient will be zero.
#'
#' @author David Gerard
#'
grad_wrapp_all <- function(parvec, ocounts, osize, weight_vec, ploidy, p1geno, p2geno,
                           p1counts = NULL, p1size = NULL, p1weight = NULL,
                           p2counts = NULL, p2size = NULL, p2weight = NULL,
                           bound_bias = TRUE,
                           update_bias_val = TRUE,
                           update_seq_error = TRUE,
                           update_od_param = TRUE) {
  gout <- grad_offspring_weights(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                                 ploidy = ploidy, p1geno = p1geno, p2geno = p2geno, s = parvec[1],
                                 ell = parvec[2], r = parvec[3])

  if (!is.null(p1counts) & !is.null(p1size) & !is.null(p1weight)) {
    gout <- gout + grad_parent_reparam(pcounts = p1counts, psize = p1size,
                                       ploidy = ploidy, pgeno = p1geno,
                                       s = parvec[1], ell = parvec[2],
                                       r = parvec[3], weight = p1weight)
  }
  if (!is.null(p2counts) & !is.null(p2size) & !is.null(p2weight)) {
    gout <- gout + grad_parent_reparam(pcounts = p2counts, psize = p2size,
                                       ploidy = ploidy, pgeno = p2geno,
                                       s = parvec[1], ell = parvec[2],
                                       r = parvec[3], weight = p2weight)
  }
  if (!update_bias_val) {
    gout[1] <- 0
  }
  if (!update_seq_error) {
    gout[2] <- 0
  }
  if (!update_od_param) {
    gout[3] <- 0
  }
  return(gout)
}

#' Update the OK points
#'
#' @inheritParams obj_offspring_vec
#' @inheritParams obj_wrapp_all
#' @param bound_bias A logical. Should we bound the bias parameter \code{parvec[1]} by a somewhat arbitrary value
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_bias_val A logical. Should we update the bias parameter
#'     (the first position of \code{parvec})?
#' @param update_seq_error A logical. Should we update the sequencing error parameter
#'     (the second position of \code{parvec})?
#' @param update_od_param A logical. Should we update the overdispersion parameter
#'     (the third position of \code{parvec})?
#'
#' @author David Gerard
#'
update_good <- function(parvec, ocounts, osize, weight_vec, ploidy,
                        p1counts = NULL, p1size = NULL, p1weight = NULL,
                        p2counts = NULL, p2size = NULL, p2weight = NULL,
                        p1geno = NULL, p2geno = NULL,
                        bound_bias = TRUE,
                        update_bias_val = TRUE,
                        update_seq_error = TRUE,
                        update_od_param = TRUE) {
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
                             ploidy = ploidy, p1geno = p1geno, p2geno = p2geno,
                             p1counts = p1counts, p1size = p1size, p1weight = p1weight,
                             p2counts = p2counts, p2size = p2size, p2weight = p2weight,
                             method = "BFGS",
                             control = list(fnscale = -1, maxit = 1000),
                             bound_bias = bound_bias,
                             update_bias_val = update_bias_val,
                             update_seq_error = update_seq_error,
                             update_od_param = update_od_param)
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
                         ploidy = ploidy, p1geno = p1geno, p2geno = p2geno,
                         p1counts = p1counts, p1size = p1size, p1weight = p1weight,
                         p2counts = p2counts, p2size = p2size, p2weight = p2weight,
                         method = "BFGS",
                         control = list(fnscale = -1, maxit = 1000),
                         bound_bias = bound_bias,
                         update_bias_val = update_bias_val,
                         update_seq_error = update_seq_error,
                         update_od_param = update_od_param)
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

#' Wrapper for \code{\link{outlier_obj}}.
#'
#' @param parvec A numeric vector of length two. The first element is the outlier mean.
#'     The second element is the outlier overdispersion parameter.
#' @param ocounts The offspring counts of the reference allele.
#' @param osize The offspring counts of reads.
#' @param weight_vec The probability a point is an outlier.
#' @param min_disp The bound on the overdispersion parameter (can't be too small).
#' @param p1counts The parent 1 counts of the reference allele.
#' @param p1size The parent 1 counts of reads.
#' @param p1weight The probability parent 1 is an outlier.
#' @param p2counts The parent 2 counts of the reference allele.
#' @param p2size The parent 2 counts of reads.
#' @param p2weight The probability parent 2 is an outlier.
#'
#' @author David Gerard
#'
#'
out_obj_wrapp <- function(parvec, ocounts, osize, weight_vec,
                          p1counts = NULL, p1size = NULL, p1weight = NULL,
                          p2counts = NULL, p2size = NULL, p2weight = NULL,
                          min_disp = 0) {
  if ((parvec[2] < min_disp) | (parvec[2] > 1 - 10 ^ -8) | (parvec[1] < 10 ^ -8) | (parvec[1] > 1 - 10 ^ -8)) {
    return(-Inf)
  } else {
    obj <- outlier_obj(ocounts = ocounts, osize = osize, weight_vec = weight_vec, out_mean = parvec[1], out_disp = parvec[2])
  }

  ## add parent data if we have it.
  if (!is.null(p1counts) & !is.null(p1size) & !is.null(p1weight)) {
    obj <- obj + p1weight * dbetabinom_mu_rho_cpp_double(x = p1counts, size = p1size, mu = parvec[1], rho = parvec[2], return_log = TRUE)
  }
  if (!is.null(p2counts) & !is.null(p2size) & !is.null(p2weight)) {
    obj <- obj + p2weight * dbetabinom_mu_rho_cpp_double(x = p2counts, size = p2size, mu = parvec[1], rho = parvec[2], return_log = TRUE)
  }
  return(obj)
}

#' Wrapper for \code{\link{outlier_grad}}.
#'
#' @inheritParams out_obj_wrapp
#'
#' @author David Gerard
#'
out_grad_wrapp <- function(parvec, ocounts, osize, weight_vec,
                           p1counts = NULL, p1size = NULL, p1weight = NULL,
                           p2counts = NULL, p2size = NULL, p2weight = NULL,
                           min_disp = 0) {
  grad <- outlier_grad(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                       out_mean = parvec[1], out_disp = parvec[2])
  if (!is.null(p1counts) & !is.null(p1size) & !is.null(p1weight)) { ## update grad if have parent 1 info
    grad <- grad + out_grad_parent(pcounts = p1counts, psize = p1size,
                                   out_mean = parvec[1], out_disp = parvec[2],
                                   weight = p1weight)
  }
  if (!is.null(p2counts) & !is.null(p2size) & !is.null(p2weight)) { ## update grad if have parent 2 info
    grad <- grad + out_grad_parent(pcounts = p2counts, psize = p2size,
                                   out_mean = parvec[1], out_disp = parvec[2],
                                   weight = p2weight)
  }
  return(grad)
}



#' This is the main optimization function for updog.
#'
#'
#' @inheritParams updog
#' @param print_val A logical. Should we print the updates (\code{TRUE}) or not (\code{FALSE})?
#' @param tol The stopping criterion
#' @param maxiter The maximum number of iterations
#' @param commit_num The number of consecutive iterations where the parental
#'     genotypes do not change before we commit to those parental genotypes.
#' @param min_disp The minimum value for the over-dispersion parameter of the
#'     outlier distribution. If this gets too small, then we can be overly confident in
#'     points being outliers.
#' @param update_outmean A logical. Should we update \code{out_mean}
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_outdisp A logical. Should we update \code{out_mean}
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_outprop A logical. Should we update \code{out_prop}
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_bias_val A logical. Should we update \code{bias_val}
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_seq_error A logical. Should we update \code{seq_error}
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_od_param A logical. Should we update \code{od_param}
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_pgeno A logical. Should we update \code{p1geno} and \code{p1geno}
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param p1geno The initial value of the first parental genotype.
#' @param p2geno The initial value of the second parental genotype.
#' @param seq_error The initial value of the sequencing error rate.
#' @param od_param The initial value of the overdispersion parameter.
#' @param bias_val The initial value of the bias parameter.
#' @param out_prop The initial value of the proportion of points that are outliers.
#' @param out_mean The initial value of the mean of the outlier distribution.
#' @param out_disp The initial value of the over-dispersion parameter of the outlier distribution.
#' @param non_mono_max The maximum number of iterations to allow non-monotonicity of likelihood.
#' @param bound_bias A logical. Should we bound \code{bias_val} by a somewhat arbitrary value
#'     (\code{TRUE}) or not (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @export
#'
updog_update_all <- function(ocounts, osize, ploidy,
                             p1counts = NULL,
                             p1size = NULL,
                             p2counts = NULL,
                             p2size = NULL,
                             print_val = TRUE,
                             tol = 10 ^ -4,
                             maxiter = 500,
                             commit_num = 4,
                             min_disp = 0,
                             update_outmean = FALSE,
                             update_outdisp = FALSE,
                             update_outprop = TRUE,
                             update_bias_val = TRUE,
                             update_seq_error = TRUE,
                             update_od_param = TRUE,
                             update_pgeno = TRUE,
                             p1geno = 3,
                             p2geno = 3,
                             seq_error = 0.01,
                             od_param = 0.01,
                             bias_val = 1,
                             out_prop = 0.001,
                             out_mean = 0.5,
                             out_disp = 1/3,
                             non_mono_max = 2,
                             bound_bias = TRUE) {
  ## Check input ----------------------------------------------------
  if (!is.null(p1counts) & !is.null(p1size)) {
    assertthat::are_equal(length(p1counts), length(p1size), 1)
  }
  if (!is.null(p2counts) & !is.null(p2size)) {
    assertthat::are_equal(length(p2counts), length(p2size), 1)
  }


  if (!update_pgeno) {
    commit_num <- -1 ## parent_count is initialized at 0, so we just don't update parental genotypes if we set commit_num to be less than or equal to 0.
  }

  ## starting values ------------------------------------------
  weight_vec <- rep(out_prop, length = length(ocounts)) ## probability each point is an outlier.
  p1weight   <- out_prop ## probability parent 1 is an outlier
  p2weight   <- out_prop ## probability parent 2 is an outlier.

  non_mono   <- 0 ## counts number of times likelihood does not increase.

  llike_new <- -Inf

  index <- 1
  err <- tol + 1
  parental_count <- 0 ## counts number of consecutive times parental genotypes do not change.
  while ((index <= maxiter) & (err > tol)) {
    llike_old <- llike_new

    ## E-step ----------------------------------------------------------------------------
    ## Skip first iteration to get good estimates of other paramters before getting outlier weights.
    if (index > 1) {
      weight_vec <- get_out_prop(ocounts = ocounts, osize = osize,
                                 ploidy = ploidy, p1geno = p1geno,
                                 p2geno = p2geno, d = bias_val,
                                 eps = seq_error, tau = od_param,
                                 out_prop = out_prop, out_mean = out_mean,
                                 out_disp = out_disp)
      if (!is.null(p1counts) & !is.null(p1size)) {
        p1weight <- get_parent_outprop(pcounts = p1counts, psize = p1size,
                                       ploidy = ploidy, pgeno = p1geno,
                                       d = bias_val, eps = seq_error,
                                       tau = od_param, out_prop = out_prop,
                                       out_mean = out_mean, out_disp = out_disp)
      }
      if (!is.null(p2counts) & !is.null(p2size)) {
        p2weight <- get_parent_outprop(pcounts = p2counts, psize = p2size,
                                       ploidy = ploidy, pgeno = p2geno,
                                       d = bias_val, eps = seq_error,
                                       tau = od_param, out_prop = out_prop,
                                       out_mean = out_mean, out_disp = out_disp)
      }
    }


    ## Update out_prop -----------------------------------------------------------------
    if (update_outprop) {
      weight_sum  <- sum(weight_vec)
      weight_size <- length(weight_vec)
      if (!is.null(p1counts) & !is.null(p1size)) {
        weight_sum <- weight_sum + p1weight
        weight_size <- weight_size + 1
      }
      if (!is.null(p2counts) & !is.null(p2size)) {
        weight_sum <- weight_sum + p2weight
        weight_size <- weight_size + 1
      }
      out_prop <- weight_sum / weight_size
    }

    ## reparameterization --------------------------------------------------------------
    s   <- log(bias_val)
    ell <- log(seq_error / (1 - seq_error))
    r   <- log((1 - od_param) / od_param)
    parvec <- c(s, ell, r)

    ## update good ---------------------------------------------------------------------
    if (parental_count >= commit_num) {
      gout <- update_good(parvec = parvec, ocounts = ocounts, osize = osize,
                          weight_vec = 1 - weight_vec, ploidy = ploidy,
                          p1geno = p1geno, p2geno = p2geno,
                          p1counts = p1counts, p1size = p1size, p1weight = 1 - p1weight,
                          p2counts = p2counts, p2size = p2size, p2weight = 1 - p2weight,
                          bound_bias = bound_bias,
                          update_bias_val = update_bias_val,
                          update_seq_error = update_seq_error,
                          update_od_param = update_od_param)
    } else {
      gout <- update_good(parvec = parvec, ocounts = ocounts, osize = osize,
                          weight_vec = 1 - weight_vec, ploidy = ploidy,
                          p1counts = p1counts, p1size = p1size, p1weight = 1 - p1weight,
                          p2counts = p2counts, p2size = p2size, p2weight = 1 - p2weight,
                          bound_bias = bound_bias,
                          update_bias_val = update_bias_val,
                          update_seq_error = update_seq_error,
                          update_od_param = update_od_param)
    }

    bias_val  <- gout$bias
    seq_error <- gout$seq_error
    od_param  <- gout$od


    ## update parental_count if no change ---------------------------------------------
    if (p1geno == gout$p1geno & p2geno == gout$p2geno) {
      parental_count <- parental_count + 1
    } else {
      parental_count <- 0
    }
    p1geno <- gout$p1geno
    p2geno <- gout$p2geno

    ## update bad -------- ------------------------------------------------------------
    # if (out_disp < min_disp) {
    #   start_disp <- min_disp + 10 ^ -3
    # } else if (out_disp > 1 - 10 ^ -8) {
    #   start_disp <- 1 - 10 ^ -6
    # } else {
       start_disp <- out_disp
    # }
    # if (out_mean < 10 ^ -8) {
    #   out_mean <- 10 ^ -6
    # } else if (out_mean > 1 - 10 ^ -8) {
    #   out_mean <- 1 - 10 ^ -6
    # }
    if (update_outmean & update_outdisp) {
      oout <- stats::optim(par = c(out_mean, start_disp),
                           fn = out_obj_wrapp, gr = out_grad_wrapp,
                           ocounts = ocounts, osize = osize,
                           weight_vec = weight_vec, min_disp = min_disp,
                           p1counts = p1counts, p1size = p1size, p1weight = p1weight,
                           p2counts = p2counts, p2size = p2size, p2weight = p2weight,
                           method = "BFGS",
                           control = list(fnscale = -1))
      out_mean <- oout$par[1]
      out_disp <- oout$par[2]
    } else if (update_outdisp) {
      oout <- stats::optim(par = start_disp, fn = outlier_obj,
                           method = "Brent", control = list(fnscale = -1),
                           lower = min_disp, upper = 1 - 10 ^ -3,
                           p1counts = p1counts, p1size = p1size, p1weight = p1weight,
                           p2counts = p2counts, p2size = p2size, p2weight = p2weight,
                           ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                           out_mean = out_mean)
      out_disp <- oout$par
    } else if (update_outmean) {
      oout <- stats::optim(par = out_mean, fn = outlier_obj,
                           method = "Brent", control = list(fnscale = -1),
                           lower = 0, upper = 1,
                           p1counts = p1counts, p1size = p1size, p1weight = p1weight,
                           p2counts = p2counts, p2size = p2size, p2weight = p2weight,
                           ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                           out_disp = out_disp)
      out_mean <- oout$par
    }


    ## Calculate log-likelihood and update err and index -------
    llike_new <- obj_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy,
                               p1geno = p1geno, p2geno = p2geno, bias_val = bias_val,
                               seq_error = seq_error, od_param = od_param, outlier = TRUE,
                               out_prop = out_prop, out_mean = out_mean, out_disp = out_disp)
    if (!is.null(p1counts) & !is.null(p1size)) { ## add parent 1 if available
      llike_new <- llike_new + obj_parent(pcounts = p1counts, psize = p1size, ploidy = ploidy,
                                          pgeno = p1geno, bias_val = bias_val, seq_error = seq_error,
                                          od_param = od_param, outlier = TRUE, out_prop = out_prop,
                                          out_mean = out_mean, out_disp = out_disp)
    }
    if (!is.null(p2counts) & !is.null(p2size)) { ## add parent 2 if available
      llike_new <- llike_new + obj_parent(pcounts = p2counts, psize = p2size, ploidy = ploidy,
                                          pgeno = p2geno, bias_val = bias_val, seq_error = seq_error,
                                          od_param = od_param, outlier = TRUE, out_prop = out_prop,
                                          out_mean = out_mean, out_disp = out_disp)
    }

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
        warning("Non-monotone likelihood.")
        non_mono <- non_mono + 1
      }
      if (non_mono >= non_mono_max) {
        stop("Your likelihood is not increasing.\n")
      }
    }

    index <- index + 1
  }

  if (abs(out_disp - min_disp) < 10 ^ -6) {
    warning("out_disp estimated at minimum. Be wary of prob_ok results.\n")
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
#' @inheritParams updog_update_all
#'
#' @author David Gerard
#'
#' @export
#'
updog_vanilla <- function(ocounts, osize, ploidy,
                          p1counts = NULL,
                          p1size = NULL,
                          p2counts = NULL,
                          p2size = NULL,
                          commit_num = 4,
                          min_disp = 0,
                          print_val = FALSE,
                          update_outmean = FALSE,
                          update_outdisp = FALSE,
                          update_outprop = TRUE,
                          update_bias_val = TRUE,
                          update_seq_error = TRUE,
                          update_od_param = TRUE,
                          update_pgeno = TRUE,
                          p1geno = 3,
                          p2geno = 3,
                          seq_error = 0.01,
                          od_param = 0.01,
                          bias_val = 1,
                          out_prop = 0.001,
                          out_mean = 0.5,
                          out_disp = 1/3,
                          non_mono_max = 2,
                          bound_bias = TRUE) {
  ## Check input -----------------------------------------------------------------------
  assertthat::assert_that(is.logical(print_val))
  assertthat::assert_that(is.logical(update_outmean))
  assertthat::assert_that(is.logical(update_outdisp))
  assertthat::assert_that(is.logical(update_pgeno))
  assertthat::are_equal(length(ocounts), length(osize))
  assertthat::assert_that(all(ocounts <= osize))
  assertthat::assert_that(all(ocounts >= 0))
  assertthat::are_equal(ploidy %% 2, 0)
  assertthat::assert_that(ploidy > 0)
  assertthat::are_equal(length(ploidy), length(print_val), length(update_outdisp), length(update_outmean), length(min_disp), 1)
  assertthat::assert_that(min_disp >= 0, min_disp < 1)

  if (!is.null(p1counts) & !is.null(p1size)) {
    assertthat::are_equal(length(p1counts), length(p1size))
    if (length(p1counts) > 1) {
      message("Aggregating parent 1 counts into one sample.\nIncorporating variability between samples of the same parent is currently not supported.")
      p1counts <- sum(p1counts)
      p1size   <- sum(p1size)
    }
    assertthat::assert_that(p1counts >= 0, p1counts <= p1size)
  } else if (!is.null(p1counts) | !is.null(p1size)) {
    warning("You can't just specify one of p1counts and p1size. Ignoring parent 1 data.")
    p1counts <- NULL
    p1size <- NULL
  }

  if (!is.null(p2counts) & !is.null(p2size)) {
    assertthat::are_equal(length(p2counts), length(p2size))
    if (length(p2counts) > 1) {
      message("Aggregating parent 2 counts into one sample.\nIncorporating variability between samples of the same parent is currently not supported.")
      p2counts <- sum(p2counts)
      p2size   <- sum(p2size)
    }
    assertthat::assert_that(p2counts >= 0, p2counts <= p2size)
  } else if (!is.null(p2counts) | !is.null(p2size)) {
    warning("You can't just specify one of p2counts and p2size. Ignoring parent 2 data.")
    p2counts <- NULL
    p2size <- NULL
  }

  ## Get the best parameters
  parout <- updog_update_all(ocounts, osize, ploidy,
                             p1counts = p1counts, p1size = p1size,
                             p2counts = p2counts, p2size = p2size,
                             print_val = print_val, commit_num = commit_num,
                             min_disp = min_disp, update_outprop = update_outprop,
                             update_outmean = update_outmean, update_outdisp = update_outdisp,
                             update_pgeno = update_pgeno,
                             p1geno = p1geno, p2geno = p2geno, seq_error = seq_error,
                             od_param = od_param, bias_val = bias_val, out_prop = out_prop,
                             out_mean = out_mean, out_disp = out_disp, non_mono_max = non_mono_max,
                             update_bias_val = update_bias_val,
                             update_seq_error = update_seq_error,
                             update_od_param = update_od_param)

  if (parout$od_param < 10 ^ -13) { ## fix for getting some weird postmat's
    parout$od_param <- 0
  }

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


