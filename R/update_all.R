
#' wrapper for \code{\link{obj_offspring_weights_reparam}}
#' @inheritParams obj_offspring_vec
#' @param parvec A vector of three elements, s, ell, and r. We have s = log(bias_val) = log(d),
#' ell = logit(seq_error) = logit(eps), and r = - logit(od_param) = - logit(tau).
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
#' @param seq_error_mean The mean of the logit-normal prior on the sequencing error rate, which corresponds
#'     to \code{parvec[2]}. Set \code{seq_error_sd = Inf} to have no penalty on the sequencing error rate.
#'     The default is -4.7, which roughly corresponds to a mean sequencing error value of 0.009.
#' @param seq_error_sd The standard deviation of the logit-normal prior on the sequencing error rate, which
#'     corresponds to \code{parvec[2]}. The default is 1, which at three standard deviations is about a sequencing
#'     error rate of 0.15.
#' @param bias_val_mean The prior mean on the log of \code{bias_val} (corresponding to \code{parvec[1]}).
#'     Set \code{bias_val_sd = Inf} to have no penalty on the bias parameter.
#' @param bias_val_sd The prior standard deviation on the log of \code{bias_val}
#'     (corresponding to \code{parvec[1]}).
#'     Set \code{bias_val_sd = Inf} to have no penalty on the bias parameter.
#' @param p1geno The genotype of parent 1 if \code{model = "f1"} or \code{model = "s1"}.
#' @param p2geno The genotype of parent 2 if \code{model = "f1"}. This needs to be null if \code{model = "s1"}
#' @param model The model for the genotype distribution. Should we assume an F1 population (\code{"f1"}),
#'     an S1 population (\code{"s1"}),
#'     Hardy-Weinberg equilibrium (\code{"hw"}), or a uniform distribution (\code{"uniform"})?
#' @param allele_freq The allele-frequency if \code{model = "hw"}
#' @param is_seq_ninf Hacky way to get around fact that optim won't allow \code{-Inf} in
#'     \code{par}, even if \code{gr} will always be zero for that parameter.
#'     How dare they not anticipate my unique situation!
#'
#' @author David Gerard
#'
obj_wrapp_all <- function(parvec, ocounts, osize, weight_vec,
                          ploidy, p1geno, p2geno, allele_freq = 0.5,
                          p1counts = NULL, p1size = NULL, p1weight = NULL,
                          p2counts = NULL, p2size = NULL, p2weight = NULL,
                          bound_bias = FALSE,
                          bound_od = FALSE,
                          update_bias_val = TRUE,
                          update_seq_error = TRUE,
                          update_od_param = TRUE,
                          seq_error_mean = -4.7,
                          seq_error_sd = 1,
                          bias_val_mean = 0,
                          bias_val_sd = 0.7,
                          model = c("f1", "s1", "hw", "uniform"),
                          is_seq_ninf = FALSE) {

  assertthat::assert_that(is.logical(is_seq_ninf))
  if (is_seq_ninf) {
    ell <- -Inf
    if (update_seq_error) {
      stop("can't have seq_error == 0 and update_seq_error == TRUE")
    }
  } else {
    ell <- parvec[2]
  }

  model <- match.arg(model)
  eps <- expit(ell)
  if (model == "s1") {
    model <- "f1"
    assertthat::are_equal(p1geno, p2geno)
  }

  ## get genotype frequencies --------------------------------------------------------
  prob_geno <- get_prob_geno(ploidy = ploidy, model = model, p1geno = p1geno, p2geno = p2geno, allele_freq = allele_freq)

  ## Normal value ------------------------------------------------------------------------
  val <- obj_offspring_weights_reparam(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                                        ploidy = ploidy, prob_geno = prob_geno,
                                        s = parvec[1], ell = ell, r = parvec[3])

  if (!is.null(p1counts) & !is.null(p1size) & !is.null(p1weight)) { ## add contribution from parent 1
    val <- val + obj_parent_reparam(pcounts = p1counts, psize = p1size, ploidy = ploidy,
                                    pgeno = p1geno, s = parvec[1], ell = ell,
                                    r = parvec[3], weight = p1weight)
  }
  if (!is.null(p2counts) & !is.null(p2size) & !is.null(p2weight)) { ## add contribution from parent 2
    val <- val + obj_parent_reparam(pcounts = p2counts, psize = p2size, ploidy = ploidy,
                                    pgeno = p2geno, s = parvec[1], ell = ell,
                                    r = parvec[3], weight = p2weight)
  }

  ## seq_error penalty ---
  if (seq_error_mean != -Inf) { ## to specify a zero sequencing error rate
    if (ell == -Inf) {
      stop("seq_error_mean != -Inf but seq_error = 0, this is not allowed")
    }
    val <- val - (ell - seq_error_mean) ^ 2 / (2 * seq_error_sd ^ 2)
  }

  ## bias_val penalty ---
  val <- val - (parvec[1] - bias_val_mean) ^ 2 / (2 * bias_val_sd ^ 2)

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
grad_wrapp_all <- function(parvec, ocounts, osize, weight_vec, ploidy, p1geno, p2geno, allele_freq = 0.5,
                           p1counts = NULL, p1size = NULL, p1weight = NULL,
                           p2counts = NULL, p2size = NULL, p2weight = NULL,
                           bound_bias = FALSE,
                           update_bias_val = TRUE,
                           update_seq_error = TRUE,
                           update_od_param = TRUE,
                           seq_error_mean = -4.7,
                           seq_error_sd = 1,
                           bias_val_mean = 0,
                           bias_val_sd = 0.7,
                           model = c("f1", "s1", "hw", "uniform"),
                           is_seq_ninf = FALSE) {

  assertthat::assert_that(is.logical(is_seq_ninf))
  if (is_seq_ninf) {
    ell <- -Inf
    if (update_seq_error) {
      stop("can't have seq_error == 0 and update_seq_error == TRUE")
    }
  } else {
    ell <- parvec[2]
  }


  model <- match.arg(model)
  if (model == "s1") {
    model <- "f1"
    assertthat::are_equal(p1geno, p2geno)
  }

  ## get genotype frequencies --------------------------------------------------------
  prob_geno <- get_prob_geno(ploidy = ploidy, model = model, p1geno = p1geno, p2geno = p2geno, allele_freq = allele_freq)

  gout <- grad_offspring_weights(ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                                 ploidy = ploidy, prob_geno = prob_geno, s = parvec[1],
                                 ell = ell, r = parvec[3])

  if (!is.null(p1counts) & !is.null(p1size) & !is.null(p1weight)) {
    gout <- gout + grad_parent_reparam(pcounts = p1counts, psize = p1size,
                                       ploidy = ploidy, pgeno = p1geno,
                                       s = parvec[1], ell = ell,
                                       r = parvec[3], weight = p1weight)
  }
  if (!is.null(p2counts) & !is.null(p2size) & !is.null(p2weight)) {
    gout <- gout + grad_parent_reparam(pcounts = p2counts, psize = p2size,
                                       ploidy = ploidy, pgeno = p2geno,
                                       s = parvec[1], ell = ell,
                                       r = parvec[3], weight = p2weight)
  }

  ## sequencing error penalty -----------------------------
  if (seq_error_mean != -Inf) {
    gout[2] <- gout[2] - (ell - seq_error_mean) / (seq_error_sd ^ 2)
  } else {
    if (ell != -Inf & seq_error_sd != Inf) {
      gout[2] <- NA
    }
  }


  ## Bias parameter penalty -------------------------------
  gout[1] <- gout[1] - (parvec[1] - bias_val_mean) / (bias_val_sd ^ 2)

  ## set to 0 if not to be updated -------------------------
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
#' @param verbose A logical. Should we write more output \code{TRUE} or not \code{FALSE}?
#'
#' @author David Gerard
#'
update_good <- function(parvec, ocounts, osize, weight_vec, ploidy,
                        p1counts = NULL, p1size = NULL, p1weight = NULL,
                        p2counts = NULL, p2size = NULL, p2weight = NULL,
                        p1geno = NULL, p2geno = NULL, allele_freq = 0.5,
                        bound_bias = FALSE,
                        update_bias_val = TRUE,
                        update_seq_error = TRUE,
                        update_od_param = TRUE,
                        seq_error_mean = -4.7,
                        seq_error_sd = 1,
                        bias_val_mean = 0,
                        bias_val_sd = 0.7,
                        model = c("f1", "s1", "hw", "uniform"),
                        verbose = FALSE) {

  ## deal with ell = -Inf
  if (parvec[2] == -Inf) {
    is_seq_ninf <- TRUE
    if (update_seq_error) {
      stop("cannot have both update_seq_error = TRUE and seq_error = 0")
    }
    if (seq_error_mean != -Inf) {
      stop("cannot have both seq_error_mean != -Inf and seq_error = 0")
    }
  } else {
    is_seq_ninf <- FALSE
  }


  model <- match.arg(model)
  best_par <- parvec
  if ((is.null(p1geno) | is.null(p2geno)) & (model == "f1" | model == "s1")) { ## try all p1geno/p2geno combos --
    best_p1 <- 0
    best_p2 <- 0
    best_llike <- -Inf

    for (p1geno in 0:ploidy) {
      if (model == "f1") {
        if (is.null(p1counts) | is.null(p2counts) | is.null(p1size) | is.null(p2size)) {
          possible_p2geno_vec <- 0:p1geno ## to deal with identifiability, constain p2geno to be less than or equal to p1geno
        } else {
          possible_p2geno_vec <- 0:ploidy
        }
      } else if (model == "s1") {
        possible_p2geno_vec <- p1geno
      }
      for (p2geno in possible_p2geno_vec) {

        if (is_seq_ninf & ((p1geno == 0 & p2geno == 0) | (p1geno == ploidy & p2geno == ploidy))) {
          ## skip iteration
        } else {
          ## get genotype frequencies --------------------------------------------------------
          prob_geno <- get_prob_geno(ploidy = ploidy, model = model, p1geno = p1geno, p2geno = p2geno, allele_freq = allele_freq)

          ## If start position has -Inf, start somewhere else
          val <- obj_offspring_weights_reparam(ocounts = ocounts, osize = osize,
                                               weight_vec = weight_vec,
                                               ploidy = ploidy, prob_geno = prob_geno,
                                               s = parvec[1], ell = parvec[2],
                                               r = parvec[3])
          if (val == -Inf) {
            start_vec <- c(0, -4.5, 4.5) ## about (1, 0.01, 0.01) for (bias_val, seq_error, od_param)
          } else {
            start_vec <- parvec
          }

          if (is_seq_ninf) { ## deals with ell = -Inf
            start_vec[2] <- 0
          }
          errout <- try({ ## sometimes L-BFGS-B fails but NM works, but is slower.
            oout <- stats::optim(par = start_vec, fn = obj_wrapp_all, gr = grad_wrapp_all,
                                 hessian = TRUE,
                                 upper = c(Inf, 0, Inf),
                                 method = "L-BFGS-B",
                                 control = list(fnscale = -1, maxit = 1000),
                                 ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                                 ploidy = ploidy, p1geno = p1geno, p2geno = p2geno,
                                 p1counts = p1counts, p1size = p1size, p1weight = p1weight,
                                 p2counts = p2counts, p2size = p2size, p2weight = p2weight,
                                 bound_bias = bound_bias,
                                 update_bias_val = update_bias_val,
                                 update_seq_error = update_seq_error,
                                 update_od_param = update_od_param,
                                 seq_error_mean = seq_error_mean,
                                 seq_error_sd = seq_error_sd,
                                 bias_val_mean = bias_val_mean,
                                 bias_val_sd = bias_val_sd,
                                 model = model,
                                 is_seq_ninf = is_seq_ninf)
          }, TRUE)
          if ("try-error" %in% class(errout)) {
            oout <- stats::optim(par = start_vec, fn = obj_wrapp_all, gr = grad_wrapp_all,
                                 hessian = TRUE,
                                 method = "Nelder-Mead",
                                 control = list(fnscale = -1, maxit = 1000),
                                 ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                                 ploidy = ploidy, p1geno = p1geno, p2geno = p2geno,
                                 p1counts = p1counts, p1size = p1size, p1weight = p1weight,
                                 p2counts = p2counts, p2size = p2size, p2weight = p2weight,
                                 bound_bias = bound_bias,
                                 update_bias_val = update_bias_val,
                                 update_seq_error = update_seq_error,
                                 update_od_param = update_od_param,
                                 seq_error_mean = seq_error_mean,
                                 seq_error_sd = seq_error_sd,
                                 bias_val_mean = bias_val_mean,
                                 bias_val_sd = bias_val_sd,
                                 model = model,
                                 is_seq_ninf = is_seq_ninf)
          }

          if (is_seq_ninf) {## deals with ell = -Inf
            oout$par[2] <- -Inf
          }

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
    }
    allele_freq <- allele_freq ## not updated
  } else if (model == "f1" | model == "s1" | model == "uniform") { ## run if we know prob_geno ahead of time ----
    ## get genotype frequencies --------------------------------------------------------
    prob_geno <- get_prob_geno(ploidy = ploidy, model = model, p1geno = p1geno, p2geno = p2geno, allele_freq = allele_freq)

    ## If start position has -Inf, start somewhere else
    val <- obj_offspring_weights_reparam(ocounts = ocounts, osize = osize,
                                         weight_vec = weight_vec,
                                         ploidy = ploidy, prob_geno = prob_geno,
                                         s = parvec[1], ell = parvec[2],
                                         r = parvec[3])
    if (val == -Inf) {
      start_vec <- c(0, -4.5, 4.5) ## about (1, 0.01, 0.01) for (bias_val, seq_error, od_param)
    } else {
      start_vec <- parvec
    }


    if (is_seq_ninf) { ## deals with ell = -Inf
      start_vec[2] <- 0
    }
    oout <- stats::optim(par = start_vec, fn = obj_wrapp_all, gr = grad_wrapp_all,
                         hessian = TRUE,
                         ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                         ploidy = ploidy, p1geno = p1geno, p2geno = p2geno,
                         p1counts = p1counts, p1size = p1size, p1weight = p1weight,
                         p2counts = p2counts, p2size = p2size, p2weight = p2weight,
                         method = "BFGS",
                         control = list(fnscale = -1, maxit = 1000),
                         bound_bias = bound_bias,
                         update_bias_val = update_bias_val,
                         update_seq_error = update_seq_error,
                         update_od_param = update_od_param,
                         seq_error_mean = seq_error_mean,
                         seq_error_sd = seq_error_sd,
                         bias_val_mean = bias_val_mean,
                         bias_val_sd = bias_val_sd,
                         model = model,
                         is_seq_ninf = is_seq_ninf)
    if (is_seq_ninf) { ## deals with ell = -Inf
      oout$par[2] <- -Inf
    }

    best_p1 <- p1geno
    best_p2 <- p2geno
    best_llike <- oout$value
    best_par   <- oout$par
    allele_freq <- allele_freq ## not updated
  } else if (model == "hw") { # iteratively estimate allele frequency and updog parameters ---
    start_vec <- parvec
    ## update a couple times
    if (verbose) {
      cat("updating allele freq and params\n\n")
    }
    obj_oaf_par <- -Inf
    for (index in 1:3) {
      oaf_out <- stats::optim(par = allele_freq, fn = obj_wrapp_all,
                              method = "Brent", lower = 10 ^ -6,
                              upper = 1 - 10 ^ -6,
                              control = list(fnscale = -1),
                              parvec = start_vec, ocounts = ocounts,
                              osize = osize, weight_vec = weight_vec,
                              ploidy = ploidy, p1geno = 0, p2geno = 0,
                              p1counts = NULL, p1size = NULL, p1weight = NULL,
                              p2counts = NULL, p2size = NULL, p2weight = NULL,
                              update_bias_val = update_bias_val,
                              update_seq_error = update_seq_error,
                              update_od_param = update_od_param,
                              seq_error_mean = seq_error_mean,
                              seq_error_sd = seq_error_sd,
                              bias_val_mean = bias_val_mean,
                              bias_val_sd = bias_val_sd,
                              model = model,
                              is_seq_ninf = is_seq_ninf)

      if (oaf_out$value < obj_oaf_par - 10^-10) {
        stop("likelihood not increasing.")
      } else {
        obj_oaf_par <- oaf_out$value
      }

      allele_freq <- oaf_out$par ## new allele_freq

      if (verbose) {
        cat("New allele_freq:", allele_freq, "\n")
        cat("like:", oaf_out$value, "\n\n")
      }

      if (is_seq_ninf) { ## deals with ell = -Inf
        start_vec[2] <- 0
      }

      errout <- try({ ## sometimes L-BFGS-B fails but NM works, but is slower.
        oout <- stats::optim(par = start_vec, fn = obj_wrapp_all, gr = grad_wrapp_all,
                             hessian = TRUE,
                             method = "L-BFGS-B",
                             upper = c(Inf, 0, Inf), ## seq error can be at most 50%
                             control = list(fnscale = -1),
                             ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                             ploidy = ploidy, p1geno = 0, p2geno = 0,
                             allele_freq = allele_freq,
                             p1counts = NULL, p1size = NULL, p1weight = NULL,
                             p2counts = NULL, p2size = NULL, p2weight = NULL,
                             bound_bias = bound_bias,
                             update_bias_val = update_bias_val,
                             update_seq_error = update_seq_error,
                             update_od_param = update_od_param,
                             seq_error_mean = seq_error_mean,
                             seq_error_sd = seq_error_sd,
                             bias_val_mean = bias_val_mean,
                             bias_val_sd = bias_val_sd,
                             model = model,
                             is_seq_ninf = is_seq_ninf)
      }, silent = TRUE)
      if ("try-error" %in% class(errout)) {
        oout <- stats::optim(par = start_vec, fn = obj_wrapp_all, gr = grad_wrapp_all,
                             hessian = TRUE,
                             method = "Nelder-Mead",
                             control = list(fnscale = -1),
                             ocounts = ocounts, osize = osize, weight_vec = weight_vec,
                             ploidy = ploidy, p1geno = 0, p2geno = 0,
                             allele_freq = allele_freq,
                             p1counts = NULL, p1size = NULL, p1weight = NULL,
                             p2counts = NULL, p2size = NULL, p2weight = NULL,
                             bound_bias = bound_bias,
                             update_bias_val = update_bias_val,
                             update_seq_error = update_seq_error,
                             update_od_param = update_od_param,
                             seq_error_mean = seq_error_mean,
                             seq_error_sd = seq_error_sd,
                             bias_val_mean = bias_val_mean,
                             bias_val_sd = bias_val_sd,
                             model = model,
                             is_seq_ninf = is_seq_ninf)
      }
      if (is_seq_ninf) { ## deals with ell = -Inf
        oout$par[2] <- -Inf
      }

      if (oout$value < obj_oaf_par - 10 ^ -10) {
        stop("likelihood not increasing")
      } else {
        obj_oaf_par <- oout$value
      }

      if (verbose) {
        cat("New (s, ell, r):", oout$par, "\n")
        cat("like:", oout$value, "\n\n")
      }

      start_vec <- oout$par
    }
    best_llike <- oout$value
    best_par   <- oout$par
    best_p1 <- 0
    best_p2 <- 0

    if (verbose) {
      cat("done updating allele freq and params\n\n")
    }
  } else {
    stop("check update_good because corner case observed.")
  }

  return_list             <- list()
  return_list$p1geno      <- best_p1
  return_list$p2geno      <- best_p2
  return_list$bias        <- exp(best_par[1])
  return_list$seq_error   <- expit(best_par[2])
  return_list$od          <- 1 / (1 + exp(best_par[3]))
  return_list$allele_freq <- allele_freq
  return_list$llike       <- best_llike
  return_list$hessian     <- oout$hessian
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
#' @inheritParams updog_old
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
#' @param seq_error_mean The mean of the logit-normal prior on the sequencing error rate, which corresponds
#'     to \code{parvec[2]}. Set \code{seq_error_sd = Inf} to have no penalty on the sequencing error rate.
#'     The default is -4.7, which roughly corresponds to a mean sequencing error value of 0.009. If you want
#'     to constain \code{seq_error} to be zero, you need to set \code{update_seq_error = FALSE},
#'     \code{seq_error = 0}, and\code{seq_error_mean = -Inf}.
#' @param seq_error_sd The standard deviation of the logit-normal prior on the sequencing error rate, which
#'     corresponds to \code{parvec[2]}. The default is 1, which at three standard deviations is about a sequencing
#'     error rate of 0.15. Set \code{seq_error_sd = Inf} to have no penalty on the sequencing error rate.
#' @param bias_val_mean The prior mean on the log of \code{bias_val} (corresponding to \code{parvec[1]}).
#'     Set \code{bias_val_sd = Inf} to have no penalty on the bias parameter.
#' @param bias_val_sd The prior standard deviation on the log of \code{bias_val}
#'     (corresponding to \code{parvec[1]}).
#'     Set \code{bias_val_sd = Inf} to have no penalty on the bias parameter.
#' @param allele_freq The allele frequency if \code{model = "hw"}.
#' @param model The model for the genotype distribution. Do we assume an
#'    F1 population (\code{"f1"}), an S1 population (\code{"s1"}), Hardy-Weinberg equilibrium (\code{"hw"}),
#'    or a uniform distribution (\code{"uniform"}).
#' @param verbose A logical. Should we output a lot more (\code{TRUE}) or not (\code{FALSE})?
#'
#' @return A list of the following elements
#'     \itemize{
#'         \item{\code{bias_val}:}{ The estimated bias parameter.}
#'         \item{\code{seq_error}:}{ The estimated sequencing error rate.}
#'         \item{\code{od_param}:}{ The estimated overdispersion parameter.}
#'         \item{\code{p1geno}:}{ The estimated genotype of one parent.}
#'         \item{\code{p1geno}:}{ The estimated genotype of the other parent.}
#'         \item{\code{out_prop}:}{ The estimated proportion of points that are outliers.}
#'         \item{\code{out_mean}:}{ The estimated mean of the outlier distribution.}
#'         \item{\code{out_disp}:}{ The estimated overdispersion parameter of the outlier distribution.}
#'         \item{\code{prob_out}:}{ A vector. Each element of which is the posterior probability that a point is an outlier.}
#'         \item{\code{allele_freq}:}{ The estimated allele-frequency of the reference allele.}
#'         \item{\code{p1_prob_out}:}{ The posterior probability that parent 1 is an outlier.}
#'         \item{\code{p2_prob_out}:}{ The posterior probability that parent 2 is an outlier.}
#'         \item{\code{num_iter}:}{ The number of iterations the optimization program was run.}
#'         \item{\code{convergence}:}{ 1 is we reached \code{maxiter} and 0 otherwise.}
#'         \item{\code{llike}:}{ The final log-likelihood of the estimates.}
#'         \item{\code{hessian}:}{ The Fisher information under the parameterization (s, ell, r), where s = log(bias_val) = log(d),
#' ell = logit(seq_error) = logit(eps), and r = - logit(od_param) =  - logit(tau).}
#'     }
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
                             p1geno = round(ploidy / 2),
                             p2geno = round(ploidy / 2),
                             seq_error = 0.01,
                             od_param = 0.01,
                             bias_val = 1,
                             out_prop = 0.001,
                             out_mean = 0.5,
                             out_disp = 1/3,
                             non_mono_max = 2,
                             bound_bias = FALSE,
                             seq_error_mean = -4.7,
                             seq_error_sd = 1,
                             bias_val_mean = 0,
                             bias_val_sd = 0.7,
                             allele_freq = 0.5,
                             model = c("f1", "s1", "hw", "uniform"),
                             verbose = FALSE) {
  ## Check input ----------------------------------------------------
  if (!is.null(p1counts) & !is.null(p1size)) {
    assertthat::are_equal(length(p1counts), length(p1size), 1)
  }
  if (!is.null(p2counts) & !is.null(p2size)) {
    assertthat::are_equal(length(p2counts), length(p2size), 1)
  }

  model <- match.arg(model)
  assertthat::assert_that(allele_freq > 0, allele_freq < 1)

  if (model == "hw" | model == "uniform") {
    if (!is.null(p1counts) | !is.null(p1size) | !is.null(p2counts) | !is.null(p2size)) {
      warning("p1counts, p1size, p2counts, p2size being ignored because model != 'f1'")
      p1counts <- NULL
      p2counts <- NULL
      p1size   <- NULL
      p2size   <- NULL
    }
  }
  if (model == "s1" & (!is.null(p2counts) | !is.null(p2size))) {
    stop("when model = 's1', `p2counts` and `p2size` must be NULL.\nYou can place parental info into `p1counts` and `p1size`")
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
      ## get genotype frequencies --------------------------------------------------------
      prob_geno <- get_prob_geno(ploidy = ploidy, model = model, p1geno = p1geno, p2geno = p2geno, allele_freq = allele_freq)
      weight_vec <- get_out_prop(ocounts = ocounts, osize = osize,
                                 ploidy = ploidy, prob_geno = prob_geno,
                                 d = bias_val,
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
    if (parental_count >= commit_num | (model != "f1" & model != "s1")) {
      gout <- update_good(parvec = parvec, ocounts = ocounts, osize = osize,
                          weight_vec = 1 - weight_vec, ploidy = ploidy,
                          p1geno = p1geno, p2geno = p2geno,
                          p1counts = p1counts, p1size = p1size, p1weight = 1 - p1weight,
                          p2counts = p2counts, p2size = p2size, p2weight = 1 - p2weight,
                          bound_bias = bound_bias,
                          update_bias_val = update_bias_val,
                          update_seq_error = update_seq_error,
                          update_od_param = update_od_param,
                          seq_error_mean = seq_error_mean,
                          seq_error_sd = seq_error_sd,
                          bias_val_mean = bias_val_mean,
                          bias_val_sd = bias_val_sd,
                          allele_freq = allele_freq,
                          model = model, verbose = verbose)
    } else {
      gout <- update_good(parvec = parvec, ocounts = ocounts, osize = osize,
                          weight_vec = 1 - weight_vec, ploidy = ploidy,
                          p1counts = p1counts, p1size = p1size, p1weight = 1 - p1weight,
                          p2counts = p2counts, p2size = p2size, p2weight = 1 - p2weight,
                          bound_bias = bound_bias,
                          update_bias_val = update_bias_val,
                          update_seq_error = update_seq_error,
                          update_od_param = update_od_param,
                          seq_error_mean = seq_error_mean,
                          seq_error_sd = seq_error_sd,
                          bias_val_mean = bias_val_mean,
                          bias_val_sd = bias_val_sd,
                          model = model, verbose = verbose)
    }

    bias_val    <- gout$bias
    seq_error   <- gout$seq_error
    od_param    <- gout$od
    allele_freq <- gout$allele_freq


    ## update parental_count if no change ---------------------------------------------
    if (model == "s1") {
      stopifnot(gout$p1geno == gout$p2geno)
    }
    if (p1geno == gout$p1geno & p2geno == gout$p2geno) {
      parental_count <- parental_count + 1
    } else {
      parental_count <- 0
    }
    p1geno <- gout$p1geno
    p2geno <- gout$p2geno

    ## update bad -------- ------------------------------------------------------------
    start_disp <- out_disp
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
    ## get genotype frequencies --------------------------------------------------------
    prob_geno <- get_prob_geno(ploidy = ploidy, model = model, p1geno = p1geno, p2geno = p2geno, allele_freq = allele_freq)
    llike_new <- obj_offspring(ocounts = ocounts, osize = osize, ploidy = ploidy,
                               prob_geno = prob_geno, bias_val = bias_val,
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

    ## Add sequencing error penalty ---
    if (seq_error != 0 & seq_error_mean != -Inf) {
      llike_new <- llike_new - (log(seq_error / (1 - seq_error)) - seq_error_mean) ^ 2 / (2 * seq_error_sd ^ 2)
    }

    ## add bias_val penalty ---
    llike_new <- llike_new - (log(bias_val) - bias_val_mean) ^ 2 / (2 * bias_val_sd ^ 2)

    ## Calculate error ----
    err <- abs(llike_new - llike_old)

    if (print_val) {
      cat("    Log-Likelihood:", llike_new, "\n")
      cat("Parental Genotypes:", gout$p1geno, gout$p2geno, "\n")
      cat("  Allele Frequency:", gout$allele_freq, "\n")
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
  return_list$allele_freq <- allele_freq
  if (!is.null(p1counts) & !is.null(p1size)) {
    return_list$p1_prob_out <- p1weight
  }
  if (!is.null(p2counts) & !is.null(p2size)) {
    return_list$p2_prob_out <- p2weight
  }
  return_list$num_iter    <- index
  return_list$convergence <- (index >= maxiter) * 1
  return_list$llike       <- llike_new
  return_list$hessian     <- gout$hessian
  return_list$log_bias           <- log(return_list$bias_val)
  return_list$logit_seq_error    <- log(return_list$seq_error / (1 - return_list$seq_error))
  return_list$neg_logit_od_param <- log((1 - return_list$od_param) / return_list$od_param)

  ## get covariance matrix of updated parameters
  if (update_bias_val | update_seq_error | update_od_param) {
    keep_vec <- c(update_bias_val, update_seq_error, update_od_param)
    return_list$covmat <- matrix(NA, nrow = 3, ncol = 3)
    ## Check to see if any parameters are excessively close to the boundary.
    which_boundary <- c(return_list$bias_val < .Machine$double.eps,
                        (return_list$seq_error < .Machine$double.eps) | (return_list$seq_error > 1 - .Machine$double.eps),
                        (return_list$od_param < .Machine$double.eps) | (return_list$od_param > 1 - .Machine$double.eps))
    try({
      return_list$covmat[keep_vec & !which_boundary, keep_vec & !which_boundary] <- -1 * solve(return_list$hessian[keep_vec & !which_boundary, keep_vec & !which_boundary])
    }, TRUE) ## covmat is NA otherwise.
  } else {
    return_list$covmat <- NULL
  }

  return(return_list)
}

#' Calculate the negative inverse hessian of the parameters from
#' \code{\link{obj_offspring_reparam}}.
#'
#' This is the covariance matrix if evaluated at the MLE's of these parameters.
#'
#' @inheritParams obj_offspring_vec
#' @inheritParams updog_old
#' @inheritParams updog_update_all
#'
#' @author David Gerard
get_cov_mle <- function(ocounts, osize, ploidy, prob_geno, bias_val,
                        seq_error, od_param, out_prop,
                        p1counts = NULL, p1size = NULL, p1geno = NULL,
                        p2counts = NULL, p2size = NULL, p2geno = NULL,
                        bias_val_mean = 0, bias_val_sd = 0.7,
                        seq_error_mean = -4.7, seq_error_sd = 1,
                        model = c("f1", "s1", "hw", "uniform")) {
  model    <- match.arg(model)
  s        <- log(bias_val)
  ell      <- log(seq_error / (1 - seq_error))
  r        <- log((1 - od_param) / od_param)
  logitout <- log(out_prop / (1 - out_prop))
  parvec <- c(s, ell, r, logitout)
  hessout <- stats::optimHess(par = parvec, fn = fn_cov_mle,
                       ocounts = ocounts, osize = osize, ploidy = ploidy,
                       prob_geno = prob_geno, model = model,
                       p1counts = p1counts, p1size = p1size, p1geno = p1geno,
                       p2counts = p2counts, p2size = p2size, p2geno = p2geno,
                       bias_val_mean = bias_val_mean, bias_val_sd = bias_val_sd,
                       seq_error_mean = seq_error_mean, seq_error_sd = seq_error_sd)
  fisherinfo <- -1 * solve(hessout)
  return(list(coefficients = parvec, cov = fisherinfo))
}

#' Wrapper for \code{\link{obj_offspring_reparam}}.
#'
#' This exists to calculate the Hessian of the parameters at the end of
#' \code{\link{updog_update_all}}.
#'
#' @inheritParams obj_offspring_vec
#' @inheritParams updog_old
#' @inheritParams updog_update_all
#' @param par A numeric vector of length 4. The elements are \code{s},
#'     \code{ell}, \code{r}, and the logit of \code{out_prop} from
#'     \code{\link{obj_offspring_reparam}}.
#'
#' @author David Gerard
#'
fn_cov_mle <- function(par, ocounts, osize, ploidy, prob_geno, p1counts = NULL,
                       p1size = NULL, p2counts = NULL, p2size = NULL, p1geno = NULL,
                       p2geno = NULL,
                       bias_val_mean = 0, bias_val_sd = 0.7,
                       seq_error_mean = -4.7, seq_error_sd = 1,
                       model = c("f1", "s1", "hw", "uniform")) {
  model <- match.arg(model)
  out_prop <- expit(par[4])
  llike_new <- obj_offspring_reparam(ocounts = ocounts, osize = osize, ploidy = ploidy,
                                     prob_geno = prob_geno, s = par[1], ell = par[2], r = par[3],
                                     outlier = TRUE, out_prop = out_prop)
  if (!is.null(p1counts) & !is.null(p1size) & !is.null(p1geno) &
      (model == "f1" | model == "s1")) {
    llike_new <- llike_new + obj_parent_reparam(pcounts = p1counts, psize = p1size,
                                                ploidy = ploidy, pgeno = p1geno,
                                                s = par[1], ell = par[2],
                                                r = par[3], weight = 1,
                                                outlier = TRUE, out_prop = out_prop)
  }
  if (!is.null(p2counts) & !is.null(p2size) & !is.null(p2geno) & model == "f1") {
    llike_new <- llike_new + obj_parent_reparam(pcounts = p2counts, psize = p2size,
                                                ploidy = ploidy, pgeno = p2geno,
                                                s = par[1], ell = par[2],
                                                r = par[3], weight = 1,
                                                outlier = TRUE, out_prop = out_prop)
  }

  ## Add sequencing error penalty ---
  llike_new <- llike_new - (par[2] - seq_error_mean) ^ 2 / (2 * seq_error_sd ^ 2)

  ## add bias_val penalty ---
  llike_new <- llike_new - (par[1] - bias_val_mean) ^ 2 / (2 * bias_val_sd ^ 2)

  return(llike_new)
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


#' Using Parental Data for Offspring Genotyping.
#'
#' This function fits a hierarchical model to sequence counts from a
#' collection of siblings --- or a population of individuals
#' in Hardy-Weinberg equilibrium --- and returns genotyped information. The
#' hierarchy comes from either the fact that they share the same parents or they come
#' from a population in Hardy-Weinberg equilibrium. If
#' you also have parental sequencing data, then you can include this
#' to improve estimates.
#'
#' The key improvements in \code{updog} are its abilities to account for common features in
#' GBS data: sequencing error rate, read-mapping bias, overdispersion, and outlying points.
#'
#' @seealso \code{\link{plot.updog}} For plotting the output of \code{\link{updog_vanilla}}.
#'
#' @inheritParams updog_old
#' @inheritParams updog_update_all
#'
#' @return A list of class \code{updog} with some or all of the following elements:
#'     \describe{
#'         \item{\code{ogeno}}{A vector. Each element of which is the maximum a posteriori estimate of each individual's genotype.}
#'         \item{\code{maxpostprob}}{A vector. The maximum posterior probability of a genotype for each individual.}
#'         \item{\code{postmean}}{A vector. The posterior mean genotype for each individual.}
#'         \item{\code{bias_val}}{The estimated bias parameter. This is a value greater than 0 which is the ratio of the probability of correctly mapping a read containing the alternative allele to the probability of correctly mapping a read containing the reference allele. A value of 1 indicates no bias. A value less than one indicates bias towards the reference alelle. A value greater than 1 indiciates bias towards the alternative allele.}
#'         \item{\code{seq_error}}{The estimated sequencing error rate. This is between 0 and 1.}
#'         \item{\code{od_param}}{The estimated overdispersion parameter. Also known as the "intra-class correlation", this is the overdispersion parameter in the underlying beta of the beta-binomial distribution of the counts. Between 0 and 1, a value closer to 0 indicates less overdispersion and a value greater than 1 indicates greater overdispersion. In real data, we we typically see estimates between 0 and 0.01.}
#'         \item{\code{p1geno}}{The estimated genotype of one parent. The number of copies of the reference allele one of the parents has.}
#'         \item{\code{p1geno}}{The estimated genotype of the other parent. The number of copies of the reference allele the other parent has.}
#'         \item{\code{allele_freq}}{The estimated allele-frequency of the reference allele. This is the binomial proportion. Between 0 and 1, a value closer to 1 indicates a larger amount of reference alleles in the population.}
#'         \item{\code{out_prop}}{The estimated proportion of points that are outliers.}
#'         \item{\code{out_mean}}{The estimated mean of the outlier distribution. The outlier distribution is beta-binomial.}
#'         \item{\code{out_disp}}{The estimated overdispersion parameter of the outlier distribution. This is the "intra-class correlation" parameter of the beta-binomial outlier distribution.}
#'         \item{\code{prob_out}}{A vector. Each element of which is the posterior probability that a point is an outlier.}
#'         \item{\code{prob_ok}}{The posterior probability that a point is a non-outlier.}
#'         \item{\code{p1_prob_out}}{The posterior probability that parent 1 is an outlier.}
#'         \item{\code{p2_prob_out}}{The posterior probability that parent 2 is an outlier.}
#'         \item{\code{num_iter}}{The number of iterations the optimization program was run.}
#'         \item{\code{convergence}}{\code{1} if we reached \code{maxiter} and \code{0} otherwise.}
#'         \item{\code{llike}}{The final log-likelihood of the estimates.}
#'         \item{\code{hessian}}{The negative-Fisher information under the parameterization (s, ell, r), where s = log(bias_val) = log(d), ell = logit(seq_error) = logit(eps), and r = - logit(od_param) = - logit(tau). If you want standard errors for these parameters (in the described parameterization), simply take the negative inverse of the hessian.}
#'         \item{\code{input}}{A list with the input counts (\code{ocounts}), the input sizes (\code{osize}), input parental counts (\code{p1counts} and \code{p2counts}), input parental sizes (\code{p2size} and \code{p1size}), the ploidy (\code{ploidy}) and the model (\code{model}).}
#'         \item{\code{log_bias}}{The log of \code{bias_val}}
#'         \item{\code{logit_seq_error}}{The logit of \code{seq_error}. I.e. \code{log(seq_error / (1 - seq_error))}}
#'         \item{\code{neg_logit_od_param}}{The negative logit of \code{od_param}. I.e. \code{log((1 - od_param) / od_param)}}
#'         \item{\code{covmat}}{The observed Fisher information matrix of \code{c(log_bias, logit_seq_error, neg_logit_od_param)}. This can be used as the covariance matrix of these estimates.}
#'     }
#' @author David Gerard
#'
#' @export
#'
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
                          bound_bias = FALSE,
                          seq_error_mean = -4.7,
                          seq_error_sd = 1,
                          bias_val_mean = 0,
                          bias_val_sd = 0.7,
                          allele_freq = 0.5,
                          model = c("f1", "s1", "hw", "uniform"),
                          verbose = FALSE) {
  ## Deal with missing data -----------------------------------------------------------
  which_na <- is.na(ocounts) | is.na(osize)
  ocounts  <- ocounts[!which_na]
  osize    <- osize[!which_na]

  ## Hacky way so that I don't have to recode everything for od_param == 0 ------------
  if (od_param < 10 ^ -100) {
    od_param <- 10 ^ -100
  }

  # if (seq_error < 10 ^ -12) {
  #   stop("We don't support seq_error = 0 right now.")
  # }


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
  assertthat::assert_that(allele_freq > 0, allele_freq < 1)

  model <- match.arg(model)

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
                             update_od_param = update_od_param,
                             seq_error_mean = seq_error_mean,
                             seq_error_sd = seq_error_sd,
                             bias_val_mean = bias_val_mean,
                             bias_val_sd = bias_val_sd,
                             allele_freq = allele_freq,
                             model = model,
                             verbose = verbose)

  if (parout$od_param < 10 ^ -13) { ## fix for getting some weird postmat's
    parout$od_param <- 0
  }

  ## get genotype frequencies --------------------------------------------------------
  prob_geno <- get_prob_geno(ploidy = ploidy, model = model, p1geno = parout$p1geno, p2geno = parout$p2geno, allele_freq = parout$allele_freq)

  parout$postmat <- bbpost_tot(ocounts = ocounts, osize = osize,
                               ploidy = ploidy, prob_geno = prob_geno,
                               seq_error = parout$seq_error,
                               bias_val = parout$bias_val,
                               od_param = parout$od_param,
                               outlier = TRUE, out_prop = parout$out_prop,
                               out_mean = parout$out_mean,
                               out_disp = parout$out_disp)
  parout$ogeno <- apply(parout$postmat, 1, which.max) - 1
  parout$ogeno[abs(parout$prob_out - 1) < 10 ^ -3] <- NA   ## Put NA for ogeno when prob_out is almost 1 ----
  parout$maxpostprob   <- parout$postmat[cbind(1:nrow(parout$postmat), parout$ogeno + 1)]
  parout$postmean      <- c(parout$postmat %*% 0:ploidy)
  parout$input         <- list()
  parout$input$ocounts <- ocounts
  parout$input$osize   <- osize
  if (!is.null(p1counts) & !is.null(p1size)) {
    parout$input$p1counts <- p1counts
    parout$input$p1size <- p1size
  }
  if (!is.null(p2counts) & !is.null(p2size)) {
    parout$input$p2counts <- p2counts
    parout$input$p2size <- p2size
  }
  parout$input$ploidy <- ploidy
  parout$input$model  <- model

  ## Fix output according to model ----------------------------------------
  if (model == "hw") {
    parout$p1geno <- -1
    parout$p2geno <- -1
  } else if (model == "s1") {
    parout$p2geno <- -1
    parout$allele_freq <- -1
  } else if (model == "f1") {
    parout$allele_freq <- -1
  } else if (model == "uniform") {
    parout$allele_freq <- -1
    parout$p1geno <- -1
    parout$p2geno <- -1
  }

  ## deal with missingness again -----------------------------------------
  temp <- rep(NA, length = length(which_na))
  temp[!which_na] <- parout$prob_out
  parout$prob_out <- temp
  parout$prob_ok  <- 1 - parout$prob_out

  temp <- matrix(NA, nrow = length(which_na), ncol = ncol(parout$postmat))
  temp[!which_na, ] <- parout$postmat
  parout$postmat <- temp

  temp <- rep(NA, length = length(which_na))
  temp[!which_na] <- parout$ogeno
  parout$ogeno <- temp

  temp <- rep(NA, length = length(which_na))
  temp[!which_na] <- parout$maxpostprob
  parout$maxpostprob <- temp

  temp <- rep(NA, length = length(which_na))
  temp[!which_na] <- parout$postmean
  parout$postmean <- temp

  temp <- rep(NA, length = length(which_na))
  temp[!which_na] <- parout$input$ocounts
  parout$input$ocounts <- temp

  temp <- rep(NA, length = length(which_na))
  temp[!which_na] <- parout$input$osize
  parout$input$osize <- temp

  class(parout) <- "updog"
  return(parout)
}

