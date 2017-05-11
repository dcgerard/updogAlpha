
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
#'
#'
#' @author David Gerard
#'
#'
#'
updog_update_all <- function(ocounts, osize, ploidy, seq_error) {

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
      cat(oout$value, "\n")
      cat(oout$par, "\n\n")
      if (best_llike < oout$value) {
        best_par    <- oout$par
        best_llike <- oout$value
        best_p1    <- p1geno
        best_p2    <- p2geno
      }
    }
  }

  best_bias <- exp(best_par[1])
  best_seq  <- expit(best_par[2])
  best_od   <- 1 / (1 + exp(best_par[3]))
}


