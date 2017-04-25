

#' Using Parental Data for Offspring Genotyping.
#'
#' This function fits a hierarchical model to sequence counts from a
#' collection of siblings and returns genotyped information. The
#' hierarchy comes from the fact that they share the same parents. If
#' you also have parental sequencing data, then you can include this
#' to improve estimates.
#'
#' If you have a lot of parental sequencing data, then it could
#' suffice to run \code{updog} with \code{update_geno} set to \code{FALSE}, which
#' would save a lot of time.
#' Otherwise, you will probably want to borrow strength
#' between the offspring by setting \code{update_geno} to \code{TRUE}.
#'
#' @param ocounts A vector of non-negative integers. The ith element
#'     is the number of reads of the reference allele the ith child.
#' @param osize A vector of positive integers. The ith element is the
#'     total number of reads for the ith child.
#' @param p1counts A vector of non-negative integers. The ith element
#'     is the number of reads of the reference allele in the ith
#'     sample of parent 1.  If \code{NULL} then the prior
#'     probabilities on parent 1's genotype will default to uniform.
#' @param p1size A vector of positive integers. The ith element is the
#'     total number of reads in the ith sample of parent 1.  If
#'     \code{NULL} then the prior probabilities on parent 1's genotype
#'     will default to uniform.
#' @param p2counts A vector of non-negative integers. The ith element
#'     is the number of reads of the reference allele in the ith
#'     sample of parent 2.  If \code{NULL} then the prior
#'     probabilities on parent 2's genotype will default to uniform.
#' @param p2size A vector of positive integers. The ith element is the
#'     total number of reads in the ith sample of parent 2.  If
#'     \code{NULL} then the prior probabilities on parent 2's genotype
#'     will default to uniform.
#' @param ploidy A positive integer. The number of copies of the
#'     genome in the species.
#' @param seq_error A non-negative numeric. This is the known
#'     sequencing error rate. The default is to estimate this
#'     using data that are all approximately the reference allele.
#' @param integrate A logical. Should we integrate over our
#'     uncertainty in the parental genotypes (\code{TRUE}) or not
#'     (\code{FALSE}). The default is \code{FALSE} because we usually
#'     know the parental genotypes with near certainty so it's not
#'     important to integrate over our uncertainty in them. This is only
#'     implemented if \code{do_eb = FALSE}
#' @param update_geno A logical. Update the parental genotypes? If \code{FALSE} and if you have parental data, then
#'     we fix the parental genotypes to be the maximum a posteriori values. If you do not have parental data,
#'     then this should not be set to \code{FALSE}.
#' @param update_pi A logical. Update the mixing proporiton? If \code{FALSE}, then
#'     1\% of the observations are assumed to be outliers.
#' @param update_outlier A logical. Update the outlier distribution? I \code{FALSE}, then
#'     the outlier distribution is assumed to just be a uniform from 0 to 1.
#' @param update_rho A logical. Update the overdispersion parameter?
#' @param do_eb Should we do empirical Bayes (\code{TRUE}) or not (\code{FALSE})?
#'     You should have a lot of parental data to be able to set this to \code{FALSE}.
#' @param overdispersion A logical. Should we fit a beta-binomial model
#'     to account for overdispersion (\code{TRUE}) or not (\code{FALSE})? If
#'     \code{overdispersion = TRUE} then we start the overdispersion parameter,
#'     \code{rho} at 0.001, a very small value. If parental information is
#'     provided, then we use that data as the starting values for \code{rho}. If
#'     \code{update_rho} is \code{FALSE}, then these values are fixed throughout the
#'     estimation procedure.
#'
#'
#' @return A list with some or all of the following elements:
#'
#'   \code{opostprob}: A matrix of proportions whose (i, j)th
#'     element is the posterior probability that child j has i - 1
#'     copies of the reference allele.  That is, the rows index the
#'     genotype and the columns index the offspring.
#'
#'   \code{p1postprob}: A vector of proportions whose ith element is
#'   the posterior probability that parent 1 has i - 1 copies of the
#'   reference allele.  These are derived ONLY from parent 1's
#'   sequence data and not jointly with all of the data.
#'
#'   \code{p2postprob}: A vector of proportions whose ith element is
#'   the posterior probability that parent 2 has i - 1 copies of the
#'   reference allele.  These are derived ONLY from parent 2's
#'   sequence data and not jointly with all of the data.
#'
#'   \code{pival}: The estimated proportion of observations that are not outliers.
#'
#'   \code{rho}: The overdispersion parameter.
#'
#'   \code{out_mu}: The outlier distribution mean.
#'
#'   \code{out_rho}: The outlier distribution overdispersion parameter.
#'
#'   \code{p1geno}: The estimated genopype of parent 1.
#'
#'   \code{p2geno}: The estimated genotype of parent 2.
#'
#'   \code{prob_ok}: A vector of proportions. The ith element is the posterior
#'       probability that the ith element is not an outlier.
#'
#'   \code{ogeno}: A vector of integers. The ith element is the estimated genotype
#'       of the ith offspring.
#'
#'   \code{alpha}: The outlier distributions's shape 1 parameter.
#'
#'   \code{beta}: The outlier distributions's shape 2 parameter.
#'
#'   \code{seq_error}: The sequencing error rate used during the updog iterates.
#'
#'
#' @author David Gerard
#'
#' @export
#'
updog <- function(ocounts, osize,  ploidy, p1counts = NULL,
                  p1size = NULL, p2counts = NULL, p2size = NULL,
                  seq_error = NULL, integrate = FALSE, do_eb = TRUE, overdispersion = TRUE,
                  update_geno = TRUE, update_pi = TRUE, update_outlier = TRUE,
                  update_rho = TRUE) {

  ## check input -------------------------------------------------------------
  assertthat::assert_that(all(ocounts >= 0))
  assertthat::assert_that(all(osize >= 1))
  assertthat::assert_that(all(ocounts <= osize))
  assertthat::assert_that(is.logical(overdispersion))

  ## estimate seq_error ------------------------------------------------------
  if (is.null(seq_error)) {
    seq_error <- est_seq_error(ncounts = ocounts, ssize = osize, ploidy = ploidy)
  }

  if (seq_error == 0 & overdispersion) {
    stop("It is not allowed that both `overdisperion = TRUE` and `seq_error = 0`.\nThis is because the beta-binomial isn't defined with means 0 or 1.")
  }

  ## get priors on parental genotypes ----------------------------------------
  rho1 <- NULL
  rho2 <- NULL
  if (is.null(p1counts) | is.null(p1size)) {
    r1vec <- rep(1 / (ploidy + 1), times = ploidy + 1)
  } else {
    assertthat::assert_that(all(p1counts >= 0))
    assertthat::assert_that(all(p1size >= 1))
    assertthat::assert_that(all(p1counts <= p1size))
    r1vec <- bin_post(ncounts = p1counts, ssize = p1size, prior = ploidy,
                      seq_error = seq_error)
  }

  # else if (overdispersion) {
  #   assertthat::assert_that(all(p1counts >= 0))
  #   assertthat::assert_that(all(p1size >= 1))
  #   assertthat::assert_that(all(p1counts <= p1size))
  #   bbout <- bb_post(ncounts = p1counts, ssize = p1size, prior = ploidy, seq_error = seq_error)
  #   r1vec <- bbout$prob
  #   rho1 <- bbout$rho
  # }

  if (is.null(p2counts) | is.null(p2size)) {
    r2vec <- rep(1 / (ploidy + 1), times = ploidy + 1)
  } else {
    assertthat::assert_that(all(p2counts >= 0))
    assertthat::assert_that(all(p2size >= 1))
    assertthat::assert_that(all(p2counts <= p2size))
    r2vec <- bin_post(ncounts = p2counts, ssize = p2size, prior = ploidy,
                      seq_error = seq_error)
  }

  # else if (overdispersion){
  #   assertthat::assert_that(all(p2counts >= 0))
  #   assertthat::assert_that(all(p2size >= 1))
  #   assertthat::assert_that(all(p2counts <= p2size))
  #   bbout <- bb_post(ncounts = p2counts, ssize = p2size, prior = ploidy, seq_error = seq_error)
  #   r2vec <- bbout$prob
  #   rho2  <- bbout$rho
  # }

  ## initialize overdispersion parameter
  if (is.null(rho1) & !is.null(rho2)) {
    rho <- rho2
  } else if (!is.null(rho1) & is.null(rho2)) {
    rho <- rho1
  } else if (!is.null(rho1) & !is.null(rho2)) {
    rho <- (rho1 + rho2) / 2
  } else {
    rho <- 0.001
  }

  assertthat::assert_that(ploidy >= 1)

  ## derive offspring genotype probabilities given parental genotypes.
  qarray <- get_q_array(ploidy = ploidy)

  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)


  if (do_eb & !overdispersion) {
    umout <- updog_maximize(ocounts = ocounts, osize = osize, qarray = qarray, r1vec = r1vec, r2vec = r2vec,
                            pk = pk, pival = 0.99, alpha = 1, beta = 1, est_fudge = TRUE,
                            tol = 10 ^ -4, itermax = 1000,
                            update_geno = update_geno, update_pi = update_pi, update_beta = update_outlier)
  } else if (do_eb & overdispersion) {
    umout <- up_max_bb(ocounts, osize, qarray, r1vec, r2vec, pk, pival = 0.99, out_mu = 1/2,
                       out_rho = 1/3, rho = rho, est_fudge = TRUE, tol = 10 ^ -4, itermax = 1000,
                       update_geno = update_geno, update_pival = update_pi, update_outlier = update_outlier,
                       update_rho = update_rho)
  }

  if (do_eb) {
    umout$seq_error <- seq_error
    umout$input <- list()
    umout$input$ocounts <- ocounts
    umout$input$osize <- osize
    umout$input$p1counts <- p1counts
    umout$input$p1size <- p1size
    umout$input$p2counts <- p2counts
    umout$input$p2size <- p2size
    umout$input$ploidy <- ploidy
    class(umout) <- "updog"
    return(umout)
  }

  if (overdispersion) {
    warning("overdispersion not yet implemented when do_eb = FALSE")
  }

  ## naive estimator -------------------------------------------------------------
  if (integrate) {
    harray <- sweep(qarray, MARGIN = 1, STATS = r1vec, FUN = `*`)
    harray <- sweep(harray, MARGIN = 2, STATS = r2vec, FUN = `*`)
    hl <- apply(harray, 3, sum)
  } else {
    genop1 <- which.max(r1vec) - 1
    genop2 <- which.max(r2vec) - 1
    hl <- qarray[genop1 + 1, genop2 + 1, ]
  }

  ## get posterior probabilities of offspring genotypes.
  postprob <- mapply(FUN = bin_post, ocounts, osize,
                     MoreArgs = list(prior = c(hl), seq_error = seq_error))

  return_list <- list()
  return_list$opostprob  <- postprob
  return_list$p1postprob <- r1vec
  return_list$p2postprob <- r2vec
  return_list$seq_error  <- seq_error

  return(return_list)
}


#' Estimate parental genotypes using offspring genotypes.
#'
#' @inheritParams updog
#' @param qarray An three-way array of proportions. The (i, j, k)th
#'     element is the probability of an offspring having k - 1
#'     reference alleles given that parent 1 has i - 1 refrerence
#'     alleles and parent 2 has j - 1 reference alleles. Each
#'     dimension of the array is \code{ploidy + 1}.
#' @param r1vec A vector of prior probabilities whose ith element is
#'     the prior probability that parent 1 has i -1 copies of allele
#'     A.
#' @param r2vec A vector of prior probabilities whose ith element is
#'     the prior probability that parent 2 has i -1 copies of allele
#'     A.
#' @param pival A fudge factor. The probability of not begin a
#'     mistake. Estimated if \code{est_fudge = TRUE}. Set to 0.99. Only used if \code{est_fudge = TRUE}.
#' @param est_fudge A logical. Should we estimate the fudge factor
#'     \code{TRUE} or not \code{FALSE}?
#' @param pk The parameter space. Should all be between 0 and 1.
#' @param tol A positive numeric. The tolerance for the stopping criterion.
#' @param itermax A positive integer. The maximum number of iterations to run in the EM.
#' @param alpha The alpha parameter in the beta distribution of the outliers. Only used if \code{est_fudge = TRUE}.
#' @param beta The beta parameter in the beta distribution of the outliers. Only used if \code{est_fudge = TRUE}.
#' @param update_geno A logical. Update the parental genotypes?
#' @param update_pi A logical. Update the mixing proporiton?
#' @param update_beta A logical. Update the beta distribution?
#'
#' @author David Gerard
#'
updog_maximize <- function(ocounts, osize, qarray, r1vec, r2vec, pk, pival = 0.99,
                           alpha = 0.1, beta = 0.1,
                           est_fudge = TRUE,
                           tol = 10 ^ -4, itermax = 1000,
                           update_geno = FALSE, update_pi = FALSE, update_beta = FALSE) {

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

    dbinommat <- mapply(FUN = stats::dbinom, ocounts, osize, MoreArgs = list(prob = pk))
    ldbinommat <- mapply(FUN = stats::dbinom, ocounts, osize, MoreArgs = list(prob = pk, log = TRUE))

    if (update_geno) {
      marg_lik_mat <- matrix(NA, nrow = ploidy + 1, ncol = ploidy + 1)
      for (ell1 in 0:ploidy) {
        for (ell2 in ell1:ploidy) {
          ## avec <- qarray[ell1 + 1, ell2 + 1, ]
          ## marg_lik_mat[ell1 + 1, ell2 + 1] <- sum(log(colSums(avec * dbinommat))) +
          ##    log(r1vec[ell1 + 1]) + log(r2vec[ell2 + 1])

          ## this is more numerically stable.
          lavec <- log(qarray[ell1 + 1, ell2 + 1, ])
          summands <- lavec + ldbinommat
          which_not_inf <- !is.infinite(lavec)
          colmax <- apply(summands[which_not_inf, , drop = FALSE], 2, max)
          llval <- sum(log(colSums(exp(sweep(x = summands[which_not_inf, , drop = FALSE],
                                             MARGIN = 2, STATS = colmax, FUN = `-`)))) +
                         colmax) + log(r1vec[ell1 + 1]) + log(r2vec[ell2 + 1])
          marg_lik_mat[ell1 + 1, ell2 + 1] <- llval
        }
      }

      which_best <- which(marg_lik_mat == max(marg_lik_mat, na.rm = TRUE), arr.ind = TRUE)

      p1geno <- which_best[1] - 1
      p2geno <- which_best[2] - 1
    } else {
      p1geno <- which.max(r1vec) - 1
      p2geno <- which.max(r2vec) - 1
    }

    if (est_fudge) {
      objval <- up_obj(pival = pival, p1geno = p1geno, p2geno = p2geno,
                       alpha = alpha, beta = beta, ocounts = ocounts, osize = osize,
                        dbinommat = dbinommat, qarray = qarray, r1vec = r1vec, r2vec = r2vec)

      err <- tol + 1
      iterindex <- 1
      while (err > tol & iterindex < itermax) {
        objold <- objval
        fixout <- up_fix(pival = pival, p1geno = p1geno, p2geno = p2geno, alpha = alpha,
                         beta = beta, ocounts = ocounts, osize = osize,
                         ldbinommat = ldbinommat, qarray = qarray, r1vec = r1vec, r2vec = r2vec,
                         update_pi = update_pi, update_geno = update_geno, update_beta = update_beta)

        p1geno <- fixout$p1geno
        p2geno <- fixout$p2geno
        alpha  <- fixout$alpha
        beta   <- fixout$beta
        pival  <- fixout$pival
        theta  <- fixout$theta

        objval <- up_obj(pival = pival, p1geno = p1geno, p2geno = p2geno,
                         alpha = alpha, beta = beta, ocounts = ocounts, osize = osize,
                         dbinommat = dbinommat, qarray = qarray, r1vec = r1vec, r2vec = r2vec)

        err <- abs(objval / objold - 1)

        assertthat::assert_that(objval - objold > -10 ^ -14)
      }
    }

    avec_final <- qarray[p1geno + 1, p2geno + 1, ]

    bbvec <- dbetabinom(x = ocounts, size = osize, alpha = alpha, beta = beta, log = FALSE)
    probmat <- sweep(x = pival * dbinommat, MARGIN = 2, STATS = (1 - pival) * bbvec, FUN = `+`) * avec_final

    postprob <- sweep(x = probmat, MARGIN = 2, STATS = colSums(probmat), FUN = `/`)

    ogeno <- apply(postprob, 2, which.max) - 1

    return_list           <- list()
    return_list$p1geno    <- p1geno
    return_list$p2geno    <- p2geno
    return_list$ogeno     <- ogeno
    return_list$prob_ok   <- theta
    return_list$pival     <- pival
    return_list$alpha     <- alpha
    return_list$beta      <- beta
    return_list$opostprob <- postprob

    return(return_list)
}

#' Fixed-point iteration for updog EM.
#'
#' @inheritParams updog
#' @inheritParams updog_maximize
#' @inheritParams up_obj
#' @param update_pi A logical. Should we update the fudge factor
#'     (\code{TRUE}) or not (\code{FALSE})?
#' @param update_geno A logical. Should we update the parental
#'     genotypes (\code{TRUE}) or not (\code{FALSE})?
#' @param update_beta A logical. Should we update the outlier model
#'     \code{TRUE} or not (\code{FALSE})?
#' @param ldbinommat The log binomial densities of the counts. The rows index the
#'     parental genotypes and the columns index the individuals.
#'
#' @author Davod Gerard
#'
#'
up_fix <- function(pival, p1geno, p2geno, alpha, beta, ocounts, osize,
                   ldbinommat, qarray, r1vec, r2vec, update_pi = TRUE,
                   update_geno = TRUE, update_beta = TRUE) {

    ploidy <- length(r1vec) - 1
    assertthat::assert_that(p1geno %in% 0:ploidy)
    assertthat::assert_that(p2geno %in% 0:ploidy)
    assertthat::are_equal(ploidy, length(r2vec) - 1)
    assertthat::assert_that(all(dim(qarray) == ploidy + 1))
    assertthat::are_equal(ploidy, nrow(ldbinommat) - 1)
    assertthat::assert_that(pival >= 0)
    assertthat::assert_that(pival <= 1)
    assertthat::assert_that(all(ocounts >= 0))
    assertthat::assert_that(all(osize >= ocounts))
    assertthat::are_equal(sum(r1vec), 1)
    assertthat::are_equal(sum(r2vec), 1)
    assertthat::assert_that(alpha > 0)
    assertthat::assert_that(beta > 0)

    avec <- qarray[p1geno + 1, p2geno + 1, ]

    bbvec <- dbetabinom(x = ocounts, size = osize, alpha = alpha, beta = beta, log = FALSE)
    num <- pival * colSums(avec * exp(ldbinommat))
    denom <- num + (1 - pival) * bbvec

    theta <- num / denom

    ## update parental genotypes ----------------------------------------------------
    if (update_geno) {
        cat(p1geno, p2geno, "\n")
        marg_lik_mat <- matrix(NA, nrow = ploidy + 1, ncol = ploidy + 1)
        for (ell1 in 0:ploidy) {
            for (ell2 in ell1:ploidy) {
                # avec <- qarray[ell1 + 1, ell2 + 1, ]
                # marg_lik_mat[ell1 + 1, ell2 + 1] <- sum(theta * log(colSums(avec * dbinommat))) +
                #    log(r1vec[ell1 + 1]) + log(r2vec[ell2 + 1])

                ## this is more numerically stable
                lavec <- log(qarray[ell1 + 1, ell2 + 1, ])
                summands <- lavec + ldbinommat
                which_not_inf <- !is.infinite(lavec)
                colmax <- apply(summands[which_not_inf, , drop = FALSE], 2, max)
                llval <- sum(theta * log(colSums(exp(sweep(x = summands[which_not_inf, , drop = FALSE],
                                                   MARGIN = 2, STATS = colmax, FUN = `-`)))) +
                               colmax) + log(r1vec[ell1 + 1]) + log(r2vec[ell2 + 1])
                marg_lik_mat[ell1 + 1, ell2 + 1] <- llval
            }
        }
        which_best <- which(marg_lik_mat == max(marg_lik_mat, na.rm = TRUE), arr.ind = TRUE)
        p1geno <- which_best[1] - 1
        p2geno <- which_best[2] - 1
    }

    ## update pival ------------------------------------------------------------------
    if (update_pi) {
        pival <- mean(theta)
    }


    ## Update alpha and beta ---------------------------------------------------------
    if (update_beta) {
      abeps <- 0.0001
        oout <- stats::optim(par = c(alpha, beta), fn = dbbwrapper, lower = c(abeps, abeps), upper = c(1, 1),
                             method = "L-BFGS-B",
                             control = list(fnscale = -1),
                             x = ocounts, size = osize, theta = theta, log = TRUE)
        alpha <- oout$par[1]
        beta  <- oout$par[2]
    }

    return(list(p1geno = p1geno, p2geno = p2geno, theta = theta,
                pival = pival, beta = beta, alpha = alpha))
}

dbbwrapper <- function(par, x, size, theta, log = TRUE) {
  sum((1 - theta) * dbetabinom(x = x, size = size, alpha = par[1], beta = par[2], log = TRUE))
}

#' The objective function
#'
#' @inheritParams updog_maximize
#' @inheritParams updog
#' @param p1geno The genotype of parent 1
#' @param p2geno The genotype of parent 2
#' @param pival A proportion. The probability of not being a mistake.
#' @param dbinommat A matrix. The rows index the genotype and the
#'     columns index the offspring.  These are the binomial
#'     probabilities of the count given the genotype.
#' @param alpha The alpha parameter in \code{\link{dbetabinom}}.
#' @param beta The beta parameter in \code{\link{dbetabinom}}.
#'
#' @author David Gerard
#'
up_obj <- function(pival, p1geno, p2geno, alpha, beta, ocounts, osize,
                   dbinommat, qarray, r1vec, r2vec) {



    ploidy <- length(r1vec) - 1
    assertthat::assert_that(p1geno %in% 0:ploidy)
    assertthat::assert_that(p2geno %in% 0:ploidy)
    assertthat::are_equal(ploidy, length(r2vec) - 1)
    assertthat::assert_that(all(dim(qarray) == ploidy + 1))
    assertthat::are_equal(ploidy, nrow(dbinommat) - 1)
    assertthat::assert_that(pival >= 0)
    assertthat::assert_that(pival <= 1)
    assertthat::assert_that(all(ocounts >= 0))
    assertthat::assert_that(all(osize >= ocounts))
    assertthat::are_equal(sum(r1vec), 1)
    assertthat::are_equal(sum(r2vec), 1)
    assertthat::assert_that(alpha > 0)
    assertthat::assert_that(beta > 0)

    avec <- qarray[p1geno + 1, p2geno + 1, ]


    bbout <- dbetabinom(x = ocounts, size = osize, alpha = alpha, beta = beta)
    binout <- colSums(avec * dbinommat)

    obj <- sum(log(pival * binout + (1 - pival) * bbout)) + log(r1vec[p1geno + 1]) + log(r2vec[p2geno + 1])

    return(obj)
}



#' Derive the probabilities of an offspring's genotype given its
#' parental genotypes for all possible combinations of parental and
#' offspring genotypes.
#'
#' @param ploidy A positive integer. The ploidy of the species.
#'
#' @author David Gerard
#'
#' @return An three-way array of proportions. The (i, j, k)th element
#'     is the probability of an offspring having k - 1 reference
#'     alleles given that parent 1 has i - 1 refrerence alleles and
#'     parent 2 has j - 1 reference alleles. Each dimension of the
#'     array is \code{ploidy + 1}. In the dimension names, "A" stands
#'     for the reference allele and "a" stands for any other allele.
#'
#' @export
#'
get_q_array <- function(ploidy) {
  assertthat::assert_that(ploidy > 0)
  assertthat::are_equal(ploidy %% 2, 0)

  qarray <- array(0, dim = rep(ploidy + 1, 3))

  for(oindex in 0:ploidy) {
    for (p1index in 0:ploidy) {
      for (p2index in 0:ploidy) {
        if (p1index + p2index < oindex) {
          qarray[p1index + 1, p2index + 1, oindex + 1] <- 0
        } else {
          minval <- max(0, oindex - p2index)
          maxval <- min(ploidy / 2, p1index)
          aseq <- minval:maxval

          p1prob <- stats::dhyper(x = aseq, m = p1index, n = ploidy - p1index, k = ploidy / 2)
          p2prob <- stats::dhyper(x = oindex - aseq, m = p2index, n = ploidy - p2index, k = ploidy / 2)
          qarray[p1index + 1, p2index + 1, oindex + 1] <- sum(p1prob * p2prob)
        }
      }
    }
  }

  ## get dimnames

  dimvec <- get_dimname(ploidy)
  dimnames(qarray) <- list(parent1 = dimvec, parent2 = dimvec, offspring = dimvec)


  assertthat::assert_that(all(abs(apply(qarray, c(1, 2), sum) - 1) < 10 ^ -14))

  return(qarray)
}


#' Get the transition probabilities from going from one ploidy to
#' another assuming a constant error rate.
#'
#' @param ploidy The ploidy of the species.
#' @param eps The probability of one allele switching.
#'
#' @author David Gerard
#'
#' @export
#'
get_transition_mat <- function(ploidy, eps) {
  bmat <- matrix(NA, nrow = ploidy + 1, ncol = ploidy + 1)
  ## Transition from i to j
  for (i in 0:ploidy) {
    for (j in 0:ploidy) {
      bmat[i + 1, j + 1] <- sum(stats::dbinom(x = 0:j, size = i, prob = 1 - eps) *
                                  stats::dbinom(x = j:0, size = ploidy - i, prob = eps))
    }
  }
  assertthat::assert_that(all(abs(rowSums(bmat) - 1) < 10 ^ -14))
  return(bmat)
}

#' Make a genotype plot.
#'
#' The x-axis will be the counts of the non-reference allele,
#' and the y-axis will be the counts of the reference allele.
#' Transparency is controlled by the \code{prob_ok} vector,
#' size is controlled by the \code{maxpostprob} vector.
#'
#' If parental genotypes are provided (\code{p1geno} and \code{p2geno}) then
#' the will be colored the same as the offspring. Since they are often hard to see,
#' a small black dot will also indicate their position.
#'
#' @param ocounts A vector of non-negative integers. The number of
#'     reference alleles observed in the offspring.
#' @param osize A vector of positive integers. The total number of
#'     reads in the offspring.
#' @param p1counts A vector of non-negative integers. The number of
#'     reference alleles observed in parent 1.
#' @param p1size A vector of positive integers. The total number of
#'     reads in parent 1.
#' @param p2counts A vector of non-negative integers. The number of
#'     reference alleles observed in parent 2.
#' @param p2size A vector of positive integers. The total number of
#'     reads in parent 2.
#' @param ploidy A non-negative integer. The ploidy of the species.
#' @param ogeno The child genotypes
#' @param seq_error The average sequencing error rate.
#' @param prob_ok A vector of posterior probabilities of not being a mistake.
#' @param maxpostprob A vector of the posterior probabilities of begin at the modal probability.
#' @param p1geno Parent 1's genotype.
#' @param p2geno Parent 2's genotype.
#'
#' @export
#'
#' @author David Gerard
#'
plot_geno <- function(ocounts, osize, ploidy, p1counts = NULL, p1size = NULL, p2counts = NULL,
                      p2size = NULL, ogeno = NULL, seq_error = 0.01, prob_ok = NULL, maxpostprob = NULL,
                      p1geno = NULL, p2geno = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 must be installed to use this function")
  }

  assertthat::assert_that(all(ocounts >= 0, na.rm = TRUE))
  assertthat::assert_that(all(osize >= ocounts, na.rm = TRUE))
  assertthat::assert_that(ploidy >= 1)
  assertthat::assert_that(seq_error >= 0)

  ## get probabilities
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)

  dfdat <- data.frame(A = ocounts, a = osize - ocounts)
  maxcount <- max(max(dfdat$A, na.rm = TRUE), max(dfdat$a, na.rm = TRUE))
  if (!is.null(ogeno)) {
    assertthat::are_equal(length(ogeno), length(ocounts))
    dfdat$genotype <- factor(ogeno, levels = 0:ploidy)
  }

  if (!is.null(prob_ok)) {
    assertthat::assert_that(all(prob_ok >= 0))
    assertthat::assert_that(all(prob_ok <= 1))
    dfdat$prob_ok <- prob_ok
  }

  if (!is.null(maxpostprob)) {
    assertthat::assert_that(all(maxpostprob >= 0))
    assertthat::assert_that(all(maxpostprob <= 1))
    dfdat$maxpostprob <- maxpostprob
  }

  slopevec <- pk / (1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1),
                         xend = xend, yend = yend)

  ## Plot children
  if (is.null(prob_ok) & is.null(maxpostprob)) {
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(y = "A", x = "a"))
  } else if (!is.null(prob_ok) & is.null(maxpostprob)) {
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(y = "A", x = "a", alpha = "prob_ok"))
  } else if (is.null(prob_ok) & !is.null(maxpostprob)) {
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(y = "A", x = "a", size = "maxpostprob"))
  } else if (!is.null(prob_ok) & !is.null(maxpostprob)) {
    pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(y = "A", x = "a", alpha = "prob_ok", size = "maxpostprob"))
  }

  ## add offspring genotypes ------------------------------------------------
  if (!is.null(ogeno)) {
    pl <- pl + ggplot2::geom_point(ggplot2::aes_string(colour = "genotype"))
  } else {
    pl <- pl + ggplot2::geom_point()
  }

  pl <- pl + ggplot2::scale_size(range = c(0.7, 3))

  pl <- pl + ggplot2::theme_bw() +
    ggplot2::xlim(0, maxcount) +
    ggplot2::ylim(0, maxcount) +
    ggplot2::ylab("Counts A") +
    ggplot2::xlab("Counts a")  +
    ggplot2::geom_segment(data = df_lines, mapping = ggplot2::aes_string(x = "x", y = "y", xend = "xend", yend = "yend"),
                          lty = 2, alpha = 1/2, color = "black", size = 0.5)

   ## add parents if we have them
   if (!is.null(p1size) & !is.null(p1counts)) {
     assertthat::assert_that(all(p1counts >= 0, na.rm = TRUE))
     assertthat::assert_that(all(p1size >= p1counts, na.rm = TRUE))
     p1dat <- data.frame(A = p1counts, a = p1size - p1counts)
     if (!is.null(p1geno)) {
       p1dat$genotype <- factor(p1geno, levels = 0:ploidy)
       pl <- pl + ggplot2::geom_point(data = p1dat, mapping = ggplot2::aes_string(color = "genotype"),
                                      size = 3, pch = 3, alpha = 1, show.legend = FALSE)
       pl <- pl + ggplot2::geom_point(data = p1dat, size = 1, color = "black", pch = 16, alpha = 1)
     } else {
       pl <- pl + ggplot2::geom_point(data = p1dat, size = 3, color = "black", pch = 3, alpha = 1)
     }
   }
   if (!is.null(p2size) & !is.null(p2counts)) {
     assertthat::assert_that(all(p2counts >= 0, na.rm = TRUE))
     assertthat::assert_that(all(p2size >= p2counts, na.rm = TRUE))
     p2dat <- data.frame(A = p2counts, a = p2size - p2counts)
     if (!is.null(p2geno)) {
       p2dat$genotype <- factor(p2geno, levels = 0:ploidy)
       pl <- pl + ggplot2::geom_point(data = p2dat, mapping = ggplot2::aes_string(color = "genotype"),
                                      size = 3, pch = 4, alpha = 1, show.legend = FALSE)
       pl <- pl + ggplot2::geom_point(data = p2dat, size = 1, color = "black", pch = 16, alpha = 1)
     } else {
       pl <- pl + ggplot2::geom_point(data = p2dat, size = 3, color = "black", pch = 4, alpha = 1)
     }
   }

  if (!is.null(ogeno) | !is.null(p1geno) | !is.null(p2geno)) {
    pl <- pl + ggplot2::scale_color_hue(drop = FALSE)
  }

   return(pl)
}

#' Returns a vector character strings that are all of the possible
#' combinations of the reference allele and the non-reference allele.
#'
#' @param ploidy The ploidy of the species.
#'
#' @return For example, if \code{ploidy = 3} then this will return
#'     c("aaa", "Aaa", "AAa", "AAA")
#'
#' @author David Gerard
#'
#'
get_dimname <- function(ploidy) {
  dimvec <- sapply(mapply(FUN = c, lapply(X = 0:ploidy, FUN = rep.int, x = "A"),
                          lapply(X = ploidy:0, FUN = rep.int, x = "a"),
                          SIMPLIFY = FALSE),
                   FUN = paste, collapse = "")
  return(dimvec)
}

#' Calculates posterior probabilities of a genotype given just the
#' sequence counts.
#'
#' The function assumes that all data come from independent samples on
#' the same individual.
#'
#' @param ncounts A vector of non-negative integers. The ith element
#'     is the number of counts of the ith sample.
#' @param ssize A vector of positive integers. The ith element is the
#'     total number of counts of the ith sample.
#' @param prior A vector of non-negative numerics that sum to one. The
#'     prior probability on the genotype. The first element is the
#'     prior probability of zero reference alleles, the second element
#'     is the prior probability of one reference allele, etc. The
#'     length of \code{prior} is one more than the ploidy of the
#'     species. You can alternatively specify \code{prior} as the
#'     ploidy of the individual, for which it will set a uniform prior
#'     on the genotype. For example, setting \code{prior = 3} will
#'     result in using \code{(1/4, 1/4, 1/4, 1/4)} as the prior
#'     probability for the genotypes (Aaaa, AAaa, AAAa, AAAA) where
#'     "A" is the reference allele in a 4-ploid individual.
#' @param seq_error A non-negative numeric. This is the known
#'     sequencing error rate. This is a rough high-ball error rate
#'     given by Li et. al. (2011).
#'
#' @return A vector of probabilities. The ith element is the posterior
#'     probability that the individual has i - 1 copies of the
#'     reference allele.
#'
#' @author David Gerard
#'
#' @export
#'
#' @references Li, Yun, Carlo Sidore, Hyun Min Kang, Michael Boehnke,
#'     and Gonçalo R. Abecasis.
#'     \href{https://www.ncbi.nlm.nih.gov/pubmed/21460063}{"Low-coverage sequencing: implications for design of complex trait association studies."}
#'     Genome research (2011).
#'
#'
bin_post <- function(ncounts, ssize, prior, seq_error = 0.01) {

  ## check input -------------------------------------------------------------
  assertthat::assert_that(all(ncounts >= 0))
  assertthat::assert_that(all(ssize >= 1))
  assertthat::assert_that(all(prior >= 0))
  assertthat::assert_that(all(ssize >= ncounts))
  assertthat::assert_that(seq_error >= 0)
  assertthat::assert_that(seq_error <= 1)


  if (abs(sum(prior) - 1) < 10 ^ -14) {
    ploidy <- length(prior) - 1
  } else {
    ploidy <- prior
    prior  <- rep(1 / (ploidy + 1), length = ploidy + 1)
  }

  assertthat::are_equal(sum(prior), 1)

  ## calculate log-probabilities ---------------------------------------------
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities

  ## deal with error rate
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)



  nsum <- sum(ncounts)
  tsum <- sum(ssize)

  logprob <- nsum * log(pk) + (tsum - nsum) * log(1 - pk) + log(prior)
  ## temp <- stats::dbinom(x = nsum, size = tsum, prob = pk, log = TRUE) + log(prior)
  ## assertthat::are_equal(exp(logprob) / sum(exp(logprob)), exp(temp) / sum(exp(temp)))


  ## Exponentiate and deal with special cases --------------------------------
  if (seq_error < 10 ^ -14) {
    logprob_b <- logprob
    if (nsum == 0) {
      logprob_b[1] <- log(prior[1])
      logprob_b[1:ploidy] <- logprob_b[1:ploidy] - max(logprob_b[1:ploidy])
      postprob <- exp(logprob_b) / sum(exp(logprob_b[1:ploidy]))
    } else if (nsum == tsum) {
      logprob_b[ploidy + 1] <- log(prior[ploidy + 1])
      logprob_b[2:(ploidy + 1)] <- logprob_b[2:(ploidy + 1)] - max(logprob_b[2:(ploidy + 1)])
      postprob <- exp(logprob_b) / sum(exp(logprob_b[2:(ploidy + 1)]))
    } else {
      logprob_b[2:ploidy] <- logprob_b[2:ploidy] - max(logprob_b[2:ploidy])
      postprob <- exp(logprob_b) / sum(exp(logprob_b[2:ploidy]))
    }
  } else {
    logprob <- logprob - max(logprob)
    postprob <- exp(logprob) / sum(exp(logprob))
  }

  ## test when third case only
  ## assertthat::are_equal(postprob, exp(logprob) / sum(exp(logprob)))

  names(postprob) <- get_dimname(ploidy)

  return(postprob)
}




# segreg_poly = function to compute Mendelian segregation in polyploid
# <m> ploidy
# <dP> parent 1 genotype - Ex: {AAAAAA=6},{AAAAAa=5}, ... , {aaaaaa=0}
# <dQ> parent 2 genotype
segreg_poly <- function(m, dP, dQ) {
 if (m%%2 != 0)
   stop("m must be an even number")
 p.dose <- numeric((m + 1))
 p.names <- character((m + 1))
 seg.p1 <- stats::dhyper(x = c(0:(m + 1)), m = dP, n = (m - dP), k = m/2)
 seg.p2 <- stats::dhyper(x = c(0:(m + 1)), m = dQ, n = (m - dQ), k = m/2)
 M <- tcrossprod(seg.p1, seg.p2)
 for (i in 1:nrow(M)) {
   for (j in 1:ncol(M)) {
     p.dose[i + j - 1] <- p.dose[i + j - 1] + M[i, j]
   }
 }
 p.dose <- p.dose[!is.na(p.dose)]
 for (i in 0:m) p.names[i + 1] <- paste(paste(rep("A", i), collapse = ""), paste(rep("a", (m - i)), collapse = ""), sep = "")
 names(p.dose) <- p.names
 return(p.dose)
}
