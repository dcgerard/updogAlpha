

#' Using Parental Data for Offspring Genotyping.
#'
#' This function fits a hierarchical model to sequence counts from a collection
#' of siblings and return genotyped information. The hierarchy comes from the fact
#' that they share the same parents. If you also have parental sequencing data,
#' then you can include this to improve estimates.
#'
#' If you have a lot of parental sequencing data, then it could suffice to run
#' \code{updog} with \code{do_mcmc} set to \code{FALSE}. Otherwise, you will
#' probably want to borrow strength between the offspring by setting \code{do_mcmc}
#' to \code{TRUE}.
#'
#' @param ocounts A vector of non-negative integers. The ith element is the
#'   number of reads of the reference allele the ith child.
#' @param osize A vector of positive integers. The ith element is the total
#'   number of reads for the ith child.
#' @param p1counts A vector of non-negative integers. The ith element is the
#'   number of reads of the reference allele in the ith sample of parent 1.
#'   If \code{NULL} then the prior probabilities on parent 1's genotype will default
#'   to uniform.
#' @param p1size A vector of positive integers. The ith element is the total
#'   number of reads in the ith sample of parent 1.
#'   If \code{NULL} then the prior probabilities on parent 1's genotype will default
#'   to uniform.
#' @param p2counts A vector of non-negative integers. The ith element is the
#'   number of reads of the reference allele in the ith sample of parent 2.
#'   If \code{NULL} then the prior probabilities on parent 2's genotype will default
#'   to uniform.
#' @param p2size A vector of positive integers. The ith element is the total
#'   number of reads in the ith sample of parent 2.
#'   If \code{NULL} then the prior probabilities on parent 2's genotype will default
#'   to uniform.
#' @param ploidy A positive integer. The number of copies of the genome in the
#'   species.
#' @param seq_error A non-negative numeric. This is the known sequencing error
#'   rate. This is a rough high-ball error rate given by Li et. al. (2011).
#' @param do_mcmc A logical. Should we also run a Gibbs sampler to jointly estimate
#'   the parental and child genotypes (\code{TRUE}) or not (\code{FALSE})? If \code{TRUE},
#'   the total number of iterations run is \code{burnin + iteramx}.
#' @param burnin A non-negative integer. The number of iterations to ignore
#'   in the Gibbs sampler if \code{do_mcmc = TRUE}.
#' @param itermax A positive integer. The number of iterations to collect
#'   in the Gibbs sampler if \code{do_mcmc = TRUE}.
#' @param iterate A logical. Should we perform the iterative posterior approximation
#'   scheme (\code{TRUE}) or note (\code{FALSE})?
#'
#'
#' @return A list with some or all of the following elements:
#'   \code{opostprob}: A matrix of proportions whose (i, j)th element
#'   is the posterior probability that child j has i - 1 copies of the reference allele.
#'   That is, the rows index the genotype and the columns index the offspring.
#'   These are derived by the one-offspring-at-a-time procedure.
#'
#'   \code{p1postprob}: A vector of proportions whose ith element is the posterior
#'   probability that parent 1 has i - 1 copies of the reference allele.
#'   These are derived ONLY from parent 1's sequence data and not jointly with
#'   all of the data.
#'
#'   \code{p2postprob}: A vector of proportions whose ith element is the posterior
#'   probability that parent 2 has i - 1 copies of the reference allele.
#'   These are derived ONLY from parent 2's sequence data and not jointly with
#'   all of the data.
#'
#'   \code{m_opostprob}: A matrix of proportions whose (i, j)th element
#'   is the posterior probability that child j has i - 1 copies of the reference allele.
#'   These are derived from the joint analysis with MCMC.
#'
#'   \code{m_p1postprob}: A vector of proportions whose ith element is the posterior
#'   probability that parent 1 has i - 1 copies of the reference allele.
#'   These are derived from the joint analysis with MCMC.
#'
#'   \code{m_p2postprob}: A vector of proportions whose ith element is the posterior
#'   probability that parent 2 has i - 1 copies of the reference allele.
#'   These are derived from the joint analysis with MCMC.
#'
#' @author David Gerard
#'
#' @export
#'
#' @references Li, Yun, Carlo Sidore, Hyun Min Kang, Michael Boehnke, and
#'   Gonçalo R. Abecasis.
#'   \href{https://www.ncbi.nlm.nih.gov/pubmed/21460063}{"Low-coverage sequencing: implications for design of complex trait association studies."}
#'   Genome research (2011).
#'
updog <- function(ocounts, osize,  ploidy, p1counts = NULL,
                  p1size = NULL, p2counts = NULL, p2size = NULL,
                  seq_error = 0.1, do_mcmc = FALSE, iterate = FALSE,
                  burnin = 250, itermax = 1000) {

  ## check input -------------------------------------------------------------
  assertthat::assert_that(all(ocounts >= 0))
  assertthat::assert_that(all(osize >= 1))
  assertthat::assert_that(all(ocounts <= osize))

  ## get priors on parental genotypes
  if (is.null(p1counts) | is.null(p1size)) {
    r1vec <- rep(1 / (ploidy + 1), times = ploidy + 1)
  } else {
    assertthat::assert_that(all(p1counts >= 0))
    assertthat::assert_that(all(p1size >= 1))
    assertthat::assert_that(all(p1counts <= p1size))
    r1vec <- bin_post(ncounts = p1counts, ssize = p1size, prior = ploidy,
                      seq_error = seq_error)
  }

  if (is.null(p2counts) | is.null(p2size)) {
    r2vec <- rep(1 / (ploidy + 1), times = ploidy + 1)
  } else {
    assertthat::assert_that(all(p2counts >= 0))
    assertthat::assert_that(all(p2size >= 1))
    assertthat::assert_that(all(p2counts <= p2size))
    r2vec <- bin_post(ncounts = p2counts, ssize = p2size, prior = ploidy,
                      seq_error = seq_error)
  }

  assertthat::assert_that(ploidy >= 1)
  assertthat::assert_that(is.logical(do_mcmc))
  assertthat::assert_that(burnin >= 0)
  assertthat::assert_that(burnin < itermax)

  ## derive offspring genotype probabilities given parental genotypes.
  qarray <- get_q_array(ploidy = ploidy)

  ## iterate to get r1vec and r2vec
  if (iterate) {
    itout <- updog_iterate(ocounts = ocounts, osize = osize,
                           qarray = qarray, r1vec = r1vec,
                           r2vec = r2vec, seq_error = seq_error)
    phi_vec <- itout$r1vec
    psi_vec <- itout$r2vec
  } else {
    ## Derive prior probabilities on offspring genotypes
    phi_vec <- r1vec ## posterior prob of p1 genotype
    psi_vec <- r2vec ## posterior prob of p2 genotype
  }

  harray <- sweep(qarray, MARGIN = 1, STATS = phi_vec, FUN = `*`)
  harray <- sweep(harray, MARGIN = 2, STATS = psi_vec, FUN = `*`)
  hl <- apply(harray, 3, sum)

  ## get posterior probabilities of offspring genotypes.
  postprob <- mapply(FUN = bin_post, ocounts, osize,
                     MoreArgs = list(prior = c(hl), seq_error = seq_error))

  return_list <- list()
  return_list$opostprob  <- postprob
  return_list$p1postprob <- phi_vec
  return_list$p2postprob <- psi_vec

  if (do_mcmc) {
    mout <- updog_mcmc(ocounts = ocounts, osize = osize, qarray = qarray, r1vec = r1vec,
                       r2vec = r2vec, seq_error = seq_error, itermax = itermax, burnin = burnin)
    return_list$m_opostprob  <- mout$opostprob
    return_list$m_p1postprob <- mout$p1postprob
    return_list$m_p2postprob <- mout$p2postprob
  }

  return(return_list)
}


#' Iterate between indiividual level estimates and joint averaging.
#'
#' @inheritParams updog
#' @inheritParams updog_mcmc
#' @param itermax The number of times to iterate between the two approaches.
#' @param tol A positive numeric. The tolderance for the stopping criteria.
#'
#' @author David Gerard
#'
updog_iterate <- function(ocounts, osize, qarray, r1vec, r2vec, seq_error = 0.01, itermax = 1000,
                          tol = 10 ^ -4) {
  ploidy <- length(r1vec) - 1
  assertthat::are_equal(length(r2vec), ploidy + 1)
  # assertthat::are_equal(length(ocounts), length(osize))
  assertthat::assert_that(all(ocounts <= osize))
  assertthat::assert_that(all(dim(qarray) == ploidy + 1))
  assertthat::assert_that(tol > 0)

  ## calculate log-probabilities ---------------------------------------------
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  ## deal with error rate ----------------------------------------------------
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)

  ## binom density where the rows indexs the genotypes and the columns index the individuals
  dbinommat <- mapply(FUN = stats::dbinom, ocounts, osize, MoreArgs = list(prob = pk, log = FALSE))
  dimnames(dbinommat) = list(genotype = get_dimname(ploidy), offspring = 1:ncol(dbinommat))
  assertthat::assert_that(all(abs(dbinommat[, 1] -
                                    stats::dbinom(x = ocounts[1], size = osize[1], prob = pk, log = FALSE)) < 10 ^ -14))

  ## update parent 1 and parent
  err <- tol + 1
  while (index < itermax & err > tol) {
    r1old <- r1vec
    harray <- sweep(qarray, MARGIN = 1, STATS = r1vec, FUN = `*`)
    harray <- sweep(harray, MARGIN = 2, STATS = r2vec, FUN = `*`)
    h2array <- apply(dbinommat, 2, function(x) {sweep(harray, MARGIN = 3, STATS = x, FUN = `*`) })
    h3array <- array(h2array, dim = c(rep(ploidy + 1, 3), length(ocounts)))
    dimvec <- get_dimname(ploidy)
    dimnames(h3array) = list(parent1 = dimvec, parent2 = dimvec, offspring = dimvec, individual = 1:length(ocounts))
    up1prob <- apply(h3array, c(1, 4), sum)
    up2prob <- apply(h3array, c(2, 4), sum)

    p1prob <- sweep(up1prob, MARGIN = 2, STATS = colSums(up1prob), FUN = `/`)
    p2prob <- sweep(up2prob, MARGIN = 2, STATS = colSums(up2prob), FUN = `/`)

    r1vec <- rowMeans(p1prob)

    index <- index + 1
    err <- sum(abs(r1vec - r1old))

    cat(err, "\n")
  }
  return(list(r1vec = r1vec, r2vec = r2vec))
}

#' The updog Gibbs sampler to jointly estimate the parental and offspring genotypes.
#'
#' @inheritParams updog
#' @param qarray An array of proportions. Each dimension size is equal to the ploidy plus 1.
#' @param r1vec The posterior probabilities of the genotypes of parent 1 conditional only
#'   on the sequence data from parent 1.
#' @param r2vec The posterior probabilities of the genotypes of parent 2 conditional only
#'   on the sequence data from parent 2.
#'
#' @author David Gerard
#'
updog_mcmc <- function(ocounts, osize, qarray, r1vec, r2vec, seq_error = 0.01, itermax = 1000,
                       burnin = 250) {

  ploidy <- length(r1vec) - 1
  assertthat::are_equal(length(r2vec), ploidy + 1)
  assertthat::are_equal(length(ocounts), length(osize))
  assertthat::assert_that(all(ocounts <= osize))
  assertthat::assert_that(all(dim(qarray) == ploidy + 1))

  ## calculate log-probabilities ---------------------------------------------
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  ## deal with error rate ----------------------------------------------------
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)

  ## the total number of samples for each genotype
  tot_p1 <- rep(0, length = ploidy + 1)
  tot_p2 <- rep(0, length = ploidy + 1)
  tot_o  <- matrix(0, nrow = ploidy + 1, ncol = length(osize))

  current_p1 <- sample(0:ploidy, size = 1, prob = r1vec)
  current_p2 <- sample(0:ploidy, size = 1, prob = r2vec)

  ## dbinom matrix for offspring
  dbinommat <- mapply(FUN = stats::dbinom, ocounts, osize, MoreArgs = list(prob = pk, log = TRUE))
  assertthat::are_equal(dbinommat[, 1], stats::dbinom(x = ocounts[1], size = osize[1], prob = pk, log = TRUE))


  for (iterindex in 1:(itermax + burnin)) {

    ## Update offspring ------------------------------------------------------
    qvec <- log(qarray[current_p1 + 1, current_p2 + 1, ])
    uprob_log<- qvec + dbinommat
    max_prob <- apply(uprob_log[qvec != -Inf, ], 2, max)
    uprobexp <- exp((sweep(x = uprob_log, MARGIN = 2, STATS = max_prob, FUN = "-")))
    opmat <- sweep(x = uprobexp, MARGIN = 2, STATS = colSums(uprobexp), FUN = "/")
    ## assertthat::are_equal(opmat, sweep(exp(uprob_log), MARGIN = 2, STATS = colSums(exp(uprob_log)), FUN = "/"))
    current_o  <- apply(opmat, 2, function (x) { sample(0:ploidy, size = 1, prob = x) })
    if (iterindex > burnin) {
      tot_o[cbind(current_o + 1, 1:length(ocounts))] <- tot_o[cbind(current_o + 1, 1:length(ocounts))] + 1
    }

    ## Update parent1 --------------------------------------------------------
    qmat <- log(qarray[, current_p2 + 1, ])
    ek <- table(c(current_o, 0:ploidy)) - 1

    tmat1 <- sweep(x = qmat, MARGIN = 2, STATS = ek, FUN = `*`)
    tmat1[is.nan(tmat1)] <- 0
    uproblog <- rowSums(tmat1) + r1vec
    uproblog <- uproblog - max(uproblog[uproblog != -Inf])
    p1pvec <- exp(uproblog) / sum(exp(uproblog))

    current_p1 <- sample(0:ploidy, size = 1, prob = p1pvec)

    if (iterindex > burnin) {
      tot_p1[current_p1 + 1] <- tot_p1[current_p1 + 1] + 1
    }

    ## update parent2 --------------------------------------------------------
    qmat <- log(qarray[current_p1 + 1, , ])
    tmat2 <- sweep(x = qmat, MARGIN = 2, STATS = ek, FUN = `*`)
    tmat2[is.nan(tmat2)] <- 0
    uproblog <- rowSums(tmat2) + r2vec
    uproblog <- uproblog - max(uproblog[uproblog != -Inf])
    p2pvec <- exp(uproblog) / sum(exp(uproblog))

    current_p2 <- sample(0:ploidy, size = 1, prob = p2pvec)

    if (iterindex > burnin) {
      tot_p2[current_p2 + 1] <- tot_p2[current_p2 + 1] + 1
    }
  }


  p1postprob <- tot_p1 / itermax
  p2postprob <- tot_p2 / itermax
  opostprob  <- tot_o / itermax
  return_list <- list()
  return_list$p1postprob <- p1postprob
  return_list$p2postprob <- p2postprob
  return_list$opostprob  <- opostprob
  return(return_list)
}

#' Derive the probabilities of an offspring's genotype given its parental
#' genotypes for all possible combinations of parental and offspring genotypes.
#'
#' @param ploidy A positive integer. The ploidy of the species.
#'
#' @author David Gerard
#'
#' @return An three-way array of proportions. The (i, j, k)th element is the probability
#'   of an offspring having k - 1 reference alleles given that parent 1 has i - 1
#'   refrerence alleles and parent 2 has j - 1 reference alleles. Each dimension of
#'   the array is \code{ploidy + 1}. In the dimension names, "A" stands for the
#'   reference allele and "a" stands for any other allele.
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

#' Make those awesome plots that Felipe showed me, but using ggplot2.
#'
#' @param ocounts A vector of non-negative integers. The number of reference alleles observed in the offspring.
#' @param osize A vector of positive integers. The total number of reads in the offspring.
#' @param p1counts A vector of non-negative integers. The number of reference alleles observed in parent 1.
#' @param p1size A vector of positive integers. The total number of reads in parent 1.
#' @param p2counts A vector of non-negative integers. The number of reference alleles observed in parent 2.
#' @param p2size A vector of positive integers. The total number of reads in parent 2.
#' @param ploidy A non-negative integer. The ploidy of the species.
#' @param col The color labels.
#' @param seq_error The average sequencing error rate.
#'
#' @export
#'
#' @author David Gerard
#'
plot_geno <- function(ocounts, osize, ploidy, p1counts = NULL, p1size = NULL, p2counts = NULL,
                      p2size = NULL, col = NULL, seq_error = 0.01) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 must be installed to use this function")
  }

  assertthat::assert_that(all(ocounts >= 0))
  assertthat::assert_that(all(osize >= ocounts))
  assertthat::assert_that(ploidy >= 1)
  assertthat::assert_that(seq_error>= 0)


  ## get probabilities
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)

  dfdat <- data.frame(A = ocounts, a = osize - ocounts)
  maxcount <- max(max(dfdat$A), max(dfdat$a))
  if (!is.null(col)) {
    assertthat::are_equal(length(col), length(ocounts))
    dfdat$genotype <- col
  }

  slopevec <- pk / (1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1),
                         xend = xend, yend = yend)

  ## Plot children
   pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes(y = A, x = a)) +
     ggplot2::geom_point() +
     ggplot2::theme_bw() +
     ggplot2::xlim(0, maxcount) +
     ggplot2::ylim(0, maxcount) +
     ggplot2::ylab("Counts A") +
     ggplot2::xlab("Counts a")  +
     ggplot2::geom_segment(data = df_lines, mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                           lty = 2, alpha = 1/2)

   ## add parents if we have them
   if (!is.null(p1size) & !is.null(p1counts)) {
     assertthat::assert_that(all(p1counts >= 0))
     assertthat::assert_that(all(p1size >= p1counts))
     p1dat <- data.frame(A = p1counts, a = p1size - p1counts)
     pl <- pl + ggplot2::geom_point(data = p1dat, size = 3, color = "red", pch = 3)
   }
   if (!is.null(p2size) & !is.null(p2counts)) {
     assertthat::assert_that(all(p2counts >= 0))
     assertthat::assert_that(all(p2size >= p2counts))
     p2dat <- data.frame(A = p2counts, a = p2size - p2counts)
     pl <- pl + ggplot2::geom_point(data = p2dat, size = 3, color = "blue", pch = 4)
   }

   print(pl)
   return(pl)
}

#' Returns a vector character strings that are all of the possible combinations of the
#' reference allele and the non-reference allele.
#'
#' @param ploidy The ploidy of the species.
#'
#' @return For example, if \code{ploidy = 3} then this will return c("aaa", "Aaa", "AAa", "AAA")
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

#' Calculates posterior probabilities of a genotype given just the sequence
#' counts.
#'
#' The function assumes that all data come from independent samples on the same
#' individual.
#'
#' @param ncounts A vector of non-negative integers. The ith element is the
#'   number of counts of the ith sample.
#' @param ssize A vector of positive integers. The ith element is the total
#'   number of counts of the ith sample.
#' @param prior A vector of non-negative numerics that sum to one. The prior
#'   probability on the genotype. The first element is the prior probability of
#'   zero reference alleles, the second element is the prior probability of one
#'   reference allele, etc. The length of \code{prior} is one more than the
#'   ploidy of the species. You can alternatively specify \code{prior} as the
#'   ploidy of the individual, for which it will set a uniform prior on the
#'   genotype. For example, setting \code{prior = 3} will result in using
#'   \code{(1/4, 1/4, 1/4, 1/4)} as the prior probability for the genotypes
#'   (Aaaa, AAaa, AAAa, AAAA) where "A" is the reference allele in a 4-ploid
#'   individual.
#' @param seq_error A non-negative numeric. This is the known sequencing error
#'   rate. This is a rough high-ball error rate given by Li et. al. (2011).
#' @param divide_by_three A logical. Should we use asymmetry in the error rates
#'   (\code{TRUE}) or not (\code{FALSE})? Defaults to \code{FALSE}.
#'
#' @return A vector of probabilities. The ith element is the posterior probability
#'   that the individual has i - 1 copies of the reference allele.
#'
#' @author David Gerard
#'
#' @export
#'
#' @references Li, Yun, Carlo Sidore, Hyun Min Kang, Michael Boehnke, and
#'   Gonçalo R. Abecasis.
#'   \href{https://www.ncbi.nlm.nih.gov/pubmed/21460063}{"Low-coverage sequencing: implications for design of complex trait association studies."}
#'   Genome research (2011).
#'
#'
bin_post <- function(ncounts, ssize, prior, seq_error = 0.01, divide_by_three = FALSE) {

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
  if (!divide_by_three) {
    pk <- (1 - seq_error) * pk + seq_error * (1 - pk)
  } else {
    pk <- (1 - seq_error) * pk + seq_error * (1 - pk) / 3
  }



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
