

#' Using Parental DNA for Offspring Genotyping.
#'
#' This function takes parental and offspring sequence counts and returns
#' posterior probabilities on the genotype of all offspring. This isn't fully
#' Bayesian because I only update the parental probabilities with just their
#' sequencing data.
#'
#' What is updog? Nothing much, what's up with you dawg!
#'
#' @param ocounts A vector of non-negative integers. The ith element is the
#'   number of reads of the reference allele the ith child.
#' @param osize A vector of positive integers. The ith element is the total
#'   number of reads for the ith child.
#' @param p1counts A vector of non-negative integers. The ith element is the
#'   number of reads of the reference allele in the ith sample of parent 1.
#' @param p1size A vector of positive integers. The ith element is the total
#'   number of reads in the ith sample of parent 1.
#' @param p2counts A vector of non-negative integers. The ith element is the
#'   number of reads of the reference allele in the ith sample of parent 2.
#' @param p2size A vector of positive integers. The ith element is the total
#'   number of reads in the ith sample of parent 2.
#' @param ploidy A positive integer. The number of copies of the genome in the
#'   species.
#' @param seq_error A proportion. The known sequencing error rate.
#'
#' @author David Gerard
#'
#' @export
#'
updog <- function(ocounts, osize, p1counts, p1size, p2counts, p2size, ploidy,
                  seq_error = 0.1) {

  ## check input -------------------------------------------------------------
  assertthat::assert_that(all(ocounts >= 0))
  assertthat::assert_that(all(p1counts >= 0))
  assertthat::assert_that(all(p2counts >= 0))
  assertthat::assert_that(all(osize >= 1))
  assertthat::assert_that(all(p1size >= 1))
  assertthat::assert_that(all(p2size >= 1))
  assertthat::assert_that(all(ocounts <= osize))
  assertthat::assert_that(all(p1counts <= p1size))
  assertthat::assert_that(all(p2counts <= p2size))
  assertthat::assert_that(ploidy >= 1)

  ## derive posteriors of parental genotypes given just parental sequence data.
  r1vec <- bin_post(ncounts = p1counts, ssize = p1size, prior = ploidy)
  r2vec <- bin_post(ncounts = p2counts, ssize = p2size, prior = ploidy)

  ## derive offspring genotype probabilities given parental genotypes.
  qarray <- get_q_array(ploidy = ploidy)

  ## Derive prior probabilities on offspring genotypes
  phi_vec <- r1vec ## posterior prob of p1 genotype
  psi_vec <- r2vec ## posterior prob of p2 genotype

  harray <- sweep(qarray, MARGIN = 1, STATS = phi_vec, FUN = `*`)
  harray <- sweep(harray, MARGIN = 2, STATS = psi_vec, FUN = `*`)
  hl <- apply(harray, 3, sum)
  postprob <- mapply(FUN = bin_post, ocounts, osize,
                     MoreArgs = list(prior = c(hl), seq_error = seq_error))

  return_list <- list()
  return_list$opostprob <- postprob
  return_list$p1postprob <- phi_vec
  return_list$p2postprob <- psi_vec

  return(return_list)
}

#' Derive the probabilities of an offspring's genotype given its parental
#' genotypes for all possible combinations of parental and offspring genotypes.
#'
#' @param ploidy A positive integer. The ploidy of the species.
#'
#' @author David Gerard
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

  assertthat::assert_that(all(abs(apply(qarray, c(1, 2), sum) - 1) < 10 ^ -14))

  return(qarray)
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
#'   rate. This is a course high-ball error rate given by Li et. al. (2011).
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
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk) / 3

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
