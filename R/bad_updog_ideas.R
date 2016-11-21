
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
  index <- 1
  while (index < itermax & err > tol) {
    r1old <- r1vec
    r2old <- r2vec
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
    r2vec <- rowMeans(p2prob)

    index <- index + 1
    err <- sum(abs(r1vec - r1old)) + sum(abs(r2vec - r2old))

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
#' @param itermax The maximum number of iterations to run the Gibbs sampler.
#' @param burnin The number of iterations to skip.
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
