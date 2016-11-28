## estimate sequencing error

#' Estimate sequencing error from counts.
#'
#' Assuming that the reference allele is in great enough quantities so that there are enough points
#'
#' @inheritParams bin_post
#' @param eps To estimate the sequencing erorr, we will only use points such that \code{ncounts / ssize >= eps}.
#'     If \code{ploidy} is specified and \code{eps = NULL}, then \code{eps} will be chosen by a heuristic.
#' @param ploidy The ploidy of the species. Cannot be \code{NULL} if \code{eps} is \code{NULL}.
#'
#' @author David Gerard
#'
#' @export
#'
est_seq_error <- function(ncounts, ssize, ploidy = NULL, eps = NULL) {

  assertthat::assert_that(all(ncounts >= 0))
  assertthat::assert_that(all(ssize >= ncounts))

  if (is.null(ploidy) & is.null(eps)) {
    stop("ploidy and eps cannot both be NULL.")
  } else if (is.null(eps)) {
    eps <- 0.3 / ploidy ## 20 percent between the theoretical no a and one a.
  } else if (!is.null(ploidy) & !is.null(eps)) {
    if (eps > 0.5 / ploidy) {
      warning("eps is on the same order as 1/ploidy. You should probably make eps smaller.")
    }
    message("Using eps even though ploidy is specified.")
  }

  assertthat::assert_that(eps >= 0)

  phat <- ncounts / ssize

  ## first way
  which_A <- phat >= (1 - eps)

  if (sum(which_A) == 0) {
    stop("Not enough points to estimate sequencing error.")
  } else if (sum(which_A) <= 5) {
    warning("Very few points to estimate sequencing error.")
  }

  seq_error <- 1 - sum(ncounts[which_A]) / sum(ssize[which_A])

  ## second way
  # varphat <- phat * (1 - phat) / ssize
  # lower_limit <- phat - 1.96 * sqrt(varphat)
  # lower_limit >= (1 - eps) * (1 - 1 / ploidy) + eps * (1 / ploidy)

  return(seq_error)
}
