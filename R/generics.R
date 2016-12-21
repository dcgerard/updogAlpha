### summary and plot generics.

#' Draw a genotype plot from the output of \code{\link{updog}}.
#'
#' Iff ggplot2 is installed, then this function will use it to make the genotype plot.
#' Otherwise, "classic" R base graphics will be used. The ggplot2 version is made using the function
#' \code{\link{plot_geno}}.
#'
#' @param x The output from \code{\link{updog}}.
#' @param gg Should we use ggplot2 to plot the genotypes (\code{TRUE}), or not (\code{FALSE}).
#'     If ggplot2 is present, then this defaults to \code{TRUE}. If it is not present, then it
#'     defaults to \code{FALSE}.
#' @param plot_beta A logical. If true, then we'll also plot the estimated beta density of the outlier model.
#' @param ... Not used.
#'
#' @return A plot object.
#'
#' @author David Gerard
#'
#' @export
#'
plot.updog <- function(x, gg = requireNamespace("ggplot2", quietly = TRUE), plot_beta = TRUE, ...) {
  assertthat::assert_that(is.updog(x))
  assertthat::assert_that(is.logical(plot_beta))
  assertthat::assert_that(is.logical(gg))

  if (!is.null(x$opostprob)) {
    maxpostprob <- apply(x$opostprob, 2, max)
  } else {
    maxpostprob <- NULL
  }

  if (requireNamespace("ggplot2", quietly = TRUE) & gg) {
    pl <- plot_geno(ocounts = x$input$ocounts, osize = x$input$osize,
                    p1counts = x$input$p1counts, p1size = x$input$p1size,
                    p2counts = x$input$p2counts, p2size = x$input$p2size,
                    ploidy = x$input$ploidy, ogeno = x$ogeno, seq_error = x$seq_error,
                    prob_ok = x$prob_ok, maxpostprob = maxpostprob,
                    p1geno = x$p1geno, p2geno = x$p2geno)
    print(pl)
  } else if (!gg) {
    plot_geno_base(ocounts = x$input$ocounts, osize = x$input$osize,
                   p1counts = x$input$p1counts, p1size = x$input$p1size,
                   p2counts = x$input$p2counts, p2size = x$input$p2size,
                   ploidy = x$input$ploidy, ogeno = x$ogeno, seq_error = x$seq_error,
                   prob_ok = x$prob_ok, maxpostprob = maxpostprob,
                   p1geno = x$p1geno, p2geno = x$p2geno)
  } else {
    stop("gg = TRUE but ggplot2 not installed")
  }

  if (plot_beta) {
    cat ("Press [enter] to continue")
    line <- readline()
    if (!is.null(x$out_mu) & !is.null(x$out_rho)) {
      plot_beta_dist(mu = x$out_mu, rho = x$out_rho)
    } else if (!is.null(x$alpha) & !is.null(x$beta)) {
      plot_beta_dist(alpha = x$alpha, beta = x$beta)
    } else {
      message("No outlier distribution to plot.")
    }
  }

}

#' Summary method for class "\code{updog}".
#'
#' @param object An \code{updog} object. Usually what has been returned from \code{\link{updog}}.
#' @param ... Not used.
#'
#' @return A list with three elements:
#'
#'     \code{prop_ok}: The estimated proportion of observations that are outliers.
#'
#'     \code{genotypes}: Counts for how many of each estimated genotype are observed.
#'
#'     \code{summ_prob}: A data frame with summaries on two variables: The maximum posterior
#'     probability of a genotype that each observation has (\code{maxpostprob}), and the posterior probability
#'     of not being an outlier (\code{prob_ok})
#'
#' @author David Gerard
#'
#' @export
#'
summary.updog <- function(object, ...) {
  genotypes <- table(c(object$ogeno, 0:object$input$ploidy)) - 1
  maxpostprob_all <- apply(object$opostprob, 2, max)
  maxpostprob <- summary(maxpostprob_all)
  prob_ok <- summary(object$prob_ok)

  summ_prob <- cbind(maxpostprob, prob_ok)

  return_list <- list()
  return_list$prop_ok <- object$pival
  return_list$genotypes    <- genotypes
  return_list$summ_prob    <- summ_prob
  return(return_list)
}

#' Plot the beta distribution.
#'
#' This function will plot the beta distribution under two different
#' parameterizations. The first one, with \code{alpha} and \code{beta},
#' is the usual parameterization of the beta. The second one parameterizes
#' the beta in terms of its mean and overdispersion. Thus,
#' either both \code{alpha} and \code{beta} must be specified,
#' or both \code{mu} and \code{rho} must be specified.
#'
#' We have the relationships
#' \deqn{\mu = \alpha / (\alpha + \beta)}
#' and
#' \deqn{\rho = 1 / (1 + \alpha + \beta).}
#' Thus, the inverse relationships are
#' \deqn{\alpha = \mu(1 -\rho) / \rho}
#' and
#' \deqn{\beta = (1 - \mu)(1 - \rho) / \rho.}
#'
#' @param alpha The shape1 parameter.
#' @param beta The shape2 parameter.
#' @param mu The mean parameter.
#' @param rho The overdispersion parameter.
#' @param ... Anything else you want to pass to \code{\link[graphics]{plot.default}}.
#'
#' @author David Gerard
#'
#' @export
#'
plot_beta_dist <- function(alpha = NULL, beta = NULL, mu = NULL, rho = NULL, ...) {
  if(!((!is.null(alpha) & !is.null(beta)) | (!is.null(mu) & !is.null(rho)))) {
    stop("you need to specify either both alpha and beta or both mu and rho")
  }
  if(!is.null(alpha) & !is.null(beta) & !is.null(mu) & !is.null(rho)) {
    stop("you can't specify all alpha and beta and mu and rho")
  }

  if (!is.null(mu) & !is.null(rho)) {
    assertthat::assert_that(mu >= 0)
    assertthat::assert_that(mu <= 1)
    assertthat::assert_that(rho > 0)
    assertthat::assert_that(rho < 1)

    alpha <- mu * (1 - rho) / rho
    beta  <- (1 - mu) * (1 - rho) / rho
  } else {
    assertthat::assert_that(alpha > 0)
    assertthat::assert_that(beta > 0)
    mu  <- alpha / (alpha + beta)
    rho <- 1 / (1 + alpha + beta)
  }

  ## Find x boundaries that cover 98% of area
  xlower <- stats::qbeta(0.01, shape1 = alpha, shape2 = beta)
  xupper <- stats::qbeta(0.99, shape1 = alpha, shape2 = beta)

  x <- seq(xlower, xupper, length = 500)
  y <- stats::dbeta(x = x, shape1 = alpha, shape2 = beta)

  old_options <- graphics::par(no.readonly = TRUE) ## save current options
  graphics::par(mar = c(3, 3, 0.5, 0.5), mfrow = c(1, 1), mgp = c(2, 1, 0))
  graphics::plot(x, y, type = "l", xlab = "quantile", ylab = "density", lwd = 2, col = "gray25", ...)
  if (mu <= 1/2 & (alpha >= 1 & beta >= 1)) {
    graphics::legend("topright", bty = "n",
                     legend = c(paste0("mu=", format(mu, digits = 2)), paste0("rho=", format(rho, digits = 2))))
  } else if (mu > 1/2 & (alpha >= 1 & beta >= 1)) {
    graphics::legend("topleft", bty = "n",
                     legend = c(paste0("mu=", format(mu, digits = 2)), paste0("rho=", format(rho, digits = 2))))
  } else {
    graphics::legend("top", bty = "n",
                     legend = c(paste0("mu=", format(mu, digits = 2)), paste0("rho=", format(rho, digits = 2))))
  }
  on.exit(graphics::par(old_options), add = TRUE)
}

#' The base R graphics version of \code{\link{plot_geno}}.
#'
#' @inheritParams plot_geno
#'
#' @export
#'
#' @seealso \code{\link{plot_geno}}
#'
#' @author David Gerard
#'
plot_geno_base <- function(ocounts, osize, ploidy, p1counts = NULL, p1size = NULL, p2counts = NULL,
                           p2size = NULL, ogeno = NULL, seq_error = 0.01, prob_ok = NULL, maxpostprob = NULL,
                           p1geno = NULL, p2geno = NULL) {

  ## copied from plot_geno. Probably a better way to do this
  assertthat::assert_that(all(ocounts >= 0))
  assertthat::assert_that(all(osize >= ocounts))
  assertthat::assert_that(ploidy >= 1)
  assertthat::assert_that(seq_error>= 0)

  ## get probabilities
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)

  dfdat <- data.frame(A = ocounts, a = osize - ocounts)
  maxcount <- max(max(dfdat$A), max(dfdat$a))
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


  col_vec <- rep(NA, length = length(ogeno))
  ## Deal with colors and transparancies
  possible_colors <- grDevices::rainbow(ploidy + 1, alpha = 1)
  if (!is.null(ogeno)) {
    col_vec <- possible_colors[ogeno + 1]
  } else {
    col_vec <- rep("#000000FF", length = length(ocounts))
  }

  if (!is.null(prob_ok)) {
    possible_transparancies <- seq(0, max(prob_ok), length = 5)
    hexmode_trans <- as.hexmode(round(255 * possible_transparancies))
    point_trans <- as.character(hexmode_trans[.bincode(prob_ok, breaks = possible_transparancies)])

    col_vec <- mapply(FUN = sub, x = col_vec, replace = point_trans, MoreArgs = list(pattern = "FF$"))
    names(col_vec) <- NULL
  }

  ## Plot and deal with sizes -----------------------------------------------
  old_options <- graphics::par(no.readonly = TRUE) ## save current options

  graphics::par(mar = c(3, 3, 0.5, 0.5), mfrow = c(1, 2), mgp = c(2, 1, 0))


  graphics::plot(NULL, xlim = c(0, maxcount), ylim = c(0, maxcount), xlab = "Counts a", ylab = "Counts A")
  graphics::arrows(x0 = df_lines$x, y0 = df_lines$y, x1 = df_lines$xend, y1 = df_lines$yend,
                   length = 0, col = "grey50", lty = 2)

  cex_val <- 1.5
  if (is.null(maxpostprob)) {
    graphics::points(x = dfdat$a, y = dfdat$A, pch = 16, col = col_vec)
  } else if (!is.null(maxpostprob)) {
    graphics::points(x = dfdat$a, y = dfdat$A, pch = 16, cex = dfdat$maxpostprob * cex_val, col = col_vec)
  }

  ## Parental data
  if (!is.null(p1counts) & !is.null(p1size)) {
    assertthat::assert_that(all(p1counts >= 0))
    assertthat::assert_that(all(p1size >= p1counts))
    if (!is.null(p1geno)) {
      p1colvec <- rep(possible_colors[p1geno + 1], length = length(p1counts))
      graphics::points(x = p1size - p1counts, y = p1counts, pch = 3, col = p1colvec)
      graphics::points(x = p1size - p1counts, y = p1counts, pch = 16, col = "black", cex = 0.3)
    } else {
      graphics::points(x = p1size - p1counts, y = p1counts, pch = 3)
    }
  }

  if (!is.null(p2counts) & !is.null(p2size)) {
    assertthat::assert_that(all(p2counts >= 0))
    assertthat::assert_that(all(p2size >= p2counts))
    if (!is.null(p2geno)) {
      p2colvec <- rep(possible_colors[p2geno + 1], length = length(p2counts))
      graphics::points(x = p2size - p2counts, y = p2counts, pch = 3, col = p2colvec)
      graphics::points(x = p2size - p2counts, y = p2counts, pch = 16, col = "black", cex = 0.3)
    } else {
      graphics::points(x = p2size - p2counts, y = p2counts, pch = 3)
    }
  }

  ## Legends depending on what was provided.
  graphics::plot.new()
  if (!is.null(ogeno) | !is.null(p1geno) | !is.null(p2geno)) {
    graphics::legend("topright", pch = 16, col = possible_colors, legend = 0:ploidy, title = "Genotype",
                     bty = "n")
  }
  if (!is.null(prob_ok)) {
    graphics::legend("bottomleft", pch = 16, col = paste0("#000000", hexmode_trans), title = "prob_ok",
                     legend = format(possible_transparancies, digits = 2), bty = "n")
  }
  if (!is.null(maxpostprob)) {
    maxpostprob_legend <- seq(0.1, max(maxpostprob), length = 4)
    maxpostprob_val <- maxpostprob_legend * cex_val
    graphics::legend("topleft", col = "black", pch = 16, pt.cex = maxpostprob_val,
                     legend = format(maxpostprob_legend, digits = 2), title = "maxpostprob", bty = "n")
  }

  ## reset old options before exiting
  on.exit(graphics::par(old_options), add = TRUE)
}


#' Tests if its argument is an updog.
#'
#' @param x Anything.
#'
#' @return A logical. \code{TRUE} if \code{x} is an updog, and \code{FALSE} otherwise.
#'
#' @author David Gerard
#'
#' @export
#'
is.updog <- function(x) inherits(x, "updog")
