### summary and plot generics.

#' Draw a genotype plot from the output of \code{\link{updog}}.
#'
#' Iff ggplot2 is installed, then this function will use it to make
#' the genotype plot.  Otherwise, "classic" R base graphics will be
#' used. The ggplot2 version is made using the function
#' \code{\link{plot_geno}}.
#'
#' @param x The output from \code{\link{updog}}.
#' @param gg Should we use ggplot2 to plot the genotypes
#' (\code{TRUE}), or not (\code{FALSE}).  If ggplot2 is present, then
#' this defaults to \code{TRUE}. If it is not present, then it
#' defaults to \code{FALSE}.
#' @param plot_beta A logical. If true, then we'll also plot the
#' estimated beta density of the outlier model, followed by the
#' estimated beta distributions of the overdispersion models.
#' @param ask A logical. Should we ask before continuing on to the
#'     next plot (\code{TRUE}) or not (\code{FALSE})?
#' @param use_colorblind A logical. Should we use a colorblind safe palette (\code{TRUE}),
#'     or not (\code{FALSE})?
#' @param show_maxpostprob A logical. Should we scale the sizes by \code{x$maxpostprob} (\code{TRUE}),
#'     or not (\code{FALSE})?
#' @param show_ogeno A logical. Should we color code by \code{x$ogeno} (\code{TRUE}),
#'     or not (\code{FALSE})?
#' @param show_outlier A logical. Should we scale the transparency by \code{x$prob_ok} (\code{TRUE}),
#'     or not (\code{FALSE})?
#' @param ... Not used.
#'
#' @return A plot object.
#'
#' @author David Gerard
#'
#' @export
#'
plot.updog <- function(x, gg = requireNamespace("ggplot2", quietly = TRUE),
                       plot_beta = TRUE, ask = TRUE,
                       use_colorblind = FALSE,
                       show_maxpostprob = FALSE,
                       show_ogeno = TRUE,
                       show_outlier = TRUE, ...) {

  assertthat::assert_that(is.updog(x))
  assertthat::assert_that(is.logical(plot_beta))
  assertthat::assert_that(is.logical(gg))
  assertthat::assert_that(is.logical(ask))
  assertthat::assert_that(is.logical(use_colorblind))
  assertthat::assert_that(is.logical(show_maxpostprob))
  assertthat::assert_that(is.logical(show_ogeno))
  assertthat::assert_that(is.logical(show_outlier))

  if (!show_maxpostprob) {
    x$maxpostprob <- NULL
  }
  if (!show_ogeno) {
    x$ogeno <- NULL
  }
  if (!show_outlier) {
    x$prob_ok <- NULL
  }

  if (requireNamespace("ggplot2", quietly = TRUE) & gg) {
    pl <- plot_geno(ocounts = x$input$ocounts, osize = x$input$osize,
                    p1counts = x$input$p1counts, p1size = x$input$p1size,
                    p2counts = x$input$p2counts, p2size = x$input$p2size,
                    ploidy = x$input$ploidy, ogeno = x$ogeno, seq_error = x$seq_error,
                    bias_val = x$bias_val,
                    prob_ok = x$prob_ok, maxpostprob = x$maxpostprob,
                    p1geno = x$p1geno, p2geno = x$p2geno, use_colorblind = use_colorblind)
    print(pl)
  } else if (!gg) {
    plot_geno_base(ocounts = x$input$ocounts, osize = x$input$osize,
                   p1counts = x$input$p1counts, p1size = x$input$p1size,
                   p2counts = x$input$p2counts, p2size = x$input$p2size,
                   ploidy = x$input$ploidy, ogeno = x$ogeno, seq_error = x$seq_error,
                   bias_val = x$bias_val,
                   prob_ok = x$prob_ok, maxpostprob = x$maxpostprob,
                   p1geno = x$p1geno, p2geno = x$p2geno, use_colorblind = use_colorblind)
  } else {
    stop("gg = TRUE but ggplot2 not installed")
  }

  ## rename od_param, out_mean and out_disp because of legacy code.
  if (!is.null(x$out_mean)) {
    x$out_mu <- x$out_mean
  }
  if (!is.null(x$out_disp)) {
    x$out_rho <- x$out_disp
  }
  if (!is.null(x$od_param)) {
    x$rho <- x$od_param
  }

  if (plot_beta) {
    if (ask) {
      cat ("Press [enter] to continue")
      line <- readline()
    }
    if (!is.null(x$out_mu) & !is.null(x$out_rho)) {
      if (gg) {
        plot_beta_dist_gg(mu = x$out_mu, rho = x$out_rho)
      } else {
        plot_beta_dist(mu = x$out_mu, rho = x$out_rho)
      }
    } else if (!is.null(x$alpha) & !is.null(x$beta)) {
      if (gg) {
        plot_beta_dist_gg(alpha = x$alpha, beta = x$beta)
      } else {
        plot_beta_dist(alpha = x$alpha, beta = x$beta)
      }
    } else {
      message("No outlier distribution to plot.")
    }

    if (!is.null(x$rho)) {
      if (ask) {
        cat ("Press [enter] to continue")
        line <- readline()
      }
        pk <- seq(0, x$input$ploidy) / x$input$ploidy ## the possible probabilities
        pk <- (1 - x$seq_error) * pk + x$seq_error * (1 - pk)
        if (gg) {
          plot_beta_dist_gg(mu = pk, rho = rep(x$rho, x$input$ploidy + 1))
        } else {
          plot_beta_dist(mu = pk, rho = rep(x$rho, x$input$ploidy + 1))
        }
    } else {
        message("No overdispersion distribution to plot.")
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

  ## beta density if provided and txtplot installed --------------------------
  if (!is.null(object$out_mu) & !is.null(object$out_rho) & requireNamespace("txtplot", quietly = TRUE)) {
    alpha <- object$out_mu * (1 - object$out_rho) / object$out_rho
    beta  <- (1 - object$out_mu) * (1 - object$out_rho) / object$out_rho
    x <- seq(0.015, 0.985, length = 100)
    y <- stats::dbeta(x = x, shape1 = alpha, shape2 = beta)
    cat("Outlier Distribution:\n")
    txtplot::txtplot(x = x, y = y)
  }


  return(return_list)
}


#' Plot the beta distribution using ggplot2.
#'
#' This is the same as \code{\link{plot_beta_dist}}, except we use \code{ggplot2} to do the plotting.
#' See \code{\link{plot_beta_dist}} for details.
#'
#' @inheritParams plot_beta_dist
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso \code{\link{plot_beta_dist}} for the R Base graphics version of this function.
#'
plot_beta_dist_gg <- function(alpha = NULL, beta = NULL, mu = NULL, rho = NULL) {
  if(!((!is.null(alpha) & !is.null(beta)) | (!is.null(mu) & !is.null(rho)))) {
    stop("you need to specify either both alpha and beta or both mu and rho")
  }
  if(!is.null(alpha) & !is.null(beta) & !is.null(mu) & !is.null(rho)) {
    stop("you can't specify all alpha and beta and mu and rho")
  }

  if (!is.null(mu) & !is.null(rho)) {
    assertthat::assert_that(all(mu >= 0))
    assertthat::assert_that(all(mu <= 1))
    assertthat::assert_that(all(rho > 0))
    assertthat::assert_that(all(rho < 1))
    assertthat::are_equal(length(mu), length(rho))

    alpha <- mu * (1 - rho) / rho
    beta  <- (1 - mu) * (1 - rho) / rho
  } else {
    assertthat::assert_that(all(alpha > 0))
    assertthat::assert_that(all(beta > 0))
    assertthat::are_equal(length(alpha), length(beta))
    mu  <- alpha / (alpha + beta)
    rho <- 1 / (1 + alpha + beta)
  }

  ## Find x boundaries
  xlower <- 0.015
  xupper <- 0.985

  for (index in 1:length(mu)) {
    x <- seq(xlower, xupper, length = 500)
    y <- stats::dbeta(x = x, shape1 = alpha[index], shape2 = beta[index])
    dfdat_temp <- data.frame(quantile = x, density = y)
    dfdat_temp$mu  <- mu[index]
    dfdat_temp$rho <- rho[index]
    if (index == 1) {
      dfdat <- dfdat_temp
    } else {
      dfdat <- rbind(dfdat, dfdat_temp)
    }
  }

  dfdat$murho <- paste0("mu=", format(dfdat$mu, digits = 2), ", rho=", format(dfdat$rho, digits = 2))

  pl <- ggplot2::ggplot(data = dfdat, mapping = ggplot2::aes_string(x = "quantile", y = "density")) +
    ggplot2::facet_wrap(~murho) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
    ggplot2::geom_vline(mapping = ggplot2::aes_string(xintercept = "mu"), lty = 2, col = "gray50")
  print(pl)

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
    assertthat::assert_that(all(mu >= 0))
    assertthat::assert_that(all(mu <= 1))
    assertthat::assert_that(all(rho > 0))
    assertthat::assert_that(all(rho < 1))
    assertthat::are_equal(length(mu), length(rho))

    alpha <- mu * (1 - rho) / rho
    beta  <- (1 - mu) * (1 - rho) / rho
  } else {
    assertthat::assert_that(all(alpha > 0))
    assertthat::assert_that(all(beta > 0))
    assertthat::are_equal(length(alpha), length(beta))
    mu  <- alpha / (alpha + beta)
    rho <- 1 / (1 + alpha + beta)
  }

  ## Find x boundaries
  xlower <- 0.015
  xupper <- 0.985

  old_options <- graphics::par(no.readonly = TRUE) ## save current options

  ncol <- round(sqrt(length(mu)))
  nrow <- ceiling(length(mu) / ncol)
  graphics::par(mfrow = c(nrow, ncol))

  ymax <- NA
  for (index in 1:length(mu)) {
    x <- seq(xlower, xupper, length = 500)
    y <- stats::dbeta(x = x, shape1 = alpha[index], shape2 = beta[index])
    ymax <- max(c(ymax, y), na.rm = TRUE)
  }

  for (index in 1:length(mu)) {
      x <- seq(xlower, xupper, length = 500)
      y <- stats::dbeta(x = x, shape1 = alpha[index], shape2 = beta[index])
      graphics::par(mar = c(3, 3, 0.5, 0.5), mgp = c(2, 1, 0))
      graphics::plot(x, y, type = "l", xlab = "quantile", ylab = "density",
                     lwd = 2, col = "gray25", ylim = c(0, ymax), ...)
      graphics::abline(v = mu[index], col = "gray25", lty = 2)
      if (mu[index] <= 1/2 & (alpha[index] >= 1 & beta[index] >= 1)) {
          graphics::legend("topright", bty = "n",
                           legend = c(paste0("mu=", format(mu[index], digits = 2)),
                           paste0("rho=", format(rho[index], digits = 2))))
      } else if (mu[index] > 1/2 & (alpha[index] >= 1 & beta[index] >= 1)) {
          graphics::legend("topleft", bty = "n",
                           legend = c(paste0("mu=", format(mu[index], digits = 2)),
                           paste0("rho=", format(rho[index], digits = 2))))
      } else {
          graphics::legend("top", bty = "n",
                           legend = c(paste0("mu=", format(mu[index], digits = 2)),
                           paste0("rho=", format(rho[index], digits = 2))))
      }
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
                           p2size = NULL, ogeno = NULL, seq_error = 0, bias_val = 1, prob_ok = NULL, maxpostprob = NULL,
                           p1geno = NULL, p2geno = NULL, use_colorblind = TRUE) {

  if (ploidy > 6 & use_colorblind) {
    warning("use_colorblind is only supported when ploidy <= 6")
    use_colorblind <- FALSE
  }

  if (is.null(seq_error)) {
    seq_error <- 0
  }
  if (is.null(bias_val)) {
    bias_val <- 1
  }


  ## copied from plot_geno. Probably a better way to do this
  assertthat::assert_that(all(ocounts >= 0))
  assertthat::assert_that(all(osize >= ocounts))
  assertthat::assert_that(ploidy >= 1)
  assertthat::assert_that(seq_error>= 0)

  ## get probabilities
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)
  pk <- pk / (bias_val * (1 - pk) + pk)

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
    assertthat::assert_that(all(maxpostprob >= 0, na.rm = TRUE))
    assertthat::assert_that(all(maxpostprob <= 1, na.rm = TRUE))
    dfdat$maxpostprob <- maxpostprob
  }

  slopevec <- pk / (1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1),
                         xend = xend, yend = yend)



  ## Deal with colors and transparancies -------------------------------------
  if (use_colorblind & requireNamespace("ggthemes", quietly = TRUE)) {
    possible_colors <- paste0(ggthemes::colorblind_pal()(ploidy + 1), "FF")
  } else if (use_colorblind) {
    possible_colors <- grDevices::rainbow(ploidy + 1, alpha = 1)
    warning("ggthemes needs to be installed to set use_colorblind = TRUE.")
  } else {
    possible_colors <- grDevices::rainbow(ploidy + 1, alpha = 1)
  }

  col_vec <- rep(NA, length = length(ogeno))
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
    maxpostprob_legend <- c(0.25, 0.5, 0.75, 1)
    maxpostprob_val <- maxpostprob_legend * cex_val
    graphics::legend("topleft", col = "black", pch = 16, pt.cex = maxpostprob_val,
                     legend = format(maxpostprob_legend, digits = 2), title = "maxpostprob", bty = "n")
  }

  ## reset old options before exiting
  on.exit(graphics::par(old_options), add = TRUE)
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
#' @param bias_val The bias parameter.
#' @param prob_ok A vector of posterior probabilities of not being a mistake.
#' @param maxpostprob A vector of the posterior probabilities of begin at the modal probability.
#' @param p1geno Parent 1's genotype.
#' @param p2geno Parent 2's genotype.
#' @param use_colorblind A logical. Should we use a colorblind safe palette (\code{TRUE}),
#'     or not (\code{FALSE})?
#'
#' @export
#'
#' @author David Gerard
#'
plot_geno <- function(ocounts, osize, ploidy, p1counts = NULL, p1size = NULL, p2counts = NULL,
                      p2size = NULL, ogeno = NULL, seq_error = 0, bias_val = 1,
                      prob_ok = NULL, maxpostprob = NULL,
                      p1geno = NULL, p2geno = NULL, use_colorblind = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 must be installed to use this function")
  }

  if (is.null(seq_error)) {
    seq_error <- 0
  }
  if (is.null(bias_val)) {
    bias_val <- 1
  }

  if (ploidy > 6 & use_colorblind) {
    warning("use_colorblind is only supported when ploidy <= 6")
    use_colorblind <- FALSE
  }

  assertthat::assert_that(all(ocounts >= 0, na.rm = TRUE))
  assertthat::assert_that(all(osize >= ocounts, na.rm = TRUE))
  assertthat::assert_that(ploidy >= 1)
  assertthat::assert_that(seq_error >= 0)

  ## get probabilities
  pk <- seq(0, ploidy) / ploidy ## the possible probabilities
  pk <- (1 - seq_error) * pk + seq_error * (1 - pk)
  pk <- pk / (bias_val * (1 - pk) + pk)

  dfdat <- data.frame(A = ocounts, a = osize - ocounts)
  maxcount <- max(max(dfdat$A, na.rm = TRUE), max(dfdat$a, na.rm = TRUE)) + 1
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
    assertthat::assert_that(all(maxpostprob >= 0, na.rm = TRUE))
    assertthat::assert_that(all(maxpostprob <= 1, na.rm = TRUE))
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
    pl <- pl + ggplot2::geom_point(ggplot2::aes_string(color = "genotype"))
  } else {
    pl <- pl + ggplot2::geom_point()
  }



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
    if (any(p1counts > maxcount | p1size - p1counts > maxcount, na.rm = TRUE)) {
      bad_parent_points <- p1counts > maxcount | p1size - p1counts > maxcount
      bad_parent_points[is.na(bad_parent_points)] <- FALSE
      parent_mult <- (maxcount - 1) / pmax(p1counts[bad_parent_points], p1size[bad_parent_points] - p1counts[bad_parent_points])
      p1counts[bad_parent_points] <- parent_mult * p1counts[bad_parent_points]
      p1size[bad_parent_points]   <- parent_mult * p1size[bad_parent_points]
      pl <- pl + ggplot2::annotate("text", x = p1size[bad_parent_points] - p1counts[bad_parent_points], y = p1counts[bad_parent_points], label = "(scaled)", alpha = 0.3)
    }
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
    if (any(p2counts > maxcount | p2size - p2counts > maxcount, na.rm = TRUE)) {
      bad_parent_points <- p2counts > maxcount | p2size - p2counts > maxcount
      bad_parent_points[is.na(bad_parent_points)] <- FALSE
      parent_mult <- (maxcount - 1) / pmax(p2counts[bad_parent_points], p2size[bad_parent_points] - p2counts[bad_parent_points])
      p2counts[bad_parent_points] <- parent_mult * p2counts[bad_parent_points]
      p2size[bad_parent_points]   <- parent_mult * p2size[bad_parent_points]
      pl <- pl + ggplot2::annotate("text", x = p2size[bad_parent_points] - p2counts[bad_parent_points], y = p2counts[bad_parent_points], label = "(scaled)", alpha = 0.3)
    }
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


  ## Set color scale based on use_colorblind --------------------------------------
  if (!is.null(ogeno) | !is.null(p1geno) | !is.null(p2geno)) {
    if (use_colorblind & requireNamespace("ggthemes", quietly = TRUE)) {
      pl <- pl + ggthemes::scale_color_colorblind(drop = FALSE)
    } else if (use_colorblind) {
      pl <- pl + ggplot2::scale_color_hue(drop = FALSE)
      warning("ggthemes needs to be installed to set use_colorblind = TRUE.")
    } else {
      pl <- pl + ggplot2::scale_color_hue(drop = FALSE)
    }
  }

  ## Set transparency scale --------------------------------------------------------
  if (!is.null(prob_ok)) {
    pl <- pl + ggplot2::scale_alpha_continuous(breaks = c(-0.01, 0.25, 0.5, 0.75, 1.01))
  }

  ## Set size scale -----------------------------------------------------------------
  if (!is.null(maxpostprob)) {
    pl <- pl + ggplot2::scale_size_continuous(breaks = c(-0.01, 0.25, 0.5, 0.75, 1.01),
                                              range = c(0.5, 3))
  }


  return(pl)
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
