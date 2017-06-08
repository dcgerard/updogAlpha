

#' Plots the genotype lines for a given sequencing error rate, bias parameter, and ploidy.
#'
#' @param seq_error The sequencing error rate.
#' @param bias_val The bias parameter.
#' @param ploidy The ploidy.
#'
#' @author David Gerard
#'
#' @export
#'
plot_problines <- function(seq_error, bias_val, ploidy) {
  if (!requireNamespace(ggplot2, quietly = TRUE)) {
    stop("ggplot2 needs to be installed to run plot_problines.")
  }
  porig <- as.factor(paste0(0:ploidy, "/", ploidy))
  pvec <- updog::get_pvec(ploidy = ploidy, bias_val = bias_val,
                          seq_error = seq_error)
  slopevec <- pvec / (1 - pvec)
  xend <- pmin(rep(1, ploidy + 1), 1 / slopevec)
  yend <- pmin(rep(1, ploidy + 1), 1 * slopevec)
  df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1),
                         xend = xend, yend = yend, porig = porig)
  df_lines$bias <- bias_val
  df_lines$seq_error <- seq_error

  ## now plot ----
  pl <- ggplot2::ggplot(data = df_lines, mapping = ggplot2::aes_string(x = "x", y = "y",
                                                                       xend = "xend", yend = "yend",
                                                                       color = "porig")) +
    ggplot2::geom_segment() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::facet_grid(seq_error ~ bias) +
    ggthemes::scale_color_colorblind(name = "Original\nProbabilities")
  return(pl)
}
