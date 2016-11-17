context("test plot")

test_that("plot works", {
  ocounts <- c(10, 20, 41, 55)
  osize   <- c(11, 26, 50, 72)

  p1counts <- c(11, 12)
  p1size   <- c(20, 40)

  p2counts <- c(7, 13)
  p2size   <- c(12, 13)

  ploidy <- 6
  seq_error <- 0.01

  plot_geno(ocounts = ocounts, osize = osize, ploidy = 6, p1counts = p1counts, p1size = p1size,
            p2counts = p2counts, p2size = p2size)

}
)
