context("NA's in plots")


test_that("plot_geno accepts NA's", {
  ocounts <- c(10, 20, 41, NA)
  osize   <- c(11, 26, 50, 72)

  ogeno   <- c(NA, 1, 1, 2)

  p1counts <- c(11, NA)
  p1size   <- c(20, 40)

  p2counts <- c(NA, 13)
  p2size   <- c(12, 13)

  ploidy <- 6
  seq_error <- 0.01

  plot_geno(ocounts = ocounts, osize = osize, ploidy = ploidy, ogeno = ogeno,
            p1counts = p1counts, p1size = p1size, p2counts = p2counts, p2size = p2size)

}
)
