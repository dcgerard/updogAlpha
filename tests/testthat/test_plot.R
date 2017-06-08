context("test plot")

test_that("plot works", {
  ocounts <- c(10, 20, 41, 55)
  osize   <- c(11, 26, 50, 72)

  p1counts <- c(11, 12)
  p1size   <- c(1000, 40)

  p2counts <- c(7, 900)
  p2size   <- c(12, 1000)

  ploidy <- 6
  seq_error <- 0.01

  maxpostprob <- rbeta(n = 4, shape1 = 1, shape2 = 1)
  prob_ok <- rbeta(n = 4, shape1 = 1, shape2 = 1)
  ogeno <- rep(3, length = length(ocounts))

  plot_geno(ocounts = ocounts, osize = osize, ploidy = 6, p1counts = p1counts, p1size = p1size,
            p2counts = p2counts, p2size = p2size, ogeno = ogeno,
            maxpostprob = maxpostprob, prob_ok = prob_ok)

}
)

test_that("plot.updog works", {
  umout <- readRDS(file = "updog_class.RDS")
  umout$maxpostprob <- runif(length(umout$ogeno))
  umout$ogeno[1] <- NA
  plot(umout, plot_beta = FALSE, show_maxpostprob = FALSE)

  plot(umout, gg = TRUE, plot_beta = TRUE, ask = FALSE)
  plot(umout, gg = FALSE, plot_beta = TRUE, ask = FALSE)
  plot(umout, gg = TRUE, plot_beta = FALSE, show_ogeno = FALSE)
  plot(umout, gg = TRUE, plot_beta = FALSE, show_outlier = FALSE)
  plot(umout, gg = FALSE, plot_beta = FALSE, show_ogeno = FALSE)
  plot(umout, gg = FALSE, plot_beta = FALSE, show_outlier = FALSE)
  plot(umout, gg = TRUE, show_maxpostprob = TRUE, plot_beta = FALSE)
  plot(umout, gg = FALSE, show_maxpostprob = TRUE, plot_beta = FALSE)

  umout$prob_ok <- NULL
  plot(umout, gg = TRUE, plot_beta = FALSE)
  plot(umout, gg = FALSE, plot_beta = FALSE)


  umout$ogeno <- NULL
  plot(umout, gg = TRUE, plot_beta = FALSE)
  plot(umout, gg = FALSE, plot_beta = FALSE)

  umout$opostprob <- NULL
  plot(umout, gg = TRUE, plot_beta = FALSE)
  plot(umout, gg = FALSE, plot_beta = FALSE)
})
