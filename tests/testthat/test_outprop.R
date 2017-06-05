context("test outlier")

test_that("get_parent_outprop works ok", {
  set.seed(346)
  pcounts <- 10
  psize <- 21
  ploidy <- 4
  pgeno <- 1
  d <- 0.9
  eps <- 0.01
  tau <- 0.01
  out_prop <- 0.2
  out_mean <- 0.5
  out_disp <- 0.3

  theta_cpp <- get_parent_outprop(pcounts = pcounts, psize = psize, ploidy = ploidy,
                                  pgeno = pgeno, d = d, eps = eps, tau = tau,
                                  out_prop = out_prop, out_mean = out_mean,
                                  out_disp = out_disp)
  xi <- pbias_double(prob = pgeno / ploidy, bias = d, seq_error = eps)
  pgood <- dbetabinom_mu_rho_cpp_double(x = pcounts, size = psize, mu = xi, rho = tau, return_log = FALSE)
  pbad  <- dbetabinom_mu_rho_cpp_double(x = pcounts, size = psize, mu = out_mean, rho = out_disp, return_log = FALSE)

  theta_manual <- out_prop * pbad / (out_prop * pbad + (1 - out_prop) * pgood)

  expect_equal(theta_manual, theta_cpp)
}
)
