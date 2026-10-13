test_that("rho2_reml recovers a known rho2_in when p > n", {
  set.seed(11)
  n <- 150; p <- 300; r2 <- 0.3
  est <- replicate(30, {
    X <- matrix(stats::rnorm(n * p), n, p)
    s <- as.vector(X %*% stats::rnorm(p)); s <- s / stats::sd(s) * sqrt(r2)
    y <- s + stats::rnorm(n, sd = sqrt(1 - r2))
    rho2_reml(X, y)$estimate
  })
  expect_lt(abs(mean(est) - r2), 0.04)
})

test_that("rho2_reml agrees with the adjusted R2 when p + 1 < n and the prior holds", {
  set.seed(12)
  n <- 400; p <- 40; r2 <- 0.2
  d <- replicate(30, {
    X <- matrix(stats::rnorm(n * p), n, p)
    s <- as.vector(X %*% stats::rnorm(p)); s <- s / stats::sd(s) * sqrt(r2)
    y <- s + stats::rnorm(n, sd = sqrt(1 - r2))
    rho2_reml(X, y)$estimate - rho2_from_data(X, y)$estimate
  })
  expect_lt(abs(mean(d)), 0.01)
})

test_that("the REML interval contains the estimate and lies in [0, 1)", {
  set.seed(13)
  X <- matrix(stats::rnorm(100 * 200), 100, 200)
  y <- as.vector(X[, 1:20] %*% stats::rnorm(20)) + stats::rnorm(100, sd = 3)
  s <- rho2_reml(X, y)
  expect_identical(s$method, "REML")
  expect_true(s$lower <= s$estimate && s$estimate <= s$upper)
  expect_true(s$lower >= 0 && s$upper < 1)
})

test_that("rho2_from_data and pred_power fall back to REML when p + 1 >= n", {
  set.seed(14)
  X <- matrix(stats::rnorm(80 * 120), 80, 120)
  y <- as.vector(X %*% stats::rnorm(120)) + stats::rnorm(80, sd = 8)
  expect_message(s <- rho2_from_data(X, y), "REML")
  expect_identical(s$method, "REML")
  expect_message(a <- pred_power(X = X, y = y), "only ridge")
  expect_identical(a$table$estimator, "ridge")
  expect_equal(a$rho2_in, s$estimate)
})

test_that("breakeven_rho2(X = ) gives the threshold pred_power reports", {
  set.seed(15)
  X <- matrix(stats::rlnorm(200 * 30), 200, 30)
  a <- pred_power(X = X, rho2_in = 0.2)
  expect_equal(breakeven_rho2(X = X, estimator = "ridge"),
               a$table$requires[a$table$estimator == "ridge"], tolerance = 1e-10)
  expect_equal(breakeven_rho2(X = X, estimator = "ols"),
               a$table$requires[a$table$estimator == "ols"], tolerance = 1e-10)
})
