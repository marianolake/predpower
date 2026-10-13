test_that("the interval reproduces the application values", {
  # Paper, Section 5 (Table 1): the raw designs of p = 25 and 100 features.
  for (pp in c(25, 100)) {
    d <- subset(phonestudy$designs, p == pp & scale == "raw")
    s <- rho2_from_R2(d$R2_in, d$n, d$p)
    expect_equal(s$estimate, d$rho2_in, tolerance = 1e-8)
    expect_equal(s$lower, d$rho2_lower, tolerance = 1e-6)
    expect_equal(s$upper, d$rho2_upper, tolerance = 1e-6)
  }
  d <- subset(phonestudy$designs, p == 100 & scale == "raw")
  s <- rho2_from_R2(d$R2_in, d$n, d$p)
  expect_equal(round(c(s$estimate, s$lower, s$upper), 3), c(0.134, 0.073, 0.199))
})

test_that("a small raw R2 gives a negative estimate and a non-degenerate interval", {
  s <- rho2_from_R2(0.03, 620, 25)
  expect_lt(s$estimate, 0)
  expect_equal(s$lower, 0)
  expect_gt(s$upper, 0)
})

test_that("an F below the central lower quantile gives an NA upper limit, not [0, 0]", {
  for (R2 in c(0, -0.01, 0.005)) {
    expect_warning(s <- rho2_from_R2(R2, 620, 25), "upper limit is NA")
    expect_equal(s$lower, 0)
    expect_true(is.na(s$upper))
  }
})

test_that("rho2_from_data uses the numerical rank for the degrees of freedom", {
  set.seed(9)
  X <- matrix(stats::rnorm(100 * 5), 100, 5)
  y <- X[, 1] + stats::rnorm(100)
  Xd <- cbind(X, X[, 1] + X[, 2])                 # one redundant column
  expect_warning(s <- rho2_from_data(Xd, y), "rank")
  expect_equal(s$p, 5)
  expect_equal(s$estimate, rho2_from_data(X, y)$estimate, tolerance = 1e-10)
})

test_that("the interval widens with p at fixed signal and n", {
  w <- vapply(c(10, 50, 200, 400), function(p) {
    R2 <- 1 - (1 - 0.05) * (620 - p - 1) / (620 - 1)
    s <- rho2_from_R2(R2, 620, p); s$upper - s$lower
  }, numeric(1))
  expect_true(all(diff(w) > 0))
})

test_that("the interval covers the adjusted estimate", {
  s <- rho2_from_R2(0.30, 200, 10)
  expect_lte(s$lower, s$estimate)
  expect_gte(s$upper, s$estimate)
})

test_that("no residual degrees of freedom is an error", {
  expect_error(rho2_from_R2(0.5, 30, 29), "residual degrees")
})
