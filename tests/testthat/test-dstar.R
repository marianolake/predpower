test_that("a balanced orthogonal design has D* very close to 1", {
  set.seed(4)
  X <- matrix(stats::rnorm(2000 * 5), 2000, 5)
  expect_gt(design_efficiency(X)$Dstar, 0.98)
})

test_that("skew lowers D* and rescaling recovers it", {
  set.seed(4)
  X <- matrix(stats::rlnorm(300 * 8, sdlog = 1.6), 300, 8)
  cs <- compare_scalings(X)
  expect_lt(cs$Dstar[cs$scale == "raw"], 0.9)
  expect_gt(cs$Dstar[cs$scale == "rank-normal"],
            cs$Dstar[cs$scale == "raw"])
})

test_that("D* lies in (0, 1]", {
  set.seed(4)
  for (i in 1:5) {
    X <- matrix(stats::rlnorm(200 * 6), 200, 6)
    d <- design_efficiency(X)$Dstar
    expect_gt(d, 0); expect_lte(d, 1)
  }
})

test_that("ridge leverage is flatter, so ridge D* is the larger", {
  set.seed(4)
  X <- matrix(stats::rlnorm(200 * 10), 200, 10)
  expect_gt(design_efficiency(X, "ridge", rho2_in = 0.2)$Dstar,
            design_efficiency(X)$Dstar)
})

test_that("bad input is caught", {
  expect_error(design_efficiency(matrix(NA_real_, 10, 2)), "missing")
  expect_error(design_efficiency(matrix(stats::rnorm(20), 4, 5)),
               "not less than n")
})

test_that("factors are dummy-coded and p counts the dummy columns", {
  set.seed(4)
  df <- data.frame(x = stats::rnorm(120),
                   g = factor(rep(c("a", "b", "c", "d"), each = 30)))
  de <- design_efficiency(df)
  expect_equal(de$p, 4)
  expect_equal(de$Dstar,
               design_efficiency(stats::model.matrix(~ ., df)[, -1])$Dstar)
  df$g <- as.character(df$g)
  expect_equal(design_efficiency(df)$p, 4)
  df$g[1] <- NA
  expect_error(design_efficiency(df), "missing")
})

test_that("leverage of 1 makes least-squares D* NA, with a warning", {
  set.seed(4)
  df <- data.frame(x = stats::rnorm(50),
                   g = factor(c("solo", rep(c("a", "b"), length.out = 49))))
  expect_warning(de <- design_efficiency(df), "leverage 1")
  expect_true(is.na(de$Dstar))
  expect_true(is.finite(design_efficiency(df, "ridge", rho2_in = 0.2)$Dstar))
})
