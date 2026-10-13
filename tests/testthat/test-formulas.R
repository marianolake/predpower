test_that("break-even values match the planning example", {
  # At n = 600, p = 100 and D* = 0.80, least squares needs rho2_in >= 0.203
  # while ridge needs 0.056, with planning p_eff about 25. The ridge value agrees with
  # rho2_min_ridge() in the replication code's R/design_formulas.R.
  expect_equal(round(breakeven_rho2(600, 100, 0.8), 3), 0.203)
  expect_equal(round(breakeven_rho2(600, 100, 0.8, "ridge"), 3), 0.056)
  expect_equal(round(p_eff(600, 100, 0.056)), 25)
})

test_that("R2oos is zero at the break-even signal", {
  for (est in c("ols", "ridge")) {
    for (d in c(0.4, 0.8, 1)) {
      be <- breakeven_rho2(500, 40, d, est)
      expect_equal(R2oos(500, 40, be, d, est), 0, tolerance = 1e-8,
                   info = paste(est, d))
    }
  }
})

test_that("D* = 1 is the optimum and lower D* always costs accuracy", {
  d <- seq(0.2, 1, by = 0.1)
  v <- R2oos(400, 30, 0.3, d)
  expect_true(all(diff(v) > 0))
  expect_equal(max(v), R2oos(400, 30, 0.3, 1))
})

test_that("accuracy rises with n and falls with p", {
  expect_true(all(diff(R2oos(seq(100, 1000, 50), 20, 0.2, 0.9)) > 0))
  expect_true(all(diff(R2oos(600, seq(5, 200, 5), 0.2, 0.9)) < 0))
})

test_that("roos exceeds R2oos, as correlation ignores miscalibration", {
  r2 <- R2oos(500, 20, 0.3, 0.9)
  rr <- roos(500, 20, 0.3, 0.9)
  expect_gt(rr^2, r2)
})

test_that("least squares is NA where it is undefined", {
  expect_true(is.na(R2oos(50, 60, 0.3, 0.9)))
  expect_true(is.na(breakeven_rho2(50, 60, 0.9)))
  expect_false(is.na(R2oos(50, 60, 0.3, 0.9, "ridge")))
})

test_that("ridge needs less signal than least squares", {
  expect_lt(breakeven_rho2(600, 100, 0.8, "ridge"),
            breakeven_rho2(600, 100, 0.8, "ols"))
})

test_that("invalid inputs are rejected", {
  expect_error(R2oos(100, 10, 0, 0.5), "rho2_in")
  expect_error(R2oos(100, 10, 0.3, 1.5), "Dstar")
  expect_error(R2oos(2, 10, 0.3, 0.5), "n must exceed")
})
