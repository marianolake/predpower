test_that("supplying X reproduces the manual mode exactly for least squares", {
  set.seed(7)
  X <- matrix(stats::rnorm(400 * 15), 400, 15)
  d <- design_efficiency(X)$Dstar
  a_data <- pred_power(X = X, rho2_in = 0.2)
  a_manual <- pred_power(n = 400, p = 15, rho2_in = 0.2, Dstar = d, estimator = "ols")
  ols <- function(a) a$table$R2oos[a$table$estimator == "ols"]
  expect_equal(ols(a_data), ols(a_manual), tolerance = 1e-10)
})

test_that("without X the ridge line uses the planning p_eff", {
  # With X the eigenvalues are the design's own, so p_eff is exact; without it
  # the Marchenko-Pastur planning form is used and the two differ slightly.
  # Passing the exact p_eff through recovers agreement.
  set.seed(7)
  X <- matrix(stats::rnorm(400 * 15), 400, 15)
  dr <- design_efficiency(X, "ridge", 0.2)$Dstar
  a_data <- pred_power(X = X, rho2_in = 0.2, estimator = "ridge")
  approx <- pred_power(n = 400, p = 15, rho2_in = 0.2, Dstar = dr, Dstar_ridge = dr,
                         estimator = "ridge")
  exact <- pred_power(n = 400, p = 15, rho2_in = 0.2, Dstar = dr, Dstar_ridge = dr,
                        peff = p_eff(400, 15, 0.2, X = X), estimator = "ridge")
  expect_equal(a_data$table$R2oos, exact$table$R2oos, tolerance = 1e-10)
  expect_equal(a_data$table$R2oos, approx$table$R2oos, tolerance = 1e-2)
})

test_that("supplying R2_in produces an interval and values, not labels", {
  d <- subset(phonestudy$designs, p == 100 & scale == "raw")
  a <- pred_power(n = d$n, p = d$p, R2_in = d$R2_in,
                    Dstar = d$Dstar, Dstar_ridge = d$Dstar_ridge)
  expect_false(is.null(a$signal))
  expect_equal(a$rho2_in, d$rho2_in, tolerance = 1e-8)
  expect_equal(a$table$requires[a$table$estimator == "ols"], d$requires_ols,
               tolerance = 1e-8)
  expect_true(all(c("R2oos_lower", "R2oos_upper", "requires") %in% names(a$table)))
  expect_false("verdict" %in% names(a$table))
})

test_that("without X, ridge needs Dstar_ridge to be stated", {
  expect_error(pred_power(n = 600, p = 100, rho2_in = 0.05, Dstar = 0.8), "Dstar_ridge")
  expect_silent(pred_power(n = 600, p = 100, rho2_in = 0.05, Dstar = 0.8, estimator = "ols"))
})

test_that("the ridge break-even is found wherever a root exists", {
  # p > n with an uneven design
  be <- breakeven_rho2(250, 600, 0.6, "ridge")
  expect_true(is.finite(be))
  expect_equal(R2oos(250, 600, be, 0.6, "ridge"), 0, tolerance = 1e-6)
})

test_that("Dstar is required when it cannot be computed", {
  expect_error(pred_power(n = 100, p = 10, rho2_in = 0.2), "Dstar")
})

test_that("a design with no estimable signal is refused rather than fudged", {
  expect_error(suppressWarnings(pred_power(n = 100, p = 90, R2_in = 0.05, Dstar = 0.9)),
               "cannot be analyzed")
})

test_that("solve_n returns the smallest n meeting the target", {
  n <- solve_n(50, 0.15, 0.8, target = 0.05)
  expect_gte(R2oos(n, 50, 0.15, 0.8), 0.05)
  expect_lt(R2oos(n - 1, 50, 0.15, 0.8), 0.05)
})

test_that("solve_p returns the largest admissible p", {
  p <- solve_p(600, 0.2, 0.8)
  expect_gte(R2oos(600, p, 0.2, 0.8), 0)
  expect_lt(R2oos(600, p + 1, 0.2, 0.8), 0)
})

test_that("unreachable targets return NA rather than a wrong number", {
  expect_true(is.na(solve_dstar(620, 100, 0.05)))
  expect_true(is.na(solve_n(50, 0.02, 0.5, target = 0.9, max_n = 1e4)))
})

test_that("solve_dstar inverts R2oos", {
  d <- solve_dstar(500, 30, 0.25, target = 0.05)
  expect_equal(R2oos(500, 30, 0.25, d), 0.05, tolerance = 1e-8)
})

test_that("sensitivity covers the requested grid", {
  s <- sensitivity(n = c(200, 600), p = 50, rho2_in = c(0.05, 0.2),
                   Dstar = 0.8, Dstar_ridge = 0.95)
  expect_equal(nrow(s), 2 * 1 * 2 * 2)   # n x p x rho2 x estimator
  expect_true(all(c("R2oos", "breakeven") %in% names(s)))
  expect_equal(unique(s$Dstar[s$estimator == "ols"]), 0.8)
  expect_equal(unique(s$Dstar[s$estimator == "ridge"]), 0.95)
  expect_error(sensitivity(n = 200, p = 50, rho2_in = 0.1, Dstar = 0.8), "Dstar_ridge")
  expect_silent(sensitivity(n = 200, p = 50, rho2_in = 0.1, Dstar_ridge = 0.9,
                            estimator = "ridge"))
})

test_that("with X the ridge break-even and interval use the design itself", {
  # The ridge threshold for a realized design is the root of the ridge R2oos
  # computed from that design's eigenvalues and leverage at each rho2_in.
  set.seed(11)
  X <- matrix(stats::rexp(300 * 40), 300, 40)
  a <- pred_power(X = X, R2_in = 0.15, estimator = "ridge")
  be <- a$table$requires
  expect_true(is.finite(be))
  eff <- design_efficiency(X, "ridge", be)
  expect_equal(R2oos(300, 40, be, eff$Dstar, "ridge", eff$p_eff), 0, tolerance = 1e-6)
  lo <- max(a$signal$lower, 1e-8)
  effl <- design_efficiency(X, "ridge", lo)
  expect_equal(a$table$R2oos_lower,
               R2oos(300, 40, lo, effl$Dstar, "ridge", effl$p_eff), tolerance = 1e-10)
})

test_that("Dstar has no default anywhere", {
  expect_error(R2oos(600, 100, 0.1), "no safe default")
  expect_error(roos(600, 100, 0.1), "no safe default")
  expect_error(breakeven_rho2(600, 100), "no safe default")
  expect_error(solve_n(100, 0.1), "no safe default")
  expect_error(solve_p(600, 0.1), "no safe default")
  expect_error(sensitivity(600, 100, 0.1, Dstar_ridge = 0.9), "no safe default")
  expect_error(plot_breakeven(p = 10), "no safe default")
})

test_that("solve_p reports Inf when ridge meets the target at the search cap", {
  expect_message(v <- solve_p(200, 0.9, 1, estimator = "ridge"), "every p up to")
  expect_identical(v, Inf)
})

test_that("solve_dstar returns 0 when any D* meets the target", {
  expect_message(v <- solve_dstar(1000, 2, 0.5, target = -1e6), "any D\\*")
  expect_identical(v, 0)
})

test_that("with X, the ridge row is computed from X even if Dstar_ridge is given", {
  set.seed(3)
  X <- matrix(stats::rexp(300 * 20), 300, 20)
  ref <- pred_power(X = X, R2_in = 0.2)
  expect_warning(a <- pred_power(X = X, R2_in = 0.2, Dstar = 0.5,
                                   Dstar_ridge = 0.5, peff = 3),
                 "ignored")
  expect_equal(a$table, ref$table)
  eff <- design_efficiency(X, "ridge", a$rho2_in)
  r <- a$table[a$table$estimator == "ridge", ]
  expect_equal(r$Dstar, eff$Dstar, tolerance = 1e-10)
  expect_equal(r$p_eff, eff$p_eff, tolerance = 1e-10)
  expect_equal(r$R2oos, R2oos(300, 20, a$rho2_in, eff$Dstar, "ridge", eff$p_eff),
               tolerance = 1e-10)
})

test_that("with X and p + 1 >= n, least squares is dropped, not faked", {
  set.seed(5)
  X <- matrix(stats::rnorm(60 * 80), 60, 80)
  expect_message(a <- pred_power(X = X, rho2_in = 0.3), "least squares is undefined")
  expect_identical(a$table$estimator, "ridge")
  expect_true(is.na(a$Dstar))
  out <- paste(utils::capture.output(print(a)), collapse = "\n")
  expect_false(grepl("D\\* = 1.000", out))
  expect_error(pred_power(X = X, rho2_in = 0.3, estimator = "ols"), "undefined")
  expect_silent(pred_power(X = X, rho2_in = 0.3, estimator = "ridge"))
})

test_that("the fine-grid ridge break-even matches a direct root on a test design", {
  set.seed(11)
  X <- matrix(stats::rexp(300 * 40), 300, 40)
  f <- function(x) {
    e <- design_efficiency(X, "ridge", x)
    R2oos(300, 40, x, e$Dstar, "ridge", e$p_eff)
  }
  direct <- stats::uniroot(f, c(0.01, 0.2), tol = 1e-12)$root
  expect_equal(predpower:::.breakeven_ridge_X(X), direct, tolerance = 1e-8)
})

test_that("plots draw without error, with and without X", {
  pdf(NULL); on.exit(grDevices::dev.off())
  d <- subset(phonestudy$designs, p == 100 & scale == "raw")
  expect_silent(plot(pred_power(n = d$n, p = d$p, R2_in = d$R2_in,
                                  Dstar = d$Dstar, Dstar_ridge = d$Dstar_ridge)))
  set.seed(2)
  expect_silent(plot(pred_power(X = matrix(stats::rlnorm(200 * 10), 200, 10),
                                  rho2_in = 0.2)))
})
