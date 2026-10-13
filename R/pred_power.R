#' Predictive power: expected out-of-sample accuracy of a design
#'
#' The main entry point. Reports what a design is expected to deliver --
#' \eqn{R^2_{oos}} and \eqn{r_{oos}} under least squares and under ridge -- and
#' what signal it would need to beat the sample mean at all.
#'
#' Three ways of supplying the design are accepted, matching three situations a
#' planner is actually in.
#'
#' \describe{
#'   \item{No prior data}{Give \code{n}, \code{p}, \code{rho2_in}, and the
#'     efficiency factor of each estimator requested: \code{Dstar} for least
#'     squares and \code{Dstar_ridge} for ridge. Neither has a default; 1 is a
#'     perfectly balanced design and is optimistic for anything measured on a
#'     skewed scale.}
#'   \item{Prior covariate data}{Give \code{X}. Then \code{n}, \code{p},
#'     \code{Dstar}, \code{Dstar_ridge} and \code{peff} are all computed from
#'     it; values supplied for them alongside \code{X} are ignored, with a
#'     warning. This is the right mode when a cohort, a pilot sample, or a
#'     planned design matrix already exists. Factor, character and logical
#'     columns are dummy-coded, and \code{p} counts the dummy columns (see
#'     \code{\link{design_efficiency}}).}
#'   \item{Prior signal estimate}{Give \code{R2_in} alongside \code{n} and
#'     \code{p} (or give \code{X} and \code{y}). The adjusted R-squared and its
#'     confidence interval are computed by inverting the noncentral \eqn{F}
#'     (\code{\link{rho2_from_R2}}), and the calculations are run at the point estimate
#'     and at both ends of the interval.}
#' }
#'
#' \strong{Least squares and ridge use different efficiency factors.}
#' \eqn{D^*} describes least-squares leverage; \eqn{D^*(\kappa^*)} describes
#' ridge leverage at the design-implied penalty \eqn{\kappa^*}, which is
#' flatter, so \eqn{D^*(\kappa^*)} is typically the higher of the two. Without
#' \code{X} both must be stated; using \code{Dstar} for both is conservative
#' for ridge.
#'
#' \strong{With \code{X}}, every ridge quantity is computed from the design
#' itself at the value of \eqn{\rho^2_{in}} it is evaluated at: \eqn{p_{eff}}
#' from the eigenvalues of the standardized covariates and \eqn{D^*(\kappa^*)}
#' from the ridge leverage. This holds for the point row, the break-even
#' \code{requires} (the root of the design's own ridge \eqn{R^2_{oos}}) and the
#' interval ends. \strong{Without \code{X}}, the ridge point row uses
#' \code{Dstar_ridge} and \code{peff} (planning \eqn{p_{eff}} by default),
#' while \code{requires} and the interval ends hold \code{Dstar_ridge} fixed
#' and use the planning \eqn{p_{eff}} at each \eqn{\rho^2_{in}}.
#'
#' \strong{When least squares is undefined} (\eqn{p + 1 \ge n}), only ridge
#' can be analyzed: \code{estimator = "ols"} is an error, and
#' \code{estimator = "both"} reports the ridge row alone, with a message.
#'
#' The results are reported as values: the expected statistics, the
#' break-even \eqn{\rho^2_{in}} the design requires, and, when a signal
#' interval is available, the expected \eqn{R^2_{oos}} at both ends of it. They
#' describe expected values; a realized \eqn{R^2_{oos}} will vary around them.
#'
#' @param n Sample size. Taken from \code{X} when that is supplied.
#' @param p Number of covariates, excluding the intercept. Taken from \code{X}.
#' @param rho2_in Assumed in-sample signal, in (0, 1). Not needed when
#'   \code{R2_in} or \code{y} is given.
#' @param Dstar Least-squares design efficiency factor \eqn{D^*}, in (0, 1].
#'   Computed from \code{X} when that is supplied; otherwise required whenever
#'   least-squares results are requested.
#' @param Dstar_ridge Design efficiency factor at the ridge penalty,
#'   \eqn{D^*(\kappa^*)}. Computed from \code{X} when that is supplied; without
#'   \code{X} it must be given whenever ridge results are requested. Ridge
#'   leverage is flatter than least-squares leverage, so \eqn{D^*(\kappa^*)} is
#'   typically higher than \code{Dstar} -- often close to 1 even when
#'   \code{Dstar} is low. Setting it equal to \code{Dstar} gives conservative
#'   ridge results.
#' @param peff Effective number of parameters at \eqn{\kappa^*}, used for the
#'   ridge point row when \code{X} is not supplied. Defaults to the planning
#'   (Marchenko-Pastur) form of \code{\link{p_eff}}, which is exact for
#'   uncorrelated covariates and conservative for correlated ones. Computed
#'   from \code{X} when that is supplied.
#' @param X Optional covariate matrix or data frame, without an intercept.
#' @param y Optional outcome, used only to estimate \code{rho2_in}.
#' @param R2_in Optional raw in-sample R-squared from prior work.
#' @param level Confidence level for the signal interval.
#' @param estimator \code{"both"}, \code{"ols"} or \code{"ridge"}.
#' @return An object of class \code{"pred_power"}. Use \code{summary()} for
#'   a sentence per estimator, \code{as.data.frame()} for the numbers, and
#'   \code{plot()} to see how the results move with sample size (or, when
#'   \code{X} was supplied, with \eqn{\rho^2_{in}}).
#' @seealso \code{\link{design_efficiency}}, \code{\link{breakeven_rho2}},
#'   \code{\link{solve_n}}, \code{\link{sensitivity}}
#' @examples
#' # 1. no prior data: state D* and D*(kappa*)
#' pred_power(n = 600, p = 100, rho2_in = 0.10, Dstar = 0.5, Dstar_ridge = 0.95)
#'
#' # 2. prior covariate data: everything but rho2_in comes from X
#' set.seed(1)
#' X <- matrix(rlnorm(300 * 20), 300, 20)
#' pred_power(X = X, rho2_in = 0.10)
#'
#' # 3. prior signal estimate, with its interval: the p = 100 PhoneStudy design
#' d <- subset(phonestudy$designs, p == 100 & scale == "raw")
#' pred_power(n = d$n, p = d$p, R2_in = d$R2_in,
#'              Dstar = d$Dstar, Dstar_ridge = d$Dstar_ridge)
#' @export
pred_power <- function(n = NULL, p = NULL, rho2_in = NULL, Dstar = NULL,
                         Dstar_ridge = NULL, peff = NULL,
                         X = NULL, y = NULL, R2_in = NULL, level = 0.95,
                         estimator = c("both", "ols", "ridge")) {
  estimator <- match.arg(estimator)
  sig <- NULL
  basis <- NULL

  ## ---- resolve the design ------------------------------------------------
  if (!is.null(X)) {
    X <- .as_design(X)
    if (!is.null(n) && n != nrow(X))
      warning("n was supplied but X has ", nrow(X), " rows; using X.",
              call. = FALSE)
    if (!is.null(p) && p != ncol(X))
      warning("p was supplied but X has ", ncol(X), " columns (after dummy ",
              "coding and dropping constant columns); using X.", call. = FALSE)
    given <- c(Dstar = !is.null(Dstar), Dstar_ridge = !is.null(Dstar_ridge),
               peff = !is.null(peff))
    if (any(given))
      warning("X was supplied, so ", paste(names(given)[given], collapse = ", "),
              " ", if (sum(given) > 1) "are" else "is",
              " computed from it; the value", if (sum(given) > 1) "s",
              " supplied ", if (sum(given) > 1) "are" else "is", " ignored.",
              call. = FALSE)
    Dstar <- Dstar_ridge <- peff <- NULL
    n <- nrow(X); p <- ncol(X)
    if (!is.null(y)) {
      sig <- rho2_from_data(X, y, level)
      rho2_in <- sig$estimate
    }
  }
  if (!is.null(R2_in)) {
    if (is.null(n) || is.null(p))
      stop("R2_in needs n and p.", call. = FALSE)
    sig <- rho2_from_R2(R2_in, n, p, level)
    rho2_in <- sig$estimate
  }
  if (is.null(n) || is.null(p))
    stop("Supply n and p, or a covariate matrix X.", call. = FALSE)
  if (is.null(rho2_in))
    stop("Supply rho2_in, or supply y or R2_in so it can be estimated.",
         call. = FALSE)
  if (!is.finite(rho2_in) || rho2_in <= 0)
    stop("The adjusted in-sample R-squared is ", signif(rho2_in, 3),
         ". A design with no estimable signal cannot be analyzed: the ridge ",
         "penalty diverges and every quantity degenerates to the ",
         "intercept-only fit.", call. = FALSE)

  ## ---- least squares needs p + 1 < n -------------------------------------
  if (p + 1 >= n) {
    if (estimator == "ols")
      stop("p + 1 = ", p + 1, " is not less than n = ", n, ": least squares ",
           "is undefined. Use estimator = \"ridge\".", call. = FALSE)
    if (estimator == "both") {
      message("p + 1 = ", p + 1, " is not less than n = ", n, ": least ",
              "squares is undefined, so only ridge results are reported.")
      estimator <- "ridge"
    }
  }
  ests <- if (estimator == "both") c("ols", "ridge") else estimator

  ## ---- efficiency factors ------------------------------------------------
  eff <- NULL
  if ("ols" %in% ests) {
    if (!is.null(X)) {
      eff <- design_efficiency(X, "ols")
      Dstar <- eff$Dstar
    } else if (is.null(Dstar)) .stop_no_dstar()
  } else if (!is.null(X)) Dstar <- NA_real_

  if ("ridge" %in% ests) {
    if (!is.null(X)) {
      basis <- .ridge_basis(X)
      rq <- .ridge_X_at(basis, rho2_in)
      Dstar_ridge <- rq$Dstar; peff <- rq$p_eff
    } else {
      if (is.null(Dstar_ridge))
        stop("Supply Dstar_ridge, the design efficiency at the ridge penalty, ",
             "D*(kappa*), or supply X so it can be computed. It is typically ",
             "higher than Dstar, because ridge leverage is flatter; setting it ",
             "equal to Dstar gives conservative ridge results. Use ",
             "estimator = \"ols\" to skip ridge.", call. = FALSE)
      if (is.null(peff)) peff <- p_eff(n, p, rho2_in)
    }
  }

  ## ---- the four statistics ----------------------------------------------
  rows <- list()
  if ("ols" %in% ests)
    rows$ols <- .power_row("ols", n, p, rho2_in, Dstar, NULL, sig)
  if ("ridge" %in% ests)
    rows$ridge <- .power_row("ridge", n, p, rho2_in, Dstar_ridge, peff, sig, basis)

  structure(list(table = do.call(rbind, rows), n = n, p = p,
                 rho2_in = rho2_in,
                 Dstar = if (is.null(Dstar)) NA_real_ else Dstar,
                 Dstar_ridge = if (is.null(Dstar_ridge)) NA_real_ else Dstar_ridge,
                 p_eff = if (is.null(peff)) NA_real_ else peff,
                 signal = sig, efficiency = eff,
                 estimator = estimator, has_X = !is.null(X),
                 ridge_basis = basis),
            class = "pred_power")
}

.power_row <- function(est, n, p, rho2, Dstar, peff, sig, basis = NULL) {
  ## With a design in hand (basis), every ridge quantity -- point, break-even
  ## and interval ends -- comes from the design's own ridge leverage at the
  ## rho2_in in question. Without one, the planning forms are used.
  ridge_X <- est == "ridge" && !is.null(basis)
  be <- if (ridge_X) .breakeven_ridge_X(basis = basis) else
    breakeven_rho2(n, p, Dstar, est)
  out <- data.frame(
    estimator = est, n = n, p = p,
    p_eff = if (is.null(peff)) NA_real_ else peff,
    Dstar = Dstar, rho2_in = rho2,
    R2oos = R2oos(n, p, rho2, Dstar, est, peff),
    roos  = roos(n, p, rho2, Dstar, est, peff),
    requires = be,
    stringsAsFactors = FALSE)
  if (!is.null(sig)) {
    lo <- max(sig$lower, 1e-8)
    hi <- if (is.na(sig$upper)) NA_real_ else min(sig$upper, 1 - 1e-8)
    at <- function(r) {
      if (is.na(r)) return(NA_real_)
      if (ridge_X) return(.R2oos_ridge_X(rho2_in = r, basis = basis))
      R2oos(n, p, r, Dstar, est, if (est == "ridge") p_eff(n, p, r) else NULL)
    }
    out$R2oos_lower <- at(lo)
    out$R2oos_upper <- at(hi)
  }
  rownames(out) <- NULL
  out
}


#' @export
as.data.frame.pred_power <- function(x, ...) x$table

#' @export
print.pred_power <- function(x, ...) {
  f3 <- function(v) if (is.na(v)) "NA" else sprintf("%.3f", v)
  cat("Predictive power\n\n")
  cat(sprintf("  n = %d, p = %d", as.integer(x$n), as.integer(x$p)))
  if (x$has_X) cat("   (from supplied X)")
  cat("\n")
  cat(sprintf("  rho2_in = %.4f", x$rho2_in))
  if (!is.null(x$signal))
    cat(sprintf("   [%s, %s], %.0f%% interval",
                f3(x$signal$lower), f3(x$signal$upper), 100 * x$signal$level))
  cat("\n")
  ests <- x$table$estimator
  if ("ols" %in% ests)
    cat(sprintf("  least squares:  D* = %s\n", f3(x$Dstar)))
  if ("ridge" %in% ests)
    cat(sprintf("  ridge:          D*(kappa*) = %s,  p_eff = %.1f%s\n",
                f3(x$Dstar_ridge), x$p_eff,
                if (x$has_X) "   (at rho2_in, from X)" else ""))
  cat("\n")

  tb <- x$table
  num <- function(v, w) ifelse(is.na(v), formatC("--", width = w),
                               formatC(v, format = "f", digits = 3, width = w))
  fmt <- data.frame(
    estimator = tb$estimator,
    R2oos = num(tb$R2oos, 8),
    roos  = num(tb$roos, 8),
    requires = num(tb$requires, 8),
    stringsAsFactors = FALSE)
  if (!is.null(tb$R2oos_lower)) {
    fmt$R2oos_lower <- num(tb$R2oos_lower, 8)
    fmt$R2oos_upper <- num(tb$R2oos_upper, 8)
  }
  print.data.frame(fmt, row.names = FALSE)
  cat("\n  'requires' is the break-even rho2_in: the signal this design needs\n")
  cat("  to predict better than the sample mean, in expectation.\n")
  invisible(x)
}

#' @export
summary.pred_power <- function(object, ...) {
  print(object)
  tb <- object$table
  f3 <- function(v) if (is.na(v)) "NA" else sprintf("%.3f", v)
  cat("\n")
  for (i in seq_len(nrow(tb))) {
    r <- tb[i, ]
    cat(sprintf("  %-7s", paste0(toupper(r$estimator), ".")))
    if (is.na(r$R2oos)) {
      why <- if (r$p + 1 >= r$n) "p + 1 is not less than n." else
        "D* is undefined (some observation has leverage 1)."
      cat("Undefined:", why, "\n"); next
    }
    if (r$R2oos < 0) {
      cat(sprintf("Expected R2oos is %.3f, i.e. worse than predicting\n", r$R2oos))
      cat(sprintf("         the sample mean. Breaking even needs rho2_in >= %s;\n",
                  f3(r$requires)))
      cat(sprintf("         the value used is %.3f.\n", r$rho2_in))
    } else {
      cat(sprintf("Expected R2oos is %.3f, with expected correlation %s.\n",
                  r$R2oos, f3(r$roos)))
      cat(sprintf("         Breaking even needs rho2_in >= %s; the value used is %.3f.\n",
                  f3(r$requires), r$rho2_in))
    }
    if (!is.null(object$signal))
      cat(sprintf("         The %.0f%% interval for rho2_in is [%s, %s].\n",
                  100 * object$signal$level, f3(object$signal$lower),
                  f3(object$signal$upper)))
  }
  if (object$has_X && is.finite(object$Dstar) && object$Dstar < 1)
    cat(sprintf("\n  D* = %.3f. compare_scalings() shows how D* changes under a\n", object$Dstar),
        "  change of measurement scale.\n", sep = "")
  invisible(object)
}

## Ridge quantities for a realized design at given values of rho2_in:
## D*(kappa*) and p_eff, both from the design's own ridge leverage at
## kappa*(rho2_in). `basis` comes from .ridge_basis().
.ridge_X_at <- function(basis, rho2_in) {
  n <- basis$n; p <- basis$p
  res <- vapply(rho2_in, function(r) {
    g <- .leverage_ridge(basis, kappa_star(p, r))
    pe <- sum(g) - 1
    ds <- if (n - pe - 1 > 0) ((pe + 1) / (n - pe - 1)) / (sum(g / (1 - g)) / n)
          else NA_real_
    c(Dstar = ds, p_eff = pe)
  }, numeric(2))
  list(Dstar = res["Dstar", ], p_eff = res["p_eff", ])
}

## Ridge R2oos for a realized design at given values of rho2_in, with p_eff and
## D*(kappa*) evaluated at kappa*(rho2_in) from the design's own ridge leverage.
.R2oos_ridge_X <- function(X, rho2_in, basis = .ridge_basis(.as_design(X))) {
  rq <- .ridge_X_at(basis, rho2_in)
  R2oos(basis$n, basis$p, rho2_in, rq$Dstar, "ridge", rq$p_eff)
}

## Ridge break-even for a realized design: the root of .R2oos_ridge_X = 0 on
## (1/n, 1). The function is evaluated on a 200-point grid -- 100 points
## log-spaced from 1/n to 0.5 and 100 linearly spaced from 0.5 to 1 - 1e-12 --
## and the first interval where it turns from negative to non-negative is
## refined by uniroot() to tolerance 1e-10. NA, with a warning, if there is no
## sign change on the grid.
.breakeven_ridge_X <- function(X, basis = .ridge_basis(.as_design(X))) {
  f <- function(x) .R2oos_ridge_X(rho2_in = x, basis = basis)
  lo <- 1 / basis$n + 1e-12; hi <- 1 - 1e-12
  grid <- unique(c(exp(seq(log(lo), log(0.5), length.out = 100)),
                   seq(0.5, hi, length.out = 101)[-1]))
  v <- f(grid)
  i <- which(v[-length(v)] < 0 & v[-1] >= 0)[1]
  if (is.na(i)) {
    warning("The design's ridge R2oos does not cross zero on (1/n, 1); ",
            "the ridge break-even is NA.", call. = FALSE)
    return(NA_real_)
  }
  if (v[i + 1] == 0) return(grid[i + 1])
  stats::uniroot(f, grid[c(i, i + 1)], tol = 1e-10)$root
}
