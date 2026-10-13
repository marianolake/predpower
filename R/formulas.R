#' @keywords internal
"_PACKAGE"

## The four design quantities, in the notation of the paper:
##   n        sample size
##   p        number of covariates; the fitted model has p + 1 parameters
##   rho2_in  in-sample signal, the population R-squared of the fitted model
##   Dstar    design efficiency factor, in (0, 1]
##
## Everything in this file is a function of those four alone. No outcome data
## enters, which is what makes these usable at the planning stage.

#' Design-implied ridge penalty
#'
#' The penalty \eqn{\kappa^* = p(1-\rho^2_{in})/\rho^2_{in}}, the ratio of
#' residual variance to per-coefficient prior variance under an isotropic prior
#' on the coefficients. It requires no outcome data and no tuning, so a planner
#' can compute it in advance.
#'
#' @param p Number of covariates, excluding the intercept.
#' @param rho2_in Assumed in-sample signal, in (0, 1).
#' @return The penalty, a numeric vector.
#' @examples
#' kappa_star(p = 100, rho2_in = 0.2)
#' @export
kappa_star <- function(p, rho2_in) {
  stopifnot(all(p > 0, na.rm = TRUE))
  p * (1 - rho2_in) / rho2_in
}

#' Effective number of parameters under ridge
#'
#' \eqn{p_{eff}(\kappa^*) = \sum_j \lambda_j/(\lambda_j + \kappa^*)}, the ridge
#' analogue of \eqn{p}. Two forms are available. Given a covariate matrix, the
#' eigenvalues are those of the design and the result is exact. Without one --
#' the usual planning situation -- the eigenvalues are approximated by those of
#' uncorrelated standardized covariates, which spread out by sampling alone
#' according to the Marchenko-Pastur law. With \eqn{\gamma = p/(n-1)} and
#' \eqn{k = \kappa^*/(n-1)} this gives the closed form
#' \deqn{p_{eff} \approx p - \frac{n-1}{2}\left[\sqrt{(1-\gamma+k)^2 + 4\gamma k}
#'   - (1-\gamma+k)\right],}
#' a function of \code{n}, \code{p} and \code{rho2_in} alone, never above
#' \eqn{\min(p, n-1)}.
#'
#' The planning form is exact (to sampling error) for uncorrelated covariates
#' at any \eqn{p/n}. Correlated or heavy-tailed covariates spread the
#' eigenvalues further and lower \eqn{p_{eff}}, so for them it overstates
#' \eqn{p_{eff}} and the ridge results are conservative. Supply \code{X} when
#' you have it.
#'
#' @param n Sample size.
#' @param p Number of covariates.
#' @param rho2_in Assumed in-sample signal.
#' @param X Optional covariate matrix or data frame, without an intercept
#'   column. When supplied, \code{n} and \code{p} are taken from it (after the
#'   dummy coding and constant-column removal described in
#'   \code{\link{design_efficiency}}) and the eigenvalues are exact.
#' @return The effective number of parameters, a numeric vector.
#' @examples
#' p_eff(n = 600, p = 100, rho2_in = 0.05)
#' @export
p_eff <- function(n, p, rho2_in, X = NULL) {
  if (!is.null(X)) {
    X <- .as_design(X)
    n <- nrow(X); p <- ncol(X)
    lam <- svd(scale(X), nu = 0, nv = 0)$d^2
    k <- kappa_star(p, rho2_in)
    return(vapply(k, function(kk) sum(lam / (lam + kk)), numeric(1)))
  }
  ## Marchenko-Pastur planning form (paper, Section 3; Supplement A.16)
  g <- p / (n - 1)
  k <- kappa_star(p, rho2_in) / (n - 1)
  a <- 1 - g + k
  p - (n - 1) / 2 * (sqrt(a^2 + 4 * g * k) - a)
}

## D/T, the share of prediction variance that is estimation error rather than
## signal. `q` is p for OLS and p_eff for ridge.
.d_over_t <- function(n, q, rho2_in, Dstar) {
  out <- (1 - rho2_in) * (q + 1) / ((n - q - 1) * Dstar)
  out[n - q - 1 <= 0] <- Inf
  out
}

#' Expected out-of-sample R-squared
#'
#' The design formula for \eqn{R^2_{oos}} under leave-one-out cross-validation,
#' \deqn{R^2_{oos} \approx (n\psi - 1)/(n - 1), \qquad \psi = \rho^2_{in} - D/T,}
#' for ordinary least squares, and its ridge counterpart at \eqn{\kappa^*} with
#' \eqn{p_{eff}} in place of \eqn{p}.
#'
#' Arguments are recycled, so any one of them may be a vector; this is what the
#' plotting functions use to trace a curve.
#'
#' @param n Sample size.
#' @param p Number of covariates, excluding the intercept.
#' @param rho2_in Assumed in-sample signal, in (0, 1).
#' @param Dstar Design efficiency factor, in (0, 1]. Required: there is no
#'   safe default. For \code{estimator = "ols"} this is \eqn{D^*}, the
#'   least-squares factor; for \code{estimator = "ridge"} it is
#'   \eqn{D^*(\kappa^*)}, the factor at the ridge penalty. The two differ.
#'   Ridge leverage is flatter than least-squares leverage, so
#'   \eqn{D^*(\kappa^*)} is typically higher than \eqn{D^*} -- often close to 1
#'   when \eqn{D^*} is low -- and using \eqn{D^*} for both is conservative for
#'   ridge. \code{\link{design_efficiency}} computes either from a covariate
#'   matrix.
#' @param estimator Either \code{"ols"} or \code{"ridge"}.
#' @param peff Effective number of parameters, for \code{estimator = "ridge"}.
#'   Defaults to the planning (Marchenko-Pastur) form of \code{\link{p_eff}}.
#' @return Expected out-of-sample R-squared. Values below zero mean the design
#'   is expected to predict worse than the sample mean. \code{NA} where least
#'   squares is undefined (\eqn{p + 1 \ge n}).
#' @examples
#' # the p = 100 PhoneStudy design (see ?phonestudy)
#' R2oos(n = 620, p = 100, rho2_in = 0.134, Dstar = 0.346)
#' # for ridge, Dstar is D*(kappa*), the factor at the ridge penalty
#' R2oos(n = 620, p = 100, rho2_in = 0.134, Dstar = 0.970, estimator = "ridge")
#' @export
R2oos <- function(n, p, rho2_in, Dstar,
                  estimator = c("ols", "ridge"), peff = NULL) {
  if (missing(Dstar)) .stop_no_dstar()
  estimator <- match.arg(estimator)
  .check_inputs(n, p, rho2_in, Dstar)
  if (estimator == "ols") {
    out <- (n * (rho2_in - .d_over_t(n, p, rho2_in, Dstar)) - 1) / (n - 1)
    out[p + 1 >= n] <- NA_real_
    return(out)
  }
  q <- if (is.null(peff)) p_eff(n, p, rho2_in) else peff
  (n * (rho2_in - .d_over_t(n, q, rho2_in, Dstar)) - 1) / (n - 1)
}

#' Expected out-of-sample correlation
#'
#' The design formula for \eqn{r_{oos}}, the expected correlation between
#' observed and held-out predicted outcomes. Under ordinary least squares it is
#' built from a slope-like intermediate; under ridge it is a function of
#' \eqn{R^2_{oos}} alone.
#'
#' \eqn{r_{oos}^2} always exceeds \eqn{R^2_{oos}}, because correlation does not
#' charge for miscalibration while \eqn{R^2_{oos}} does. A design can therefore
#' show a respectable correlation while predicting worse than the sample mean.
#'
#' @inheritParams R2oos
#' @return Expected out-of-sample correlation. \code{NA} where the design is
#'   too weak for the quantity to be defined.
#' @examples
#' roos(n = 620, p = 100, rho2_in = 0.134, Dstar = 0.346)
#' roos(n = 620, p = 100, rho2_in = 0.134, Dstar = 0.970, estimator = "ridge")
#' @export
roos <- function(n, p, rho2_in, Dstar,
                 estimator = c("ols", "ridge"), peff = NULL) {
  if (missing(Dstar)) .stop_no_dstar()
  estimator <- match.arg(estimator)
  .check_inputs(n, p, rho2_in, Dstar)
  if (estimator == "ridge") {
    R2  <- R2oos(n, p, rho2_in, Dstar, "ridge", peff)
    nn  <- rep_len(n, length(R2))
    den <- nn * ((nn - 2) * R2 + 2)
    out <- rep(NA_real_, length(R2))
    ok  <- is.finite(den) & den > 0
    ## den <= 0 means R2 is far enough below zero that the correlation is not
    ## defined -- the design predicts so badly the quantity loses meaning.
    out[ok] <- ((nn - 1) * R2 + 1)[ok] / sqrt(den[ok])
    return(out)
  }
  ## S is the intermediate of the proof, not a reported quantity; it is formed
  ## here only because r_oos is assembled from it.
  len <- max(length(n), length(p), length(rho2_in), length(Dstar))
  n <- rep_len(n, len); p <- rep_len(p, len)
  rho2_in <- rep_len(rho2_in, len); Dstar <- rep_len(Dstar, len)
  q   <- rho2_in - 1 / n
  den <- q * (n - p - 1) * Dstar + (1 - rho2_in) * (p + 1)
  S   <- ifelse(is.finite(den) & den != 0, q * (n - p - 1) * Dstar / den, NA_real_)
  rad <- S * (n * rho2_in - 1) / (n - 1)
  out <- ifelse(!is.na(rad) & rad >= 0,
                sign(n * rho2_in - 1) * sqrt(abs(rad)), NA_real_)
  out[p + 1 >= n] <- NA_real_
  out
}

#' Break-even signal: the minimum in-sample R-squared a design needs
#'
#' The smallest \eqn{\rho^2_{in}} at which a design is expected to beat the
#' sample mean, that is at which \eqn{R^2_{oos}} reaches zero. Under least
#' squares this is closed form. Under ridge it is implicit, because
#' \eqn{\kappa^*} depends on \eqn{\rho^2_{in}} and so \eqn{p_{eff}} does too;
#' it is found by bracketed root-finding (\code{\link[stats]{uniroot}}) on
#' \eqn{(1/n, 1)}. Without a covariate matrix, \eqn{p_{eff}} at each candidate
#' \eqn{\rho^2_{in}} comes from the planning form of \code{\link{p_eff}},
#' while \eqn{D^*(\kappa^*)} is held at the value supplied.
#'
#' This is the single most useful planning quantity in the package. It answers
#' "how strong would the signal have to be for this design to be worth
#' running?" without requiring any guess about the signal at all.
#'
#' @param n Sample size.
#' @param p Number of covariates.
#' @param Dstar Design efficiency factor, in (0, 1]. Required. For
#'   \code{estimator = "ridge"} this is \eqn{D^*(\kappa^*)}, the factor at the
#'   ridge penalty, which is typically higher than the least-squares
#'   \eqn{D^*}; see \code{\link{R2oos}}.
#' @param estimator Either \code{"ols"} or \code{"ridge"}.
#' @param X Optional covariate matrix or data frame. With \code{X} and
#'   \code{estimator = "ridge"}, the threshold is found with the design's own
#'   ridge leverage at each candidate \eqn{\rho^2_{in}} (so \eqn{p_{eff}} and
#'   \eqn{D^*(\kappa^*)} come from \code{X}), and \code{n}, \code{p} and
#'   \code{Dstar} are taken from \code{X}; this is the threshold
#'   \code{\link{pred_power}} reports. With \code{X} and
#'   \code{estimator = "ols"}, \code{Dstar} is computed from \code{X}.
#' @return The break-even in-sample signal. \code{NA} where the design admits
#'   no solution in (0, 1) -- for least squares, when \eqn{p + 1 \ge n}.
#' @examples
#' # the p = 100 PhoneStudy design: D* = 0.346, D*(kappa*) = 0.970
#' breakeven_rho2(n = 620, p = 100, Dstar = 0.346)
#' breakeven_rho2(n = 620, p = 100, Dstar = 0.970, estimator = "ridge")
#' @export
breakeven_rho2 <- function(n, p, Dstar, estimator = c("ols", "ridge"), X = NULL) {
  estimator <- match.arg(estimator)
  if (!is.null(X)) {
    X <- .as_design(X)
    if (estimator == "ridge") return(.breakeven_ridge_X(basis = .ridge_basis(X)))
    n <- nrow(X); p <- ncol(X)
    if (p + 1 >= n) return(NA_real_)
    Dstar <- design_efficiency(X, "ols")$Dstar
  }
  if (missing(Dstar)) .stop_no_dstar()
  len <- max(length(n), length(p), length(Dstar))
  n <- rep_len(n, len); p <- rep_len(p, len); Dstar <- rep_len(Dstar, len)
  if (estimator == "ols") {
    out <- 1 / n + (1 - 1 / n) * (p + 1) / ((p + 1) + (n - p - 1) * Dstar)
    out[n - p - 1 <= 0] <- NA_real_
    return(out)
  }
  vapply(seq_len(len), function(i) .breakeven_ridge1(n[i], p[i], Dstar[i]),
         numeric(1))
}

.breakeven_ridge1 <- function(n, p, Dstar, tol = 1e-10) {
  if (!is.finite(n) || !is.finite(p) || !is.finite(Dstar) || n <= 2) return(NA_real_)
  ## rho2_in appears on both sides (p_eff depends on it through kappa*); the
  ## threshold is the root of x - rhs(x) on (1/n, 1), found by bracketing.
  rhs <- function(x) {
    pe <- p_eff(n, p, x)
    1 / n + (1 - 1 / n) * (pe + 1) / ((pe + 1) + (n - pe - 1) * Dstar)
  }
  f <- function(x) x - rhs(x)
  lo <- 1 / n + 1e-12; hi <- 1 - 1e-12
  flo <- f(lo); fhi <- f(hi)
  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) return(NA_real_)
  stats::uniroot(f, c(lo, hi), tol = tol)$root
}

## The design efficiency factor has no safe default; every function that needs
## it says so in the same words.
.stop_no_dstar <- function(arg = "Dstar") {
  stop("Supply ", arg, ". There is no safe default: Dstar = 1 assumes ",
       "perfectly even leverage, which is optimistic for most real designs. ",
       "design_efficiency() computes D* from a covariate matrix, and ",
       "D*(kappa*) with estimator = \"ridge\". Without a design matrix, state ",
       "both: D*(kappa*) is typically higher than D*, because ridge leverage ",
       "is flatter, and using D* for both is conservative for ridge.",
       call. = FALSE)
}

.check_inputs <- function(n, p, rho2_in, Dstar) {
  if (any(rho2_in <= 0 | rho2_in >= 1, na.rm = TRUE))
    stop("rho2_in must be strictly between 0 and 1.", call. = FALSE)
  if (any(Dstar <= 0 | Dstar > 1 + 1e-8, na.rm = TRUE))   # tolerance for rounding in a computed D*
    stop("Dstar must be in (0, 1]. It is a ratio of an even-leverage ",
         "penalty to the realized one, so it cannot exceed 1.", call. = FALSE)
  if (any(n <= 2, na.rm = TRUE)) stop("n must exceed 2.", call. = FALSE)
  if (any(p < 1, na.rm = TRUE)) stop("p must be at least 1.", call. = FALSE)
  invisible(TRUE)
}
