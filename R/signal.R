#' In-sample signal from a fitted R-squared, with a confidence interval
#'
#' The design formulas are exact given \eqn{\rho^2_{in}}; they do not estimate
#' it. When it has to come from data rather than from prior work, two points of
#' practice matter, and this function enforces both.
#'
#' It must be the \emph{adjusted} R-squared: the raw statistic is badly
#' upward-biased once \eqn{p} is an appreciable fraction of \eqn{n}. And it
#' must carry an interval. Inverting the noncentral \eqn{F} gives an exact
#' confidence interval for \eqn{\rho^2_{in}} from \eqn{n}, \eqn{p} and the raw
#' R-squared alone, at no extra cost and with no extra data.
#'
#' The interval is not a formality. Its width is governed by residual degrees
#' of freedom, so at large \eqn{p} a study can report a point estimate of a
#' third and still be unable to distinguish its signal from zero. Feeding the
#' bounds through \code{\link{pred_power}} shows what that uncertainty means
#' for the design.
#'
#' \strong{The point estimate} is the adjusted R-squared,
#' \eqn{1 - (1 - R^2)(n - 1)/(n - p - 1)}. It is not constrained to be
#' non-negative: when the raw R-squared is small relative to \eqn{p/(n-1)} it
#' is negative, and it is reported as it is rather than floored at zero.
#'
#' \strong{The interval} collects the values of \eqn{\rho^2_{in}} whose
#' noncentral \eqn{F} (noncentrality \eqn{\lambda = n\rho^2_{in}/(1 -
#' \rho^2_{in})}) places the observed \eqn{F} between its \eqn{\alpha/2} and
#' \eqn{1 - \alpha/2} quantiles. A limit that would fall below zero is set to
#' zero, so a lower limit of 0 means the data cannot distinguish the signal
#' from none. The same inversion is used whatever the observed R-squared (a raw
#' R-squared at or below zero is treated as \eqn{F = 0}). When the observed
#' \eqn{F} lies below the \eqn{\alpha/2} quantile of the central \eqn{F} --
#' smaller than even a signal of exactly zero would typically produce -- no
#' value of \eqn{\rho^2_{in} \ge 0} satisfies the upper condition. The upper
#' limit is then returned as \code{NA}, with a warning, rather than as 0,
#' which would wrongly suggest the signal is known to be exactly zero.
#'
#' @param R2 The RAW in-sample R-squared of the fitted model, which is what the
#'   \eqn{F} statistic is built from.
#' @param n Sample size.
#' @param p Number of covariates, excluding the intercept.
#' @param level Confidence level.
#' @return An object of class \code{"signal_estimate"} carrying the adjusted
#'   point estimate and the interval.
#' @references Steiger, J. H. and Fouladi, R. T. (1992). R2: A computer program
#'   for interval estimation, power calculation, and hypothesis testing for the
#'   squared multiple correlation.
#' @examples
#' # the p = 100 PhoneStudy design (see ?phonestudy)
#' d <- subset(phonestudy$designs, p == 100 & scale == "raw")
#' rho2_from_R2(R2 = d$R2_in, n = d$n, p = d$p)
#'
#' # a small raw R-squared: the adjusted estimate is negative, the interval
#' # starts at 0 and its upper limit is still informative
#' rho2_from_R2(R2 = 0.03, n = 620, p = 25)
#' @export
rho2_from_R2 <- function(R2, n, p, level = 0.95) {
  if (level <= 0 || level >= 1) stop("level must be in (0, 1).", call. = FALSE)
  df1 <- p; df2 <- n - p - 1
  if (df2 <= 0)
    stop("n - p - 1 = ", df2, "; the model has no residual degrees of ",
         "freedom, so no interval exists.", call. = FALSE)
  if (!is.finite(R2) || R2 >= 1)
    stop("R2 must be a finite value below 1.", call. = FALSE)
  adj <- 1 - (1 - R2) * (n - 1) / df2
  Fobs <- max(R2, 0) / df1 / ((1 - R2) / df2)
  alpha <- 1 - level
  ## Each limit is the noncentrality at which the observed F sits at the
  ## requested tail, solved on a log scale because lambda spans decades. A limit
  ## below zero is set to zero; NA marks a limit with no solution.
  solve_ncp <- function(target) {
    f <- function(lg) stats::pf(Fobs, df1, df2, ncp = exp(lg)) - target
    if (f(log(1e-8)) < 0) return(0)
    if (f(log(1e10)) > 0) return(NA_real_)
    exp(stats::uniroot(f, c(log(1e-8), log(1e10)), tol = 1e-10)$root)
  }
  l <- solve_ncp(1 - alpha / 2)
  u <- if (stats::pf(Fobs, df1, df2) < alpha / 2) {
    warning("The observed F (", signif(Fobs, 3), ") is below the ",
            100 * alpha / 2, "% point of the central F: even a signal of zero ",
            "would typically produce a larger value. No rho2_in >= 0 meets ",
            "the upper condition, so the upper limit is NA. The data give no ",
            "evidence of signal.", call. = FALSE)
    NA_real_
  } else solve_ncp(alpha / 2)
  ci <- c(l / (n + l), u / (n + u))
  structure(list(estimate = adj, raw = R2, lower = unname(ci[1]),
                 upper = unname(ci[2]), n = n, p = p, level = level,
                 method = "adjusted R2"),
            class = "signal_estimate")
}

#' @export
print.signal_estimate <- function(x, ...) {
  cat("In-sample signal\n\n")
  if (identical(x$method, "REML")) {
    cat(sprintf("  REML estimate       %.3f   <- use this as rho2_in\n", x$estimate))
    cat(sprintf("  %.0f%% interval        [%.3f, %.3f]   (profile likelihood)\n",
                100 * x$level, x$lower, x$upper))
    cat(sprintf("  n = %d, p = %d; estimated under the isotropic ridge prior\n",
                x$n, x$p))
    return(invisible(x))
  }
  cat(sprintf("  raw R-squared       %.3f\n", x$raw))
  cat(sprintf("  adjusted R-squared  %.3f   <- use this as rho2_in\n", x$estimate))
  cat(sprintf("  %.0f%% interval        [%.3f, %s]\n",
              100 * x$level, x$lower,
              if (is.na(x$upper)) "NA" else sprintf("%.3f", x$upper)))
  cat(sprintf("  residual df         %d\n", x$n - x$p - 1))
  if (x$estimate <= 0)
    cat("\n  The adjusted estimate is not positive; rho2_in cannot be analyzed\n",
        "  at it (the ridge penalty diverges).\n", sep = "")
  if (x$lower <= 0)
    cat("\n  The interval reaches zero: this design cannot distinguish its\n",
        "  signal from none. Treat the point estimate with caution.\n", sep = "")
  invisible(x)
}

#' In-sample signal from a covariate matrix and an outcome
#'
#' Convenience wrapper that fits the least-squares model and passes its raw
#' R-squared to \code{\link{rho2_from_R2}}.
#'
#' The degrees of freedom are those of the fit actually made: \code{p} is the
#' numerical rank of \code{cbind(1, X)} (from its QR decomposition) minus one.
#' When some columns are linear combinations of others, the rank is smaller
#' than the number of columns; a warning reports this and the rank is used.
#' Factor, character and logical columns of a data frame are dummy-coded as in
#' \code{\link{design_efficiency}}.
#'
#' When \eqn{p + 1 \ge n} the least-squares fit is exact and the adjusted
#' R-squared is undefined; \eqn{\rho^2_{in}} is then estimated by
#' \code{\link{rho2_reml}}, with a message.
#'
#' @param X Covariate matrix or data frame, without an intercept column.
#' @param y Outcome vector.
#' @param level Confidence level.
#' @return An object of class \code{"signal_estimate"}.
#' @export
rho2_from_data <- function(X, y, level = 0.95) {
  X <- .as_design(X)
  if (length(y) != nrow(X)) stop("length(y) must equal nrow(X).", call. = FALSE)
  if (anyNA(y)) stop("y contains missing values.", call. = FALSE)
  n <- nrow(X)
  if (ncol(X) + 1 >= n) {
    message("cbind(1, X) has ", ncol(X) + 1, " columns and n = ", n,
            ": the least-squares fit is exact, so rho2_in is estimated by ",
            "REML under the isotropic ridge prior (see ?rho2_reml).")
    return(rho2_reml(X, y, level))
  }
  Xi <- cbind(1, scale(X))
  rk <- qr(Xi)$rank
  if (rk < ncol(Xi))
    warning("cbind(1, X) has ", ncol(Xi), " columns but numerical rank ", rk,
            "; the degrees of freedom use the rank, p = ", rk - 1, ".",
            call. = FALSE)
  p <- rk - 1
  if (p + 1 >= n) {
    message("The rank of cbind(1, X), ", p + 1, ", is not less than n = ", n,
            ": the least-squares fit is exact, so rho2_in is estimated by ",
            "REML under the isotropic ridge prior (see ?rho2_reml).")
    return(rho2_reml(X, y, level))
  }
  fit <- stats::lm.fit(Xi, y)
  R2 <- 1 - sum(fit$residuals^2) / sum((y - mean(y))^2)
  rho2_from_R2(R2, n, p, level)
}

#' In-sample signal by restricted maximum likelihood, for any p
#'
#' Estimates \eqn{\rho^2_{in}} when the adjusted R-squared cannot -- in
#' particular when \eqn{p + 1 \ge n} -- under the same isotropic prior on the
#' coefficients from which the ridge penalty \eqn{\kappa^*} is derived.
#'
#' With the covariates standardized and held fixed, the prior
#' \eqn{\beta_z \sim (0, \nu^2 I_p)} and noise \eqn{\epsilon \sim (0,
#' \sigma^2_{|x} I_n)} give the outcome covariance
#' \eqn{\nu^2 ZZ' + \sigma^2_{|x} I}, and \eqn{\rho^2_{in} = p\nu^2 /
#' (p\nu^2 + \sigma^2_{|x})}. The two variance components are separated by how
#' the outcome's variance changes across the eigen-directions of \eqn{ZZ'}:
#' the signal adds variance to each direction in proportion to its
#' eigenvalue, the noise the same amount to every direction. They are
#' estimated by restricted maximum likelihood (Patterson and Thompson, 1971),
#' which removes the intercept by using the \eqn{n - 1} contrasts orthogonal
#' to the constant, with a Gaussian working likelihood. This is the estimator
#' used for the share of trait variance explained by genome-wide markers
#' (Yang et al., 2010), where \eqn{p} also far exceeds \eqn{n}.
#'
#' \strong{Computation.} With \eqn{ZZ' = \sum_j \lambda_j u_j u_j'}, the
#' centred outcome has independent components \eqn{u_j'y} with variance
#' \eqn{\sigma^2_{|x}(\theta\lambda_j + 1)}, \eqn{\theta = \nu^2/\sigma^2_{|x}},
#' on the directions with \eqn{\lambda_j > 0}, and variance \eqn{\sigma^2_{|x}}
#' on the remaining contrasts. \eqn{\sigma^2_{|x}} is profiled out; the
#' restricted likelihood is maximized over \eqn{\rho^2_{in}} on a grid of
#' step \code{step}, and the interval is the profile-likelihood set
#' \eqn{\{\rho^2 : 2(\ell_{max} - \ell(\rho^2)) \le \chi^2_{1, level}\}}.
#'
#' \strong{What it assumes.} The estimate is only as good as the prior. On
#' realistic designs it is close to unbiased when the signal is spread across
#' many directions or concentrated in the largest ones, can overstate
#' \eqn{\rho^2_{in}} somewhat when the signal sits in a few of the largest
#' directions with \eqn{p} near \eqn{n}, and fails -- returning values near
#' zero -- when the signal lies in the directions of smallest eigenvalue. It
#' estimates the signal share expected under the prior, not the share for the
#' particular coefficients that generated the data. When \eqn{p + 1 < n},
#' the adjusted R-squared (\code{\link{rho2_from_R2}}) makes no such
#' assumption and is preferable unless \eqn{p} is close to \eqn{n}.
#'
#' @param X Covariate matrix or data frame, without an intercept column.
#' @param y Outcome vector.
#' @param level Confidence level for the profile-likelihood interval.
#' @param step Grid step for \eqn{\rho^2_{in}}.
#' @return An object of class \code{"signal_estimate"} with \code{method =
#'   "REML"}: the estimate and interval, \code{n}, \code{p} and \code{level}.
#' @references
#' Patterson, H. D. and Thompson, R. (1971). Recovery of inter-block
#' information when block sizes are unequal. \emph{Biometrika}, 58(3),
#' 545--554.
#'
#' Yang, J., Benyamin, B., McEvoy, B. P., et al. (2010). Common SNPs explain a
#' large proportion of the heritability for human height. \emph{Nature
#' Genetics}, 42(7), 565--569.
#' @seealso \code{\link{rho2_from_data}}, \code{\link{rho2_from_R2}},
#'   \code{\link{pred_power}}
#' @examples
#' set.seed(1)
#' n <- 200; p <- 400                     # more predictors than observations
#' X <- matrix(rnorm(n * p), n, p)
#' s <- as.vector(X %*% rnorm(p)); s <- s / sd(s) * sqrt(0.3)
#' y <- s + rnorm(n, sd = sqrt(0.7))      # rho2_in = 0.3
#' rho2_reml(X, y)
#' @export
rho2_reml <- function(X, y, level = 0.95, step = 0.001) {
  if (level <= 0 || level >= 1) stop("level must be in (0, 1).", call. = FALSE)
  X <- .as_design(X)
  if (length(y) != nrow(X)) stop("length(y) must equal nrow(X).", call. = FALSE)
  if (anyNA(y)) stop("y contains missing values.", call. = FALSE)
  Z <- scale(X); n <- nrow(Z); p <- ncol(Z)
  e <- eigen(tcrossprod(Z), symmetric = TRUE)
  r <- sum(e$values > 1e-8 * e$values[1])
  U <- e$vectors[, seq_len(r), drop = FALSE]; lam <- e$values[seq_len(r)]
  yc <- y - mean(y)
  ys <- drop(crossprod(U, yc)); rss0 <- max(sum(yc^2) - sum(ys^2), 0)
  m <- n - 1
  ll <- function(x) {
    w <- x / (p * (1 - x)) * lam + 1
    s2 <- (sum(ys^2 / w) + rss0) / m
    -0.5 * (sum(log(w)) + m * log(s2))
  }
  grid <- c(1e-6, seq(step, 1 - step, by = step))
  L <- vapply(grid, ll, numeric(1))
  k <- which.max(L)
  inside <- grid[L >= L[k] - stats::qchisq(level, 1) / 2]
  structure(list(estimate = grid[k], raw = NA_real_, lower = min(inside),
                 upper = max(inside), n = n, p = p, level = level,
                 method = "REML"),
            class = "signal_estimate")
}
