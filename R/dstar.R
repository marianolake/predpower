#' Design efficiency of a covariate matrix
#'
#' Computes the design efficiency factor \eqn{D^*} from a covariate matrix,
#' together with the leverage diagnostics that explain it. No outcome is
#' involved, so this can be run on pilot data, on an existing cohort, or on a
#' planned design matrix before any outcome is collected.
#'
#' \eqn{D^*} is the ratio of the prediction penalty an evenly balanced design
#' would incur to the penalty this design actually incurs. It equals 1 when
#' leverage is spread evenly and falls toward 0 as leverage concentrates.
#' \eqn{D^* = 0.5} means the design delivers half the predictive efficiency its
#' sample size and dimension would otherwise buy.
#'
#' In practice the dominant driver is skew, not imbalance: a few extreme rows
#' on a long-tailed covariate pull leverage toward themselves. When \eqn{D^*}
#' is well below 1, \code{\link{compare_scalings}} shows how much of the loss a
#' change of measurement scale recovers.
#'
#' \strong{Categorical covariates.} A data frame may contain factor, character
#' or logical columns. These are expanded to treatment (dummy) indicators with
#' \code{model.matrix(~ ., data = X)}, dropping the intercept column, so a
#' factor with \eqn{k} levels contributes \eqn{k - 1} columns and \code{p}
#' counts those dummy columns, not the original variables. Numeric matrices are
#' used as they are.
#'
#' \strong{When leave-one-out is undefined.} If some observation has
#' least-squares leverage of 1 (to within \eqn{10^{-10}}) -- for example the
#' only member of a category -- its leave-one-out prediction does not exist,
#' and neither does \eqn{D^*}. \code{design_efficiency()} then returns
#' \code{Dstar = NA} with a warning naming how many observations are affected.
#' Ridge leverage at \eqn{\kappa^* > 0} is always below 1, so the ridge factor
#' remains defined.
#'
#' @param X Covariate matrix or data frame, one row per observation, WITHOUT an
#'   intercept column. Factor, character and logical columns of a data frame
#'   are dummy-coded (see Details). Constant columns are dropped with a
#'   warning. Missing values are an error.
#' @param estimator \code{"ols"} for hat-matrix leverage, or \code{"ridge"} for
#'   ridge leverage at the design-implied penalty, which needs
#'   \code{rho2_in}.
#' @param rho2_in Assumed in-sample signal; required when
#'   \code{estimator = "ridge"}.
#' @return An object of class \code{"design_efficiency"}: a list with
#'   \code{Dstar} (for \code{estimator = "ridge"} this is
#'   \eqn{D^*(\kappa^*)}), the leverage vector, \code{max_leverage}, \code{n},
#'   \code{p} (the number of columns after dummy coding), the effective
#'   parameter count \code{p_eff} for ridge, and per-column skewness.
#' @examples
#' set.seed(1)
#' X <- matrix(rlnorm(200 * 5), 200, 5)   # skewed covariates
#' design_efficiency(X)
#'
#' # a data frame with a 4-level factor: p counts its 3 dummy columns
#' df <- data.frame(age = rnorm(200), site = factor(sample(letters[1:4], 200, TRUE)))
#' design_efficiency(df)$p
#' @export
design_efficiency <- function(X, estimator = c("ols", "ridge"), rho2_in = NULL) {
  estimator <- match.arg(estimator)
  X <- .as_design(X)
  n <- nrow(X); p <- ncol(X)
  if (estimator == "ols" && p + 1 >= n)
    stop("p + 1 = ", p + 1, " is not less than n = ", n,
         ". Least squares leverage is undefined here; use ",
         'estimator = "ridge".', call. = FALSE)

  if (estimator == "ols") {
    h <- .leverage_ols(cbind(1, scale(X)))
    peff <- NULL
  } else {
    if (is.null(rho2_in))
      stop('estimator = "ridge" needs rho2_in, since the penalty ',
           "kappa* depends on it.", call. = FALSE)
    h <- .leverage_ridge(.ridge_basis(X), kappa_star(p, rho2_in))
    peff <- sum(h) - 1
  }
  q <- if (estimator == "ols") p else peff
  saturated <- estimator == "ols" && any(h >= 1 - 1e-10)
  if (saturated) {
    warning(sum(h >= 1 - 1e-10), " observation(s) have least-squares ",
            "leverage 1, so their leave-one-out predictions -- and D* -- are ",
            "undefined (this happens, for example, when a category has a ",
            "single member). D* is returned as NA.", call. = FALSE)
    Dstar <- NA_real_
  } else {
    lev_sum <- sum(h / (1 - h)) / n
    Dstar <- if (n - q - 1 > 0) ((q + 1) / (n - q - 1)) / lev_sum else NA_real_
  }

  skew <- apply(X, 2, function(z) {
    s <- stats::sd(z); if (!is.finite(s) || s == 0) return(0)
    mean((z - mean(z))^3) / s^3
  })
  structure(list(Dstar = Dstar, leverage = h, max_leverage = max(h),
                 n = n, p = p, p_eff = peff, estimator = estimator,
                 rho2_in = rho2_in, skewness = skew,
                 even_leverage = (q + 1) / n),
            class = "design_efficiency")
}

## Converts a covariate matrix or data frame to the numeric design the formulas
## use: factor, character and logical columns of a data frame become treatment
## dummies (intercept column dropped), constant columns are removed, and
## missing values are refused before any expansion so that model.matrix() never
## drops rows silently.
.as_design <- function(X) {
  if (is.data.frame(X)) {
    if (anyNA(X)) stop("X contains missing values; leverage is undefined. ",
                       "Drop or impute them first.", call. = FALSE)
    if (ncol(X) == 0) stop("X has no columns.", call. = FALSE)
    if (all(vapply(X, is.numeric, logical(1)))) {
      X <- data.matrix(X)
    } else {
      X <- stats::model.matrix(~ ., data = X)[, -1, drop = FALSE]
      attr(X, "assign") <- NULL; attr(X, "contrasts") <- NULL
    }
  }
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.", call. = FALSE)
  if (anyNA(X)) stop("X contains missing values; leverage is undefined. ",
                     "Drop or impute them first.", call. = FALSE)
  keep <- apply(X, 2, function(z) stats::sd(z) > 0)
  if (any(!keep)) {
    warning(sum(!keep), " constant column(s) dropped.", call. = FALSE)
    X <- X[, keep, drop = FALSE]
  }
  if (ncol(X) == 0) stop("X has no non-constant columns.", call. = FALSE)
  X
}

## Least-squares leverage of a design that already carries its intercept column.
## Not capped: a leverage of 1 makes D* undefined, and design_efficiency()
## reports that rather than replacing it with a large finite number.
.leverage_ols <- function(Xi) {
  XtXi <- tryCatch(solve(crossprod(Xi)), error = function(e) NULL)
  if (is.null(XtXi))
    stop("The cross-product matrix of X is singular: some columns are ",
         "linear combinations of others. Remove the redundant columns.",
         call. = FALSE)
  pmin(rowSums((Xi %*% XtXi) * Xi), 1)
}

## Ridge leverage with an unpenalised intercept, from the singular value
## decomposition of the centred and scaled covariates Z = U diag(d) V'. Because
## the columns of Z are orthogonal to the intercept,
##   g_i(kappa) = 1/n + sum_j U_ij^2 d_j^2 / (d_j^2 + kappa),
## which is the diagonal of [1, Z] ([1, Z]'[1, Z] + diag(0, kappa I))^{-1} [1, Z]'.
## The decomposition is computed once per design, so the leverage at many
## penalties (as the break-even search needs) costs one matrix-vector product
## each. Capped just below 1, as in the replication code.
.ridge_basis <- function(X) {
  Z <- scale(X)
  s <- svd(Z, nv = 0)
  list(U2 = s$u^2, d2 = s$d^2, n = nrow(Z), p = ncol(Z))
}

.leverage_ridge <- function(basis, kappa) {
  w <- basis$d2 / (basis$d2 + kappa)
  pmin(1 / basis$n + drop(basis$U2 %*% w), 1 - 1e-12)
}

#' @export
print.design_efficiency <- function(x, ...) {
  cat("Design efficiency\n\n")
  cat(sprintf("  n = %d, p = %d", x$n, x$p))
  if (!is.null(x$p_eff)) cat(sprintf(", p_eff = %.1f", x$p_eff))
  cat(sprintf("   (%s)\n", x$estimator))
  cat(sprintf("  D*            %.3f\n", x$Dstar))
  cat(sprintf("  max leverage  %.3f   (even-leverage value %.3f)\n",
              x$max_leverage, x$even_leverage))
  hi <- mean(abs(x$skewness) > 2)
  cat(sprintf("  |skew| > 2    %.0f%% of columns (median %.1f)\n",
              100 * hi, stats::median(abs(x$skewness))))
  if (is.finite(x$Dstar)) {
    cat(sprintf("\n  This design delivers %.0f%% of the predictive efficiency of an even\n", 100 * x$Dstar))
    cat("  design of the same size.")
    if (x$Dstar < 1) cat(" compare_scalings() shows whether a change of\n  measurement scale raises it.")
    cat("\n")
  } else cat("\n  D* could not be computed.\n")
  invisible(x)
}


#' Rescale a covariate matrix
#'
#' Applies a monotone transformation to every column. Because \eqn{D^*} depends
#' on the measurement scale but the underlying participants do not, rescaling
#' can raise design efficiency without collecting anything.
#'
#' @param X Covariate matrix or data frame.
#' @param method \code{"log"} shifts each column to be positive and applies
#'   \code{log1p}; \code{"rank"} maps each column to normal scores;
#'   \code{"none"} returns \code{X} unchanged.
#' @return A numeric matrix.
#' @export
rescale_design <- function(X, method = c("log", "rank", "none")) {
  method <- match.arg(method)
  X <- .as_design(X)
  switch(method,
    none = X,
    log  = apply(X, 2, function(v) log1p(v - min(v))),
    rank = apply(X, 2, function(v) stats::qnorm((rank(v) - 0.5) / length(v))))
}

#' Compare design efficiency across measurement scales
#'
#' Runs \code{\link{design_efficiency}} on the raw covariates and on log and
#' rank-normal rescalings of them, and reports what each scale would deliver.
#' The predictor set and the sample are held fixed, so any movement in
#' \eqn{D^*} is attributable to the transformation alone.
#'
#' @inheritParams design_efficiency
#' @param rho2_in Assumed in-sample signal. Optional; when given, the expected
#'   \eqn{R^2_{oos}} under each scaling is reported too.
#' @return A data frame with one row per scaling, of class
#'   \code{"scaling_comparison"}.
#' @examples
#' set.seed(1)
#' X <- matrix(rlnorm(200 * 5), 200, 5)
#' compare_scalings(X, rho2_in = 0.3)
#' @export
compare_scalings <- function(X, rho2_in = NULL, estimator = c("ols", "ridge")) {
  estimator <- match.arg(estimator)
  X <- .as_design(X)
  scales <- c("raw", "log", "rank-normal")
  rows <- lapply(seq_along(scales), function(i) {
    Z <- rescale_design(X, c("none", "log", "rank")[i])
    de <- design_efficiency(Z, estimator, rho2_in)
    data.frame(scale = scales[i], Dstar = de$Dstar,
               max_leverage = de$max_leverage,
               R2oos = if (is.null(rho2_in)) NA_real_ else
                 R2oos(de$n, de$p, rho2_in, de$Dstar, estimator, de$p_eff),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(rho2_in)) out$R2oos <- NULL
  structure(out, class = c("scaling_comparison", "data.frame"))
}

#' @export
print.scaling_comparison <- function(x, ...) {
  cat("Design efficiency by measurement scale\n\n")
  print.data.frame(format(as.data.frame(x), digits = 3), row.names = FALSE)
  best <- x$scale[which.max(x$Dstar)]
  if (best != "raw")
    cat(sprintf("\n  %s scaling raises D* from %.3f to %.3f.\n",
                best, x$Dstar[x$scale == "raw"], max(x$Dstar)))
  invisible(x)
}
