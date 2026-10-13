## Solvers. Each inverts one of the design formulas for one input, holding the
## others fixed. All are monotone in the quantity being solved for, which is
## what makes a bracketed root reliable:
##   R2oos increases in n, in rho2_in and in Dstar, and decreases in p.

#' Solve for the sample size a design would need
#'
#' The prediction analogue of a power calculation: the smallest \eqn{n} at
#' which the design reaches a target \eqn{R^2_{oos}}. With the default target
#' of 0 this is the sample size at which the design begins to beat the sample
#' mean at all.
#'
#' \strong{Search.} The search starts at \eqn{n = p + 3}, the smallest sample
#' size at which the formulas are defined. If the target is already met there,
#' that is the answer. Otherwise \eqn{n} is doubled until the target is met,
#' and the last doubling interval is then bisected on the integers, which
#' returns the smallest qualifying \eqn{n} because expected \eqn{R^2_{oos}}
#' increases with \eqn{n}. If the target is not met at \code{max_n}, the
#' result is \code{NA}. \code{Dstar} (and, for ridge, the planning
#' \eqn{p_{eff}}) is held fixed as \eqn{n} changes.
#'
#' @param p Number of covariates.
#' @param rho2_in Assumed in-sample signal.
#' @param Dstar Design efficiency factor. Required. For
#'   \code{estimator = "ridge"} this is \eqn{D^*(\kappa^*)}, the factor at the
#'   ridge penalty, which is typically higher than the least-squares
#'   \eqn{D^*}; see \code{\link{R2oos}}.
#' @param target Target expected \eqn{R^2_{oos}}.
#' @param estimator \code{"ols"} or \code{"ridge"}.
#' @param max_n Largest sample size to consider.
#' @return The required sample size, rounded up. \code{NA} if the target is
#'   unreachable at any \eqn{n \le} \code{max_n}, which happens when the
#'   assumed signal is simply too weak.
#' @examples
#' solve_n(p = 100, rho2_in = 0.134, Dstar = 0.346)
#' solve_n(p = 100, rho2_in = 0.134, Dstar = 0.970, estimator = "ridge")
#' solve_n(p = 20, rho2_in = 0.20, Dstar = 0.9, target = 0.10)
#' @export
solve_n <- function(p, rho2_in, Dstar, target = 0,
                    estimator = c("ols", "ridge"), max_n = 1e6) {
  if (missing(Dstar)) .stop_no_dstar()
  estimator <- match.arg(estimator)
  f <- function(nn) R2oos(nn, p, rho2_in, Dstar, estimator) - target
  lo <- p + 3
  if (!is.finite(f(max_n)) || f(max_n) < 0) return(NA_real_)
  if (is.finite(f(lo)) && f(lo) >= 0) return(lo)
  ## step out to a bracket, then bisect on the integers
  while (lo < max_n && (!is.finite(f(lo)) || f(lo) < 0)) lo <- lo * 2
  hi <- lo; lo <- max(p + 3, lo %/% 2)
  while (hi - lo > 1) {
    mid <- (lo + hi) %/% 2
    if (is.finite(f(mid)) && f(mid) >= 0) hi <- mid else lo <- mid
  }
  hi
}

#' Solve for the largest predictor set a design can carry
#'
#' The number of covariates at which a design of given size and signal still
#' reaches the target. Adding predictors beyond this is expected to make
#' out-of-sample prediction worse, not better.
#'
#' \strong{Search.} The largest qualifying \eqn{p} is found by bisection on
#' the integers between 1 and a cap: \eqn{n - 2} for least squares, the
#' largest \eqn{p} at which it is defined, and \eqn{10n} for ridge, which is
#' defined at any \eqn{p}. For least squares a target met at \eqn{n - 2} is
#' reported as \eqn{n - 2}. For ridge a target still met at the cap is
#' reported as \code{Inf}, with a message: the target is met for every
#' \eqn{p} up to \eqn{10n}, and the search does not establish a limit.
#'
#' @param n Sample size.
#' @param rho2_in Assumed in-sample signal. Held FIXED as \code{p} varies,
#'   which is conservative: in practice adding predictors usually raises the
#'   in-sample signal somewhat, and the honest comparison holds constant what
#'   the extra predictors are actually expected to contribute.
#' @param Dstar Design efficiency factor, held fixed as \code{p} varies.
#'   Required. For \code{estimator = "ridge"} this is \eqn{D^*(\kappa^*)}; see
#'   \code{\link{R2oos}}.
#' @param target Target expected \eqn{R^2_{oos}}.
#' @param estimator \code{"ols"} or \code{"ridge"}.
#' @return The largest admissible number of covariates; \code{NA} if even one
#'   covariate fails the target; \code{Inf} (ridge only) if the target is met
#'   at every \eqn{p} up to the search cap of \eqn{10n}.
#' @examples
#' solve_p(n = 620, rho2_in = 0.134, Dstar = 0.346)
#' solve_p(n = 620, rho2_in = 0.134, Dstar = 0.970, estimator = "ridge")
#' @export
solve_p <- function(n, rho2_in, Dstar, target = 0,
                    estimator = c("ols", "ridge")) {
  if (missing(Dstar)) .stop_no_dstar()
  estimator <- match.arg(estimator)
  f <- function(pp) {
    v <- R2oos(n, pp, rho2_in, Dstar, estimator)
    if (!is.finite(v)) -Inf else v - target
  }
  if (f(1) < 0) return(NA_integer_)
  hi <- if (estimator == "ols") n - 2 else 10 * n
  if (f(hi) >= 0) {
    if (estimator == "ols") return(as.integer(hi))
    message("The target is met for every p up to the search cap of 10n = ",
            hi, "; returning Inf.")
    return(Inf)
  }
  lo <- 1
  while (hi - lo > 1) {
    mid <- (lo + hi) %/% 2
    if (f(mid) >= 0) lo <- mid else hi <- mid
  }
  as.integer(lo)
}

#' Solve for the design efficiency a study would need
#'
#' How evenly leverage would have to be spread for the design to reach the
#' target. A value above 1 is impossible, and returns \code{NA}: no
#' rearrangement of the covariates can rescue that design, and only more data
#' or a smaller predictor set will.
#'
#' For \code{estimator = "ridge"} the quantity solved for is
#' \eqn{D^*(\kappa^*)}, the efficiency factor at the ridge penalty, not the
#' least-squares \eqn{D^*}.
#'
#' \strong{Search.} The root is found by \code{\link[stats]{uniroot}} on
#' \eqn{[10^{-6}, 1]}. If the target is already met at \eqn{10^{-6}}, the
#' result is 0, with a message that any \eqn{D^*} in (0, 1] meets it.
#'
#' @param n Sample size.
#' @param p Number of covariates.
#' @param rho2_in Assumed in-sample signal.
#' @param target Target expected \eqn{R^2_{oos}}.
#' @param estimator \code{"ols"} or \code{"ridge"}.
#' @return The required \eqn{D^*} (or \eqn{D^*(\kappa^*)} for ridge);
#'   \code{NA} when the target is unreachable even at 1; 0 when any value in
#'   (0, 1] reaches it.
#' @examples
#' solve_dstar(n = 620, p = 100, rho2_in = 0.134)
#' solve_dstar(n = 620, p = 100, rho2_in = 0.134, estimator = "ridge")
#' @export
solve_dstar <- function(n, p, rho2_in, target = 0,
                        estimator = c("ols", "ridge")) {
  estimator <- match.arg(estimator)
  f <- function(d) R2oos(n, p, rho2_in, d, estimator) - target
  if (!is.finite(f(1)) || f(1) < 0) return(NA_real_)
  lo <- 1e-6
  if (f(lo) >= 0) {
    message("The target is met at D* = ", lo, ", so any D* in (0, 1] meets ",
            "it; returning 0.")
    return(0)
  }
  stats::uniroot(f, c(lo, 1), tol = 1e-10)$root
}

#' Sensitivity table over a grid of assumptions
#'
#' Because \eqn{\rho^2_{in}} is assumed rather than known, a single analysis is
#' rarely enough. This crosses the design inputs and reports the expected
#' \eqn{R^2_{oos}} (or \eqn{r_{oos}}) and the break-even signal at each
#' combination, so the conclusion can be stated as a range rather than a point.
#'
#' Least squares and ridge use different efficiency factors: \code{Dstar}
#' (\eqn{D^*}) for least squares and \code{Dstar_ridge}
#' (\eqn{D^*(\kappa^*)}) for ridge. Each is required when its estimator is
#' requested. Without a design matrix, \eqn{D^*(\kappa^*)} is typically higher
#' than \eqn{D^*}, because ridge leverage is flatter; setting
#' \code{Dstar_ridge = Dstar} is conservative for ridge. Ridge rows use the
#' planning form of \code{\link{p_eff}}.
#'
#' @param n,p,rho2_in Design inputs. Supply vectors for those you want crossed
#'   and scalars for the rest.
#' @param Dstar Least-squares design efficiency factor \eqn{D^*}; a vector is
#'   crossed with the other inputs. Required when least-squares results are
#'   requested.
#' @param Dstar_ridge Design efficiency factor at the ridge penalty,
#'   \eqn{D^*(\kappa^*)}; a vector is crossed with the other inputs. Required
#'   when ridge results are requested.
#' @param estimator \code{"ols"}, \code{"ridge"} or \code{"both"}.
#' @param statistic \code{"R2oos"} or \code{"roos"}.
#' @return A data frame of class \code{"sensitivity"} with one row per
#'   combination. Its \code{Dstar} column holds the efficiency factor used in
#'   that row: \eqn{D^*} on least-squares rows and \eqn{D^*(\kappa^*)} on ridge
#'   rows.
#' @examples
#' sensitivity(n = c(300, 620, 1200), p = 100,
#'             rho2_in = c(0.07, 0.13, 0.20),
#'             Dstar = 0.346, Dstar_ridge = 0.970)
#' @export
sensitivity <- function(n, p, rho2_in, Dstar, Dstar_ridge,
                        estimator = c("both", "ols", "ridge"),
                        statistic = c("R2oos", "roos")) {
  estimator <- match.arg(estimator)
  statistic <- match.arg(statistic)
  ests <- if (estimator == "both") c("ols", "ridge") else estimator
  if ("ols" %in% ests && missing(Dstar)) .stop_no_dstar("Dstar")
  if ("ridge" %in% ests && missing(Dstar_ridge))
    stop("Supply Dstar_ridge, the design efficiency at the ridge penalty, ",
         "D*(kappa*). It is typically higher than Dstar, because ridge ",
         "leverage is flatter; setting it equal to Dstar gives conservative ",
         "ridge results. Use estimator = \"ols\" to skip ridge.", call. = FALSE)
  g <- do.call(rbind, lapply(ests, function(e)
    expand.grid(n = n, p = p, rho2_in = rho2_in,
                Dstar = if (e == "ols") Dstar else Dstar_ridge,
                estimator = e, stringsAsFactors = FALSE)))
  fn <- if (statistic == "R2oos") R2oos else roos
  g$value <- vapply(seq_len(nrow(g)), function(i) {
    pe <- if (g$estimator[i] == "ridge")
      p_eff(g$n[i], g$p[i], g$rho2_in[i]) else NULL
    fn(g$n[i], g$p[i], g$rho2_in[i], g$Dstar[i], g$estimator[i], pe)
  }, numeric(1))
  g$breakeven <- vapply(seq_len(nrow(g)), function(i)
    breakeven_rho2(g$n[i], g$p[i], g$Dstar[i], g$estimator[i]), numeric(1))
  names(g)[names(g) == "value"] <- statistic
  structure(g, class = c("sensitivity", "data.frame"), statistic = statistic)
}

#' @export
print.sensitivity <- function(x, ...) {
  cat("Sensitivity:", attr(x, "statistic"), "\n\n")
  print.data.frame(format(as.data.frame(x), digits = 3), row.names = FALSE)
  invisible(x)
}
