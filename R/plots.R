## Shared visual conventions. Colour is never the only channel: every ordered
## series carries its position redundantly in line type, so the figures survive
## a black-and-white printout or any form of colour blindness.

.pal <- function(k) {
  base <- c("#08343A", "#12626C", "#2E97A2", "#63B0B7", "#9ACFD4")
  if (k <= 5) base[round(seq(1, 5, length.out = k))]
  else grDevices::colorRampPalette(base)(k)
}
.lty <- function(k) c(1, 2, 4, 3, 5, 6)[((seq_len(k) - 1) %% 6) + 1]
.ink <- "#101619"; .rule <- "#C8D0D3"; .gy <- "#7A8890"

.frame <- function(xlim, ylim, xlab, ylab, main = NULL) {
  graphics::plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
  graphics::box(col = .rule)
  graphics::axis(1, col = .rule, col.axis = .gy, cex.axis = 0.8, tck = -0.02)
  graphics::axis(2, col = .rule, col.axis = .gy, cex.axis = 0.8, tck = -0.02,
                 las = 1)
  graphics::mtext(xlab, 1, line = 2.2, cex = 0.85, col = .ink)
  graphics::mtext(ylab, 2, line = 2.8, cex = 0.85, col = .ink)
  if (!is.null(main))
    graphics::mtext(main, 3, line = 0.6, cex = 0.9, col = .ink, adj = 0)
}

#' Plot a predictive power analysis
#'
#' Two panels, whose horizontal axis depends on how the analysis was made.
#'
#' \strong{Without a covariate matrix}, the panels trace the design against
#' sample size. The left shows expected \eqn{R^2_{oos}} against \eqn{n},
#' holding \eqn{p}, \eqn{\rho^2_{in}} and the efficiency factors fixed, and
#' marks where the analyzed design sits. The right traces the break-even signal
#' against \eqn{n} and overlays the assumed \eqn{\rho^2_{in}}, with its
#' confidence interval as a band when one is available. Ridge curves use the
#' planning form of \code{\link{p_eff}} at every \eqn{n}.
#'
#' \strong{With a covariate matrix} (\code{pred_power(X = ...)}), the sample
#' size is fixed by the design, so the panels trace it against
#' \eqn{\rho^2_{in}} instead. The left shows expected \eqn{R^2_{oos}} and the
#' right expected \eqn{r_{oos}}. Ridge curves are computed from the design
#' itself: at each \eqn{\rho^2_{in}}, \eqn{p_{eff}} and \eqn{D^*(\kappa^*)}
#' come from the design's ridge leverage at \eqn{\kappa^*(\rho^2_{in})}, as in
#' the analysis. Open circles mark the analyzed \eqn{\rho^2_{in}}; filled circles
#' mark the break-even values (\code{requires}), where each \eqn{R^2_{oos}}
#' curve crosses zero. The interval for \eqn{\rho^2_{in}}, when available, is
#' shaded.
#'
#' @param x A \code{pred_power} object.
#' @param n_range Range of sample sizes to trace when the analysis has no
#'   covariate matrix. Defaults to a window around the analyzed \code{n}.
#'   Ignored when it has one.
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for the plot.
#' @examples
#' d <- subset(phonestudy$designs, p == 100 & scale == "raw")
#' a <- pred_power(n = d$n, p = d$p, R2_in = d$R2_in,
#'                   Dstar = d$Dstar, Dstar_ridge = d$Dstar_ridge)
#' plot(a)
#'
#' set.seed(1)
#' plot(pred_power(X = matrix(rlnorm(300 * 20), 300, 20), rho2_in = 0.10))
#' @export
plot.pred_power <- function(x, n_range = NULL, ...) {
  if (x$has_X) return(.plot_power_X(x))
  if (is.null(n_range)) n_range <- c(max(10, x$p + 3), max(4 * x$n, x$p * 8))
  nn <- seq(n_range[1], n_range[2], length.out = 300)
  ests <- x$table$estimator
  cols <- .pal(max(2, length(ests))); ltys <- .lty(length(ests))
  dfor <- function(e) if (e == "ridge") x$Dstar_ridge else x$Dstar

  op <- graphics::par(mfrow = c(1, 2), mar = c(3.6, 4.2, 2.2, 0.8),
                      mgp = c(2, 0.5, 0))
  on.exit(graphics::par(op), add = TRUE)

  ## ---- panel (a): expected accuracy against n ----
  curves <- lapply(ests, function(e) R2oos(nn, x$p, x$rho2_in, dfor(e), e))
  fin <- unlist(curves); fin <- fin[is.finite(fin)]
  yl <- if (length(fin)) range(c(fin, 0)) else c(-1, 1)
  .frame(range(nn), yl, "sample size n", expression(expected~italic(R)[oos]^2),
         "(a) accuracy against sample size")
  graphics::abline(h = 0, col = .ink, lty = 3)
  for (i in seq_along(ests))
    graphics::lines(nn, curves[[i]], col = cols[i], lwd = 2, lty = ltys[i])
  for (i in seq_along(ests))
    graphics::points(x$n, R2oos(x$n, x$p, x$rho2_in, dfor(ests[i]), ests[i]),
                     pch = 21, bg = "white", col = cols[i], lwd = 1.6, cex = 1.1)
  graphics::legend("bottomright", legend = ests, col = cols[seq_along(ests)],
                   lty = ltys, lwd = 2, bty = "n", cex = 0.8, seg.len = 2.2)

  ## ---- panel (b): required signal against n ----
  be <- lapply(ests, function(e) breakeven_rho2(nn, x$p, dfor(e), e))
  fb <- unlist(be); fb <- fb[is.finite(fb)]
  hi_sig <- if (!is.null(x$signal) && is.finite(x$signal$upper)) x$signal$upper else 0
  yl2 <- c(0, max(c(fb, x$rho2_in, hi_sig), na.rm = TRUE))
  .frame(range(nn), yl2, "sample size n",
         expression("break-even " * rho[i*n]^2), "(b) signal required")
  if (!is.null(x$signal)) {
    top <- if (is.finite(x$signal$upper)) x$signal$upper else yl2[2]
    graphics::rect(min(nn), x$signal$lower, max(nn), top,
                   col = grDevices::adjustcolor(.gy, 0.15), border = NA)
  }
  graphics::abline(h = x$rho2_in, col = .ink, lty = 2)
  graphics::text(max(nn), x$rho2_in, "assumed signal", adj = c(1, -0.5),
                 cex = 0.72, col = .ink)
  for (i in seq_along(ests))
    graphics::lines(nn, be[[i]], col = cols[i], lwd = 2, lty = ltys[i])
  graphics::points(rep(x$n, length(ests)), x$table$requires,
                   pch = 21, bg = "white", col = cols[seq_along(ests)],
                   lwd = 1.6, cex = 1.1)
  invisible(x)
}

## The analysis of a realized design, traced against rho2_in at the design's own n.
.plot_power_X <- function(x) {
  ests <- x$table$estimator
  cols <- .pal(max(2, length(ests))); ltys <- .lty(length(ests))
  top <- max(c(x$rho2_in, x$table$requires,
               if (!is.null(x$signal)) x$signal$upper), na.rm = TRUE)
  rr <- seq(max(1 / x$n, 1e-4), min(0.99, 1.5 * top), length.out = 120)
  curve_at <- function(e, r, stat) {
    if (e == "ridge") {
      rq <- .ridge_X_at(x$ridge_basis, r)
      d <- rq$Dstar; pe <- rq$p_eff
    } else { d <- x$Dstar; pe <- NULL }
    if (stat == "R2oos") R2oos(x$n, x$p, r, d, e, pe)
    else vapply(seq_along(r), function(j)
      roos(x$n, x$p, r[j], if (length(d) > 1) d[j] else d, e,
           if (length(pe) > 1) pe[j] else pe), numeric(1))
  }
  op <- graphics::par(mfrow = c(1, 2), mar = c(3.6, 4.2, 2.2, 0.8),
                      mgp = c(2, 0.5, 0))
  on.exit(graphics::par(op), add = TRUE)

  for (stat in c("R2oos", "roos")) {
    cv <- lapply(ests, function(e) curve_at(e, rr, stat))
    pt <- vapply(ests, function(e) curve_at(e, x$rho2_in, stat), numeric(1))
    fin <- c(unlist(cv), pt); fin <- fin[is.finite(fin)]
    yl <- if (length(fin)) c(max(min(c(fin, 0)), -1), max(c(fin, 0.05))) else c(-1, 1)
    if (stat == "R2oos")
      .frame(range(rr), yl, expression("assumed " * rho[i*n]^2),
             expression(expected~italic(R)[oos]^2), "(a) accuracy against signal")
    else
      .frame(range(rr), yl, expression("assumed " * rho[i*n]^2),
             expression(expected~italic(r)[oos]), "(b) correlation against signal")
    if (!is.null(x$signal)) {
      up <- if (is.finite(x$signal$upper)) x$signal$upper else max(rr)
      graphics::rect(x$signal$lower, yl[1] - 1, up, yl[2] + 1,
                     col = grDevices::adjustcolor(.gy, 0.15), border = NA)
    }
    graphics::abline(h = 0, col = .ink, lty = 3)
    graphics::abline(v = x$rho2_in, col = .ink, lty = 2)
    for (i in seq_along(ests))
      graphics::lines(rr, cv[[i]], col = cols[i], lwd = 2, lty = ltys[i])
    graphics::points(rep(x$rho2_in, length(ests)), pt, pch = 21, bg = "white",
                     col = cols[seq_along(ests)], lwd = 1.6, cex = 1.1)
    if (stat == "R2oos") {
      ok <- is.finite(x$table$requires)
      graphics::points(x$table$requires[ok], rep(0, sum(ok)), pch = 21,
                       bg = cols[seq_along(ests)][ok], col = cols[seq_along(ests)][ok],
                       cex = 1.1)
      graphics::legend("bottomright", legend = ests, col = cols[seq_along(ests)],
                       lty = ltys, lwd = 2, bty = "n", cex = 0.8, seg.len = 2.2)
    }
  }
  invisible(x)
}

#' Plot the break-even signal against sample size
#'
#' The planning curve: how much in-sample signal a design would need to beat
#' the sample mean, as a function of sample size, for one or more predictor-set
#' sizes. Reading across at an assumed signal gives the sample size required;
#' reading up at a planned sample size gives the signal required.
#'
#' @param p Number of covariates. A vector draws one curve per value.
#' @param Dstar Design efficiency factor, held fixed across the curves.
#'   Required. For \code{estimator = "ridge"} this is \eqn{D^*(\kappa^*)}, the
#'   factor at the ridge penalty, which is typically higher than the
#'   least-squares \eqn{D^*}; see \code{\link{R2oos}}.
#' @param estimator \code{"ols"} or \code{"ridge"}.
#' @param n_range Range of sample sizes.
#' @param rho2_in Optional assumed signal, drawn as a horizontal reference.
#' @param mark_n Optional sample size, drawn as a vertical reference.
#' @return \code{NULL}, invisibly. Called for the plot.
#' @examples
#' plot_breakeven(p = c(25, 100, 200), Dstar = 0.346, rho2_in = 0.134,
#'                mark_n = 620)
#' plot_breakeven(p = c(25, 100, 200), Dstar = 0.970, estimator = "ridge",
#'                rho2_in = 0.134, mark_n = 620)
#' @export
plot_breakeven <- function(p = c(10, 25, 50, 100), Dstar,
                           estimator = c("ols", "ridge"),
                           n_range = c(30, 2000), rho2_in = NULL,
                           mark_n = NULL) {
  if (missing(Dstar)) .stop_no_dstar()
  estimator <- match.arg(estimator)
  nn <- seq(n_range[1], n_range[2], length.out = 400)
  cols <- .pal(length(p)); ltys <- .lty(length(p))
  vals <- lapply(p, function(pp) breakeven_rho2(nn, pp, Dstar, estimator))
  fin <- unlist(vals); fin <- fin[is.finite(fin)]
  .frame(range(nn), c(0, max(c(fin, rho2_in), na.rm = TRUE)),
         "sample size n", expression("break-even " * rho[i*n]^2),
         sprintf("%s,  %s = %.2f", estimator,
                 if (estimator == "ridge") "D*(kappa*)" else "D*", Dstar))
  if (!is.null(mark_n)) graphics::abline(v = mark_n, col = .ink, lty = 3)
  if (!is.null(rho2_in)) {
    graphics::abline(h = rho2_in, col = .ink, lty = 2)
    graphics::text(max(nn), rho2_in, "assumed signal", adj = c(1, -0.5),
                   cex = 0.72, col = .ink)
  }
  for (i in seq_along(p))
    graphics::lines(nn, vals[[i]], col = cols[i], lwd = 2, lty = ltys[i])
  graphics::legend("topright", legend = paste("p =", p), col = cols,
                   lty = ltys, lwd = 2, bty = "n", cex = 0.8, seg.len = 2.2)
  invisible(NULL)
}

#' Plot leverage concentration
#'
#' The sorted leverage of a design, with the even-leverage value marked. A
#' design with a few points far above the reference line is one whose
#' efficiency is being spent on a handful of observations.
#'
#' @param x A \code{design_efficiency} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for the plot.
#' @examples
#' set.seed(1)
#' plot(design_efficiency(matrix(rlnorm(200 * 5), 200, 5)))
#' @export
plot.design_efficiency <- function(x, ...) {
  h <- sort(x$leverage)
  op <- graphics::par(mar = c(3.6, 4.2, 2.2, 0.8), mgp = c(2, 0.5, 0))
  on.exit(graphics::par(op), add = TRUE)
  .frame(c(1, length(h)), c(0, max(h)), "observations, ordered by leverage",
         "leverage", sprintf("D* = %.3f", x$Dstar))
  graphics::lines(h, col = .pal(2)[1], lwd = 2)
  graphics::abline(h = x$even_leverage, col = .ink, lty = 3)
  graphics::text(1, x$even_leverage, "even leverage", adj = c(0, -0.5),
                 cex = 0.72, col = .ink)
  invisible(x)
}

#' Plot a sensitivity grid
#'
#' One curve per assumed signal, traced against whichever input was varied.
#'
#' @param x A \code{sensitivity} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for the plot.
#' @export
plot.sensitivity <- function(x, ...) {
  stat <- attr(x, "statistic")
  d <- as.data.frame(x)
  d <- d[d$estimator == d$estimator[1], ]
  grp <- unique(d$rho2_in); cols <- .pal(length(grp)); ltys <- .lty(length(grp))
  op <- graphics::par(mar = c(3.6, 4.2, 2.2, 0.8), mgp = c(2, 0.5, 0))
  on.exit(graphics::par(op), add = TRUE)
  .frame(range(d$n), range(c(d[[stat]], 0), na.rm = TRUE), "sample size n",
         stat, paste0(d$estimator[1], ",  p = ", d$p[1]))
  graphics::abline(h = 0, col = .ink, lty = 3)
  for (i in seq_along(grp)) {
    s <- d[d$rho2_in == grp[i], ]
    s <- s[order(s$n), ]
    graphics::lines(s$n, s[[stat]], col = cols[i], lwd = 2, lty = ltys[i])
  }
  graphics::legend("bottomright", legend = paste("rho2_in =", grp), col = cols,
                   lty = ltys, lwd = 2, bty = "n", cex = 0.8, seg.len = 2.2)
  invisible(x)
}
