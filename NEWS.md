# predpower 0.2.0

* The package is renamed from `designaudit` to `predpower`, and its main
  function `design_audit()` to `pred_power()`, which returns an object of class
  `"pred_power"` (with `print()`, `summary()`, `plot()` and `as.data.frame()`
  methods). The other functions keep their names.

# designaudit 0.1.3

* New `rho2_reml()`: estimates rho2_in by restricted maximum likelihood under
  the isotropic ridge prior (Patterson and Thompson, 1971; as in Yang et al.,
  2010), with a profile-likelihood interval. It works for any p, including
  p + 1 >= n, where the adjusted R-squared is undefined.

* `rho2_from_data()`, and so `design_audit(X = , y = )`, use `rho2_reml()`,
  with a message, when p + 1 >= n, instead of stopping.

* `breakeven_rho2()` gains an `X` argument: with a covariate matrix and
  `estimator = "ridge"` it returns the threshold computed from the design's
  own ridge leverage, the value `design_audit()` reports.

# designaudit 0.1.2

* `Dstar` no longer has a default in `R2oos()`, `roos()`, `breakeven_rho2()`,
  `solve_n()`, `solve_p()`, `sensitivity()` or `plot_breakeven()` (which
  defaulted to 1, or 0.8 in `plot_breakeven()`). Omitting it is an error: no
  value is safe without a design to compute it from. For a single-estimator
  function with `estimator = "ridge"`, `Dstar` is D*(kappa*), and the
  documentation now says so.

* `sensitivity()` gains a `Dstar_ridge` argument, required when ridge results
  are requested; ridge rows use it in place of `Dstar`. `design_audit()`
  without `X` requires only the efficiency factor of each estimator
  requested.

* Data frames with factor, character or logical columns are dummy-coded with
  `model.matrix()` in `design_efficiency()`, `design_audit()`,
  `rho2_from_data()` and `p_eff()`, and `p` counts the dummy columns. They were
  previously converted to integer codes.

* `design_efficiency()` returns `Dstar = NA`, with a warning, when an
  observation has least-squares leverage 1 (for example the sole member of a
  category). It previously capped the leverage and returned a D* near 0.

* `design_audit()` with `p + 1 >= n` reports ridge alone (with a message) for
  `estimator = "both"`, and is an error for `estimator = "ols"`, instead of
  failing or printing a least-squares D*.

* `design_audit()` with `X`: values passed for `Dstar`, `Dstar_ridge` or `peff`
  are ignored with a warning, and every ridge quantity -- the point row,
  `requires` and the interval ends -- is computed from `X`. Previously a
  supplied `Dstar_ridge` or `peff` was used for the point row while the other
  columns used `X`.

* `plot.design_audit()` for an audit made with `X` traces the design against
  `rho2_in` using its own ridge leverage, so the plotted points lie on their
  curves. The audit object stores what this needs.

* The ridge break-even of a realized design is found on a 200-point grid over
  (1/n, 1) refined by `uniroot()`, with a warning and `NA` if there is no
  root. The previous 14-point grid gives the same root wherever it found one.
  Ridge leverage is computed from one singular value decomposition per design.

* `rho2_from_R2()` no longer returns the interval [0, 0] when the raw R-squared
  is at or below zero. The lower limit is 0; the upper limit comes from the
  same noncentral-F inversion as for any other R-squared, and is `NA`, with a
  warning, when the observed F lies below the lower alpha/2 point of the
  central F, where no non-negative signal satisfies it.

* `rho2_from_data()` uses the numerical rank of `cbind(1, X)` for the degrees
  of freedom, with a warning when it is below the number of columns.

* `solve_p(estimator = "ridge")` returns `Inf`, with a message, when the target
  is met at the search cap of 10n, instead of returning the cap.
  `solve_dstar()` returns 0, with a message, when any D* in (0, 1] meets the
  target.

* `phonestudy` is rebuilt from the analysis in the paper: the sociability facet
  (E2.Sociableness) as outcome, all 1,821 features with median imputation,
  designs of p = 25, 100, 600 and 1,821 on the raw and rank-normal scales, the
  REML estimate of rho2_in for p = 1,821, and D*(kappa*) and p_eff at the ridge
  break-even. The previous version held an analysis of the Extraversion domain
  score on 1,049 complete-case features. Column names changed; see
  `?phonestudy`.

* Documentation, README and vignette examples use the current `phonestudy`
  values. Corrected statements about the ridge break-even search (bracketed
  root-finding, not fixed-point iteration) and about the efficiency recovered
  by rescaling in the README and vignette examples.

# designaudit 0.1.1

* Results are reported as values only. `design_audit()` no longer labels a
  design "ruled out", "undecided" or "supported", and `design_efficiency()` no
  longer describes D* in categories; it reports the share of an even design's
  efficiency the design delivers.

* Without `X`, `design_audit()` now requires `Dstar_ridge` (D*(kappa*)) whenever
  ridge results are requested, instead of silently reusing `Dstar`. D*(kappa*)
  is typically higher than D*; setting them equal is conservative.

* `breakeven_rho2(estimator = "ridge")` now finds the threshold by bracketing.
  The fixed-point iteration it used could return NA where a threshold existed
  (for example n = 250, p = 600, Dstar = 0.6).

* When `X` is supplied, `design_audit()` now computes the ridge break-even
  (`requires`) and the ridge interval ends from the design itself -- p_eff from
  its eigenvalues and D*(kappa*) from its ridge leverage, both at the relevant
  rho2_in. Previously these used the planning approximations even when the
  design was available.

* Without a covariate matrix, `p_eff()` -- and so every ridge result in
  `R2oos()`, `roos()`, `design_audit()`, `breakeven_rho2()` and the solvers --
  now approximates the eigenvalues by those of uncorrelated covariates under
  the Marchenko-Pastur law, instead of setting every eigenvalue to n - 1. The
  old form ignored the sampling spread of the eigenvalues: it overstated
  `p_eff` even for uncorrelated covariates (by a median 33% at (p+1)/n = 0.5 and
  123% at 2 in the paper's simulation grid) and could exceed n - 1 when p > n,
  leaving the formulas undefined. The new form is exact for uncorrelated
  covariates at any p/n and never exceeds min(p, n - 1). The documentation
  previously described the old form as "mildly optimistic"; it was in fact
  conservative.

# designaudit 0.1.0

First release.

* `design_audit()` accepts a design three ways: stated directly, from a
  covariate matrix, or from a prior in-sample R-squared with its interval.
* `R2oos()` and `roos()` for least squares and for ridge at the
  design-implied penalty.
* `design_efficiency()` computes D* from covariates; `compare_scalings()`
  shows what a change of measurement scale would recover.
* `rho2_from_R2()` gives an exact interval for the in-sample signal by
  inverting the noncentral F.
* `solve_n()`, `solve_p()`, `solve_dstar()` and `breakeven_rho2()` invert the
  formulas for each input in turn.
* `sensitivity()` crosses design inputs; plot methods for every object.
