# predpower

Design-stage assessment of out-of-sample prediction accuracy — a power
calculation for prediction studies.

Power analysis lets you ask, before collecting anything, whether a planned
study can detect an effect. `predpower` asks the same question about
prediction: given a sample size, a predictor count, an assumed in-sample
signal and a design efficiency factor, what out-of-sample accuracy can this
design support in expectation?

The answer is often "less than you would guess", and it is available before
any outcome exists.

## Installation

```r
# install.packages("remotes")
remotes::install_github("marianolake/predpower")
```

No dependencies beyond base R. Building the vignette additionally needs
`knitr`, `rmarkdown` and pandoc.

```r
vignette("predpower")
```

## The four inputs and the four statistics

A design is described by four quantities: the sample size `n`, the number of
predictors `p`, the in-sample signal `rho2_in` (the population R² of the
fitted model), and a design efficiency factor. For ordinary least squares and
for ridge at the design-implied penalty κ\* = p(1−ρ²ᵢₙ)/ρ²ᵢₙ, the package
returns the expected out-of-sample R² and the expected correlation between
observed and predicted outcomes.

The two estimators use different efficiency factors. D\* describes
least-squares leverage; D\*(κ\*) describes ridge leverage at κ\*, which is
flatter, so D\*(κ\*) is typically the higher of the two. Neither has a
default: D\* = 1 means perfectly even leverage, which is optimistic for most
real designs. Without a covariate matrix, state both; using D\* for both is
conservative for ridge.

```r
library(predpower)

pred_power(n = 600, p = 100, rho2_in = 0.10, Dstar = 0.5, Dstar_ridge = 0.95)
#> Predictive power
#>
#>   n = 600, p = 100
#>   rho2_in = 0.1000
#>   least squares:  D* = 0.500
#>   ridge:          D*(kappa*) = 0.950,  p_eff = 38.4
#>
#>  estimator    R2oos     roos requires
#>        ols   -0.266    0.145    0.289
#>      ridge    0.032    0.179    0.031
```

Under least squares, expected R²ₒₒₛ is **negative**: this design is expected
to predict worse than the sample mean, and `requires` — the break-even
ρ²ᵢₙ — is 0.289, almost three times the assumed signal. Under ridge, expected
R²ₒₒₛ is 0.032 and the break-even is 0.031, just below the assumed signal.

These are expected values. A realized R²ₒₒₛ varies around them, so a design
that only just meets `requires` beats the sample mean in roughly half of
samples.

Note also that `roos` is 0.145 while R²ₒₒₛ is −0.266. Correlation does not
charge for miscalibration; R² does. A design can show a respectable
correlation while predicting worse than a constant.

## Three ways to describe a design

**No prior data.** Supply the four quantities directly, as above.

**Prior covariate data.** Supply `X`. Then n, p, D\*, D\*(κ\*) and the
effective number of ridge parameters p_eff are all computed from it, the
ridge quantities at whatever ρ²ᵢₙ they are needed. Factor, character and
logical columns are dummy-coded, and p counts the dummy columns.

```r
set.seed(1)
X <- matrix(rlnorm(300 * 20), 300, 20)
pred_power(X = X, rho2_in = 0.10)
#>   n = 300, p = 20   (from supplied X)
#>   rho2_in = 0.1000
#>   least squares:  D* = 0.862
#>   ridge:          D*(kappa*) = 0.938,  p_eff = 12.3   (at rho2_in, from X)
#>
#>  estimator    R2oos     roos requires
#>        ols    0.018    0.231    0.083
#>      ridge    0.052    0.229    0.026
```

**Prior signal estimate.** Supply a raw in-sample R² from earlier work. The
adjusted estimate and an exact confidence interval follow by inverting the
noncentral F, and the analysis runs at the point estimate and both ends of the
interval. The shipped `phonestudy` data carry the values for a published
study; here is its 100-feature design:

```r
d <- subset(phonestudy$designs, p == 100 & scale == "raw")
pred_power(n = d$n, p = d$p, R2_in = d$R2_in,
             Dstar = d$Dstar, Dstar_ridge = d$Dstar_ridge)
#>   n = 620, p = 100
#>   rho2_in = 0.1339   [0.073, 0.199], 95% interval
#>   least squares:  D* = 0.346
#>   ridge:          D*(kappa*) = 0.970,  p_eff = 46.9
#>
#>  estimator    R2oos     roos requires R2oos_lower R2oos_upper
#>        ols   -0.356    0.168    0.361      -0.452      -0.253
#>      ridge    0.058    0.240    0.027       0.018       0.111
```

Set the break-even `requires` beside the interval for ρ²ᵢₙ, [0.073, 0.199].
Least squares needs 0.361, above the whole interval. Ridge needs 0.027 with
the planning approximation used here, and 0.025 when p_eff and D\*(κ\*) come
from the design's own eigenvalues and leverage (`d$requires_ridge`); both lie
below the lower limit, so ridge is expected to beat the sample mean.

## Solving for what you don't know

For ridge, the `Dstar` argument of these functions is D\*(κ\*).

```r
solve_n(p = 100, rho2_in = 0.134, Dstar = 0.346)                       # 1995
solve_n(p = 100, rho2_in = 0.134, Dstar = 0.970, estimator = "ridge")  # 127
solve_p(n = 620, rho2_in = 0.134, Dstar = 0.346)                       # 30
solve_dstar(n = 620, p = 100, rho2_in = 0.134)                         # NA
breakeven_rho2(n = 620, p = 100, Dstar = 0.346)                        # 0.361
breakeven_rho2(n = 620, p = 100, Dstar = 0.970, estimator = "ridge")   # 0.027
```

`solve_dstar` returning `NA` is an answer, not a failure: no rearrangement of
the covariates can make that least-squares design break even, because even
perfectly even leverage leaves it short. Only more data or fewer predictors
will.

Because ρ²ᵢₙ is assumed rather than known, report a range:

```r
sensitivity(n = c(300, 620, 1200), p = 100,
            rho2_in = c(0.07, 0.13, 0.20),
            Dstar = 0.346, Dstar_ridge = 0.970)
```

## A worked example

`phonestudy` ships the design summaries from a published study predicting
the sociability facet of Extraversion from smartphone behavior — 620
participants, all 1,821 behavioral features (median-imputed), and nested
predictor sets of 25, 100 and 600 of them, each on its raw scale and after a
rank-normal transform. These are derived quantities (leverage, D\*,
D\*(κ\*), p_eff, break-even values, the signal estimate); the underlying
features are not redistributed.

```r
phonestudy$designs[, c("p", "scale", "rho2_in", "Dstar",
                       "requires_ols", "requires_ridge")]
```

The vignette walks through the full analysis on it.

## Design efficiency from covariates

D\* is the ratio of the prediction penalty an evenly balanced design would
incur to the one this design actually incurs. It is 1 when leverage is spread
evenly and falls toward 0 as it concentrates.

```r
set.seed(1)
X <- matrix(rlnorm(200 * 5), 200, 5)
design_efficiency(X)
#>   D*            0.711
#>   max leverage  0.741   (even-leverage value 0.030)
#>   |skew| > 2    100% of columns (median 3.2)
#>
#>   This design delivers 71% of the predictive efficiency of an even
#>   design of the same size. compare_scalings() shows whether a change of
#>   measurement scale raises it.
```

The dominant driver is usually skew, not imbalance — and skew is a property of
the measurement scale, not of the participants. So it is often recoverable:

```r
compare_scalings(X, rho2_in = 0.3)
#>        scale Dstar max_leverage R2oos
#>          raw 0.711       0.7414 0.266
#>          log 0.984       0.1822 0.274
#>  rank-normal 0.991       0.0841 0.275
```

Same data, same predictors, same sample: the rank-normal transform raises D\*
from 0.711 to 0.991, recovering 97% of the efficiency the raw scale lost.

## Plots

```r
d <- subset(phonestudy$designs, p == 100 & scale == "raw")
plot(pred_power(n = d$n, p = d$p, R2_in = d$R2_in,
                  Dstar = d$Dstar, Dstar_ridge = d$Dstar_ridge))
plot(pred_power(X = X, rho2_in = 0.3))
plot_breakeven(p = c(25, 100, 200), Dstar = 0.346, rho2_in = 0.134, mark_n = 620)
plot(design_efficiency(X))
```

Without `X`, the plot traces the design against sample size; with `X`,
whose sample size is fixed, against ρ²ᵢₙ. All figures encode every ordered
series in line type as well as colour, so they survive a black-and-white
printout.

## Scope and caveats

- The target is **leave-one-out** cross-validation. Refit-and-replicate gives a
  different and simpler answer, with no leverage term at all.
- ρ²ᵢₙ is **assumed, not estimated**. Like the effect size in a power analysis,
  the calculation is only as good as the assumption; use `sensitivity()`.
- Without `X`, ridge approximates p_eff by its value for **uncorrelated**
  covariates (the Marchenko-Pastur form). That is exact for uncorrelated
  covariates at any p/n; correlated or heavy-tailed covariates lower p_eff, so
  for them the ridge results are conservative. Supply `X` when you have it.
- The ridge results assume an isotropic prior on the coefficients. Simulation
  shows sparsity costs nothing measurable but collinearity does bite.
- Covariates are assumed measured without error.

## Reference

Tipton, E. Designing for Prediction: Design Formulas for
Expected Out-of-Sample R^2 and r.
