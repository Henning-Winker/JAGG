# JAGG

**ggplot2 graphics and diagnostics for [JABBA](https://github.com/jabbamodel/JABBA) and [JARA](https://github.com/henning-winker/JARA)**

JAGG draws ggplot2 versions of the JABBA and JARA plots directly from their
fitted objects, hindcasts and lists of scenario runs. JABBA and JARA are not
changed: their base-graphics plots stay as they are, and JAGG is an optional
companion.

All conversion from fit objects to tidy data frames lives in JAGG
(`jagg_tidy()`), so the same plotting code serves both models.

> Status: early development (0.1.0). Function names and arguments may still change.

## Installation

```r
# install.packages("remotes")
remotes::install_github("henning-winker/JAGG")
```

Fitting models needs [JAGS](https://mcmc-jags.sourceforge.io/) plus JABBA
and/or JARA. Plotting saved fits needs only JAGG and ggplot2.

## Usage

```r
library(JABBA); library(JAGG)
jbi  <- build_jabba(catch = iccat$swos$catch, cpue = iccat$swos$cpue, se = iccat$swos$se,
                    assessment = "SWOS", model.type = "Pella", BmsyK = 0.4)
jb   <- fit_jabba(jbi, quickmcmc = TRUE)
jbhc <- hindcast_jabba(jbi, jb, peels = 1:5)

jagg_advice(jb, rel = TRUE, refpts = c(Blim = 0.4, Bpa = 0.6))
jagg_runstest(jb)
jagg_hcxval(jbhc)
```

The full demo for both models is in
[`inst/examples/jagg_demo.R`](inst/examples/jagg_demo.R).

## Functions

| Function | Plot | JABBA | JARA |
|---|---|:-:|:-:|
| `jagg_advice()` | Trajectories with 90% and 60% credible bands and colour-coded reference points, after `FLRef::plotAdvice()`; one run or several scenarios | ✓ | ✓ |
| `jagg_fits()` | Index fits with credible and posterior predictive intervals | ✓ | ✓ |
| `jagg_residuals()` | Index residuals with yearly boxplots, loess smoother and RMSE | ✓ | ✓ |
| `jagg_runstest()` | Runs test with three-sigma control limits | ✓ | ✓ |
| `jagg_kobe()` | Kobe plot with posterior quadrant probabilities | ✓ | |
| `jagg_retro()` | Retrospective analysis with Mohn's rho | ✓ | ✓ |
| `jagg_hcxval()` | Hindcast cross-validation with MASE, after `ggplotFL::plotXval()` | ✓ | ✓ |
| `jagg_mase()` | MASE table | ✓ | ✓ |
| `jagg_tidy()` | Fit, hindcast or list of runs to tidy data frames | ✓ | ✓ |

### Checks against JABBA

On the demo fits:

- **Runs test:** `jagg_runs()` gives the same p-values and limits as JABBA's
  `jbruns_sig3()` (200 random series).
- **Mohn's rho:** matches `jbplot_retro()`.
- **MASE:** matches `jbmase()`.
  - `MASE.adj` follows the ss3diags definition, which floors the mean naive
    error at 0.1.
  - JABBA's `jbmase()` floors each naive residual instead, so it gives a
    lower `MASE.adj`.
- **Hindcast MASE with gaps:** `jagg_mase()` pairs prediction and naive
  residuals by year. Indices with gaps, or ending before the terminal year,
  are therefore handled correctly.

### Colours

- `r4col()` is a subset of `rich.colors()` by Arni Magnusson (gregmisc, now
  gplots), a.k.a. `r4ss::rich.colors.short()`.
- `jagg_cols()` follows the default colour rule of
  `r4ss::SSplotComparisons()`, as `sscol()` in ss3diags and `ss3col()` in
  JABBA do.

## Gallery

JABBA, South Atlantic swordfish (ICCAT example data, quick MCMC):

![](man/figures/01_jabba_advice.png)
![](man/figures/02_jabba_advice_rel.png)
![](man/figures/14_jabba_advice_scenarios.png)
![](man/figures/03_jabba_fits.png)
![](man/figures/04_jabba_residuals.png)
![](man/figures/05_jabba_runstest.png)
![](man/figures/06_jabba_kobe.png)
![](man/figures/07_jabba_retro.png)
![](man/figures/08_jabba_hcxval.png)

JARA, smoothhound shark:

![](man/figures/09_jara_advice.png)
![](man/figures/10_jara_fits.png)
![](man/figures/11_jara_runstest.png)
![](man/figures/12_jara_retro.png)
![](man/figures/13_jara_hcxval.png)
