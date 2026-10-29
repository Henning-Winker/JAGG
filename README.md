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

#### Status and trajectories

| Function | Plot | Replaces | JABBA | JARA |
|---|---|---|:-:|:-:|
| `jagg_advice()` | Trajectories with 90% and 60% bands and colour-coded reference points, after `FLRef::plotAdvice()`; one run or scenarios | `jbplot_trj`, `jbplot_summary`, `jrplot_poptrj` | ✓ | ✓ |
| `jagg_kobe()` | Kobe plot with posterior quadrant probabilities | `jbplot_kobe` | ✓ | |
| `jagg_biplot()` | F/Fmsy vs B/Bmsy with 50/80/95% density regions | `jbplot_biplot` | ✓ | |
| `jagg_prj()` | Forecast scenarios from `fw_jabba()` | `jbplot_prj` | ✓ | |
| `jagg_ensemble()` | Several runs, ensembles or forecasts | `jbplot_ensemble` | ✓ | |
| `jagg_crisk()` | Risk against projected catch | `jbplot_Crisk` | ✓ | |
| `jagg_catch()` | Observed and estimated catch | `jbplot_catch`, `jbplot_catcherror` | ✓ | |
| `jagg_procdev()` | Process error deviations | `jbplot_procdev` | ✓ | |
| `jagg_spdyn()`, `jagg_spphase()` | Surplus production curve and phase plot | `jbplot_spdyn`, `jbplot_spphase` | ✓ | |
| `jagg_trjfit()`, `jagg_poptrj()` | Trajectory with scaled indices; with generation lengths and IUCN windows | `jrplot_trjfit`, `jrplot_poptrj` | | ✓ |
| `jagg_iucn()` | % change over 3 GL against the IUCN categories | `jrplot_iucn` | | ✓ |
| `jagg_retroiucn()` | IUCN status of each hindcast peel | `jrplot_retroiucn` | | ✓ |
| `jagg_changes()` | Annual rate of change, all years and 1-3 GL | `jrplot_changes`, `jrplot_r` | | ✓ |
| `jagg_state()` | Current and projected state | `jrplot_state` | | ✓ |
| `jagg_timeblock()` | Time-block effect | `jrplot_timeblock` | | ✓ |
| `jagg_jointindex()` | Joint index | `jrplot_jointindex` | | ✓ |

#### Model diagnostics

| Function | Plot | Replaces | JABBA | JARA |
|---|---|---|:-:|:-:|
| `jagg_indices()` | Input indices with 95% CIs | `jbplot_indices`, `jrplot_indices` | ✓ | ✓ |
| `jagg_fits()`, `jagg_logfits()` | Index fits, natural and log scale | `jbplot_cpuefits`, `jbplot_logfits`, `jrplot_fits`, `jrplot_logfits` | ✓ | ✓ |
| `jagg_residuals()`, `jagg_stdresiduals()` | Residuals with yearly boxplots, loess, RMSE / SDNR | `jbplot_residuals`, `jbplot_stdresiduals`, `jrplot_residuals` | ✓ | ✓ |
| `jagg_runstest()` | Runs test with three-sigma limits | `jbplot_runstest`, `jrplot_runstest` | ✓ | ✓ |
| `jagg_qq()` | QQ plot of standardised residuals | | ✓ | ✓ |
| `jagg_ppc()` | Posterior predictive checks | `jbplot_PPC`, `jrplot_PPC` | ✓ | ✓ |
| `jagg_ppdist()`, `jagg_bprior()` | Priors and posteriors with PPMR/PPVR | `jbplot_ppdist`, `jbplot_bprior` | ✓ | |
| `jagg_mcmc()` | MCMC traces | `jbplot_mcmc` | ✓ | |
| `jagg_retro()` | Retrospective analysis with Mohn's rho | `jbplot_retro`, `jrplot_retrobias` | ✓ | ✓ |
| `jagg_hcxval()`, `jagg_mase()` | Hindcast cross-validation with MASE (per index and joint) | `jbplot_hcxval`, `jbmase` | ✓ | ✓ |
| `jagg_table()` | Summary table of runs | | ✓ | ✓ |

#### Many stocks (JARA mixed trends), after Winker et al. (2026)

| Function | Plot |
|---|---|
| `jagg_joint()` | Joint trends of indicators across stocks: ratio or proportion of stocks (Fig. 2); scenarios with `prj = TRUE` (Figs 5-6) |
| `jagg_joint_grid()` | Joint trends by group and indicator (Fig. 3) |
| `jagg_status_area()` | Proportion of stocks by Kobe or FAO status (Fig. 4a,b) |
| `jagg_kobe_trend()` | Kobe trajectory of joint trends (Fig. 4c) |
| `jagg_spaghetti()` | Stock trajectories with final-year distribution |
| `jagg_tidy_mixed()` | Joint trend quantiles as a data frame |

#### Data and helpers

`jagg_tidy()`, `jagg_model()`, `jagg_runs()`, `jagg_iucn_status()`,
`theme_jagg()`, `scale_x_year()`, `r4col()`, `jagg_cols()`, `jagg_peelcols()`,
`jagg_pal_iucn()`, `jagg_pal_indicator()`, `jagg_pal_kobe()`, `jagg_pal_fao()`.

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
- **Joint MASE:** matches `jbmase()`.
- **JARA:** IUCN percentages, P(decline) and PPC p-values match the JARA
  functions.

### Deliberate differences from the base plots

- `jagg_spdyn()` plots production in year t, B[t+1] − B[t] + C[t], against
  B[t]. `jbplot_spdyn()` plots it against B[t+1].
- `jagg_crisk()` uses `bfrac` and `bref` in the risk itself.
  `jbplot_Crisk()` always uses 0.5 Bmsy and changes only the label.
- `jagg_prj()` works on current fits through `fw_jabba()` output.
  `jbplot_prj()` reads a `$projections` slot that `fit_jabba()` no longer
  creates.
- `jagg_iucn_status()` makes LC under criterion A1 start at −50, so the
  categories sum to 100%. `jrplot_iucn()` hard-codes −30.
- JARA hindcasts contain no reference run, so `jagg_retro()`,
  `jagg_hcxval()` and `jagg_retroiucn()` take it as `ref`.
- `jagg_ppdist()` computes PPMR/PPVR from analytic prior moments, so the
  values don't change between calls.
- Offsets are deterministic instead of random jitter.

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
