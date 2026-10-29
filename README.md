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

`jagg_advice()` and `jagg_ensemble()` take the same arguments (`type`, `plotCIs`, `quantiles`, `inner`, `joint`, `bfrac`, `bref`, `reflines`, `reflabels`, `refpts`, `cols`, `alpha`, `xlim`, `ncol`, `fmax`, `prj`, `legend`), named after `jbplot_ensemble()` and `jbplot_summary()`. They differ only in defaults: nested 90%/60% bands with labelled reference points for the advice plot, and 95% bands over six panels for the ensemble. `plotCIs = FALSE` draws medians only.

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

Rendered with [`inst/examples/jagg_gallery.R`](inst/examples/jagg_gallery.R) (quick MCMC runs).

### JABBA: South Atlantic swordfish (ICCAT example data)

![](man/figures/jabba_advice.png)
![](man/figures/jabba_advice_rel.png)
![](man/figures/jabba_advice_scenarios.png)
![](man/figures/jabba_kobe.png)
![](man/figures/jabba_biplot.png)
![](man/figures/jabba_prj.png)
![](man/figures/jabba_ensemble.png)
![](man/figures/jabba_ensemble_noci.png)
![](man/figures/jabba_crisk.png)
![](man/figures/jabba_catch.png)
![](man/figures/jabba_procdev.png)
![](man/figures/jabba_spdyn.png)
![](man/figures/jabba_spphase.png)
![](man/figures/jabba_indices.png)
![](man/figures/jabba_fits.png)
![](man/figures/jabba_residuals.png)
![](man/figures/jabba_stdresiduals.png)
![](man/figures/jabba_runstest.png)
![](man/figures/jabba_qq.png)
![](man/figures/jabba_ppc.png)
![](man/figures/jabba_ppdist.png)
![](man/figures/jabba_mcmc.png)
![](man/figures/jabba_retro.png)
![](man/figures/jabba_hcxval.png)

### JARA: smoothhound shark

![](man/figures/jara_advice.png)
![](man/figures/jara_indices.png)
![](man/figures/jara_fits.png)
![](man/figures/jara_logfits.png)
![](man/figures/jara_runstest.png)
![](man/figures/jara_ppc.png)
![](man/figures/jara_retro.png)
![](man/figures/jara_retroiucn.png)
![](man/figures/jara_hcxval.png)
![](man/figures/jara_changes.png)
![](man/figures/jara_iucn.png)
![](man/figures/jara_jointindex.png)
![](man/figures/jara_poptrj.png)
![](man/figures/jara_state.png)
![](man/figures/jara_timeblock.png)
![](man/figures/jara_trjfit.png)

### JARA mixed trends: ICES Greater North Sea stocks

Data from Winker et al. (2026) ([ICES-Bmsy-FaF](https://github.com/Henning-Winker/ICES-Bmsy-FaF)); 17 stocks from 1980, by-order grid with Clupeiformes from all ecoregions.

![](man/figures/mixed_fits.png)
![](man/figures/mixed_grid.png)
![](man/figures/mixed_joint.png)
![](man/figures/mixed_joint_pr.png)
![](man/figures/mixed_kobe_trend.png)
![](man/figures/mixed_spaghetti.png)
![](man/figures/mixed_status_fao.png)
![](man/figures/mixed_status_kobe.png)

## Reference

Winker, H., Cardinale, M., Sharma, R., Kell, L. T., Mosqueira, I. and Griffiths, C. (2026). Assessing the progress of stock rebuilding in the Northeast Atlantic against levels that can produce maximum sustainable yield. *Fish and Fisheries* 27(3). https://doi.org/10.1111/faf.70066
