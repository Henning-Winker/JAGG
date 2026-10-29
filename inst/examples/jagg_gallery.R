# Render the JAGG gallery (README figures) from saved fits.
# Fits are produced by jagg_demo.R (JABBA swordfish, JARA smoothhound) and
# by fitting mixed-trend JARA models to ICES stocks (see jagg_demo.R).
# Usage: Rscript jagg_gallery.R <directory with .rds fits> <output directory>

args <- commandArgs(trailingOnly = TRUE)
D <- args[1]; O <- args[2]
library(JAGG)
rd <- function(f) readRDS(file.path(D, f))
sv <- function(p, f, w = 8, h = 5) {
  ggplot2::ggsave(file.path(O, f), p, width = w, height = h, dpi = 110)
  cat(f, "\n")
}

jb   <- rd("jabba_swos.rds");  s2 <- rd("jabba_swos_schaefer.rds"); jbhc <- rd("jabba_swos_hc.rds")
fwF  <- rd("fw_F.rds"); fwC <- rd("fw_C.rds"); fwE <- rd("fw_ens.rds")
jr   <- rd("jara_shsh.rds"); jrhc <- rd("jara_shsh_hc.rds")
ind  <- rd("jara_ns_indicators.rds"); ord <- rd("jara_ns_orders.rds")
jb$scenario <- "Pella"; jb$kbtrj$run <- "Pella"; s2$kbtrj$run <- "Schaefer"

# JABBA: status
sv(jagg_advice(jb, ncol = 3), "jabba_advice.png", 10, 3.6)
sv(jagg_advice(jb, rel = TRUE, refpts = c(Blim = 0.4, Bpa = 0.6), ncol = 3), "jabba_advice_rel.png", 10, 3.6)
sv(jagg_advice(list(jb, s2), rel = TRUE, ncol = 3), "jabba_advice_scenarios.png", 10, 3.8)
sv(jagg_kobe(jb), "jabba_kobe.png", 6.5, 5)
sv(jagg_biplot(jb), "jabba_biplot.png", 7, 5.5)
sv(jagg_prj(fwF), "jabba_prj.png", 10, 3.8)
sv(jagg_ensemble(list(Pella = jb, Schaefer = s2)), "jabba_ensemble.png", 10, 6)
sv(jagg_ensemble(list(Pella = jb, Schaefer = s2), plotCIs = FALSE), "jabba_ensemble_noci.png", 10, 6)
sv(jagg_crisk(fwC), "jabba_crisk.png", 10, 3.8)
sv(jagg_catch(jb), "jabba_catch.png", 7, 4)
sv(jagg_procdev(jb), "jabba_procdev.png", 7, 4)
sv(jagg_spdyn(jb), "jabba_spdyn.png", 7, 5)
sv(jagg_spphase(jb), "jabba_spphase.png", 7, 5)
# JABBA: diagnostics
sv(jagg_indices(jb), "jabba_indices.png", 9, 5)
sv(jagg_fits(jb), "jabba_fits.png", 9, 6)
sv(jagg_residuals(jb), "jabba_residuals.png", 8, 4.5)
sv(jagg_stdresiduals(jb), "jabba_stdresiduals.png", 8, 4.5)
sv(jagg_runstest(jb), "jabba_runstest.png", 9, 6)
sv(jagg_qq(jb), "jabba_qq.png", 9, 6)
sv(jagg_ppc(jb), "jabba_ppc.png", 9, 7)
sv(jagg_ppdist(jb), "jabba_ppdist.png", 9, 7)
sv(jagg_mcmc(jb), "jabba_mcmc.png", 9, 7)
sv(jagg_retro(jbhc, type = c("B", "F", "B/Bmsy", "F/Fmsy")), "jabba_retro.png", 9, 6)
sv(jagg_hcxval(jbhc), "jabba_hcxval.png", 9, 6.5)
# JARA
sv(jagg_advice(jr), "jara_advice.png", 7, 4)
sv(jagg_poptrj(jr), "jara_poptrj.png", 8, 4.5)
sv(jagg_trjfit(jr), "jara_trjfit.png", 8, 4.5)
sv(jagg_iucn(jr), "jara_iucn.png", 8, 4.5)
sv(jagg_retroiucn(jrhc, ref = jr), "jara_retroiucn.png", 8, 4.5)
sv(jagg_changes(jr), "jara_changes.png", 8, 4.5)
sv(jagg_state(jr), "jara_state.png", 8, 4.5)
sv(jagg_jointindex(jr), "jara_jointindex.png", 8, 4.5)
sv(jagg_indices(jr), "jara_indices.png", 9, 5)
sv(jagg_fits(jr), "jara_fits.png", 9, 6)
sv(jagg_logfits(jr), "jara_logfits.png", 9, 6)
sv(jagg_runstest(jr), "jara_runstest.png", 9, 6)
sv(jagg_ppc(jr), "jara_ppc.png", 9, 6)
sv(jagg_retro(jrhc, ref = jr), "jara_retro.png", 7, 4)
sv(jagg_hcxval(jrhc, ref = jr), "jara_hcxval.png", 9, 7)
# Mixed trends: ICES Greater North Sea stocks
i4 <- list(`F/Fmsy` = ind$FFmsy, `B/Bmsy` = ind$BBmsy, `B/Btrigger` = ind$BBtri, ABImsy = ind$ABImsy)
sv(jagg_joint(i4), "mixed_joint.png", 8, 4)
sv(jagg_joint(i4, type = "pr", threshold.id = 2), "mixed_joint_pr.png", 8, 4)
sv(jagg_joint_grid(lapply(ord, function(o) list(`F/Fmsy` = o$FFmsy, `B/Bmsy` = o$BBmsy))), "mixed_grid.png", 9, 5)
sv(jagg_status_area(ind$BBmsy, ind$FFmsy, backcast.na = 2000), "mixed_status_kobe.png", 8, 4)
sv(jagg_status_area(ind$BBmsy, scheme = "fao", backcast.na = 2000), "mixed_status_fao.png", 8, 4)
sv(jagg_kobe_trend(ind$BBmsy, ind$FFmsy), "mixed_kobe_trend.png", 6.5, 5)
sv(jagg_spaghetti(ind$BBmsy, ymax = 4), "mixed_spaghetti.png", 9, 4.5)
sv(jagg_fits(ind$BBmsy, ncol = 6, refline = 1, ppd = FALSE), "mixed_fits.png", 12, 7)
sv(jagg_timeblock(ind$BBmsy), "jara_timeblock.png", 7, 4)
