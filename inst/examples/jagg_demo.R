# JAGG demo: fit the JABBA and JARA examples and draw the JAGG plots
# Runs in about 10 minutes (quick MCMC). Requires JAGS, JABBA and JARA:
#   remotes::install_github(c("jabbamodel/JABBA", "henning-winker/JARA"))

library(JABBA)
library(JARA)
library(JAGG)

#-----------------------------------------------------------
# JABBA: South Atlantic swordfish (ICCAT), Pella-Tomlinson
#-----------------------------------------------------------
set.seed(123)
jbi <- build_jabba(catch = iccat$swos$catch, cpue = iccat$swos$cpue, se = iccat$swos$se,
                   assessment = "SWOS", scenario = "Pella (BmsyK=0.4)", model.type = "Pella",
                   BmsyK = 0.4, sigma.est = TRUE, fixed.obsE = 0.1, verbose = FALSE)
jb   <- fit_jabba(jbi, quickmcmc = TRUE, do.ppc = TRUE, verbose = FALSE)
jbhc <- hindcast_jabba(jbi, jb, peels = 1:5, verbose = FALSE)

# Schaefer scenario for comparison
jbi2 <- build_jabba(catch = iccat$swos$catch, cpue = iccat$swos$cpue, se = iccat$swos$se,
                    assessment = "SWOS", scenario = "Schaefer", model.type = "Schaefer",
                    sigma.est = TRUE, fixed.obsE = 0.1, verbose = FALSE)
jb2  <- fit_jabba(jbi2, quickmcmc = TRUE, verbose = FALSE)

jagg_advice(jb, ncol = 3)                                        # B, F, catch
jagg_advice(jb, rel = TRUE, refpts = c(Blim = 0.4, Bpa = 0.6), ncol = 3)
jagg_advice(list(jb, jb2), rel = TRUE, ncol = 3)                 # scenarios
jagg_fits(jb)
jagg_residuals(jb)
jagg_runstest(jb)
jagg_kobe(jb)
jagg_retro(jbhc, type = c("B", "F", "B/Bmsy", "F/Fmsy"))
jagg_hcxval(jbhc)
jagg_mase(jbhc)

#-----------------------------------------------------------
# JARA: smoothhound shark, relative abundance indices
#-----------------------------------------------------------
data(jaradata)
jri  <- build_jara(I = jrdat$SmoothhoundShark$I, se = jrdat$SmoothhoundShark$SE,
                   model.type = "relative", assessment = "Smoothhound", GL = 13.7)
jr   <- fit_jara(jri, quickmcmc = TRUE, do.ppc = TRUE)
jrhc <- jara_hindcast(jri, peels = 1:5)

jagg_advice(jr)
jagg_fits(jr)
jagg_runstest(jr)
jagg_retro(jrhc, ref = jr)      # jara_hindcast() has no reference run: pass it
jagg_hcxval(jrhc, ref = jr)
jagg_mase(jrhc, ref = jr)

#-----------------------------------------------------------
# JARA mixed trends: ICES Greater North Sea stocks
# Data: ices.ts.cat1.csv from github.com/Henning-Winker/ICES-Bmsy-FaF
# (Winker et al. 2026, Fish and Fisheries, doi:10.1111/faf.70066)
#-----------------------------------------------------------
dat <- read.csv("ices.ts.cat1.csv")
dat <- dat[dat$Year >= 1980 & dat$ecoregion == "Greater North Sea", ]
wide <- function(v) {
  w <- reshape(dat[, c("Year", "stock", v)], idvar = "Year", timevar = "stock", direction = "wide")
  names(w) <- sub(paste0("^", v, "\\."), "", names(w)); w <- w[order(w$Year), ]; w[w == 0] <- NA
  w[, c(TRUE, colSums(!is.na(w[, -1])) >= 5)]
}
mixfit <- function(v, obsE = 0.05, ...) {
  b <- build_jara(I = wide(v), model.type = "mixed.trends", assessment = "ices_ns", scenario = v,
                  fixed.obsE = obsE, timeblock = 2015, backcast.na = 2000, ...)
  fit_jara(b, ni = 5000, nb = 1000, nt = 4, nc = 2, threshold = c(0.8, 1))
}
ind <- list(`F/Fmsy` = mixfit("FFmsy"), `B/Bmsy` = mixfit("BBmsy"),
            `B/Btrigger` = mixfit("BBtri"), ABImsy = mixfit("ABImsy", obsE = 0.2, proc.pen = 0.5))

jagg_joint(ind)                                   # geometric mean ratios
jagg_joint(ind, type = "pr", threshold.id = 2)    # proportion of stocks above 1 (F below 1)
jagg_status_area(ind$`B/Bmsy`, ind$`F/Fmsy`, backcast.na = 2000)
jagg_status_area(ind$`B/Bmsy`, scheme = "fao", backcast.na = 2000)
jagg_kobe_trend(ind$`B/Bmsy`, ind$`F/Fmsy`)
jagg_spaghetti(ind$`B/Bmsy`, ymax = 4)
jagg_fits(ind$`B/Bmsy`, ncol = 6, refline = 1, ppd = FALSE)
jagg_timeblock(ind$`B/Bmsy`)
