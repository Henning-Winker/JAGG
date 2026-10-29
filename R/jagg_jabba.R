# JABBA-specific plots

.need_jabba <- function(x, what) {
  if (jagg_model(x) != "jabba") stop(what, " needs a JABBA fit.", call. = FALSE)
}

# first, 70% and last years highlighted in trajectory plots (as JABBA)
.hl_years <- function(yrs) {
  n <- length(yrs)
  data.frame(year = yrs[c(1, round(stats::quantile(seq_len(n), 0.7)), n)], shape = c(22, 21, 24))
}

#' Catch
#'
#' Observed catch as a grey area and, when the fit estimated catch with
#' error (`add.catch.CV = TRUE`), the estimated catch with its 95% credible
#' interval (as `jbplot_catch()` and `jbplot_catcherror()`).
#' @param x a JABBA fit
#' @export
jagg_catch <- function(x) {
  .need_jabba(x, "jagg_catch()")
  d <- data.frame(year = x$yr, catch = x$catch)
  lab <- paste("Catch", x$settings$catch.metric)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year)) +
    ggplot2::geom_area(ggplot2::aes(y = .data$catch), fill = "grey80", colour = "grey30")
  if (is.data.frame(x$est.catch)) {
    e <- x$est.catch
    p <- p + ggplot2::geom_ribbon(data = e, ggplot2::aes(x = .data$yr, ymin = .data$lci, ymax = .data$uci),
                                  fill = "dodgerblue", alpha = 0.35) +
      ggplot2::geom_line(data = e, ggplot2::aes(x = .data$yr, y = .data$mu), colour = "blue", linewidth = 0.7) +
      ggplot2::geom_point(ggplot2::aes(y = .data$catch), shape = 21, fill = "white", size = 1.6)
  }
  p + scale_x_year() + ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = "Year", y = lab) + theme_jagg()
}

#' Process error deviations
#' @param x a JABBA fit
#' @export
jagg_procdev <- function(x) {
  .need_jabba(x, "jagg_procdev()")
  ts <- x$timeseries[, , "procB"]
  d <- data.frame(year = x$yr, mu = ts[, "mu"], lci = ts[, "lci"], uci = ts[, "uci"])
  lim <- max(0.22, abs(unlist(d[, -1])), na.rm = TRUE)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "grey75") +
    ggplot2::geom_hline(yintercept = 0, linetype = 5) +
    ggplot2::geom_line(ggplot2::aes(y = .data$mu), linewidth = 0.8) +
    ggplot2::coord_cartesian(ylim = c(-lim, lim)) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Process error deviates") + theme_jagg()
}

#' Standardised residuals
#'
#' As `jbplot_stdresiduals()`: standardised residuals of all indices with a
#' boxplot per year, a loess smoother through the yearly mean, and the SDNR
#' with its critical value sqrt(qchisq(0.95, N - 1) / (N - 1)).
#' @param x a JABBA fit
#' @export
jagg_stdresiduals <- function(x) {
  .need_jabba(x, "jagg_stdresiduals()")
  r <- x$std.residuals
  if (!is.matrix(r)) stop("No standardised residuals in this fit (catch-only model?).", call. = FALSE)
  d <- data.frame(index = rep(rownames(r), ncol(r)), year = rep(as.numeric(colnames(r)), each = nrow(r)),
                  residual = as.vector(r))
  d <- d[is.finite(d$residual), ]
  n <- nrow(d)
  sdnr <- sqrt(sum(d$residual^2) / (n - 1)); crit <- sqrt(stats::qchisq(0.95, n - 1) / (n - 1))
  ym <- stats::aggregate(residual ~ year, d, mean)
  lim <- max(1, 1.2 * max(abs(d$residual)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$year, y = .data$residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_boxplot(ggplot2::aes(group = .data$year), fill = "grey90", colour = "grey60",
                          outlier.shape = NA, width = 0.7) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$year, yend = 0, colour = .data$index),
                          position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$index), shape = 21, size = 1.8,
                        position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_smooth(data = ym, method = "loess", formula = y ~ x, se = FALSE,
                         colour = "black", linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = jagg_cols(length(unique(d$index)))) +
    ggplot2::scale_fill_manual(values = jagg_cols(length(unique(d$index)))) +
    ggplot2::coord_cartesian(ylim = c(-lim, lim)) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Standardised residuals",
                  subtitle = sprintf("SDNR = %.2f (%.2f)", sdnr, crit)) + theme_jagg()
}

#' Surplus production dynamics
#'
#' The expected production curve with the realised surplus production of each
#' year. Realised production in year t is B[t+1] - B[t] + C[t] and is plotted
#' against B[t]. (`jbplot_spdyn()` plots it against B[t+1].)
#' @param x a JABBA fit
#' @param rel plot B/K and SP/MSY
#' @export
jagg_spdyn <- function(x, rel = FALSE) {
  .need_jabba(x, "jagg_spdyn()")
  pf <- x$pfunc; K <- pf$SB0[1]
  bmsy <- x$estimates["SBmsy", 1]; msy <- x$estimates["MSY", 1]
  B <- x$timeseries[, "mu", "B"]; SP <- x$timeseries[, "mu", "SPt"]
  n <- length(B)
  tr <- data.frame(year = x$yr[-n], B = B[-n], SP = SP[-1])
  cv <- data.frame(B = pf$SB_i, SP = pf$SP)
  sx <- if (rel) K else 1; sy <- if (rel) msy else 1
  tr$B <- tr$B / sx; tr$SP <- tr$SP / sy; cv$B <- cv$B / sx; cv$SP <- cv$SP / sy
  hl <- merge(.hl_years(tr$year), tr, by = "year")
  ggplot2::ggplot(tr, ggplot2::aes(x = .data$B, y = .data$SP)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_hline(yintercept = msy / sy, linetype = 3, colour = "grey40") +
    ggplot2::geom_vline(xintercept = bmsy / sx, linetype = 3, colour = "grey40") +
    ggplot2::geom_line(data = cv, colour = "blue", linewidth = 0.9) +
    ggplot2::geom_path() + ggplot2::geom_point(size = 1) +
    ggplot2::geom_point(data = hl, ggplot2::aes(shape = factor(.data$year)), fill = "white", size = 3) +
    ggplot2::scale_shape_manual(values = stats::setNames(hl$shape, hl$year)) +
    ggplot2::expand_limits(x = 0) +
    ggplot2::labs(x = if (rel) "B/K" else "Biomass", y = if (rel) "SP/MSY" else "Surplus production",
                  shape = NULL) + theme_jagg()
}

# zone colours shared by the phase plots
.zone_fills <- c(green = "#2ca25f", yellow = "#fee391", orange = "#fd8d3c", red = "#de2d26")

#' Surplus production phase plot
#'
#' Catch against biomass over the production curve, with the zones of
#' `jbplot_spphase()`: left of Bmsy red above and yellow below the line
#' C = Fmsy B, right of Bmsy orange above and green below it.
#' @param x a JABBA fit
#' @param rel plot B/K and C/MSY
#' @export
jagg_spphase <- function(x, rel = FALSE) {
  .need_jabba(x, "jagg_spphase()")
  pf <- x$pfunc; K <- pf$SB0[1]; hmsy <- pf$Hmsy[1]
  bmsy <- x$estimates["SBmsy", 1]; msy <- x$estimates["MSY", 1]
  tr <- data.frame(year = x$yr, B = x$timeseries[, "mu", "B"], C = x$catch)
  cv <- data.frame(B = pf$SB_i, C = pf$SP)
  xmax <- max(cv$B, tr$B) * 1.02; ymax <- max(1.05 * tr$C, 1.1 * msy)
  g <- expand.grid(B = seq(0, xmax, length.out = 250), C = seq(0, ymax, length.out = 250))
  above <- g$C > hmsy * g$B
  g$zone <- ifelse(g$B < bmsy, ifelse(above, "red", "yellow"), ifelse(above, "orange", "green"))
  sx <- if (rel) K else 1; sy <- if (rel) msy else 1
  g$B <- g$B / sx; g$C <- g$C / sy; tr$B <- tr$B / sx; tr$C <- tr$C / sy; cv$B <- cv$B / sx; cv$C <- cv$C / sy
  hl <- merge(.hl_years(tr$year), tr, by = "year")
  ggplot2::ggplot(tr, ggplot2::aes(x = .data$B, y = .data$C)) +
    ggplot2::geom_raster(data = g, ggplot2::aes(fill = .data$zone), alpha = 0.85, interpolate = TRUE) +
    ggplot2::scale_fill_manual(values = .zone_fills, guide = "none") +
    ggplot2::geom_line(data = cv, colour = "blue", linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = msy / sy, colour = "blue", linetype = 5) +
    ggplot2::geom_vline(xintercept = bmsy / sx, colour = "blue", linetype = 2) +
    ggplot2::geom_path() + ggplot2::geom_point(size = 1) +
    ggplot2::geom_point(data = hl, ggplot2::aes(shape = factor(.data$year)), fill = "white", size = 3) +
    ggplot2::scale_shape_manual(values = stats::setNames(hl$shape, hl$year)) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(x = if (rel) "B/K" else "Biomass", y = if (rel) "Catch/MSY" else "Catch", shape = NULL) +
    theme_jagg()
}

# highest-density contour levels of a 2-D kernel density for given masses
.hdr_levels <- function(z, mass) {
  zs <- sort(as.vector(z), decreasing = TRUE)
  cm <- cumsum(zs) / sum(zs)
  vapply(mass, function(m) zs[which(cm >= m)[1]], numeric(1))
}

#' Biplot: F/Fmsy against B/Bmsy
#'
#' As `jbplot_biplot()`: zones defined by a target F of 0.8 Fmsy and a
#' biomass threshold of 0.2, the 50%, 80% and 95% highest-density regions of
#' the final-year posterior, the median trajectory, and the percentage of the
#' final-year posterior in each zone.
#' @param x a JABBA fit
#' @param ftarget F target relative to Fmsy (0.8)
#' @param bthreshold biomass threshold (0.2)
#' @export
jagg_biplot <- function(x, ftarget = 0.8, bthreshold = 0.2) {
  .need_jabba(x, "jagg_biplot()")
  zone <- function(f, b) ifelse(b > 1, ifelse(f < ftarget, "green", ifelse(f > 1, "red", "yellow")),
                                ifelse(b > bthreshold + f, "green",
                                       ifelse(b < bthreshold + ftarget * f, "red", "yellow")))
  trj <- data.frame(year = x$yr, f = x$timeseries[, "mu", "FFmsy"], b = x$timeseries[, "mu", "BBmsy"])
  kb <- x$kobe
  pct <- round(100 * prop.table(table(factor(zone(kb$harvest, kb$stock), names(.zone_fills)[c(1, 2, 4)]))), 1)
  xmax <- max(2, trj$f, stats::quantile(kb$harvest, 0.85)) * 1.02
  ymax <- max(x$refpts$k[1] / x$refpts$bmsy[1], trj$b) + 0.05
  g <- expand.grid(f = seq(0, xmax, length.out = 250), b = seq(0, ymax, length.out = 250))
  g$zone <- factor(zone(g$f, g$b), names(pct))
  kd <- MASS::kde2d(kb$harvest, kb$stock, n = 151,
                    lims = c(range(kb$harvest) + c(-1, 1) * 0.2 * diff(range(kb$harvest)),
                             range(kb$stock) + c(-1, 1) * 0.2 * diff(range(kb$stock))))
  lv <- .hdr_levels(kd$z, c(0.95, 0.8, 0.5))
  hdr <- do.call(rbind, lapply(seq_along(lv), function(i) {
    cl <- grDevices::contourLines(kd$x, kd$y, kd$z, levels = lv[i])
    do.call(rbind, lapply(seq_along(cl), function(j)
      data.frame(f = cl[[j]]$x, b = cl[[j]]$y, region = c("95%", "80%", "50%")[i],
                 grp = paste(i, j))))
  }))
  hl <- merge(.hl_years(trj$year), trj, by = "year")
  ggplot2::ggplot(trj, ggplot2::aes(x = .data$f, y = .data$b)) +
    ggplot2::geom_raster(data = g, ggplot2::aes(fill = .data$zone), alpha = 0.85, interpolate = TRUE) +
    ggplot2::scale_fill_manual(values = .zone_fills[names(pct)], labels = sprintf("%.1f%%", pct),
                               name = paste("Posterior", max(x$yr))) +
    ggplot2::geom_polygon(data = hdr, ggplot2::aes(group = .data$grp),
                          fill = c(`95%` = "cornsilk4", `80%` = "grey", `50%` = "cornsilk2")[hdr$region],
                          alpha = 0.9) +
    ggplot2::geom_vline(xintercept = 1, linetype = 3) + ggplot2::geom_hline(yintercept = 1, linetype = 3) +
    ggplot2::geom_path() +
    ggplot2::geom_point(data = hl, shape = hl$shape, fill = "white", size = 3) +
    ggplot2::geom_text(data = hl, ggplot2::aes(label = .data$year), vjust = -1, hjust = "inward", size = 3) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
    ggplot2::labs(x = expression(F/F[MSY]), y = expression(B/B[MSY]),
                  caption = "Shaded: 50%, 80% and 95% highest-density regions") +
    theme_jagg() + ggplot2::theme(legend.title = ggplot2::element_text(size = 8))
}

# prior density on a grid and analytic prior mean and CV
.prior <- function(par, s, xr) {
  ln <- function(p) list(d = function(v) stats::dlnorm(v, log(p[1]), p[2]),
                         mean = exp(log(p[1]) + p[2]^2 / 2), cv = sqrt(exp(p[2]^2) - 1),
                         q = function(pr) stats::qlnorm(pr, log(p[1]), p[2]))
  if (par == "K") return(ln(s$K.pr))
  if (par == "r") return(ln(s$r.pr))
  if (par == "m" && s$model.id == 4) return(ln(c(s$mu.m, s$m.CV)))
  if (par == "psi") {
    if (identical(s$psi.dist, "beta")) {
      a <- s$psi.pr[1]; b <- s$psi.pr[2]
      return(list(d = function(v) stats::dbeta(v, a, b), mean = a / (a + b),
                  cv = sqrt(a * b / ((a + b)^2 * (a + b + 1))) / (a / (a + b)),
                  q = function(pr) stats::qbeta(pr, a, b)))
    }
    return(ln(s$psi.pr))
  }
  if (par == "sigma2" && isTRUE(s$sigma.proc) && s$igamma[1] > 0.9) {
    a <- s$igamma[1]; b <- s$igamma[2]   # 1/sigma2 ~ Gamma(a, rate b)
    return(list(d = function(v) exp(a * log(b) - lgamma(a) - (a + 1) * log(v) - b / v),
                mean = if (a > 1) b / (a - 1) else NA, cv = if (a > 2) 1 / sqrt(a - 2) else NA,
                q = function(pr) 1 / stats::qgamma(1 - pr, a, rate = b)))
  }
  NULL
}

#' Prior and posterior distributions
#'
#' As `jbplot_ppdist()`: one panel per parameter with the prior (dark grey)
#' and the posterior (light), and the posterior-to-prior ratios of means
#' (PPMR) and of squared CVs (PPVR). Prior moments are computed analytically,
#' so the ratios do not change between calls. Fixed parameters are shown as
#' a vertical line.
#' @param x a JABBA fit
#' @param addPP print PPMR and PPVR
#' @param ncol number of facet columns
#' @export
jagg_ppdist <- function(x, addPP = TRUE, ncol = 3) {
  .need_jabba(x, "jagg_ppdist()")
  pp <- x$pars_posterior; s <- x$settings
  dens <- list(); fixed <- list(); lab <- list()
  for (par in names(pp)) {
    v <- pp[[par]]
    if (stats::sd(v) == 0) { fixed[[par]] <- data.frame(par = par, value = v[1]); next }
    pd <- stats::density(v, adjust = 2)
    pr <- .prior(par, s)
    xr <- range(pd$x)
    if (!is.null(pr)) xr <- range(xr, pr$q(c(1e-4, 0.95)))
    if (par %in% c("psi")) xr[1] <- max(0, min(xr[1], 0.5))
    xr[1] <- max(xr[1], 0)
    dens[[length(dens) + 1]] <- data.frame(par = par, type = "Posterior", x = pd$x, y = pd$y)
    if (!is.null(pr)) {
      gx <- seq(max(xr[1], 1e-12), xr[2], length.out = 400)
      dens[[length(dens) + 1]] <- data.frame(par = par, type = "Prior", x = gx, y = pr$d(gx))
      if (addPP) lab[[par]] <- data.frame(par = par, label = sprintf("PPMR = %.3f\nPPVR = %.3f",
        mean(v) / pr$mean, (stats::sd(v) / mean(v))^2 / pr$cv^2))
    }
  }
  d <- do.call(rbind, dens); fx <- do.call(rbind, fixed); lb <- do.call(rbind, lab)
  lev <- names(pp)
  d$par <- factor(d$par, lev); d$type <- factor(d$type, c("Prior", "Posterior"))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$type)) +
    ggplot2::geom_area(position = "identity", alpha = 0.7, colour = "grey30", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = c(Prior = "grey40", Posterior = "grey80"))
  if (!is.null(fx)) { fx$par <- factor(fx$par, lev)
    p <- p + ggplot2::geom_vline(data = fx, ggplot2::aes(xintercept = .data$value), linewidth = 0.9) }
  if (!is.null(lb)) { lb$par <- factor(lb$par, lev)
    p <- p + ggplot2::geom_text(data = lb, ggplot2::aes(x = Inf, y = Inf, label = .data$label),
                                inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, size = 2.6) }
  p + ggplot2::facet_wrap(~par, scales = "free", ncol = ncol, drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Density") + theme_jagg() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "bottom")
}

#' MCMC traces
#' @param x a JABBA fit
#' @param ncol number of facet columns
#' @export
jagg_mcmc <- function(x, ncol = 3) {
  .need_jabba(x, "jagg_mcmc()")
  pp <- x$pars_posterior
  pp <- pp[, vapply(pp, stats::sd, numeric(1)) > 0, drop = FALSE]
  d <- data.frame(par = factor(rep(names(pp), each = nrow(pp)), names(pp)),
                  iter = rep(seq_len(nrow(pp)), ncol(pp)), value = unlist(pp, use.names = FALSE))
  mu <- stats::aggregate(value ~ par, d, mean)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$iter, y = .data$value)) +
    ggplot2::geom_line(colour = "blue", linewidth = 0.15) +
    ggplot2::geom_hline(data = mu, ggplot2::aes(yintercept = .data$value), colour = "red", linewidth = 0.8) +
    ggplot2::facet_wrap(~par, scales = "free_y", ncol = ncol) +
    ggplot2::labs(x = "Iteration", y = NULL) + theme_jagg()
}

#' Biomass prior against its posterior
#'
#' As `jbplot_bprior()`, for fits with an informative prior on B/K, B/Bmsy
#' or F/Fmsy in one year (`b.prior` in `build_jabba()`).
#' @param x a JABBA fit
#' @export
jagg_bprior <- function(x) {
  .need_jabba(x, "jagg_bprior()")
  bp <- x$settings$b.pr
  if (bp[3] == 0 || !is.numeric(x$bppd)) stop("No biomass prior used in this fit.", call. = FALSE)
  pd <- stats::density(x$bppd, adjust = 2)
  xr <- range(pd$x, stats::qlnorm(c(1e-4, 0.95), log(bp[1]), bp[2])); xr[1] <- max(0, xr[1])
  gx <- seq(max(xr[1], 1e-9), xr[2], length.out = 400)
  d <- rbind(data.frame(type = "Prior", x = gx, y = stats::dlnorm(gx, log(bp[1]), bp[2])),
             data.frame(type = "Posterior", x = pd$x, y = pd$y))
  d$type <- factor(d$type, c("Prior", "Posterior"))
  ml <- c("B/K", "B/Bmsy", "F/Fmsy")[bp[4] + 1]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$type)) +
    ggplot2::geom_area(position = "identity", alpha = 0.7, colour = "grey30", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = c(Prior = "grey40", Posterior = "grey80")) +
    ggplot2::labs(x = paste(ml, "in", bp[3]), y = "Density") + theme_jagg()
}

#' Ensemble or multi-run trajectories
#'
#' As `jbplot_ensemble()`: medians and credible bands of several runs,
#' scenarios or forecasts, coloured with [jagg_cols()].
#' @param x a JABBA fit, a list of fits, or `fw_jabba()` output
#' @param type quantities: "B/Bmsy", "F/Fmsy", "B", "F", "Process deviation",
#'   "Catch", "B/B0"
#' @param joint pool all runs into one ("Joint")
#' @param level credible interval level
#' @param fmax cap for F/Fmsy in projections (as `jbplot_ensemble()`)
#' @param ncol number of facet columns
#' @export
jagg_ensemble <- function(x, type = c("B/Bmsy", "F/Fmsy", "B", "F", "Process deviation", "Catch"),
                          joint = FALSE, level = 0.95, fmax = 5, ncol = 3) {
  m <- jagg_model(x)
  kb <- switch(m, jabba = x$kbtrj, jabba_kb = x,
               jabba_list = do.call(rbind, lapply(seq_along(x), function(i) {
                 k <- x[[i]]$kbtrj
                 k$run <- if (!is.null(names(x)) && nzchar(names(x)[i])) names(x)[i] else x[[i]]$scenario
                 k })),
               stop("x must be a JABBA fit, a list of fits or fw_jabba() output.", call. = FALSE))
  kb$run <- as.character(kb$run)
  if (joint) kb$run <- "Joint"
  if ("harvest" %in% names(kb)) kb$harvest[kb$type == "prj"] <- pmin(kb$harvest[kb$type == "prj"], fmax)
  it <- .iters(kb)
  if ("Bdev" %in% names(kb))
    it <- rbind(it, data.frame(run = kb$run, year = kb$year, type = kb$type, iter = kb$iter,
                               quantity = "Process deviation", value = kb$Bdev))
  it <- it[it$quantity %in% type, ]
  a <- (1 - level) / 2
  q <- .qbands(it, c(a, 0.5, 0.5, 0.5, 1 - a))
  runs <- unique(kb$run)
  q$run <- factor(q$run, runs); q$quantity <- factor(q$quantity, type)
  cols <- stats::setNames(jagg_cols(length(runs)), runs)
  p <- ggplot2::ggplot(q, ggplot2::aes(x = .data$year, y = .data$med, colour = .data$run, fill = .data$run)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo2, ymax = .data$hi2), colour = NA, alpha = 0.25) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = ncol) +
    scale_x_year() + ggplot2::labs(x = "Year", y = NULL) + theme_jagg()
  ref <- data.frame(quantity = factor(intersect(type, c("B/Bmsy", "F/Fmsy", "Process deviation")), type))
  ref$y <- ifelse(ref$quantity == "Process deviation", 0, 1)
  if (nrow(ref)) p <- p + ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = .data$y), linetype = 2)
  if (any(type != "Process deviation")) p <- p + ggplot2::expand_limits(y = 0)
  p
}

#' Risk of projected catches
#'
#' For `fw_jabba(quant = "Catch")` projections: the probability that B is
#' below Bmsy, F above Fmsy, or B below a fraction of Bmsy or B0, against
#' the catch, one line per projection year (as `jbplot_Crisk()`, but using
#' `bfrac` and `bref` in the risk itself).
#' @param x `fw_jabba()` output with catch scenarios
#' @param type any of "B<Bmsy", "F>Fmsy", "B<Bfrac"
#' @param bfrac fraction of the biomass reference for "B<Bfrac"
#' @param bref "bmsy" or "b0"
#' @param years projection years to show (default: up to five, evenly spaced)
#' @export
jagg_crisk <- function(x, type = c("B<Bmsy", "F>Fmsy", "B<Bfrac"), bfrac = 0.5,
                       bref = c("bmsy", "b0"), years = NULL) {
  bref <- match.arg(bref)
  if (jagg_model(x) != "jabba_kb") stop("x must be fw_jabba() output.", call. = FALSE)
  k <- x[x$type == "prj", ]
  if (is.null(years)) {
    # skip implementation-lag years, where all scenarios still have the same catch
    mc <- tapply(k$Catch, list(k$year, k$run), stats::median)
    y <- as.numeric(rownames(mc))[apply(mc, 1, function(v) diff(range(v, na.rm = TRUE)) > 0)]
    years <- unique(y[round(seq(1, length(y), length.out = min(5, length(y))))])
  }
  k <- k[k$year %in% years, ]
  bf <- if (bref == "bmsy") k$stock / bfrac else k$BB0 / bfrac
  key <- paste(k$run, k$year)
  agg <- function(v, nm) data.frame(key = names(v), risk = as.vector(v), type = nm)
  r <- rbind(agg(tapply(k$stock < 1, key, mean), "B<Bmsy"),
             agg(tapply(k$harvest > 1, key, mean), "F>Fmsy"),
             agg(tapply(bf < 1, key, mean), "B<Bfrac"))
  cat <- tapply(k$Catch, key, stats::median)
  r$catch <- cat[r$key]
  r$year <- factor(sub(".* ", "", r$key))
  r <- r[r$type %in% type, ]
  lb <- c("B<Bmsy" = "P(B < Bmsy)", "F>Fmsy" = "P(F > Fmsy)",
          "B<Bfrac" = sprintf("P(B < %s %s)", bfrac, if (bref == "bmsy") "Bmsy" else "B0"))
  r$type <- factor(lb[r$type], lb[type])
  r <- r[order(r$catch), ]
  ggplot2::ggplot(r, ggplot2::aes(x = .data$catch, y = .data$risk, colour = .data$year)) +
    ggplot2::geom_hline(yintercept = c(0.2, 0.8), linetype = 3) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = jagg_cols(length(levels(r$year)))) +
    ggplot2::facet_wrap(~type) + ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = paste("Catch", "(projection scenario)"), y = "Probability") + theme_jagg()
}
