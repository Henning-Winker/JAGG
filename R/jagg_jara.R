# JARA-specific plots

.need_jara <- function(x, what) {
  if (jagg_model(x) != "jara") stop(what, " needs a JARA fit.", call. = FALSE)
}

#' IUCN Red List colours
#' @return named vector CR, EN, VU, NT, LC
#' @export
jagg_pal_iucn <- function() c(CR = "#D81E05", EN = "#FC7F3F", VU = "#F9E814", NT = "#CCE226", LC = "#60C659")

# IUCN thresholds (upper bounds of % change) under criteria A1 or A2-A4
.iucn_cuts <- function(criteria = "A2", NT_opt = FALSE) {
  cuts <- if (criteria == "A1") c(CR = -90, EN = -70, VU = -50, NT = -40) else c(CR = -80, EN = -50, VU = -30, NT = -20)
  if (!NT_opt) cuts <- cuts[c("CR", "EN", "VU")]
  cuts
}

#' IUCN category probabilities from the posterior of % change
#'
#' Categories as in `jrplot_iucn()`: CR if change <= the CR threshold, EN,
#' VU (and NT with `NT_opt`) up to their thresholds, LC above the last one.
#' Unlike `jrplot_iucn()`, LC under criterion A1 starts above -50 (the VU
#' threshold), so the categories always sum to 100%. The status is the
#' category with the largest share, and NT when that is LC with less than 50%.
#' @param change posterior draws of % change over 3 generations
#' @param criteria "A2" (A2-A4, default) or "A1"
#' @param NT_opt include Near Threatened
#' @return list with `perc` (named %), `status`, `change` (median) and
#'   `p.decline` (% of draws below 0)
#' @export
jagg_iucn_status <- function(change, criteria = c("A2", "A1"), NT_opt = FALSE) {
  criteria <- match.arg(criteria)
  cuts <- .iucn_cuts(criteria, NT_opt)
  cat <- cut(change, c(-Inf, cuts, Inf), labels = c(names(cuts), "LC"), right = TRUE)
  perc <- 100 * prop.table(table(cat))
  k <- names(perc)[which.max(perc)]
  status <- if (k == "LC" && max(perc) < 50) "NT" else k
  list(perc = round(c(perc), 1), status = status, change = stats::median(change),
       p.decline = 100 * mean(change < 0))
}

#' IUCN Red List assessment from the % change over three generations
#'
#' As `jrplot_iucn()`: the posterior of the % change in abundance over three
#' generation lengths over the IUCN category bands, with the median and the
#' percentage of the posterior in each category.
#' @param x a JARA fit
#' @param criteria "A2" (A2-A4) or "A1"
#' @param NT_opt include Near Threatened
#' @param xlim x-axis limits (default -100 to the 99% quantile, at least 30)
#' @return a ggplot object; the result of [jagg_iucn_status()] in attribute "iucn"
#' @export
jagg_iucn <- function(x, criteria = c("A2", "A1"), NT_opt = FALSE, xlim = NULL) {
  .need_jara(x, "jagg_iucn()")
  criteria <- match.arg(criteria)
  ch <- x$posteriors$pop.change
  st <- jagg_iucn_status(ch, criteria, NT_opt)
  if (is.null(xlim)) xlim <- c(-100, min(max(30, stats::quantile(ch, 0.99)), 1000))
  cuts <- .iucn_cuts(criteria, NT_opt)
  bands <- data.frame(cat = c(names(cuts), "LC"), xmin = c(-Inf, cuts), xmax = c(cuts, Inf))
  bands$cat <- factor(bands$cat, names(st$perc))
  den <- stats::density(ch, adjust = 2)
  dd <- data.frame(x = den$x, y = den$y)
  pal <- jagg_pal_iucn()[names(st$perc)]
  mp <- x$settings$mp.assess
  win <- if (!is.null(mp)) sprintf("Change between %s and %s (3-year medians), criteria %s",
                                   x$pyr[mp[1]], x$pyr[mp[2]], criteria) else NULL
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf,
                       ymax = Inf, fill = .data$cat)) +
    ggplot2::scale_fill_manual(values = pal, labels = sprintf("%s (%s%%)", names(st$perc), st$perc),
                               name = NULL) +
    ggplot2::geom_area(data = dd, ggplot2::aes(x = .data$x, y = .data$y), fill = "grey75",
                       colour = "black", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = st$change, linetype = 2, colour = "grey30") +
    ggplot2::annotate("label", x = st$change, y = max(dd$y) * 1.1, size = 3.4, fontface = "bold",
                      fill = "white", label = sprintf("Change = %s%.0f%%", ifelse(st$change > 0, "+", ""), st$change)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, max(dd$y) * 1.2), expand = FALSE) +
    ggplot2::labs(x = "Change (%)", y = "Density", caption = win,
                  subtitle = sprintf("Status: %s   P(decline) = %.1f%%", st$status, st$p.decline)) +
    theme_jagg() + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                  axis.ticks.y = ggplot2::element_blank())
  attr(p, "iucn") <- st
  p
}

#' Retrospective IUCN assessment
#'
#' As `jrplot_retroiucn()`: half-violins of the % change posterior for the
#' reference run and each peel of a JARA hindcast, coloured by IUCN
#' category, with the status of each run.
#' @param hc a JARA hindcast (`jara_hindcast()`)
#' @param ref the reference fit (`jara_hindcast()` has no reference run)
#' @param criteria "A2" or "A1"
#' @return a ggplot object; the status table in attribute "iucn"
#' @export
jagg_retroiucn <- function(hc, ref, criteria = c("A2", "A1")) {
  if (jagg_model(hc) != "jara_hc") stop("hc must be a JARA hindcast.", call. = FALSE)
  criteria <- match.arg(criteria)
  po <- rbind(data.frame(peel = 0, change = ref$posteriors$pop.change),
              data.frame(peel = as.numeric(as.character(hc$posteriors$level)), change = hc$posteriors$pop.change))
  ty <- max(ref$yr)
  peels <- sort(unique(po$peel), decreasing = TRUE)        # oldest end year first
  ymax <- min(max(30, stats::quantile(po$change, 0.99)), 1000)
  cuts <- .iucn_cuts(criteria)
  gy <- seq(-100, ymax, length.out = 300); dy <- diff(gy)[1]
  tab <- list(); sl <- list()
  for (j in seq_along(peels)) {
    v <- po$change[po$peel == peels[j]]
    st <- jagg_iucn_status(v, criteria)
    tab[[j]] <- data.frame(year = ty - peels[j], change = round(st$change, 1), t(st$perc), status = st$status)
    dn <- stats::density(v, adjust = 1, from = -100, to = ymax, n = 300)
    w <- dn$y / max(dn$y) * 0.9
    sl[[j]] <- data.frame(x = j, ymin = gy - dy / 2, ymax = gy + dy / 2, w = w,
                          cat = as.character(cut(gy, c(-Inf, cuts, Inf), labels = c(names(cuts), "LC"))))
  }
  s <- do.call(rbind, sl); tb <- do.call(rbind, tab)
  s$cat <- factor(s$cat, c("LC", rev(names(cuts))))
  lab <- ifelse(peels == 0, paste0(ty, "\n(ref)"), ty - peels)
  edge <- do.call(rbind, lapply(split(s, s$x), function(si)
    data.frame(x = c(si$x[1] - 0.45, si$x - 0.45 + si$w, si$x[1] - 0.45),
               y = c(si$ymin[1], (si$ymin + si$ymax) / 2, si$ymax[nrow(si)]), grp = si$x[1])))
  pal <- jagg_pal_iucn()[levels(s$cat)]
  ggplot2::ggplot(s) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$x - 0.45, xmax = .data$x - 0.45 + .data$w,
                       ymin = .data$ymin, ymax = .data$ymax, fill = .data$cat, colour = .data$cat),
                       linewidth = 0.3) +
    ggplot2::geom_path(data = edge, ggplot2::aes(x = .data$x, y = .data$y, group = .data$grp), linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::annotate("text", x = seq_along(peels), y = ymax, label = tb$status, vjust = 1.2, size = 3.3) +
    ggplot2::scale_x_continuous(breaks = seq_along(peels), labels = lab) +
    ggplot2::coord_cartesian(ylim = c(-100, ymax)) +
    ggplot2::labs(x = "Assessment end year", y = "Change over 3 GL (%)") + theme_jagg() -> p
  attr(p, "iucn") <- tb
  p
}

#' Annual rate of change
#'
#' As `jrplot_changes()` (and `jrplot_r()` with `scale = "r"`): posteriors of
#' the median annual rate of change over all years and the most recent one,
#' two and three generation lengths.
#' @param x a JARA fit
#' @param scale "percent" (annual % change, (exp(r) - 1) * 100) or "r"
#' @export
jagg_changes <- function(x, scale = c("percent", "r")) {
  .need_jara(x, "jagg_changes()")
  scale <- match.arg(scale)
  rs <- x$posteriors[, -1, drop = FALSE]
  nm <- c(r.All.yrs = "All years", r.1GL = "1 GL", r.2GL = "2 GL", r.3GL = "3 GL")[names(rs)]
  v <- if (scale == "percent") (exp(as.matrix(rs)) - 1) * 100 else as.matrix(rs)
  xl <- stats::quantile(v, c(0.001, 0.999))
  med <- apply(v, 2, stats::median)
  unit <- if (scale == "percent") "%" else ""
  labs <- sprintf("%s: %s%.*f%s", nm, ifelse(med > 0, "+", ""), if (scale == "percent") 2 else 3, med, unit)
  d <- do.call(rbind, lapply(seq_along(nm), function(i) {
    dn <- stats::density(v[, i], adjust = 2, from = xl[1], to = xl[2])
    data.frame(series = labs[i], x = dn$x, y = dn$y)
  }))
  d$series <- factor(d$series, labs)
  cols <- stats::setNames(c("grey45", "#0000FF", "#00CD00", "#FF0000")[seq_along(nm)], labs)
  md <- data.frame(series = factor(labs, labs), x = med)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$series)) +
    ggplot2::geom_area(position = "identity", alpha = 0.3) +
    ggplot2::geom_vline(data = md, ggplot2::aes(xintercept = .data$x, colour = .data$series), linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = if (scale == "percent") "Annual rate of change (%)" else "r", y = "Density") +
    theme_jagg() + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                  axis.ticks.y = ggplot2::element_blank())
}

#' Current and projected state relative to reference years
#'
#' As `jrplot_state()`: posteriors of abundance in the last assessment year
#' and the last projection year relative to the mean of the reference years.
#' @param x a JARA fit
#' @param type "both", "current" or "projected" (default: both if the
#'   projection extends at least 3 years)
#' @param ref.yr reference years (default: the first three)
#' @param extinction threshold of relative state for P(extinction)
#' @param level credible interval level for the table
#' @return a ggplot object; the table in attribute "state"
#' @export
jagg_state <- function(x, type = NULL, ref.yr = NULL, extinction = 0.01, level = 0.95) {
  .need_jara(x, "jagg_state()")
  pd <- as.matrix(x$pop.posterior); yr <- as.numeric(colnames(pd))
  if (is.null(ref.yr)) ref.yr <- yr[1:3]
  ey <- max(x$yr); py <- max(x$pyr)
  if (is.null(type)) type <- if (py - ey < 3) "current" else "both"
  ref <- rowMeans(pd[, yr %in% ref.yr, drop = FALSE])
  st <- list(Current = pd[, yr == ey] / ref, Projected = pd[, yr == py] / ref)
  yrs <- c(Current = ey, Projected = py)
  a <- (1 - level) / 2
  tab <- do.call(rbind, lapply(names(st), function(n) data.frame(state = n, year = yrs[[n]],
    median = stats::median(st[[n]]), lci = stats::quantile(st[[n]], a), uci = stats::quantile(st[[n]], 1 - a))))
  rownames(tab) <- NULL
  show <- switch(type, both = c("Current", "Projected"), current = "Current", projected = "Projected")
  xmax <- max(1.1, vapply(st[show], stats::quantile, numeric(1), 0.99))
  labs <- stats::setNames(sprintf("%s (%s): %.2f", show, yrs[show], tab$median[match(show, tab$state)]), show)
  d <- do.call(rbind, lapply(show, function(n) {
    dn <- stats::density(st[[n]], adjust = 2, from = 0, to = xmax)
    data.frame(state = labs[[n]], x = dn$x, y = dn$y)
  }))
  d$state <- factor(d$state, labs)
  cols <- stats::setNames(c(Current = "grey35", Projected = "#FF0000")[show], labs)
  md <- data.frame(state = factor(labs, labs), x = tab$median[match(show, tab$state)])
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$state)) +
    ggplot2::geom_area(position = "identity", alpha = 0.55) +
    ggplot2::geom_vline(data = md, ggplot2::aes(xintercept = .data$x, colour = .data$state), linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = sprintf("State relative to %s", if (length(ref.yr) > 1)
      paste0(min(ref.yr), "-", max(ref.yr)) else ref.yr), y = "Density",
      subtitle = sprintf("P(state < %s in %s) = %.3f", extinction, py, mean(st$Projected < extinction))) +
    theme_jagg()
  tab$p.extinct <- c(NA, mean(st$Projected < extinction))
  attr(p, "state") <- tab
  p
}

#' Trajectory with the scaled indices
#'
#' As `jrplot_trjfit()`: the fitted global trajectory with its credible
#' interval and the observations scaled by each index's catchability
#' (relative models; matched by name, so subsets work), or each index's
#' fitted trajectory (census models).
#' @param x a JARA fit
#' @param indices optional subset of index names
#' @export
jagg_trjfit <- function(x, indices = NULL) {
  .need_jara(x, "jagg_trjfit()")
  f <- x$fits; t <- x$trj[x$trj$estimation == "fit", ]
  if (is.null(indices)) indices <- unique(f$name)
  f <- f[f$name %in% indices, ]
  ni <- length(indices); cols <- stats::setNames(jagg_cols(ni), indices)
  if (x$settings$model.type == "relative") {
    g <- t[t$name == "global", ]
    q <- x$pars[paste0("q.", f$name), "median"]
    f$scaled <- f$obs / q
    p <- ggplot2::ggplot(g, ggplot2::aes(x = .data$yr)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "grey80") +
      ggplot2::geom_line(ggplot2::aes(y = .data$mu), linewidth = 0.9) +
      ggplot2::geom_point(data = f, ggplot2::aes(x = .data$year, y = .data$scaled, fill = .data$name),
                          shape = 21, size = 2)
  } else {
    ti <- t[t$name %in% indices, ]
    p <- ggplot2::ggplot(ti, ggplot2::aes(x = .data$yr)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci, fill = .data$name), alpha = 0.25) +
      ggplot2::geom_line(ggplot2::aes(y = .data$mu, colour = .data$name)) +
      ggplot2::geom_point(data = f, ggplot2::aes(x = .data$year, y = .data$obs, fill = .data$name),
                          shape = 21, size = 2) +
      ggplot2::scale_colour_manual(values = cols)
  }
  p + ggplot2::scale_fill_manual(values = cols) + ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Abundance") + theme_jagg()
}

#' Population trajectory with generation lengths
#'
#' As `jrplot_poptrj()`: the global trajectory over the fit and projection
#' years, the last data year, and lines one, two and three generation lengths
#' back from it. With `window = TRUE` the two three-year windows used for the
#' IUCN % change (which end in the last projection year) are shaded.
#' @param x a JARA fit
#' @param plotGL draw the generation-length lines
#' @param window shade the IUCN change windows
#' @export
jagg_poptrj <- function(x, plotGL = TRUE, window = TRUE) {
  .need_jara(x, "jagg_poptrj()")
  g <- x$trj[x$trj$name == "global", ]
  ey <- max(x$yr); GL <- x$settings$GL
  fit <- g[g$estimation == "fit", ]
  prj <- rbind(fit[fit$yr == ey, ], g[g$estimation == "prj", ])
  p <- ggplot2::ggplot(g, ggplot2::aes(x = .data$yr))
  mp <- x$settings$mp.assess
  if (window && !is.null(mp)) {
    w <- data.frame(xmin = x$pyr[mp] - 1.5, xmax = x$pyr[mp] + 1.5)
    p <- p + ggplot2::geom_rect(data = w, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
                                ymin = -Inf, ymax = Inf, inherit.aes = FALSE, fill = "#F9E814", alpha = 0.25)
  }
  p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "#9999994D") +
    ggplot2::geom_ribbon(data = fit, ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = "#B3B3B380") +
    ggplot2::geom_line(data = prj, ggplot2::aes(y = .data$mu), colour = "#DF536B", linewidth = 0.9, linetype = 5) +
    ggplot2::geom_line(data = fit, ggplot2::aes(y = .data$mu), linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = ey, linetype = 2)
  if (plotGL) {
    gl <- data.frame(x = ey - (1:3) * GL, lab = c("1 GL", "2 GL", "3 GL"), col = c("#2297E6", "#61D04F", "#DF536B"))
    gl <- gl[gl$x > min(x$yr) - 1, ]
    p <- p + ggplot2::geom_vline(data = gl, ggplot2::aes(xintercept = .data$x), colour = gl$col,
                                 linetype = 2, linewidth = 0.8) +
      ggplot2::geom_text(data = gl, ggplot2::aes(x = .data$x, y = Inf, label = .data$lab), colour = gl$col,
                         vjust = 1.5, hjust = -0.1, size = 3)
  }
  p + ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Abundance",
                  caption = if (window && !is.null(mp)) "Shaded: windows of the IUCN % change over 3 GL") +
    theme_jagg()
}

#' Time-block effect
#'
#' As `jrplot_timeblock()`: posterior of the change in the mean rate of change
#' after the time block year (`timeblock` in `build_jara()`).
#' @param x a JARA fit with a time block
#' @param type "change" (% change in the annual rate, (exp(effect) - 1) * 100)
#'   or "r" (mean rate after the time block)
#' @param level credible interval level (tails shaded lighter)
#' @export
jagg_timeblock <- function(x, type = c("change", "r"), level = 0.95) {
  .need_jara(x, "jagg_timeblock()")
  type <- match.arg(type)
  tb <- x$timeblock
  if (is.null(tb)) stop("No time block in this fit.", call. = FALSE)
  v <- if (type == "change") (exp(tb$effect) - 1) * 100 else tb$mean.r + tb$effect
  a <- (1 - level) / 2; ci <- stats::quantile(v, c(a, 1 - a))
  dn <- stats::density(v, adjust = 2)
  d <- data.frame(x = dn$x, y = dn$y, inner = dn$x >= ci[1] & dn$x <= ci[2])
  md <- stats::median(v)
  pr <- if (mean(v) > 0) sprintf("P(effect > 0) = %.3f", mean(v > 0)) else sprintf("P(effect < 0) = %.3f", mean(v < 0))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_area(fill = "grey70") +
    ggplot2::geom_area(data = d[d$inner, ], fill = "grey50") +
    ggplot2::geom_vline(xintercept = md) + ggplot2::geom_vline(xintercept = 0, colour = "blue") +
    ggplot2::annotate("text", x = md, y = Inf, vjust = 1.5, size = 3.4, fontface = "bold",
                      label = if (type == "change") sprintf("%s%.1f%%", ifelse(md > 0, "+", ""), md) else
                        sprintf("%s%.3f", ifelse(md > 0, "+", ""), md)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = sprintf("Time-block effect (%s)", tb$year[1]), y = "Density", subtitle = pr) +
    theme_jagg() + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                  axis.ticks.y = ggplot2::element_blank())
}

#' Joint index
#'
#' As `jrplot_jointindex()`: a single index combining all indices, the
#' fitted trend plus the mean residual of each year, with 95% intervals and
#' the trend, both scaled to the mean trend. Intervals are quantiles (JARA
#' uses highest-density intervals).
#' @param x a JARA fit
#' @param show.trend draw the fitted trend
#' @param level interval level
#' @return a ggplot object; the joint index in attribute "index"
#' @export
jagg_jointindex <- function(x, show.trend = TRUE, level = 0.95) {
  .need_jara(x, "jagg_jointindex()")
  n <- length(x$yr)
  trj <- as.matrix(x$pop.posterior)[, 1:n]
  mr <- colMeans(x$residuals, na.rm = TRUE)
  ip <- exp(log(trj) + matrix(mr, nrow(trj), n, byrow = TRUE))
  a <- (1 - level) / 2
  tr <- apply(trj, 2, stats::median); mu <- mean(tr)
  d <- data.frame(year = x$yr, trend = tr / mu, index = apply(ip, 2, stats::median) / mu,
                  lci = apply(ip, 2, stats::quantile, a, na.rm = TRUE) / mu,
                  uci = apply(ip, 2, stats::quantile, 1 - a, na.rm = TRUE) / mu,
                  log.se = apply(log(ip), 2, stats::sd))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year))
  if (show.trend) p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$trend), colour = "#DF536B", linetype = 2)
  p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$index), na.rm = TRUE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), width = 0.3, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(y = .data$index), shape = 21, fill = "white", size = 2, na.rm = TRUE) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Joint index") + theme_jagg()
  attr(p, "index") <- d
  p
}
