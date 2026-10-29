# Plots for JARA mixed-trend fits of many stocks (model.type = "mixed.trends"),
# after Winker et al. (2026), Fish and Fisheries 27(3), doi:10.1111/faf.70066,
# and the code in
# github.com/Henning-Winker/ICES-Bmsy-FaF

#' Indicator colours used in Winker et al. (2026)
#' @return named vector for F/Fmsy, B/Bmsy, B/Btrigger and ABImsy
#' @export
jagg_pal_indicator <- function() c(`F/Fmsy` = "#000000", `B/Bmsy` = "#006400",
                                   `B/Btrigger` = "#FF8C00", ABImsy = "#0000FF")

#' Kobe and FAO status colours used in Winker et al. (2026)
#' @return named colour vector
#' @export
jagg_pal_kobe <- function() c(Sustainable = "#00FF00", Rebuilding = "#FFFF00",
                              Overfishing = "#FFA500", Overfished = "#FF0000")

#' @rdname jagg_pal_kobe
#' @export
jagg_pal_fao <- function() c(Underfished = "#0000FF", Sustainable = "#00FF00", Overfished = "#FF0000")

.need_mixed <- function(f) {
  if (jagg_model(f) != "jara" || !isTRUE(f$settings$mixed.trends))
    stop("Needs a JARA fit with model.type = \"mixed.trends\".", call. = FALSE)
}

# draws matrix (iterations x years) of the joint ratio ("mu") or of the
# proportion of stocks above the threshold ("pr")
.mixed_draws <- function(f, type, threshold.id) {
  if (type == "mu") return(as.matrix(f$pop.posterior))
  pp <- f$prop.posterior
  if (is.null(pp)) stop("No prop.posterior in this fit.", call. = FALSE)
  if (!is.data.frame(pp)) pp <- pp[[threshold.id]]
  as.matrix(pp)
}

# threshold value used for prop.posterior[[threshold.id]]
.mixed_threshold <- function(f, threshold.id) {
  pp <- f$prop.posterior
  if (is.data.frame(pp)) return(1)
  as.numeric(sub("^t", "", names(pp)[threshold.id]))
}

#' Tidy joint trend of a mixed-trend JARA fit
#'
#' Quantiles by year of the joint ratio across stocks (`type = "mu"`, from
#' `pop.posterior`, the bootstrapped geometric mean by default) or of the
#' proportion of stocks above a threshold (`type = "pr"`, from
#' `prop.posterior`), with the observed geometric mean or observed proportion.
#' Same quantities as `JARA::mixed.trend()`, with quantiles matched to years
#' by name and labelled by their true probabilities.
#'
#' @param f a JARA mixed-trend fit
#' @param type "mu" (joint ratio) or "pr" (proportion of stocks)
#' @param threshold.id which threshold of `fit_jara(threshold = ...)` for "pr"
#' @param below for "pr": report the proportion below the threshold
#'   (e.g. F < Fmsy) instead of above
#' @param prj include projection years
#' @param run label
#' @return data frame: run, year, estimation, obs, n, q025, q25, q50, q75, q975
#' @export
jagg_tidy_mixed <- function(f, type = c("mu", "pr"), threshold.id = 1, below = FALSE, prj = FALSE,
                            run = f$scenario) {
  .need_mixed(f)
  type <- match.arg(type)
  dr <- .mixed_draws(f, type, threshold.id)
  if (type == "pr" && below) dr <- 1 - dr
  yrs <- as.numeric(colnames(dr))
  q <- t(apply(dr, 2, stats::quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE))
  ob <- f$fits[is.finite(f$fits$obs) & f$fits$obs > 0, ]
  thr <- if (type == "pr") .mixed_threshold(f, threshold.id) else NA
  o <- if (type == "mu") {
    sc <- if (is.null(f$settings$mixed.scale)) "geomean" else f$settings$mixed.scale
    fn <- switch(sc, geomean = function(v) exp(mean(log(v))), mean = mean, median = stats::median)
    tapply(ob$obs, ob$year, fn)
  } else if (below) tapply(ob$obs <= thr, ob$year, mean) else tapply(ob$obs >= thr, ob$year, mean)
  n <- tapply(ob$obs, ob$year, length)
  d <- data.frame(run = run, year = yrs, estimation = ifelse(yrs %in% f$yr, "fit", "prj"),
                  obs = as.numeric(o[as.character(yrs)]), n = as.numeric(n[as.character(yrs)]),
                  q025 = q[, 1], q25 = q[, 2], q50 = q[, 3], q75 = q[, 4], q975 = q[, 5], row.names = NULL)
  if (!prj) d <- d[d$estimation == "fit", ]
  attr(d, "threshold") <- thr
  d
}

# colours for a set of runs: the indicator palette where the names match
.run_cols <- function(runs) {
  pal <- jagg_pal_indicator()
  key <- gsub("[^A-Za-z]", "", tolower(runs)); pk <- gsub("[^A-Za-z]", "", tolower(names(pal)))
  hit <- match(key, pk)
  if (all(!is.na(hit))) stats::setNames(unname(pal[hit]), runs) else stats::setNames(jagg_cols(length(runs)), runs)
}

#' Joint trends across stocks
#'
#' After Figure 2 of Winker et al. (2026): the joint trend of one or more
#' indicators across stocks, either as the ratio (`type = "mu"`, e.g. the
#' geometric mean B/Bmsy) or as the proportion of stocks above a threshold
#' (`type = "pr"`), with 50% and 95% credible bands, the median and the
#' observed values.
#'
#' @param fits a named list of JARA mixed-trend fits (one per indicator), or a
#'   single fit
#' @param type "mu" (ratio) or "pr" (proportion of stocks)
#' @param threshold.id threshold of `fit_jara(threshold = ...)` for "pr"
#' @param below names of the fits for which "pr" is the proportion below the
#'   threshold (default: names containing "F", e.g. "F/Fmsy")
#' @param prj include projection years
#' @param obs draw the observed values
#' @param colours named colours (default [jagg_pal_indicator()] when the
#'   names match, else [jagg_cols()])
#' @param refline reference line(s) (default 1 for "mu", 0.5 for "pr")
#' @param vlines optional named years to mark, e.g. scenario start years
#' @export
jagg_joint <- function(fits, type = c("mu", "pr"), threshold.id = 1, below = NULL, prj = FALSE,
                       obs = TRUE, colours = NULL, refline = NULL, vlines = NULL) {
  type <- match.arg(type)
  if (jagg_model(fits) == "jara") fits <- stats::setNames(list(fits), fits$scenario)
  if (is.null(names(fits))) names(fits) <- vapply(fits, function(f) f$scenario, "")
  if (is.null(below)) below <- names(fits)[grepl("^F", names(fits))]
  d <- do.call(rbind, lapply(names(fits), function(nm) {
    b <- nm %in% below
    td <- jagg_tidy_mixed(fits[[nm]], type, threshold.id, below = b, prj = prj, run = nm)
    if (type == "pr") {
      thr <- attr(td, "threshold")
      op <- if (b) "<" else ">"
      td$label <- if (grepl("/", nm)) sprintf("%s %s %s", sub("/.*", "", nm), op,
                          if (thr == 1) sub(".*/", "", nm) else paste0(thr, " ", sub(".*/", "", nm))) else
        sprintf("%s %s %s", nm, op, thr)
    } else td$label <- nm
    td
  }))
  labs <- unique(d$label)
  d$label <- factor(d$label, labs)
  cols <- if (is.null(colours)) .run_cols(names(fits)) else colours
  cols <- stats::setNames(unname(cols[names(fits)]), labs)
  if (is.null(refline)) refline <- if (type == "mu") 1 else 0.5
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year, fill = .data$label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q025, ymax = .data$q975), alpha = 0.15) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q25, ymax = .data$q75), alpha = 0.2) +
    ggplot2::geom_hline(yintercept = refline, linetype = 2, colour = "grey40") +
    ggplot2::geom_line(ggplot2::aes(y = .data$q50, colour = .data$label), linewidth = 0.6)
  if (obs) p <- p + ggplot2::geom_point(ggplot2::aes(y = .data$obs, colour = .data$label), shape = 21,
                                        fill = "white", size = 1.5, na.rm = TRUE)
  if (!is.null(vlines)) p <- p + ggplot2::geom_vline(xintercept = vlines, linetype = 2,
                                                     colour = cols[seq_along(vlines)])
  p + ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::expand_limits(y = if (type == "mu") c(0, 1) else c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = function(x) seq(5 * floor(x[1] / 5), x[2], 5)) +
    ggplot2::labs(x = "Year", y = if (type == "mu") "Ratio" else "Proportion of stocks") +
    theme_jagg() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Joint trends by group
#'
#' After Figure 3 of Winker et al. (2026): joint trends with their 50% (and
#' optionally 95%) credible bands for each group (e.g. taxonomic order or
#' ecoregion, columns) and indicator (rows), with the number of stocks per
#' group in the column labels.
#'
#' @param fits a nested list `fits[[group]][[indicator]]` of JARA mixed-trend fits
#' @param type "mu" or "pr" (see [jagg_joint()])
#' @param ci90 also draw the 95% band
#' @param colours named indicator colours
#' @export
jagg_joint_grid <- function(fits, type = c("mu", "pr"), ci90 = FALSE, colours = NULL) {
  type <- match.arg(type)
  d <- do.call(rbind, lapply(names(fits), function(g) do.call(rbind, lapply(names(fits[[g]]), function(i) {
    f <- fits[[g]][[i]]
    td <- jagg_tidy_mixed(f, type, below = type == "pr" && grepl("^F", i), run = i)
    td$group <- sprintf("%s (n = %d)", g, length(f$indices))
    td
  }))))
  inds <- unique(d$run)
  d$run <- factor(d$run, inds); d$group <- factor(d$group, unique(d$group))
  cols <- if (is.null(colours)) .run_cols(inds) else colours
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year, fill = .data$run))
  if (ci90) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q025, ymax = .data$q975), alpha = 0.15)
  p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$q25, ymax = .data$q75), alpha = 0.35) +
    ggplot2::geom_hline(yintercept = if (type == "mu") 1 else 0.5, linetype = 2, colour = "grey40") +
    ggplot2::geom_line(ggplot2::aes(y = .data$q50, colour = .data$run), linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::facet_grid(run ~ group, scales = "free_y") +
    ggplot2::expand_limits(y = 0) + scale_x_year(n = 4) +
    ggplot2::labs(x = "Year", y = if (type == "mu") "Ratio" else "Proportion of stocks") +
    theme_jagg() + ggplot2::theme(legend.position = "none",
                                  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# stock-level medians (fit years) of a mixed-trend fit, masked before
# backcast.na where a stock has no observation yet (as JARA::dfidx)
.stock_medians <- function(f, backcast.na = NULL) {
  t <- f$trj[f$trj$name != "global" & f$trj$estimation == "fit", c("name", "yr", "mu")]
  if (is.null(backcast.na)) backcast.na <- f$settings$backcast.na
  if (!is.null(backcast.na)) {
    first <- tapply(f$fits$year[is.finite(f$fits$obs)], f$fits$name[is.finite(f$fits$obs)], min)
    fo <- first[as.character(t$name)]
    t$mu[t$yr < backcast.na & (is.na(fo) | t$yr < fo)] <- NA
  }
  t
}

#' Stock status classification over time
#'
#' After Figure 4a,b of Winker et al. (2026): the proportion of stocks in
#' each status class by year, from the stock-level medians of mixed-trend
#' fits. `scheme = "kobe"` needs a B/Bmsy (or B/Btrigger) fit and an F/Fmsy
#' fit and uses the Kobe quadrants (Sustainable, Rebuilding, Overfishing,
#' Overfished); `scheme = "fao"` needs only the biomass fit (Underfished
#' B/Bmsy > 1.2, Sustainable 0.8-1.2, Overfished < 0.8). Stocks are matched
#' by name.
#'
#' @param fit_b JARA mixed-trend fit of B/Bmsy (or another biomass ratio)
#' @param fit_f JARA mixed-trend fit of F/Fmsy (Kobe scheme)
#' @param scheme "kobe" or "fao"
#' @param fao_cuts biomass cut-offs for the FAO scheme
#' @param backcast.na mask stock values before this year where unobserved
#'   (default: the fit's `backcast.na`)
#' @param alpha fill transparency
#' @export
jagg_status_area <- function(fit_b, fit_f = NULL, scheme = c("kobe", "fao"), fao_cuts = c(0.8, 1.2),
                             backcast.na = NULL, alpha = 0.6) {
  scheme <- match.arg(scheme)
  b <- .stock_medians(fit_b, backcast.na)
  if (scheme == "kobe") {
    if (is.null(fit_f)) stop("The Kobe scheme needs an F/Fmsy fit.", call. = FALSE)
    fq <- .stock_medians(fit_f, backcast.na)
    d <- merge(b, fq, by = c("name", "yr"), suffixes = c(".b", ".f"))
    d <- d[is.finite(d$mu.b) & is.finite(d$mu.f), ]
    d$class <- ifelse(d$mu.b >= 1, ifelse(d$mu.f < 1, "Sustainable", "Overfishing"),
                      ifelse(d$mu.f < 1, "Rebuilding", "Overfished"))
    pal <- jagg_pal_kobe()
  } else {
    d <- b[is.finite(b$mu), ]
    d$class <- ifelse(d$mu > fao_cuts[2], "Underfished", ifelse(d$mu >= fao_cuts[1], "Sustainable", "Overfished"))
    pal <- jagg_pal_fao()
  }
  tab <- as.data.frame(prop.table(table(yr = d$yr, class = factor(d$class, names(pal))), 1))
  tab$yr <- as.numeric(as.character(tab$yr))
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$yr, y = .data$Freq, fill = .data$class)) +
    ggplot2::geom_area(alpha = alpha, colour = "grey30", linewidth = 0.2) +
    ggplot2::geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = 2, colour = "grey", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = pal, name = "Status:") +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = function(x) seq(5 * ceiling(x[1] / 5), x[2], 5)) +
    ggplot2::labs(x = "Year", y = "Proportion of stocks") + theme_jagg() +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 8),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Kobe trajectory of joint trends
#'
#' After Figure 4c of Winker et al. (2026): the path of the joint median
#' B/Bmsy against the joint median F/Fmsy of two mixed-trend fits, coloured by
#' year, over the Kobe quadrants.
#' @param fit_b JARA mixed-trend fit of B/Bmsy (or B/Btrigger)
#' @param fit_f JARA mixed-trend fit of F/Fmsy
#' @param xlab,ylab axis labels
#' @param lim axis limits
#' @export
jagg_kobe_trend <- function(fit_b, fit_f, xlab = expression(B/B[MSY]), ylab = expression(F/F[MSY]),
                            lim = c(2.05, 2)) {
  b <- jagg_tidy_mixed(fit_b); f <- jagg_tidy_mixed(fit_f)
  d <- merge(b[, c("year", "q50")], f[, c("year", "q50")], by = "year", suffixes = c(".b", ".f"))
  qd <- data.frame(xmin = c(1, 0, 1, 0), xmax = c(Inf, 1, Inf, 1), ymin = c(0, 0, 1, 1), ymax = c(1, 1, Inf, Inf),
                   class = c("Sustainable", "Rebuilding", "Overfishing", "Overfished"))
  ends <- d[c(1, nrow(d)), ]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$q50.b, y = .data$q50.f)) +
    ggplot2::geom_rect(data = qd, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
                       ymax = .data$ymax, fill = .data$class), inherit.aes = FALSE, alpha = 0.5) +
    ggplot2::scale_fill_manual(values = jagg_pal_kobe(), guide = "none") +
    ggplot2::geom_hline(yintercept = 1, colour = "grey40") + ggplot2::geom_vline(xintercept = 1, colour = "grey40") +
    ggplot2::geom_path(colour = "grey20") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$year), size = 2.2) +
    ggplot2::geom_text(data = ends, ggplot2::aes(label = .data$year), vjust = -1, size = 3) +
    ggplot2::scale_colour_viridis_c(direction = -1, name = "Year:") +
    ggplot2::coord_cartesian(xlim = c(0, lim[1]), ylim = c(0, lim[2]), expand = FALSE) +
    ggplot2::labs(x = xlab, y = ylab) + theme_jagg() +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 8))
}

#' Stock trajectories with the final-year distribution
#'
#' After the appendix figure of Winker et al. (2026): the fitted trajectory of
#' every stock of a mixed-trend fit and, with `marginal = TRUE` (needs the
#' patchwork package), a histogram of the final-year medians.
#' @param f a JARA mixed-trend fit
#' @param ymax upper y limit (values above are cut off by the axis)
#' @param refline reference line
#' @param fill histogram colour
#' @param marginal add the final-year histogram
#' @export
jagg_spaghetti <- function(f, ymax = NULL, refline = 1, fill = "darkgreen", marginal = TRUE) {
  .need_mixed(f)
  t <- f$trj[f$trj$name != "global" & f$trj$estimation == "fit", ]
  if (is.null(ymax)) ymax <- stats::quantile(t$mu, 0.99, na.rm = TRUE)
  ns <- length(unique(t$name))
  p <- ggplot2::ggplot(t, ggplot2::aes(x = .data$yr, y = .data$mu, colour = .data$name)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = refline) +
    ggplot2::scale_colour_manual(values = r4col(ns)) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) + scale_x_year() +
    ggplot2::labs(x = "Year", y = f$scenario) + theme_jagg() + ggplot2::theme(legend.position = "none")
  if (!marginal || !requireNamespace("patchwork", quietly = TRUE)) return(p)
  last <- t[t$yr == max(t$yr), ]
  h <- ggplot2::ggplot(last, ggplot2::aes(x = .data$mu)) +
    ggplot2::geom_histogram(bins = 15, fill = fill, colour = "grey50", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = refline) +
    ggplot2::coord_flip(xlim = c(0, ymax)) +
    ggplot2::labs(x = NULL, y = paste("Stocks in", max(t$yr))) + theme_jagg() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())
  patchwork::wrap_plots(p, h, widths = c(4, 1))
}
