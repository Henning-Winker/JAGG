#' Fits to abundance indices
#'
#' Observations (points), fitted median with 95% credible interval and,
#' optionally, the 95% posterior predictive interval, one facet per index.
#' @param x a JABBA or JARA fit
#' @param ppd draw the posterior predictive interval
#' @param log plot on log scale (as `jbplot_logfits()`)
#' @param colour colour of the fit
#' @export
jagg_fits <- function(x, ppd = TRUE, log = FALSE, colour = "dodgerblue") {
  d <- jagg_tidy(x, "fits")
  d <- d[d$peel == 0 | is.na(d$peel), ]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year))
  if (ppd && any(is.finite(d$ppd_lci)))
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ppd_lci, ymax = .data$ppd_uci),
                                  fill = "grey85")
  p <- p +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lci, ymax = .data$uci), fill = colour, alpha = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = .data$hat), colour = colour, linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = .data$obs), shape = 21, fill = "white", size = 1.8,
                        na.rm = TRUE) +
    ggplot2::facet_wrap(~index, scales = "free_y") + scale_x_year() +
    ggplot2::labs(x = "Year", y = "Index") + theme_jagg()
  if (log) p <- p + ggplot2::scale_y_log10() else p <- p + ggplot2::expand_limits(y = 0)
  p
}

#' Index residuals by year
#'
#' As `jbplot_residuals()`: residuals of all indices coloured by index, a
#' boxplot per year, a loess smoother through all residuals, and the RMSE.
#' @param x a JABBA or JARA fit
#' @export
jagg_residuals <- function(x) {
  d <- jagg_tidy(x, "fits")
  d <- d[(d$peel == 0 | is.na(d$peel)) & is.finite(d$residual), ]
  rmse <- sqrt(mean(d$residual^2))
  ni <- length(unique(d$index))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$year, y = .data$residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_boxplot(ggplot2::aes(group = .data$year), fill = "grey90", colour = "grey60",
                          outlier.shape = NA, width = 0.7) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$index), alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$index), size = 1.6) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "black",
                         linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = jagg_cols(ni)) + scale_x_year() +
    ggplot2::labs(x = "Year", y = "log residual",
                  subtitle = sprintf("RMSE = %.1f%%", 100 * rmse)) + theme_jagg()
}

#' Runs test of index residuals
#'
#' After `ggplotFL::plotRunstest()` and JABBA `jbplot_runstest()`: residuals as
#' segments from zero, three-sigma limits from the moving range (green when
#' the runs test passes, p >= 0.05, red when it fails), outliers in red.
#' @param x a JABBA or JARA fit
#' @param mixing alternative of the runs test, "less" by default
#' @return a ggplot object; the runs test table is in attribute "runs"
#' @export
jagg_runstest <- function(x, mixing = "less") {
  d <- jagg_tidy(x, "fits")
  d <- d[(d$peel == 0 | is.na(d$peel)) & is.finite(d$residual), ]
  rt <- do.call(rbind, lapply(split(d, d$index), function(di) {
    r <- jagg_runs(di$residual[order(di$year)], mixing)
    data.frame(index = di$index[1], p.runs = r$p.runs, lcl = r$sig3lim[1], ucl = r$sig3lim[2],
               test = ifelse(r$p.runs < 0.05, "Failed", "Passed"),
               xmin = min(di$year) - 0.5, xmax = max(di$year) + 0.5)
  }))
  d <- merge(d, rt[, c("index", "lcl", "ucl")], by = "index")
  d$outlier <- d$residual < d$lcl | d$residual > d$ucl
  lab <- stats::setNames(sprintf("%s\np = %.3f", rt$index, rt$p.runs), rt$index)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year, y = .data$residual)) +
    ggplot2::geom_rect(data = rt, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                       ymin = .data$lcl, ymax = .data$ucl, fill = .data$test),
                       inherit.aes = FALSE, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(Passed = "#cbe368", Failed = "#ef8575"), drop = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$year, yend = 0)) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$outlier), size = 1.1) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = "#d50102"), guide = "none") +
    ggplot2::facet_wrap(~index, scales = "free_x", labeller = ggplot2::as_labeller(lab)) +
    scale_x_year(n = 4) + ggplot2::labs(x = "Year", y = "log residual") + theme_jagg() +
    ggplot2::theme(legend.position = "bottom")
  attr(p, "runs") <- rt[, c("index", "p.runs", "test", "lcl", "ucl")]
  p
}

#' Kobe plot
#'
#' Median B/Bmsy and F/Fmsy trajectory over the posterior of the final year,
#' with the percentage of the posterior in each quadrant.
#' @param x a JABBA fit
#' @param year year of the posterior cloud, default the last assessment year
#' @param npoints number of posterior draws to plot
#' @export
jagg_kobe <- function(x, year = NULL, npoints = 2000) {
  if (jagg_model(x) != "jabba") stop("jagg_kobe() needs a JABBA fit.", call. = FALSE)
  kb <- x$kbtrj[x$kbtrj$type == "fit", ]
  if (is.null(year)) year <- max(kb$year)
  trj <- data.frame(year = sort(unique(kb$year)),
                    stock = tapply(kb$stock, kb$year, stats::median),
                    harvest = tapply(kb$harvest, kb$year, stats::median))
  pt <- kb[kb$year == year, ]
  zone <- ifelse(pt$stock >= 1 & pt$harvest <= 1, "green",
          ifelse(pt$stock < 1 & pt$harvest <= 1, "yellow",
          ifelse(pt$stock >= 1 & pt$harvest > 1, "orange", "red")))
  pct <- round(100 * prop.table(table(factor(zone, c("green", "yellow", "orange", "red")))), 1)
  set.seed(1); pt <- pt[sample(nrow(pt), min(npoints, nrow(pt))), ]
  xmax <- max(2, trj$stock, stats::quantile(pt$stock, 0.99)) * 1.05
  ymax <- max(2, trj$harvest, stats::quantile(pt$harvest, 0.99)) * 1.05
  qd <- data.frame(xmin = c(1, 0, 1, 0), xmax = c(Inf, 1, Inf, 1), ymin = c(0, 0, 1, 1),
                   ymax = c(1, 1, Inf, Inf), zone = c("green", "yellow", "orange", "red"))
  qd$lab <- sprintf("%.1f%%", pct[qd$zone])
  fills <- c(green = "#2ca25f", yellow = "#fee391", orange = "#fd8d3c", red = "#de2d26")
  ends <- trj[c(1, nrow(trj)), ]
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = qd, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                       ymin = .data$ymin, ymax = .data$ymax, fill = .data$zone), alpha = 0.85) +
    ggplot2::scale_fill_manual(values = fills, breaks = qd$zone, labels = qd$lab,
                               name = paste("Posterior", year)) +
    ggplot2::geom_point(data = pt, ggplot2::aes(.data$stock, .data$harvest),
                        colour = "grey30", alpha = 0.15, size = 0.6) +
    ggplot2::geom_path(data = trj, ggplot2::aes(.data$stock, .data$harvest), linewidth = 0.6) +
    ggplot2::geom_point(data = trj, ggplot2::aes(.data$stock, .data$harvest), shape = 21,
                        fill = "white", size = 1.4) +
    ggplot2::geom_point(data = ends, ggplot2::aes(.data$stock, .data$harvest),
                        shape = c(22, 24), fill = "white", size = 3) +
    ggplot2::geom_text(data = ends, ggplot2::aes(.data$stock, .data$harvest, label = .data$year),
                       vjust = -1.1, size = 3) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
    ggplot2::labs(x = expression(B/B[MSY]), y = expression(F/F[MSY])) +
    theme_jagg() + ggplot2::theme(legend.title = ggplot2::element_text(size = 8),
                                  legend.position = "right")
}
