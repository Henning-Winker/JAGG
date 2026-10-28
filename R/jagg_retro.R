# Retrospective analysis and hindcast cross-validation

# Prediction table (index, year, obs, hat, peel) from a hindcast.
# JABBA: hindcast_jabba() list; peel runs keep all observations and predict
#   the peeled years, so the fits are used directly.
# JARA: jara_hindcast() removes the peeled years from the data and labels
#   runs by peel, without a reference run; predictions are the index-level
#   trajectories (index scale) and observations come from the reference fit.
.hc_fits <- function(hc, ref = NULL) {
  m <- jagg_model(hc)
  if (m == "jabba_hc") return(jagg_tidy(hc, "fits")[, c("index", "year", "obs", "hat", "peel")])
  if (m != "jara_hc") stop("x must be a JABBA or JARA hindcast.", call. = FALSE)
  if (is.null(ref)) stop("JARA hindcasts need the reference fit: ref = fit_jara(...).", call. = FALSE)
  rf <- ref$fits
  obs <- rf[, c("name", "year", "obs")]
  t <- hc$trj[hc$trj$name != "global", ]
  pr <- data.frame(index = t$name, year = t$yr, hat = t$mu, peel = as.numeric(as.character(t$level)))
  pr <- merge(pr, obs, by.x = c("index", "year"), by.y = c("name", "year"), all.x = TRUE)
  r0 <- data.frame(index = rf$name, year = rf$year, obs = rf$obs, hat = rf$hat, peel = 0)
  rbind(r0, pr[, names(r0)])
}

# Trajectories (run, peel, quantity, year, median, lci, uci) from a hindcast
.hc_ts <- function(hc, ref = NULL) {
  m <- jagg_model(hc)
  if (m == "jabba_hc") {
    lab <- c(B = "B", F = "F", BBmsy = "B/Bmsy", FFmsy = "F/Fmsy", BB0 = "B/B0",
             procB = "Process deviation", SPt = "Surplus production")
    d <- jagg_tidy(hc, "ts")
    d$quantity <- unname(lab[d$quantity])
    return(d[d$year <= d$end_year, c("peel", "quantity", "year", "median", "lci", "uci")])
  }
  if (is.null(ref)) stop("JARA hindcasts need the reference fit: ref = fit_jara(...).", call. = FALSE)
  t <- hc$trj[hc$trj$name == "global" & hc$trj$estimation == "fit", ]
  r <- ref$trj[ref$trj$name == "global" & ref$trj$estimation == "fit", ]
  rbind(data.frame(peel = 0, quantity = "Abundance", year = r$yr, median = r$mu, lci = r$lci, uci = r$uci),
        data.frame(peel = as.numeric(as.character(t$level)), quantity = "Abundance", year = t$yr,
                   median = t$mu, lci = t$lci, uci = t$uci))
}

#' Mohn's rho
#' @param ts output of `.hc_ts()`
#' @return data frame of quantity and rho
#' @keywords internal
.mohn <- function(ts) {
  do.call(rbind, lapply(split(ts, ts$quantity), function(d) {
    peels <- setdiff(sort(unique(d$peel)), 0)
    rel <- vapply(peels, function(p) {
      ey <- max(d$year[d$peel == p])
      xp <- d$median[d$peel == p & d$year == ey]
      x0 <- d$median[d$peel == 0 & d$year == ey]
      if (length(xp) && length(x0) && x0 != 0) (xp - x0) / x0 else NA_real_
    }, numeric(1))
    data.frame(quantity = d$quantity[1], rho = mean(rel, na.rm = TRUE))
  }))
}

#' Retrospective analysis
#'
#' Reference run with its 95% credible interval, one line per peel ending in
#' its terminal year, and Mohn's rho in each facet label. Colours: reference
#' black, peels [jagg_cols()], as in `jbplot_retro()`.
#'
#' @param hc a JABBA hindcast (`hindcast_jabba()`) or JARA hindcast
#'   (`jara_hindcast()`)
#' @param type quantities to show (JABBA: "B", "F", "B/Bmsy", "F/Fmsy",
#'   "B/B0"; JARA: "Abundance")
#' @param ref JARA only: the reference fit
#' @param ncol number of facet columns
#' @return a ggplot object; Mohn's rho in attribute "rho"
#' @export
jagg_retro <- function(hc, type = NULL, ref = NULL, ncol = NULL) {
  ts <- .hc_ts(hc, ref)
  if (is.null(type)) type <- if (jagg_model(hc) == "jara_hc") "Abundance" else c("B/Bmsy", "F/Fmsy")
  ts <- ts[ts$quantity %in% type, ]
  rho <- .mohn(ts)
  lab <- stats::setNames(sprintf("%s\nMohn's rho = %.3f", rho$quantity, rho$rho), rho$quantity)
  ts$quantity <- factor(ts$quantity, levels = type)
  peels <- sort(unique(ts$peel)); np <- length(peels) - 1
  ty <- max(ts$year[ts$peel == 0])
  plab <- ifelse(peels == 0, "Reference", as.character(ty - peels))
  ts$run <- factor(plab[match(ts$peel, peels)], levels = plab)
  cols <- stats::setNames(jagg_peelcols(np), plab)
  ends <- ts[stats::ave(ts$year, ts$run, ts$quantity, FUN = max) == ts$year & ts$peel > 0, ]
  r0 <- ts[ts$peel == 0, ]
  p <- ggplot2::ggplot(ts, ggplot2::aes(x = .data$year, y = .data$median, colour = .data$run)) +
    ggplot2::geom_ribbon(data = r0, ggplot2::aes(ymin = .data$lci, ymax = .data$uci),
                         fill = "grey80", colour = NA) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(data = ends, size = 2) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = ncol,
                        labeller = ggplot2::as_labeller(lab)) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = NULL) + theme_jagg()
  if (any(type %in% c("B/Bmsy", "F/Fmsy")))
    p <- p + ggplot2::geom_hline(data = data.frame(quantity = factor(intersect(type, c("B/Bmsy", "F/Fmsy")),
                                 levels = type)), ggplot2::aes(yintercept = 1), linetype = 2, colour = "grey40")
  attr(p, "rho") <- rho
  p
}

#' Hindcast cross-validation
#'
#' After `ggplotFL::plotXval()` / `a4adiags::plotXval2()`: observations as a
#' dashed grey line with ringed points, each peel's fitted index up to its
#' one-step-ahead prediction, the observation predicted by each peel coloured
#' like that peel, and MASE (from [jagg_mase()]) in each facet label.
#'
#' @param hc a JABBA or JARA hindcast
#' @param ref JARA only: the reference fit
#' @param MAE.base.adj minimum naive error for MASE.adj (default 0.1)
#' @param adj also show MASE.adj in the labels
#' @param ncol number of facet columns
#' @param nyears number of years before the hindcast horizon to show
#' @return a ggplot object; the MASE table in attribute "mase"
#' @export
jagg_hcxval <- function(hc, ref = NULL, MAE.base.adj = 0.1, adj = FALSE, ncol = 2, nyears = 3) {
  d <- .hc_fits(hc, ref)
  ms <- suppressMessages(jagg_mase(d, MAE.base.adj = MAE.base.adj))
  peels <- sort(unique(d$peel)); np <- max(peels)
  ty <- max(d$year[d$peel == 0])
  y0 <- ty - np - nyears
  keep <- ms$Index[ms$n.eval > 0]
  d <- d[d$index %in% keep & d$year >= y0, ]
  ms <- ms[ms$Index %in% keep, ]

  plab <- ifelse(peels == 0, "Reference", as.character(ty - peels))
  cols <- stats::setNames(jagg_peelcols(np), plab)
  # fitted lines: reference to the terminal year, peel p to its prediction year
  ln <- d[is.finite(d$hat) & d$year <= ty - d$peel + ifelse(d$peel > 0, 1, 0), ]
  ln$run <- factor(plab[match(ln$peel, peels)], levels = plab)
  pr <- ln[ln$peel > 0 & ln$year == ty - ln$peel + 1, ]
  ob <- d[d$peel == 0 & is.finite(d$obs), ]
  ob$pred_by <- ifelse(ob$year > ty - np, ob$year - 1, NA)
  ob$run <- factor(ifelse(is.na(ob$pred_by), NA, as.character(ob$pred_by)), levels = plab)

  lab <- if (adj) sprintf("%s\nMASE = %.2f (adj = %.2f)", ms$Index, ms$MASE, ms$MASE.adj) else
    sprintf("%s\nMASE = %.2f", ms$Index, ms$MASE)
  lab <- stats::setNames(lab, ms$Index)

  p <- ggplot2::ggplot(ln, ggplot2::aes(x = .data$year, y = .data$hat, colour = .data$run)) +
    ggplot2::geom_line(data = ob, ggplot2::aes(y = .data$obs), linetype = 2, colour = "grey60") +
    ggplot2::geom_point(data = ob, ggplot2::aes(y = .data$obs), shape = 21, colour = "black",
                        fill = "white", size = 2.8) +
    ggplot2::geom_point(data = ob[!is.na(ob$run), ], ggplot2::aes(y = .data$obs, fill = .data$run),
                        shape = 21, colour = "black", size = 2.8, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(data = pr, size = 1.8) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::facet_wrap(~index, scales = "free_y", ncol = ncol, labeller = ggplot2::as_labeller(lab)) +
    scale_x_year() + ggplot2::labs(x = "Year", y = "Index") + theme_jagg() +
    ggplot2::theme(legend.position = "bottom")
  attr(p, "mase") <- ms
  p
}
