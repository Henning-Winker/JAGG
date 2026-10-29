#' Projections of forecast scenarios
#'
#' Plots `fw_jabba()` output: the fitted period (median and credible band,
#' black/grey) followed by the projections of each scenario (median and band,
#' coloured with [jagg_cols()]), with reference lines at 1 for B/Bmsy and
#' F/Fmsy. Works for a single fit or an ensemble (`fw_jabba(list(...))`).
#'
#' @param x output of `fw_jabba()` (or a JABBA fit's `kbtrj` with projections)
#' @inheritParams jagg_advice
#' @param type quantities: any of "B/Bmsy", "F/Fmsy", "B/B0", "B", "F", "Catch"
#' @param plotCIs draw bands for the projection scenarios (the fitted period
#'   always has its band)
#' @param quantiles probabilities of the bands
#' @param cols scenario colours (default [jagg_cols()])
#' @param alpha transparency of the scenario bands
#' @param nyears number of fitted years to show before the projections
#'   (NULL for all)
#' @return a ggplot object
#' @export
jagg_prj <- function(x, type = c("B/Bmsy", "F/Fmsy", "Catch"), plotCIs = TRUE, quantiles = c(0.05, 0.95),
                     cols = NULL, alpha = 0.12, fmax = 5, ncol = NULL, nyears = 20) {
  if (jagg_model(x) == "jabba") x <- x$kbtrj
  if (jagg_model(x) != "jabba_kb") stop("x must be fw_jabba() output.", call. = FALSE)
  it <- .traj_iters(list(model = "jabba", fits = NULL, kb = x), type, 0.5, "bmsy", fmax)
  q <- .qbands(it, c(quantiles[1], 0.5, 0.5, 0.5, quantiles[2]))
  fit <- q[q$type == "fit", ]
  # the full fit is stored under the reference run ("Ref", or "joint" for an
  # ensemble); scenario runs only repeat the last fitted year
  ny <- tapply(fit$year, fit$run, function(y) length(unique(y)))
  fit <- fit[fit$run == names(ny)[which.max(ny)], ]
  prj <- q[q$type == "prj", ]
  ly <- max(fit$year)
  if (!is.null(nyears)) fit <- fit[fit$year > ly - nyears, ]
  # join each scenario to the last fitted year
  last <- fit[fit$year == ly, ]
  prj <- rbind(prj, do.call(rbind, lapply(unique(prj$run), function(r) transform(last, run = r))))
  runs <- unique(as.character(x$run[x$type == "prj"]))
  prj$run <- factor(prj$run, levels = runs)
  fit$quantity <- factor(fit$quantity, levels = type)
  prj$quantity <- factor(prj$quantity, levels = type)
  if (is.null(cols)) cols <- jagg_cols(length(runs))
  cols <- stats::setNames(rep(cols, length.out = length(runs)), runs)

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = .data$year, y = .data$med)) +
    ggplot2::geom_ribbon(data = fit, ggplot2::aes(ymin = .data$lo2, ymax = .data$hi2),
                         fill = "grey75") +
    ggplot2::geom_line(data = fit, linewidth = 0.7)
  if (plotCIs) p <- p + ggplot2::geom_ribbon(data = prj, ggplot2::aes(ymin = .data$lo2, ymax = .data$hi2,
                                             fill = .data$run), alpha = alpha)
  p <- p + ggplot2::geom_line(data = prj, ggplot2::aes(colour = .data$run), linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = ly, linetype = 3, colour = "grey40") +
    ggplot2::scale_colour_manual(values = cols) +
    (if (plotCIs) ggplot2::scale_fill_manual(values = cols)) +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = ncol) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = NULL) + theme_jagg()
  ref <- intersect(type, c("B/Bmsy", "F/Fmsy"))
  if (length(ref)) p <- p + ggplot2::geom_hline(
    data = data.frame(quantity = factor(ref, levels = type)), ggplot2::aes(yintercept = 1),
    linetype = 2, colour = "darkgreen")
  p
}
