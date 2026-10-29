# Plots shared by JABBA and JARA: input indices, posterior predictive checks,
# QQ plots of standardised residuals, log-scale fits

# observed indices with their standard errors, original index order
.index_data <- function(x) {
  m <- jagg_model(x)
  if (m == "jabba") {
    I <- x$inputseries$cpue; se <- sqrt(x$settings$SE2)[seq_len(nrow(I)), , drop = FALSE]
  } else if (m == "jara") {
    I <- x$inputseries$I
    s <- sqrt(x$settings$SE2)[seq_len(nrow(I)), , drop = FALSE]
    se <- s; se[, x$settings$qs] <- s            # SE2 columns are in qs order
  } else stop("jagg_indices() needs a JABBA or JARA fit.", call. = FALSE)
  nm <- names(I)[-1]
  data.frame(index = rep(nm, each = nrow(I)), year = rep(I[[1]], length(nm)),
             obs = unlist(I[-1], use.names = FALSE), se = as.vector(se))
}

#' Abundance indices with 95% confidence intervals
#'
#' As `jbplot_indices()` / `jrplot_indices()`: all indices in one panel with
#' lognormal 95% intervals exp(log(I) +/- 1.96 se), using the input se plus
#' any fixed observation error. Indices are offset slightly along the year
#' axis (deterministically, instead of the random jitter of the base plots).
#' @param x a JABBA or JARA fit
#' @param normalise divide each index by its mean
#' @param facet one panel per index instead of all in one
#' @export
jagg_indices <- function(x, normalise = TRUE, facet = FALSE) {
  d <- .index_data(x)
  d <- d[is.finite(d$obs) & d$obs > 0, ]
  d$lo <- d$obs * exp(-1.96 * d$se); d$hi <- d$obs * exp(1.96 * d$se)
  if (normalise) {
    mu <- stats::ave(d$obs, d$index, FUN = mean)
    d$obs <- d$obs / mu; d$lo <- d$lo / mu; d$hi <- d$hi / mu
  }
  idx <- unique(d$index); ni <- length(idx)
  d$index <- factor(d$index, idx)
  dodge <- if (facet || ni == 1) 0 else 0.4
  pd <- ggplot2::position_dodge(width = dodge)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year, y = .data$obs, group = .data$index)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), width = 0.5,
                           colour = "grey30", position = pd) +
    ggplot2::geom_line(colour = "grey40", alpha = 0.7, position = pd) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$index), shape = 21, size = 2.2, position = pd) +
    ggplot2::scale_fill_manual(values = jagg_cols(ni)) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = if (normalise) "Normalised index" else "Abundance index") + theme_jagg()
  if (facet) p <- p + ggplot2::facet_wrap(~index, scales = "free_y") + ggplot2::theme(legend.position = "none")
  p
}

#' Fits on log scale
#'
#' As `jbplot_logfits()` / `jrplot_logfits()`: log observations with 95%
#' observation-error bars and the fitted index.
#' @param x a JABBA or JARA fit
#' @param colour colour of the fit
#' @export
jagg_logfits <- function(x, colour = "#2297E6") {
  d <- jagg_tidy(x, "fits")
  d <- d[(d$peel == 0 | is.na(d$peel)) & is.finite(d$obs) & d$obs > 0, ]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$year)) +
    ggplot2::geom_line(ggplot2::aes(y = log(.data$hat)), colour = colour, linewidth = 0.8)
  if (any(is.finite(d$obs_se)))
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = log(.data$obs) - 1.96 * .data$obs_se,
                                                 ymax = log(.data$obs) + 1.96 * .data$obs_se), width = 0.4)
  p + ggplot2::geom_point(ggplot2::aes(y = log(.data$obs)), shape = 21, fill = "white", size = 1.8) +
    ggplot2::facet_wrap(~index, scales = "free_y") + scale_x_year() +
    ggplot2::labs(x = "Year", y = "log(index)") + theme_jagg()
}

#' Posterior predictive checks
#'
#' As `jbplot_PPC()` / `jrplot_PPC()`: realised against predicted discrepancy
#' for each draw, one panel per index, with the Bayesian p-value
#' P(D_pred > D_obs). Values near 0 or 1 indicate misfit. Requires
#' `do.ppc = TRUE` in `fit_jabba()` / `fit_jara()`.
#' @param x a JABBA or JARA fit
#' @param joint pool all indices into one panel ("Combined")
#' @param npoints draws to plot per panel
#' @return a ggplot object; the p-values in attribute "p"
#' @export
jagg_ppc <- function(x, joint = FALSE, npoints = 1000) {
  m <- jagg_model(x)
  if (!m %in% c("jabba", "jara")) stop("jagg_ppc() needs a JABBA or JARA fit.", call. = FALSE)
  d <- x$PPC
  if (!is.data.frame(d)) stop("No PPC in this fit: refit with do.ppc = TRUE.", call. = FALSE)
  if (joint) d$name <- "Combined"
  bp <- do.call(rbind, lapply(split(d, d$name), function(di)
    data.frame(index = di$name[1], p = mean(di$Dy > di$Dx))))
  bp <- rbind(bp, data.frame(index = "Combined", p = mean(d$Dy > d$Dx)))[if (joint) 1 else seq_len(nrow(bp) + 1), ]
  bp <- bp[!duplicated(bp$index), ]
  set.seed(1)
  th <- do.call(rbind, lapply(split(d, d$name), function(di) di[sample(nrow(di), min(npoints, nrow(di))), ]))
  lab <- stats::setNames(sprintf("%s\np = %.3f", bp$index, bp$p), bp$index)
  th$name <- factor(th$name, unique(d$name))
  p <- ggplot2::ggplot(th, ggplot2::aes(x = .data$Dx, y = .data$Dy)) +
    ggplot2::geom_blank(ggplot2::aes(x = .data$Dy, y = .data$Dx)) +   # square panels
    ggplot2::geom_point(shape = 21, colour = "#666666B3", fill = "#B3B3B380", size = 1.2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linewidth = 0.7) +
    ggplot2::facet_wrap(~name, scales = "free", labeller = ggplot2::as_labeller(lab)) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::labs(x = expression("Realised discrepancy " * D[obs]),
                  y = expression("Predicted discrepancy " * D[pred]),
                  subtitle = if (!joint) sprintf("Combined p = %.3f", bp$p[bp$index == "Combined"])) +
    theme_jagg()
  attr(p, "p") <- bp
  p
}

#' QQ plot of standardised residuals
#'
#' After `flicc::plot_resid_qq()`: standardised index residuals against
#' standard-normal quantiles, one panel per index. JABBA: the stored
#' standardised residuals; JARA: residual / posterior observation SD.
#' @param x a JABBA or JARA fit
#' @export
jagg_qq <- function(x) {
  m <- jagg_model(x)
  if (m == "jabba") {
    r <- x$std.residuals
    d <- data.frame(index = rep(rownames(r), ncol(r)), z = as.vector(r))
  } else if (m == "jara") {
    f <- x$fits
    d <- data.frame(index = f$name, z = f$residual / f$obs.err)
  } else stop("jagg_qq() needs a JABBA or JARA fit.", call. = FALSE)
  d <- d[is.finite(d$z), ]
  d$index <- factor(d$index, unique(d$index))
  ggplot2::ggplot(d, ggplot2::aes(sample = .data$z)) +
    ggplot2::stat_qq(size = 1.2, alpha = 0.7) +
    ggplot2::stat_qq_line(colour = "#D6604D") +
    ggplot2::facet_wrap(~index) +
    ggplot2::labs(x = "Theoretical N(0,1) quantiles", y = "Standardised residual quantiles") +
    theme_jagg()
}
