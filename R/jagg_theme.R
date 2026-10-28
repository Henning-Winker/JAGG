# JAGG theme, scales and reference-point helpers

#' JAGG ggplot theme
#'
#' `theme_bw()` with the compact legend used in `FLRef::plotAdvice()`.
#' @param base_size base font size
#' @export
theme_jagg <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.key.size   = ggplot2::unit(0.3, "cm"),
      legend.key.height = ggplot2::unit(0.5, "cm"),
      legend.key.width  = ggplot2::unit(0.3, "cm"),
      legend.title      = ggplot2::element_blank(),
      legend.text       = ggplot2::element_text(size = 7),
      strip.background  = ggplot2::element_rect(fill = "grey92"),
      panel.spacing.x   = ggplot2::unit(1, "lines"),
      panel.grid.minor  = ggplot2::element_blank())
}

#' Integer year axis
#' @param n approximate number of breaks
#' @param ... passed to `scale_x_continuous()`
#' @export
scale_x_year <- function(n = 5, ...) {
  ggplot2::scale_x_continuous(breaks = function(x) unique(round(pretty(x, n = n))), ...)
}

#' Colour of a reference point from its name
#'
#' Limits red, precautionary points orange, B0/triggers/F ranges blue,
#' everything else (MSY targets) dark green, after `FLRef::plotAdvice()`.
#' @param x character vector of reference point names
#' @keywords internal
.ref_colour <- function(x) {
  red    <- c("Blim", "Bcrit", "Flim", "Fext", "Fcrash", "Fcrit")
  orange <- c("Bpa", "Bthr", "Fpa", "Fthr", "Fthresh", "Ftri", "Fp0.5")
  blue   <- c("B0", "K", "SPR0", "Btrigger", "Btri", "Flower", "Fupper",
              "Fmsy.low", "Fmsy.up")
  ifelse(x %in% red, "red", ifelse(x %in% orange, "darkorange",
         ifelse(x %in% blue, "blue", "darkgreen")))
}

#' Reference-point lines with labels
#'
#' @param refs data frame with columns `quantity` (facet), `name`, `value`
#' @param years range of years, used to stagger the labels
#' @param label.size label text size
#' @keywords internal
geom_refpts <- function(refs, years, label.size = 2.5) {
  if (is.null(refs) || !nrow(refs)) return(NULL)
  refs$colour <- .ref_colour(refs$name)
  xy <- stats::quantile(seq(min(years), max(years)),
                        c(0.2, 0.45, 0.75, 0.6, 0.3, 0.5, 0.1), names = FALSE)
  refs$x <- stats::ave(seq_len(nrow(refs)), refs$quantity,
                       FUN = function(i) xy[(seq_along(i) - 1) %% length(xy) + 1])
  list(
    ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = .data$value),
                        colour = refs$colour, linetype = 2, linewidth = 0.5),
    ggplot2::geom_label(data = refs, ggplot2::aes(x = .data$x, y = .data$value,
                        label = .data$name), colour = refs$colour, size = label.size,
                        fill = "white", label.size = 0, label.padding = ggplot2::unit(0.1, "lines"),
                        inherit.aes = FALSE))
}
