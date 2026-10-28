# JAGG colour palettes, consistent with r4ss and the JABBA/JARA base plots

#' Rich colour palette (r4ss)
#'
#' A subset of `rich.colors()` by Arni Magnusson (gregmisc package, now
#' gplots), a.k.a. `rich.colors.short()` in r4ss (authors Arni Magnusson and
#' Ian Taylor), with alpha transparency. Same function as `r4col()` in
#' FLSelex, `rc()` in ss3diags and `rc4()` in JABBA and FLRef.
#'
#' @param n number of colours
#' @param alpha transparency (0-1)
#' @return character vector of colour codes
#' @seealso `r4ss::rich.colors.short()`, `gplots::rich.colors()`
#' @export
r4col <- function(n, alpha = 1) {
  x <- seq(0, 1, length.out = n)
  r <- 1 / (1 + exp(20 - 35 * x))
  g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
  b <- stats::dnorm(x, 0.25, 0.15) / max(stats::dnorm(x, 0.25, 0.15))
  grDevices::rgb(r, g, b, alpha = alpha)
}

#' Colours for n runs, scenarios or peels (r4ss rule)
#'
#' Follows the default colour assignment of `r4ss::SSplotComparisons()`:
#' `r4col(n + 1)[-1]` for more than 3 lines, `r4col(n)` for fewer than 3,
#' and blue, red, green3 for exactly 3. Same rule as `sscol()` in ss3diags
#' and `ss3col()` in JABBA. (FLRef's `ss3col()` uses blue, green, red for 3
#' and so differs from r4ss.)
#' Transparency uses `adjustcolor()` as in r4ss (approach by Trevor Branch).
#'
#' @param n number of colours
#' @param alpha transparency (0-1)
#' @return character vector of colour codes
#' @export
jagg_cols <- function(n, alpha = 1) {
  col <- if (n > 3) r4col(n + 1)[-1] else if (n < 3) r4col(n) else c("blue", "red", "green3")
  if (alpha < 1) col <- grDevices::adjustcolor(col, alpha.f = alpha)
  col
}

#' Peel colours for retrospective and hindcast plots
#'
#' Reference run in black, peels coloured with [jagg_cols()], as in
#' `jbplot_retro()` and `jbplot_hcxval()`.
#'
#' @param npeels number of peels, excluding the reference run
#' @param alpha transparency (0-1)
#' @export
jagg_peelcols <- function(npeels, alpha = 1) {
  c(grDevices::adjustcolor("black", alpha.f = alpha), jagg_cols(npeels, alpha))
}
