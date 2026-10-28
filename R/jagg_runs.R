# Runs test and three-sigma limits, as in JABBA jbruns_sig3() and ss3diags,
# without the randtests dependency.

#' Runs test with three-sigma control limits
#'
#' Limits from the average moving range (Montgomery, eq. 6.33), after
#' removing moving ranges above 3.267 times their mean (Nelson 1982):
#' sigma = AMR / 1.128, limits = +/- 3 sigma. The runs test is the
#' Wald-Wolfowitz test on the signs of the residuals with threshold 0, using
#' the normal approximation as in `randtests::runs.test()`.
#'
#' @param x residuals
#' @param mixing "less" (default, one-sided: too few runs), "two.sided" or
#'   "greater"
#' @return list with `sig3lim` (lcl, ucl) and `p.runs`
#' @export
jagg_runs <- function(x, mixing = c("less", "two.sided", "greater")) {
  mixing <- match.arg(mixing)
  x  <- x[is.finite(x)]
  mr <- abs(diff(x))
  amr <- mean(mr)
  mr  <- mr[mr < 3.267 * amr]
  stdev <- mean(mr) / 1.128
  lim <- c(-3 * stdev, 3 * stdev)

  s  <- sign(x); s <- s[s != 0]
  n1 <- sum(s > 0); n2 <- sum(s < 0); n <- n1 + n2
  if (n1 == 0 || n2 == 0) return(list(sig3lim = lim, p.runs = 0.001))
  runs <- 1 + sum(diff(s) != 0)
  mu <- 2 * n1 * n2 / n + 1
  v  <- 2 * n1 * n2 * (2 * n1 * n2 - n) / (n^2 * (n - 1))
  z  <- (runs - mu) / sqrt(v)
  p  <- switch(mixing, less = stats::pnorm(z), greater = 1 - stats::pnorm(z),
               two.sided = 2 * min(stats::pnorm(z), 1 - stats::pnorm(z)))
  list(sig3lim = lim, p.runs = round(p, 3))
}
