# Iteration-level output and quantile bands

#' Long iteration-level trajectories
#'
#' JABBA: from `kbtrj` (saved with `save.trj = TRUE`, the default).
#' JARA: from `pop.posterior` (global abundance).
#' @param f a JABBA or JARA fit
#' @return data frame: run, year, type ("fit"/"prj"), iter, quantity, value
#' @keywords internal
.iters <- function(f) {
  m <- jagg_model(f)
  if (m %in% c("jabba", "jabba_kb")) {
    kb <- if (m == "jabba") f$kbtrj else f
    if (is.null(kb)) stop("No iteration-level output: refit with save.trj = TRUE.", call. = FALSE)
    q <- c(B = "B", H = "F", Catch = "Catch", stock = "B/Bmsy", harvest = "F/Fmsy", BB0 = "B/B0")
    q <- q[names(q) %in% names(kb)]
    run <- as.character(kb$run)
    out <- do.call(rbind, lapply(names(q), function(v) data.frame(
      run = run, year = kb$year, type = kb$type, iter = kb$iter,
      quantity = q[[v]], value = kb[[v]])))
    return(out)
  }
  if (m == "jara") {
    pp  <- as.matrix(f$pop.posterior)
    yrs <- f$pyr
    if (ncol(pp) != length(yrs)) yrs <- seq(min(f$yr), length.out = ncol(pp))
    return(data.frame(
      run = f$scenario, year = rep(yrs, each = nrow(pp)),
      type = ifelse(rep(yrs, each = nrow(pp)) %in% f$yr, "fit", "prj"),
      iter = rep(seq_len(nrow(pp)), ncol(pp)),
      quantity = "Abundance", value = as.vector(pp)))
  }
  stop("Iteration-level output is only available for single JABBA or JARA fits.", call. = FALSE)
}

#' Quantiles by run, quantity and year
#' @param it output of `.iters()`
#' @param probs five probabilities: outer low, inner low, median, inner high, outer high
#' @keywords internal
.qbands <- function(it, probs = c(0.05, 0.2, 0.5, 0.8, 0.95)) {
  key <- paste(it$run, it$quantity, it$year, it$type, sep = "\r")
  qs  <- t(vapply(split(it$value, key), stats::quantile, numeric(5),
                  probs = probs, na.rm = TRUE, names = FALSE))
  k   <- do.call(rbind, strsplit(rownames(qs), "\r", fixed = TRUE))
  data.frame(run = k[, 1], quantity = k[, 2], year = as.numeric(k[, 3]), type = k[, 4],
             lo2 = qs[, 1], lo1 = qs[, 2], med = qs[, 3], hi1 = qs[, 4], hi2 = qs[, 5],
             row.names = NULL)
}
