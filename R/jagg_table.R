#' Summary table of one or several runs
#'
#' One row per run, after `flicc::fits_summary()`: goodness of fit (as stored
#' in the fit), the number of indices passing the runs test, and terminal
#' status with credible intervals. Pass hindcasts in `hc` (same order and
#' length as the runs) to add Mohn's rho and the joint MASE.
#'
#' @param x a JABBA or JARA fit, or a list of fits
#' @param hc optional list of hindcasts, one per run (JABBA), or a single
#'   hindcast for a single run
#' @param level credible interval level
#' @param digits number of significant digits
#' @return a data frame
#' @export
jagg_table <- function(x, hc = NULL, level = 0.9, digits = 3) {
  m <- jagg_model(x)
  fits <- if (m %in% c("jabba", "jara")) list(x) else x
  if (!is.null(hc) && m %in% c("jabba", "jara")) hc <- list(hc)
  a <- (1 - level) / 2
  ci <- function(v) stats::quantile(v, c(0.5, a, 1 - a), na.rm = TRUE, names = FALSE)

  rows <- lapply(seq_along(fits), function(i) {
    f <- fits[[i]]; mi <- jagg_model(f)
    st <- f$stats
    gof <- stats::setNames(as.numeric(st[[2]]), as.character(st[[1]]))
    rt <- tryCatch(attr(jagg_runstest(f), "runs"), error = function(e) NULL)
    out <- data.frame(run = f$scenario, model = mi,
                      n_indices = length(unique(jagg_tidy(f, "fits")$index)),
                      runs_passed = if (is.null(rt)) NA else sum(rt$test == "Passed"),
                      runs_untested = if (is.null(rt)) NA else sum(rt$test == "Too few residuals"))
    for (g in intersect(c("RMSE", "SDNR", "DIC"), names(gof))) out[[g]] <- gof[[g]]
    if (mi == "jabba") {
      ty <- max(f$yr)
      kb <- f$kbtrj[f$kbtrj$year == ty & f$kbtrj$type == "fit", ]
      b <- ci(kb$stock); h <- ci(kb$harvest)
      out$year <- ty
      out$`B/Bmsy` <- b[1]; out$`B/Bmsy.lci` <- b[2]; out$`B/Bmsy.uci` <- b[3]
      out$`F/Fmsy` <- h[1]; out$`F/Fmsy.lci` <- h[2]; out$`F/Fmsy.uci` <- h[3]
      out$`P(green)` <- mean(kb$stock >= 1 & kb$harvest <= 1)
    } else {
      pc <- ci(f$posteriors$pop.change)
      out$year <- max(f$pyr)
      out$`change%` <- pc[1]; out$`change%.lci` <- pc[2]; out$`change%.uci` <- pc[3]
    }
    if (!is.null(hc)) {
      h <- hc[[i]]
      ref <- if (jagg_model(h) == "jara_hc") f else NULL
      ms <- suppressMessages(jagg_mase(h, ref = ref))
      out$MASE <- ms$MASE[ms$Index == "joint"]
      rho <- .mohn(.hc_ts(h, ref))
      for (q in intersect(c("B/Bmsy", "F/Fmsy", "Abundance"), rho$quantity))
        out[[paste0("rho.", q)]] <- rho$rho[rho$quantity == q]
    }
    out
  })
  nm <- Reduce(union, lapply(rows, names))
  out <- do.call(rbind, lapply(rows, function(r) { r[setdiff(nm, names(r))] <- NA; r[nm] }))
  num <- vapply(out, is.numeric, logical(1)) &
    !names(out) %in% c("year", "n_indices", "runs_passed", "runs_untested")
  out[num] <- lapply(out[num], signif, digits = digits)
  rownames(out) <- NULL
  out
}
