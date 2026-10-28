# JAGG: hindcast cross-validation statistics (MASE)
#
# Follows ss3diags::SSmase (Winker and Kell) for definitions:
#   - one-step-ahead prediction for year T - p + 1 from the run peeled by p
#   - prediction residual  = log(hat) - log(obs)
#   - naive residual       = log(obs_t) - log(obs_prev), obs_prev being the
#                            previous observed year within the evaluation
#                            window (terminal years of all runs)
#   - MASE     = mean|prediction residual| / mean|naive residual|
#   - MASE.adj = mean|prediction residual| / max(mean|naive residual|, MAE.base.adj)
# but pairs prediction and naive residuals by year (not by position), so gaps
# and indices ending before the terminal year are handled correctly.

#' Hindcast cross-validation statistics
#'
#' @param x a JABBA hindcast (`hindcast_jabba()`), JARA hindcast
#'   (`jara_hindcast()`), or tidy fits with columns index, year, obs, hat, peel
#' @param ref JARA only: the reference fit (`jara_hindcast()` contains only
#'   the peeled runs)
#' @param MAE.base.adj minimum denominator for MASE.adj (default 0.1), as
#'   `MAE.base.adj` in ss3diags and `naive.min` in JABBA
#' @param residuals return the prediction and naive residuals by year as well
#' @return data frame with Index, MASE, MASE.adj, MAE.PR, MAE.base, n.eval;
#'   with `residuals = TRUE` a list of `mase` and `residuals`
#' @export
jagg_mase <- function(x, ref = NULL, MAE.base.adj = 0.1, residuals = FALSE) {
  d <- if (is.data.frame(x)) x else .hc_fits(x, ref)
  d <- d[is.finite(d$hat), c("index", "year", "obs", "hat", "peel")]
  peels <- sort(unique(d$peel[d$peel > 0]))
  if (!length(peels)) stop("No peels found: x is not a hindcast.", call. = FALSE)

  ty     <- max(d$year[d$peel == 0])            # terminal year, reference run
  endyrs <- ty - c(0, peels)                     # terminal years of all runs
  eval_win <- sort(endyrs)                       # ss3diags yr.eval

  res <- lapply(split(d, d$index), function(di) {
    ref <- di[di$peel == 0 & is.finite(di$obs) & di$obs > 0, ]
    ow  <- ref[ref$year %in% eval_win, c("year", "obs")]
    ow  <- ow[order(ow$year), ]

    # prediction residuals: year t = ty - p + 1 from peel p
    pr <- do.call(rbind, lapply(peels, function(p) {
      yt  <- ty - p + 1
      hat <- di$hat[di$peel == p & di$year == yt]
      obs <- ref$obs[ref$year == yt]
      if (length(hat) == 1 && length(obs) == 1)
        data.frame(year = yt, pred.res = log(hat) - log(obs))
    }))

    # naive residuals: each observed year in the window vs the previous one
    nv <- if (nrow(ow) > 1)
      data.frame(year = ow$year[-1], naive.res = diff(log(ow$obs)),
                 gap = diff(ow$year))

    r <- if (!is.null(pr) && !is.null(nv)) merge(pr, nv, by = "year") else NULL
    if (is.null(r) || !nrow(r)) {
      out <- data.frame(Index = di$index[1], MASE = NA_real_, MASE.adj = NA_real_,
                        MAE.PR = NA_real_, MAE.base = NA_real_, n.eval = 0L)
      return(list(mase = out, res = NULL))
    }
    maepr  <- mean(abs(r$pred.res))
    scaler <- mean(abs(r$naive.res))
    out <- data.frame(Index = di$index[1], MASE = maepr / scaler,
                      MASE.adj = maepr / max(scaler, MAE.base.adj),
                      MAE.PR = maepr, MAE.base = scaler, n.eval = nrow(r))
    if (nrow(r) < length(peels) || any(r$gap > 1))
      message("Index ", di$index[1], ": MASE from ", nrow(r), " of ", length(peels),
              " prediction residuals", if (any(r$gap > 1)) " (unequal spacing of naive residuals)")
    list(mase = out, res = data.frame(Index = di$index[1], r))
  })

  mase <- do.call(rbind, lapply(res, `[[`, "mase")); rownames(mase) <- NULL
  if (!residuals) return(mase)
  rr <- do.call(rbind, lapply(res, `[[`, "res")); rownames(rr) <- NULL
  list(mase = mase, residuals = rr)
}
