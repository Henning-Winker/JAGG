# Trajectory plots of one or several runs: jagg_advice() and jagg_ensemble()
# share one argument set and one engine (.traj_plot), so that they replace
# jbplot_summary() / jbplot_trj() and jbplot_ensemble() consistently.

# runs: a fit, a list of fits, or a JABBA kobe table (fw_jabba output)
.traj_runs <- function(x) {
  m <- jagg_model(x)
  if (m == "jabba_kb") return(list(model = "jabba", fits = NULL, kb = x))
  fits <- if (m %in% c("jabba", "jara")) list(x) else if (m %in% c("jabba_list", "jara_list")) x else
    stop("x must be a JABBA or JARA fit, a list of fits, or fw_jabba() output.", call. = FALSE)
  nm <- names(fits); if (is.null(nm)) nm <- rep("", length(fits))
  nm[!nzchar(nm)] <- vapply(fits[!nzchar(nm)], function(f) as.character(f$scenario), "")
  names(fits) <- make.unique(nm)
  list(model = if (grepl("^jara", m)) "jara" else "jabba", fits = fits, kb = NULL)
}

# iteration-level values of the requested quantities
.traj_iters <- function(r, type, bfrac, bref, fmax) {
  if (r$model == "jara")
    return(do.call(rbind, lapply(names(r$fits), function(n) { i <- .iters(r$fits[[n]]); i$run <- n; i })))
  kb <- if (!is.null(r$kb)) r$kb else do.call(rbind, lapply(names(r$fits), function(n) {
    k <- r$fits[[n]]$kbtrj
    if (is.null(k)) stop("No iteration-level output in run ", n, ": refit with save.trj = TRUE.", call. = FALSE)
    k$run <- n; k }))
  kb$run <- as.character(kb$run)
  # cap projected F at fmax (as jbplot_ensemble)
  pj <- kb$type == "prj"
  if (any(pj)) {
    ratio <- stats::median(kb$H[!pj] / kb$harvest[!pj], na.rm = TRUE)
    kb$harvest[pj] <- pmin(kb$harvest[pj], fmax)
    kb$H[pj] <- pmin(kb$H[pj], fmax * ratio)
  }
  it <- .iters(kb)
  extra <- list()
  if ("Process deviation" %in% type && "Bdev" %in% names(kb)) extra$`Process deviation` <- kb$Bdev
  if ("B/Bfrac" %in% type) extra$`B/Bfrac` <- if (bref == "bmsy") kb$stock / bfrac else kb$BB0 / bfrac
  if (length(extra)) it <- rbind(it, do.call(rbind, lapply(names(extra), function(q)
    data.frame(run = kb$run, year = kb$year, type = kb$type, iter = kb$iter, quantity = q, value = extra[[q]]))))
  it[it$quantity %in% type, ]
}

# reference points: JABBA estimates, bfrac line and user points
.traj_refs <- function(r, type, refpts, bfrac, bref) {
  out <- data.frame(quantity = character(), name = character(), value = numeric(), colour = character())
  add <- function(q, n, v, col = NA) if (q %in% type) out <<- rbind(out, data.frame(quantity = q, name = n, value = v, colour = col))
  if (r$model == "jabba") {
    add("B/Bmsy", "Bmsy", 1); add("F/Fmsy", "Fmsy", 1); add("B/B0", "B0", 1); add("B/Bfrac", "Bfrac", 1, "red")
    if (!is.null(r$fits)) {
      rp <- r$fits[[1]]$refpts; rp <- rp[rp$quant == "hat", ]
      add("B", "Bmsy", rp$bmsy); add("B", "K", rp$k); add("F", "Fmsy", rp$fmsy)
      add("Catch", "MSY", rp$msy); add("B/B0", "Bmsy", rp$bmsy / rp$k)
      if (bref == "b0") add("B", paste0(bfrac, "B0"), bfrac * rp$k, "red")
    }
    if (bref == "bmsy") { add("B/Bmsy", paste0(bfrac, "Bmsy"), bfrac, "red")
      if (!is.null(r$fits)) add("B", paste0(bfrac, "Bmsy"), bfrac * r$fits[[1]]$refpts$bmsy[1], "red")
    } else add("B/B0", paste0(bfrac, "B0"), bfrac, "red")
  }
  if (length(refpts)) {
    nm <- names(refpts); isF <- substr(nm, 1, 1) == "F"
    bq <- if ("B/Bmsy" %in% type) "B/Bmsy" else "B"; fq <- if ("F/Fmsy" %in% type) "F/Fmsy" else "F"
    for (i in seq_along(refpts)) add(if (isF[i]) fq else bq, nm[i], unname(refpts[i]))
  }
  out
}

.traj_plot <- function(x, type, plotCIs, quantiles, inner, joint, bfrac, bref, reflines, reflabels,
                       refpts, cols, alpha, xlim, ncol, fmax, prj, legend, osa, label.size) {
  bref <- match.arg(bref, c("bmsy", "b0"))
  r <- .traj_runs(x)
  it <- .traj_iters(r, type, bfrac, bref, fmax)
  if (!nrow(it)) stop("None of the requested quantities is available: ", paste(type, collapse = ", "), call. = FALSE)
  if (joint) it$run <- "Joint"
  if (!prj) it <- it[it$type == "fit", ]
  if (!is.null(xlim)) it <- it[it$year >= xlim[1] & it$year <= xlim[2], ]
  pr <- if (is.null(inner)) c(quantiles[1], 0.5, 0.5, 0.5, quantiles[2]) else
    c(quantiles[1], inner[1], 0.5, inner[2], quantiles[2])
  q <- .qbands(it, pr)
  if (osa) {
    ly <- max(q$year[q$type == "fit"])
    q <- q[!(q$quantity %in% c("F", "F/Fmsy", "Catch") & q$year == ly & q$type == "fit"), ]
  }
  # join projections to the last fitted year of the same run (or of the
  # reference run, which carries the full fit in fw_jabba output)
  if (any(q$type == "prj")) {
    fit <- q[q$type == "fit", ]
    ny <- tapply(fit$year, fit$run, function(y) length(unique(y)))
    refrun <- names(ny)[which.max(ny)]
    joins <- do.call(rbind, lapply(unique(q$run[q$type == "prj"]), function(rn) {
      f <- fit[fit$run == rn, ]; if (!nrow(f)) f <- fit[fit$run == refrun, ]
      f <- f[f$year == max(f$year), ]; f$run <- rn; f$type <- "prj"; f }))
    q <- rbind(q, joins)
  }
  runs <- unique(it$run)
  q$run <- factor(q$run, runs); q$quantity <- factor(q$quantity, type[type %in% q$quantity])
  nr <- length(runs)
  if (is.null(cols)) cols <- if (nr == 1) "dodgerblue" else jagg_cols(nr)
  cols <- stats::setNames(rep(cols, length.out = nr), runs)

  p <- ggplot2::ggplot(q, ggplot2::aes(x = .data$year, group = interaction(.data$run, .data$type)))
  if (plotCIs) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo2, ymax = .data$hi2, fill = .data$run), alpha = alpha)
    if (!is.null(inner)) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo1, ymax = .data$hi1, fill = .data$run),
                                                       alpha = min(1, 2 * alpha))
  }
  p <- p + if (nr == 1) ggplot2::geom_line(ggplot2::aes(y = .data$med, linetype = .data$type), colour = "black", linewidth = 0.6) else
    ggplot2::geom_line(ggplot2::aes(y = .data$med, linetype = .data$type, colour = .data$run), linewidth = 0.7)
  if (nr > 1) p <- p + ggplot2::scale_colour_manual(values = cols)
  if (plotCIs) p <- p + ggplot2::scale_fill_manual(values = cols)
  p <- p +
    ggplot2::scale_linetype_manual(values = c(fit = 1, prj = 2), guide = "none") +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = ncol) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = NULL) + theme_jagg()
  if (!legend || nr == 1) p <- p + ggplot2::theme(legend.position = "none")
  if (reflines) {
    refs <- .traj_refs(r, type, refpts, bfrac, bref)
    refs <- refs[refs$quantity %in% levels(q$quantity), ]
    if (!reflabels) refs <- refs[!refs$name %in% c("K", "MSY") | refs$quantity %in% c("B", "Catch"), ]
    if (nrow(refs)) {
      refs$quantity <- factor(refs$quantity, levels(q$quantity))
      p <- p + geom_refpts(refs, range(q$year), label.size, labels = reflabels)
    }
    if ("Process deviation" %in% levels(q$quantity))
      p <- p + ggplot2::geom_hline(data = data.frame(quantity = factor("Process deviation", levels(q$quantity))),
                                   ggplot2::aes(yintercept = 0), linetype = 2)
  }
  p
}

#' Advice plot: stock trajectories with nested credible intervals
#'
#' After `FLRef::plotAdvice()`: one facet per quantity, a 90% band and a 60%
#' band around the median, y-axis from zero, and reference points drawn as
#' labelled lines coloured by name (limits red, precautionary points orange,
#' B0/K and triggers blue, MSY targets green). Replaces `jbplot_trj()` and,
#' with a list of fits, `jbplot_summary()`. Takes the same arguments as
#' [jagg_ensemble()]; only the defaults differ.
#'
#' @param x a JABBA or JARA fit, a list of fits (coloured by run, named by
#'   the list names or the fits' scenarios), or `fw_jabba()` output
#' @param type quantities to show, in order. JABBA: any of "B", "F",
#'   "Catch", "B/Bmsy", "F/Fmsy", "B/B0", "B/Bfrac", "Process deviation";
#'   JARA: "Abundance". Default: B, F, Catch (or B/Bmsy, F/Fmsy, Catch with
#'   `rel = TRUE`)
#' @param rel JABBA shortcut for `type = c("B/Bmsy", "F/Fmsy", "Catch")`
#' @param plotCIs draw the credible bands
#' @param quantiles probabilities of the outer band
#' @param inner probabilities of the inner band (NULL for none)
#' @param joint pool all runs into one ("Joint")
#' @param bfrac,bref biomass fraction and reference ("bmsy" or "b0") for the
#'   "B/Bfrac" panel and the red `bfrac` line, as in `jbplot_ensemble()`
#' @param reflines draw reference lines
#' @param reflabels label the reference lines (estimated Bmsy, K, Fmsy, MSY)
#' @param refpts optional named vector of extra reference points, e.g.
#'   `c(Blim = 0.4, Bpa = 0.6)`, in the units of the B/Bmsy and F/Fmsy panels
#'   when shown, else of B and F. Names starting with "F" go to the F panel.
#' @param cols run colours (default: dodgerblue for one run, [jagg_cols()]
#'   for several)
#' @param alpha transparency of the outer band (the inner band uses twice it)
#' @param xlim range of years
#' @param ncol number of facet columns
#' @param fmax cap for projected F/Fmsy, as in `jbplot_ensemble()`
#' @param prj include projection years (drawn dashed)
#' @param legend show the run legend
#' @param osa drop the last year of F and catch (no catch yet in that year)
#' @param label.size size of the reference point labels
#' @return a ggplot object
#' @export
jagg_advice <- function(x, type = NULL, rel = FALSE, plotCIs = TRUE, quantiles = c(0.05, 0.95),
                        inner = c(0.2, 0.8), joint = FALSE, bfrac = 0.5, bref = c("bmsy", "b0"),
                        reflines = TRUE, reflabels = TRUE, refpts = NULL, cols = NULL, alpha = 0.2,
                        xlim = NULL, ncol = NULL, fmax = 5, prj = TRUE, legend = TRUE, osa = FALSE,
                        label.size = 2.5) {
  if (is.null(type)) type <- if (.traj_runs(x)$model == "jara") "Abundance" else
    if (rel) c("B/Bmsy", "F/Fmsy", "Catch") else c("B", "F", "Catch")
  .traj_plot(x, type, plotCIs, quantiles, inner, joint, bfrac, bref, reflines, reflabels, refpts,
             cols, alpha, xlim, ncol, fmax, prj, legend, osa, label.size)
}

#' Ensemble or multi-run trajectories
#'
#' As `jbplot_ensemble()`: medians and 95% credible bands of several runs,
#' scenarios, an ensemble (`joint = TRUE`) or `fw_jabba()` forecasts,
#' coloured with [jagg_cols()], with reference lines at 1 (and 0 for process
#' deviations) and the `bfrac` line. Takes the same arguments as
#' [jagg_advice()]; only the defaults differ.
#'
#' @inheritParams jagg_advice
#' @param type quantities (default: base subplots 1-6 of `jbplot_ensemble()`);
#'   add "B/B0" or "B/Bfrac" for subplots 7-8
#' @return a ggplot object
#' @export
jagg_ensemble <- function(x, type = c("B/Bmsy", "F/Fmsy", "B", "F", "Process deviation", "Catch"),
                          rel = FALSE, plotCIs = TRUE, quantiles = c(0.025, 0.975), inner = NULL,
                          joint = FALSE, bfrac = 0.5, bref = c("bmsy", "b0"), reflines = TRUE,
                          reflabels = FALSE, refpts = NULL, cols = NULL, alpha = 0.25, xlim = NULL,
                          ncol = 3, fmax = 5, prj = TRUE, legend = TRUE, osa = FALSE, label.size = 2.5) {
  if (rel) type <- c("B/Bmsy", "F/Fmsy", "Catch")
  .traj_plot(x, type, plotCIs, quantiles, inner, joint, bfrac, bref, reflines, reflabels, refpts,
             cols, alpha, xlim, ncol, fmax, prj, legend, osa, label.size)
}
