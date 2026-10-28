#' Advice plot: stock trajectories with nested credible intervals
#'
#' After `FLRef::plotAdvice()`: one facet per quantity, a 90% band (alpha
#' 0.2) and a 60% band (alpha 0.4) around the median, y-axis from zero, and
#' reference points drawn as labelled lines coloured by name (limits red,
#' precautionary points orange, B0/K and triggers blue, MSY targets green).
#'
#' @param x a JABBA or JARA fit, or a list of fits of the same model
#'   (scenarios), which are then coloured by run
#' @param rel JABBA: plot B/Bmsy, F/Fmsy and catch instead of B, F and catch
#' @param panels quantities to show, in order. JABBA: any of "B", "F",
#'   "Catch", "B/Bmsy", "F/Fmsy", "B/B0"; JARA: "Abundance"
#' @param probs five probabilities: outer and inner bands and the median
#' @param refpts optional named vector of extra reference points, e.g.
#'   `c(Blim = 0.5, Bpa = 0.8)`, in relative units when `rel = TRUE`. Names
#'   starting with "B" go to the biomass panel, "F" to the F panel.
#' @param plotrefs draw reference points
#' @param colour band colour for a single run
#' @param ncol number of facet columns
#' @param osa drop the last year of F and catch (no catch yet in that year)
#' @param prj include projection years (drawn dashed)
#' @param label.size size of the reference point labels
#' @return a ggplot object
#' @export
jagg_advice <- function(x, rel = FALSE, panels = NULL,
                        probs = c(0.05, 0.2, 0.5, 0.8, 0.95), refpts = NULL,
                        plotrefs = TRUE, colour = "dodgerblue", ncol = NULL,
                        osa = FALSE, prj = TRUE, label.size = 2.5) {
  m <- jagg_model(x)
  fits <- if (m %in% c("jabba", "jara")) list(x) else if (m %in% c("jabba_list", "jara_list")) x else
    stop("jagg_advice() needs a JABBA or JARA fit or a list of fits.", call. = FALSE)
  model <- if (grepl("^jabba", m)) "jabba" else "jara"

  if (is.null(panels))
    panels <- if (model == "jara") "Abundance" else if (rel) c("B/Bmsy", "F/Fmsy", "Catch") else c("B", "F", "Catch")

  it <- do.call(rbind, lapply(fits, .iters))
  it <- it[it$quantity %in% panels, ]
  if (!prj) it <- it[it$type == "fit", ]
  q <- .qbands(it, probs)
  if (osa) {
    ly <- max(q$year[q$type == "fit"])
    q <- q[!(q$quantity %in% c("F", "F/Fmsy", "Catch") & q$year == ly & q$type == "fit"), ]
  }
  # join fit and projection lines in the last fitted year
  if (any(q$type == "prj")) {
    last <- q[q$type == "fit", ]
    last <- last[stats::ave(last$year, last$run, last$quantity, FUN = max) == last$year, ]
    last$type <- "prj"
    q <- rbind(q, last)
  }
  q$quantity <- factor(q$quantity, levels = panels)
  runs <- unique(q$run); nr <- length(runs)

  p <- ggplot2::ggplot(q, ggplot2::aes(x = .data$year))
  if (nr == 1) {
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo2, ymax = .data$hi2, group = .data$type),
                           fill = colour, alpha = 0.2) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo1, ymax = .data$hi1, group = .data$type),
                           fill = colour, alpha = 0.4) +
      ggplot2::geom_line(ggplot2::aes(y = .data$med, linetype = .data$type), linewidth = 0.6)
  } else {
    cols <- stats::setNames(jagg_cols(nr), runs)
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo2, ymax = .data$hi2, fill = .data$run,
                           group = interaction(.data$run, .data$type)), alpha = 0.15) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo1, ymax = .data$hi1, fill = .data$run,
                           group = interaction(.data$run, .data$type)), alpha = 0.25) +
      ggplot2::geom_line(ggplot2::aes(y = .data$med, colour = .data$run, linetype = .data$type,
                         group = interaction(.data$run, .data$type)), linewidth = 0.7) +
      ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols)
  }
  p <- p + ggplot2::scale_linetype_manual(values = c(fit = 1, prj = 2), guide = "none") +
    ggplot2::facet_wrap(~quantity, scales = "free_y", ncol = ncol) +
    ggplot2::expand_limits(y = 0) + scale_x_year() +
    ggplot2::labs(x = "Year", y = NULL) + theme_jagg()

  if (plotrefs) {
    refs <- .advice_refs(fits[[1]], model, rel, refpts)
    refs <- refs[refs$quantity %in% panels, ]
    refs$quantity <- factor(refs$quantity, levels = panels)
    p <- p + geom_refpts(refs, range(q$year), label.size)
  }
  p
}

# reference points for jagg_advice()
.advice_refs <- function(f, model, rel, extra) {
  if (model == "jara") return(data.frame(quantity = character(), name = character(), value = numeric()))
  rp <- f$refpts[f$refpts$quant == "hat", ]
  bmsyk <- rp$bmsy / rp$k
  r <- if (rel) data.frame(
    quantity = c("B/Bmsy", "F/Fmsy", "B/B0", "B/B0", "Catch"),
    name     = c("Bmsy", "Fmsy", "Bmsy", "B0", "MSY"),
    value    = c(1, 1, bmsyk, 1, rp$msy)) else data.frame(
    quantity = c("B", "B", "F", "Catch", "B/B0", "B/B0"),
    name     = c("Bmsy", "K", "Fmsy", "MSY", "Bmsy", "B0"),
    value    = c(rp$bmsy, rp$k, rp$fmsy, rp$msy, bmsyk, 1))
  if (length(extra)) {
    nm <- names(extra)
    bq <- if (rel) "B/Bmsy" else "B"; fq <- if (rel) "F/Fmsy" else "F"
    r <- rbind(r, data.frame(quantity = ifelse(substr(nm, 1, 1) == "F", fq, bq),
                             name = nm, value = unname(extra)))
  }
  r
}
