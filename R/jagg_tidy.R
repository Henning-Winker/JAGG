# JAGG: model detection and fit -> tidy converters (base R only)

#' Detect which model produced an object
#' @param x a JABBA or JARA fit, hindcast, list of fits, or a JABBA kobe
#'   posterior table (`fit$kbtrj`, or `fw_jabba()` projections)
#' @return one of "jabba", "jabba_hc", "jara", "jara_hc", "jabba_list",
#'   "jara_list", "jabba_kb"
#' @export
jagg_model <- function(x) {
  cl <- intersect(class(x), c("jabba", "jabba_hc", "jara", "jara_hc"))
  if (length(cl)) return(cl[1])
  if (is.data.frame(x) && all(c("year", "run", "type", "iter", "stock", "harvest") %in% names(x)))
    return("jabba_kb")

  is_jabba <- function(o) is.list(o) && all(c("timeseries", "kobe", "refpts") %in% names(o))
  is_jara  <- function(o) is.list(o) && all(c("trj", "fits", "pop.posterior") %in% names(o))

  if (is_jabba(x)) return("jabba")
  if (is_jara(x))  return("jara")
  if (is.list(x) && !is.null(x$peels) && is.data.frame(x$trj) && is.data.frame(x$fits))
    return("jara_hc")
  if (is.list(x) && length(x) > 0) {
    if (all(vapply(x, is_jabba, logical(1)))) {
      # hindcast_jabba() names elements: scenario, "-yyyy", "-yyyy", ...
      if (length(x) > 1 && !is.null(names(x)) && all(grepl("^-[0-9]{4}$", names(x)[-1])))
        return("jabba_hc")
      return("jabba_list")
    }
    if (all(vapply(x, is_jara, logical(1)))) return("jara_list")
  }
  stop("Object not recognised as a JABBA or JARA fit, hindcast or list of fits.", call. = FALSE)
}

#' Convert JABBA/JARA output to a standard tidy data frame
#' @param x fit, hindcast or list of fits
#' @param what "ts" (trajectories) or "fits" (index fits and residuals);
#'   "kobe", "ppc", "posterior" to follow
#' @export
jagg_tidy <- function(x, what = c("ts", "fits")) {
  what  <- match.arg(what)
  model <- jagg_model(x)

  switch(model,
    jabba      = .jabba_tidy(x, what),
    jara       = .jara_tidy(x, what),
    jabba_list = .bind_runs(lapply(x, .jabba_tidy, what = what)),
    jara_list  = .bind_runs(lapply(x, .jara_tidy,  what = what)),
    # peel = number of years removed (0 = reference run), as diags$retro.peels
    jabba_hc   = .bind_runs(Map(function(f, p) {
                   d <- .jabba_tidy(f, what)
                   if (is.data.frame(f$diags)) p <- f$diags$retro.peels[1]
                   d$peel <- p
                   d$end_year <- max(f$yr) - p
                   d
                 }, x, seq_along(x) - 1L)),
    jara_hc    = .jara_hc_tidy(x, what)
  )
}

# ---- JABBA ------------------------------------------------------------------

.jabba_tidy <- function(f, what) {
  run <- f$scenario
  if (what == "ts") {
    ts  <- f$timeseries                     # year x c(mu,lci,uci) x quantity
    yrs <- as.numeric(dimnames(ts)[[1]])
    out <- do.call(rbind, lapply(dimnames(ts)[[3]], function(q) data.frame(
      model = "jabba", run = run, year = yrs, quantity = q,
      estimation = ifelse(yrs %in% f$yr, "fit", "prj"),
      median = ts[, "mu", q], lci = ts[, "lci", q], uci = ts[, "uci", q],
      row.names = NULL)))
    return(out)
  }
  # what == "fits"
  d <- f$diags
  if (!is.data.frame(d))
    stop("No index fits in this JABBA fit (catch-only model?).", call. = FALSE)
  ppd <- f$cpue.ppd                         # year x c(mu,lci,uci,se,obserror) x index
  key <- cbind(as.character(d$year), as.character(d$name))
  ok  <- key[, 1] %in% dimnames(ppd)[[1]] & key[, 2] %in% dimnames(ppd)[[3]]
  pl  <- pu <- se <- rep(NA_real_, nrow(d))
  pl[ok] <- ppd[cbind(key[ok, 1], "lci", key[ok, 2])]
  pu[ok] <- ppd[cbind(key[ok, 1], "uci", key[ok, 2])]
  se[ok] <- ppd[cbind(key[ok, 1], "obserror", key[ok, 2])]
  data.frame(
    model = "jabba", run = run, index = d$name, year = d$year,
    obs = d$obs, obs_se = se, hat = d$hat, lci = d$hat.lci, uci = d$hat.uci,
    ppd_lci = pl, ppd_uci = pu, residual = d$residual,
    peel = d$retro.peels, hindcast = d$hindcast, row.names = NULL)
}

# ---- JARA -------------------------------------------------------------------

.jara_tidy <- function(f, what) {
  if (what == "ts") {
    t <- f$trj
    return(data.frame(
      model = "jara", run = t$scenario, year = t$yr,
      quantity = ifelse(t$name == "global", "global", as.character(t$name)),
      type = t$Type, estimation = t$estimation,
      median = t$mu, lci = t$lci, uci = t$uci, row.names = NULL))
  }
  d <- f$fits
  data.frame(
    model = "jara", run = d$scenario, index = d$name, year = d$year,
    obs = d$obs, obs_se = d$obs.err, hat = d$hat, lci = d$lci, uci = d$uci,
    ppd_lci = d$lpp, ppd_uci = d$upp, residual = d$residual,
    peel = d$retro.peels, hindcast = NA, row.names = NULL)
}

.jara_hc_tidy <- function(hc, what) {
  if (what == "ts") {
    t <- hc$trj                             # factor, level, name, yr, Type, ...
    return(data.frame(
      model = "jara", run = t$level, year = t$yr, quantity = t$name,
      type = t$Type, estimation = t$estimation,
      median = t$mu, lci = t$lci, uci = t$uci, row.names = NULL))
    # TODO: map run (level) -> peel via unique(hc$fits[, c("scenario", "retro.peels")])
  }
  # jara_hindcast() removes peeled years from the data: runs are labelled by
  # peel (scenario) and retro.peels stays 0; there is no reference run
  d <- .jara_tidy(list(fits = hc$fits), "fits")
  d$peel <- as.numeric(as.character(hc$fits$scenario))
  d$hindcast <- FALSE
  d
}

# ---- helpers ----------------------------------------------------------------

.bind_runs <- function(lst) {
  cols <- Reduce(union, lapply(lst, names))
  lst  <- lapply(lst, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] })
  out  <- do.call(rbind, lst)
  rownames(out) <- NULL
  out
}
