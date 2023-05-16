# Miscellaneous functions
gcorr_pd <- function(uncond_pd, exp_val, rsq, rho) {
  return(pnorm((qnorm(uncond_pd) - sqrt(rsq) * exp_val)/(sqrt(1 - rsq * rho^2))))
}

binomial_test <- function(nobs, ndefs, conf_level = 0.95) {
  if (length(nobs) != length(ndefs)) {
    stop("Lengths of nobs and ndefs must be equal.")
  }
  if (any(ndefs > nobs)) {
    stop("The number of defaults must be less than or equal to the number of observations.")
  }
  if (any((c(nobs, ndefs) %% 1) != 0)) {
    stop("Nobs and ndefs must be integers.")
  }
  if (any(c(nobs, ndefs) < 0)) {
    stop("Nobs and ndefs must be non-negative.")
  }
  if (conf_level <= 0 | conf_level >= 1) {
    stop("conf_level must be greater than 0 and less than 1.")
  }

  alph <- (1 - conf_level) / 2
  lb <- sapply(seq_along(nobs), function(x) qbeta(alph, ndefs[x], nobs[x] - ndefs[x] + 1))
  ub <- sapply(seq_along(nobs), function(x) qbeta(1 - alph, ndefs[x] + 1, nobs[x] - ndefs[x]))

  return(list(lb, ub))
}

bucket <- function(x, bins, na.rm = FALSE) {
  return(.bincode(x, quantile(x, probs = 0:bins/bins, na.rm = na.rm), right = TRUE, include.lowest = TRUE))
}

avPlots_invis <- function(mdl, ...) {
  ff <- tempfile()
  png(filename = ff)
  out <- avPlots(mdl, ...)
  dev.off()
  unlink(ff)

  return(out)
}

added_variable_plots <- function(mdl, point_colour = NA, smooth_colour = NA) {
  stopifnot(class(mdl) == "lm")

  mdls <- avPlots_invis(mdl)
  mdls <- lapply(mdls, as.data.table)

  for (i in seq_along(mdls)) {
    mdls[[i]][, Variable := names(mdls)[i]]
    colnames(mdls[[i]])[1] <- "value"
  }

  mdls <- rbindlist(mdls)
  mdls[, Variable := factor(Variable, unique(Variable))]
  ggplot(mdls, aes_string("value", colnames(mdls)[2])) +
    geom_point(colour = ifelse(!is.na(point_colour), point_colour, "black")) +
    facet_wrap(~ Variable, scales = "free") + xlab("Value | Others") +
    ylab("Variate | Others") +
    geom_smooth(
      method = "lm", se = FALSE, formula = y ~ x, colour = ifelse(!is.na(smooth_colour), smooth_colour, "#3366FF")
    ) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..rr.label.., sep = "~~~")), parse = TRUE)
}

get_column_type <- function(dt, col_type) {
  cls <- unlist(lapply(dt, class))
  return(names(cls[cls == col_type]))
}

mean_ci <- function(x, conf_level = 0.95, na.rm = FALSE) {
  if (!na.rm & any(is.na(x))) {
    stop("There are NAs in the data and na.rm is FALSE.")
  }
  y <- x[!is.na(x)]
  tst <- qt(p = 1 - (1 - conf_level) / 2, df = length(y) - 1)

  return(c(mean(y) - tst * sd(y) / sqrt(length(y)), mean(y), mean(y) + tst * sd(y) / sqrt(length(y))))
}

weighted_mean_ci <- function(x, weights, conf_level = 0.95, na.rm = FALSE) {
  if (!na.rm & any(is.na(x))) {
    stop("There are NAs in the data and na.rm is FALSE.")
  }
  wts <- weights[!is.na(x)]
  y <- x[!is.na(x)]
  nx <- length(y)
  vx <- wtd.var(y, wts, normwt = TRUE)
  tstat <- weighted.mean(y, wts) / sqrt(vx / nx)
  cint <- qt(1 - (1 - conf_level)/2, nx - 1)
  cint <- tstat + c(- cint, cint)

  return(c(cint[1] * sqrt(vx / nx), weighted.mean(y, wts), cint[2] * sqrt(vx / nx)))
}
