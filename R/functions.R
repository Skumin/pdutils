library(data.table)
library(ggplot2)

apd <- function(cpds, tenors = NULL) {
  if(any(diff(cpds) <= 0)) {
    stop('Cumulative default probabilities must be increasing.')
  }
  if(any(cpds <= 0)) {
    stop('Cumulative default probabilities must be greater than zero.')
  }
  if(is.null(tenors)) {
    tnrs <- seq_along(cpds)
  } else {
    tnrs <- tenors
  }
  return(1 - (1 - cpds) ^ (1 / tnrs))
}

cpd <- function(apds, tenors = NULL) {
  if(any(apds <= 0)) {
    stop('Annualized default probabilities must be greater than zero.')
  }
  if(is.null(tenors)) {
    tnrs <- seq_along(apds)
  } else {
    tnrs <- tenors
  }
  cpds <- 1 - (1 - apds) ^ tnrs
  if(any(diff(cpds) <= 0)) {
    stop('The annualized default probabilities result in non-increasing cumulative defaut probabilities.')
  } else {
    return(cpds)
  }
}

fpd <- function(cpds) {
  if(any(diff(cpds) <= 0)) {
    stop('Cumulative default probabilities must be increasing.')
  }
  if(any(cpds <= 0)) {
    stop('Cumulative default probabilities must be greater than zero.')
  }
  fpds <- c(cpds[1], sapply(seq_along(cpds)[-length(cpds)], function(x) 1 - (1 - cpds[x + 1])/(1 - cpds[x])))
  return(fpds)
}

bcr <- function(nobs, ndefs, rho, tau, periods, ci = 0.7, simulations = 1000) {
  if(!class(nobs) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(ndefs) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(rho) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(tau) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(periods) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(ci) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(simulations) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(any(c(length(nobs), length(ndefs), length(rho), length(tau), length(periods), length(ci), length(simulations)) != 1)) {
    stop('All arguments must be of length 1.')
  }

  if(ndefs > nobs) {
    stop('nobs must be greater than ndefs.')
  }

  if(ndefs <= 0 | nobs <= 0) {
    stop('nobs and ndefs must be greater than 0.')
  }

  if(nobs %% 1 != 0 | ndefs %% 1 != 0) {
    stop('The number of defaults and the number of observations must be integers.')
  }

  if(periods %% 1 != 0 | simulations %% 1 != 0) {
    stop('The number of periods and the numbers of simulations must be integers.')
  }

  if(periods < 1 | simulations < 1) {
    stop('periods and simulations must be greater than or equal to 1.')
  }

  if(class(rho) != 'numeric' | class(tau) != 'numeric' | class(ci) != 'numeric') {
    stop('rho, tau, and ci must be floats.')
  }

  if(!data.table::between(rho, 0., 1., incbounds = TRUE) | !data.table::between(tau, 0., 1., incbounds = TRUE)) {
    stop('rho and tau must be between 0 and 1.')
  }

  if(ci <= 0 | ci > 1) {
    stop('ci must be greater than 0 and less than or equal to 1.')
  }

  return(pt_multi_pd(nobs, ndefs, rho, tau, periods, ci, simulations))
}

gcorr_pd <- function(uncond_pd, exp_val, rsq, rho) {
  return(pnorm((qnorm(uncond_pd) - sqrt(rsq) * exp_val)/(sqrt(1 - rsq * rho^2))))
}

pluto_tasche <- function(nobs, ndefs, rho, tau, periods, ci = 0.7, simulations = 1000) {
  # From the highest risk to the lowest risk
  if(length(nobs) != length(ndefs)) {
    stop('nobs and ndefs must have the same length.')
  }

  if(length(nobs) < 2 | length(ndefs) < 2) {
    stop('nobs and ndefs must be vectors.')
  }

  if(!class(nobs) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(ndefs) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(rho) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(tau) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(periods) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(ci) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(!class(simulations) %in% c('numeric', 'integer')) {
    stop('All arguments must be numeric.')
  }

  if(any(c(length(rho), length(tau), length(periods), length(ci), length(simulations)) != 1)) {
    stop('All arguments must be of length 1.')
  }

  if(any(ndefs < 0) | any(nobs <= 0)) {
    stop('All elements of nobs must be greater than 0 and of ndefs greater than or equal to 0.')
  }

  if(sum(ndefs) < 1) {
    stop('There must be at least one default.')
  }

  if(any(nobs %% 1 != 0) | any(ndefs %% 1 != 0)) {
    stop('The numbers of defaults and the numbers of observations must be integers.')
  }

  if(periods %% 1 != 0 | simulations %% 1 != 0) {
    stop('The number of periods and the numbers of simulations must be integers.')
  }

  if(periods < 1 | simulations < 1) {
    stop('periods and simulations must be greater than or equal to 1.')
  }

  if(class(rho) != 'numeric' | class(tau) != 'numeric' | class(ci) != 'numeric') {
    stop('rho, tau, and ci must be floats.')
  }

  if(!data.table::between(rho, 0., 1., incbounds = TRUE) | !data.table::between(tau, 0., 1., incbounds = TRUE)) {
    stop('rho and tau must be between 0 and 1.')
  }

  if(ci <= 0 | ci > 1) {
    stop('ci must be greater than 0 and less than or equal to 1.')
  }

  # A basic check to see if the default rate is decreasing on average. In future better checks may be needed
  if(mean(diff(ndefs/nobs)) > 0) {
    warning('The default rate is, on average, increasing, when moving through grades. Make sure nobs and ndefs are sorted from the worst (i.e., highest risk) to the best grade.')
  }

  # Run the C++ function here, after all checks have been performed
  return(pt_multi_pd_full(nobs, ndefs, rho, tau, periods, ci, simulations))
}

binomial_test_old <- function(n, p, a) {
  return(head(which(pbinom(c(0, seq_len(n)), n, p) >= a) - 1, 1))
}

binomial_test <- function(nobs, ndefs, conf_level = 0.95) {
  if(length(nobs) != length(ndefs)) {
    stop('Lengths of nobs and ndefs must be equal.')
  }
  if(any(ndefs > nobs)) {
    stop('The number of defaults must be less than or equal to the number of observations.')
  }
  if(any((c(nobs, ndefs) %% 1) != 0)) {
    stop('Nobs and ndefs must be integers.')
  }
  if(any(c(nobs, ndefs) < 0)) {
    stop('Nobs and ndefs must be non-negative.')
  }
  if(conf_level <= 0 | conf_level >= 1) {
    stop('conf_level must be greater than 0 and less than 1.')
  }

  alph <- (1 - conf_level) / 2
  lb <- sapply(seq_along(nobs), function(x) qbeta(alph, ndefs[x], nobs[x] - ndefs[x] + 1))
  ub <- sapply(seq_along(nobs), function(x) qbeta(1 - alph, ndefs[x] + 1, nobs[x] - ndefs[x]))
  return(list(lb, ub))
}

bucket <- function(x, bins, na.rm = FALSE) {
  return(.bincode(x, quantile(x, probs = 0:bins/bins, na.rm = na.rm), right = TRUE, include.lowest = TRUE))
}

mann_whitney <- function(dat, pd_name, default_flag = 'dumdef1', na.rm = FALSE) {
  stopifnot(is.data.table(dat))
  dflt_col <- default_flag
  tmp <- copy(dat[, c(pd_name, dflt_col), with = FALSE])

  if(any(!sort(unique(unlist(tmp[, 2]))) %in% c(0, 1))) {
    stop("The default_flag column must only contains 0s and 1s.")
  }

  if(sum(complete.cases(tmp)) != nrow(tmp)) {
    if(!na.rm) {
      stop('There are NAs in the two columns and na.rm is FALSE.')
    } else {
      tmp <- tmp[complete.cases(tmp)]
    }
  }

  allobs <- nrow(tmp)
  defaults <- as.numeric(tmp[, sum(get(eval(dflt_col)))])
  tmp[, ranker := data.table::frank(get(eval(pd_name)), ties.method = 'average')]
  output <- (as.numeric(tmp[get(eval(dflt_col)) == 1, sum(ranker)]) - defaults * (defaults + 1) / 2) / defaults / (allobs - defaults)
  return(output)
}

mann_whitney_vec <- function(pds, default_flag, na.rm = FALSE) {
  if(any(!sort(unique(default_flag)) %in% c(0, 1))) {
    stop("The default_flag column must only contains 0s and 1s.")
  }

  if(length(pds) != length(default_flag)) {
    stop('pds and default_flag must have the same length.')
  }

  if(sum(is.na(pds)) > 0 | sum(is.na(default_flag)) > 0) {
    if(!na.rm) {
      stop('There are NAs in the two columns and na.rm is FALSE.')
    } else {
      ids <- union(which(is.na(pds)), which(is.na(default_flag)))
      pds <- pds[-ids]
      default_flag <- default_flag[-ids]
    }
  }

  allobs <- length(pds)
  defaults <- sum(default_flag)
  ranker <- data.table::frank(pds, ties.method = 'average')
  output <- (sum(ranker[default_flag == 1]) - defaults * (defaults + 1) / 2) / defaults / (allobs - defaults)
  return(output)
}

mann_whitney_multiclass <- function(dat, pred, resp, na.rm = FALSE) {
  stopifnot(is.data.table(dat))

  if(any(unlist(lapply(dat[, c(pred, resp), with = FALSE], class)) != 'factor')) {
    stop('The pred and resp columns must be factors.')
  }

  resp_col <- resp
  tmp <- copy(dat[, c(pred, resp_col), with = FALSE])

  if(sum(complete.cases(tmp)) != nrow(tmp)) {
    if(!na.rm) {
      stop('There are NAs in the two columns and na.rm is FALSE.')
    } else {
      tmp <- tmp[complete.cases(tmp)]
    }
  }

  actl_g <- levels(tmp[, get(eval(resp_col))])
  num_g <- data.table(lab = actl_g, lab_id = seq_along(actl_g))

  tmp <- merge(tmp, setNames(num_g, c(resp_col, 'actual_id')), by = resp_col, all.x = TRUE)
  tmp <- merge(tmp, setNames(num_g, c(pred, 'pred_id')), by = pred, all.x = TRUE)

  actuals <- tmp[, actual_id]
  predicted <- tmp[, pred_id]

  mat <- matrix(0, nrow = length(actl_g), ncol = length(actl_g))

  for(i in seq_along(actl_g)) {
    for(j in seq_along(actl_g)) {
      if(j >= i) {
        next
      } else {
        mat[i, j] <- mann_whitney_vec(predicted[actuals == i | actuals == j], as.integer(actuals[actuals == i | actuals == j] == i))
      }
    }
  }

  return(mean(mat[lower.tri(mat)]))
}

AR <- function(dff, pd_name, default_flag = "dumdef1", na.rm = FALSE) {
  return(mann_whitney(dat = dff, pd_name = pd_name, default_flag = default_flag, na.rm = na.rm) * 2 - 1)
}

AR_vec <- function(pds, default_flag, na.rm = FALSE) {
  return(mann_whitney_vec(pds = pds, default_flag = default_flag, na.rm = na.rm) * 2 - 1)
}

kolmogorov_smirnov <- function(dat, pd_name, default_flag = 'dumdef1', na.rm = FALSE) {
  stopifnot(is.data.table(dat))
  dflt_col <- default_flag
  tmp <- copy(dat[, c(pd_name, dflt_col), with = FALSE])

  if(any(sort(unique(unlist(tmp[, 2]))) != c(0, 1))) {
    stop("The default_flag column must contain 0s and 1s.")
  }

  if(sum(complete.cases(tmp)) != nrow(tmp)) {
    if(!na.rm) {
      stop('There are NAs in the two columns and na.rm is FALSE.')
    } else {
      tmp <- tmp[complete.cases(tmp)]
    }
  }

  nondefs <- tmp[get(eval(dflt_col)) == 0, get(eval(pd_name))]
  defs <- tmp[get(eval(dflt_col)) == 1, get(eval(pd_name))]
  cdf_nondefs <- ecdf(nondefs)
  cdf_defs <- ecdf(defs)

  xs <- sort(union(defs, nondefs))

  if(length(xs) <= 50) {
    warning(paste0(pd_name, ' is either discrete or the dataset is too small. The location of the maximum difference may be imprecise.'))
  }

  x_max <- xs[which.max(cdf_nondefs(xs) - cdf_defs(xs))]
  y_val <- c(cdf_nondefs(x_max), cdf_defs(x_max))

  return(data.frame(x_pos = x_max, y_min = y_val[2], y_max = y_val[1], ks_stat = y_val[1] - y_val[2], p_val = ks.test(defs, nondefs)$p.val))
}

ar_compare <- function(dat, pd_name1, pd_name2, default_flag = 'dumdef1') {
  stopifnot(is.data.table(dat))
  .dflt_col <- default_flag
  tmp <- copy(dat[, c(pd_name1, pd_name2, .dflt_col), with = FALSE])
  MW1 <- mann_whitney(tmp, pd_name1, .dflt_col)
  MW2 <- mann_whitney(tmp, pd_name2, .dflt_col)
  allobs <- nrow(tmp)
  defaults <- as.numeric(tmp[, sum(get(eval(.dflt_col)))])
  nodefaults <- allobs - defaults

  setorderv(tmp, pd_name1)
  tmp[, rank1 := .I]
  setorderv(tmp, pd_name2)
  tmp[, rank2 := .I]

  data.default <- copy(tmp[get(eval(.dflt_col)) == 1])
  data.nodefault <- copy(tmp[get(eval(.dflt_col)) == 0])

  data.default[, rank2.d := .I]
  setorderv(data.default, pd_name1)
  data.default[, rank1.d := .I]

  setorderv(data.nodefault, pd_name1)
  data.nodefault[, rank1.nd := .I]
  setorderv(data.nodefault, pd_name2)
  data.nodefault[, rank2.nd := .I]

  data.default[, V1_central := (rank1 - rank1.d)/nodefaults - MW1]
  data.default[, V2_central := (rank2 - rank2.d)/nodefaults - MW2]

  data.default[, s11 := V1_central^2]
  data.default[, s22 := V2_central^2]
  data.default[, s21 := V1_central * V2_central]

  data.nodefault[, V1_central := (rank1 - rank1.nd)/defaults - MW1]
  data.nodefault[, V2_central := (rank2 - rank2.nd)/defaults - MW2]

  data.nodefault[, s11 := V1_central^2]
  data.nodefault[, s22 := V2_central^2]
  data.nodefault[, s21 := V1_central * V2_central]

  s1d <- as.numeric(data.default[, sum(s11)])/(defaults - 1)
  s2d <- as.numeric(data.default[, sum(s22)])/(defaults - 1)
  s12d <- as.numeric(data.default[, sum(s21)])/(defaults - 1)

  s1nd <- as.numeric(data.nodefault[, sum(s11)])/(nodefaults - 1)
  s2nd <- as.numeric(data.nodefault[, sum(s22)])/(nodefaults - 1)
  s12nd <- as.numeric(data.nodefault[, sum(s21)])/(nodefaults - 1)

  Var <- (s1d/defaults + s1nd/nodefaults) + (s2d/defaults + s2nd/nodefaults) - 2 * (s12d/defaults + s12nd/nodefaults)
  p <- 1 - pchisq((MW1 - MW2)^2/Var, 1)
  AR1 <- 2 * MW1 - 1
  AR2 <- 2 * MW2 - 1
  diff <- AR1 - AR2
  output <- data.frame(Group = c("all"), N = allobs, Ndef = defaults, AR1, AR2, diff, p.value = p)
  return(output)
}

ar_ci <- function(dat, pd_name, default_flag = 'dumdef1', conf_level = 0.95, na.rm = FALSE) {
  stopifnot(is.data.table(dat))
  .dflt_col <- default_flag
  tmp <- copy(dat[, c(pd_name, .dflt_col), with = FALSE])

  if(sum(complete.cases(tmp)) != nrow(tmp)) {
    if(!na.rm) {
      stop('There are NAs in the two columns and na.rm is FALSE.')
    } else {
      tmp <- tmp[complete.cases(tmp)]
    }
  }

  if(length(unique(tmp[, get(eval(.dflt_col))])) == 1) {
    return(as.numeric(c(NA, NA)))
  }

  roc_profile <- pROC::roc(tmp[, get(eval(.dflt_col))], tmp[, get(eval(pd_name))], quiet = TRUE, direction = '<')

  if(as.numeric(roc_profile$auc) == 1) {
    return(c(1, 1))
  }

  cis <- pROC::ci.auc(roc_profile, conf.level = conf_level, quiet = TRUE)
  cis <- as.numeric(cis)[-2]

  return(pmax(pmin(c(2 * cis[1] - 1, 2 * cis[2] - 1), 1), 0))
}

yearmon <- function(dates) {
  mnths <- data.table::month(dates)
  yrs <- data.table::year(dates)
  return(as.integer(data.table::fifelse(mnths <= 9, stringi::stri_c(yrs, '0', mnths), stringi::stri_c(yrs, mnths))))
}

minimodel <- function(dat, var, numbins = 50, default_flag = 'dumdef1') {
  dflt_col <- default_flag
  dff <- copy(dat[, c(var, dflt_col), with = FALSE])
  setDT(dff)
  dff <- dff[complete.cases(dff)]
  setkeyv(dff, var)
  dff[, id := .I]
  colnames(dff) <- c('var', dflt_col, 'cumweights')
  dff[, group := bucket(cumweights, bins = numbins)]
  dff[, defrate := mean(get(eval(dflt_col))), by = group]
  setDF(dff)
  return(dff)
}

trans_var_vec <- function(var, default_flag, numbins = 50, span = 0.75, na.rm = FALSE) {
  if(any(!sort(unique(default_flag)) %in% c(0, 1))) {
    stop("The default_flag column must only contains 0s and 1s.")
  }

  if(length(var) != length(default_flag)) {
    stop('pds and default_flag must have the same length.')
  }

  tbl <- data.table(varx = var, dflt = default_flag)
  tbl[, id := .I]

  tbl[, group := bucket(varx, bins = numbins, na.rm = TRUE)]
  temp <- tbl[, .(defrate = mean(dflt, na.rm = TRUE)), by = group]
  temp <- temp[!is.na(group)]
  mdl <- loess(defrate ~ group, data = temp, span = span)
  temp[, predl := predict(mdl)]
  tbl <- merge(tbl, temp[, .(group, predl)], by = 'group', all.x = TRUE)
  setorder(tbl, id)
  return(tbl[, predl])
}

trans_var <- function(dat, var, numbins = 50, default_flag = 'dumdef1', span = 0.75, na.rm = FALSE) {
  dflt_col <- default_flag
  tbl <- copy(dat[, c(var, dflt_col), with = FALSE])
  setDT(tbl)
  colnames(tbl) <- c('varx', 'dflt')

  if(any(!sort(unique(tbl[, dflt])) %in% c(0, 1))) {
    stop("The default_flag column must only contains 0s and 1s.")
  }

  tbl[, id := .I]
  tbl[, group := bucket(varx, bins = numbins, na.rm = TRUE)]
  temp <- tbl[, .(defrate = mean(dflt, na.rm = TRUE)), by = group]
  temp <- temp[!is.na(group)]
  mdl <- loess(defrate ~ group, data = temp, span = span)
  temp[, predl := predict(mdl)]
  tbl <- merge(tbl, temp[, .(group, predl)], by = 'group', all.x = TRUE)
  setorder(tbl, id)
  return(tbl[, predl])
}

minimodel_plot <- function(dat, var, numbins = 50, span = 0.5, default_flag = 'dumdef1', perconly = FALSE, label = NULL) {
  if(is.null(label) || is.na(label)) {
    label1 <- var
  } else {
    label1 <- label
  }
  dflt_col <- default_flag
  dff <- minimodel(dat = dat, var = var, numbins = numbins, default_flag = dflt_col)
  setDT(dff)
  vrs <- dff[, median(var), group]
  dff.unique <- dff[!duplicated(group)]
  cls <- copy(colnames(dff.unique))
  dff.unique[, var := NULL]
  dff.unique <- merge(dff.unique, vrs, all.x = TRUE, by = 'group')
  colnames(dff.unique)[ncol(dff.unique)] <- 'var'
  dff.unique <- dff.unique[, cls, with = FALSE]
  lspred <- loess(defrate ~ group, data = dff.unique, span = span)
  dff[, loesspd := predict(lspred, newdata = dff)]
  dff.unique[, loesspd := predict(lspred, newdata = dff.unique)]
  ar_var <- AR(dff, 'loesspd', default_flag = dflt_col)
  plt1 <- ggplot(dff, aes_string(x = 'var')) + geom_density() + ylab('Density') + xlab(label1) + ggtitle('Ratio distribution', paste0(label1, ', AR after transform = ', round(ar_var, 4))) + theme_minimal()
  dff <- dff.unique[, c(1, 4:6)]
  dff1 <- melt(dff, id.vars = 'group')
  plt2 <- ggplot(dff1, aes(x = group, y = value, colour = variable)) + geom_point(data = dff1[variable == 'defrate']) + geom_line(data = dff1[variable == 'loesspd']) + scale_x_continuous(labels = function(x) paste0(floor(100 / numbins * x), '%')) + ggtitle('Transform in percentile space') + xlab('Percentile') + ylab('Default rate') + theme_minimal() + theme(legend.position = 'none') + scale_y_continuous(labels = scales::percent)
  dff1 <- melt(dff[, c('var', 'defrate', 'loesspd')], id.vars = 'var')
  plt3 <- ggplot(dff1, aes(x = var, y = value, colour = variable)) + geom_point(data = dff1[variable == 'defrate']) + geom_line(data = dff1[variable == 'loesspd']) + ggtitle('Transform in ratio space - Median per group') + xlab(label1) + ylab(NULL) + theme_minimal() + theme(legend.position = 'none') + scale_y_continuous(labels = scales::percent)
  lay <- rbind(c(1, 1), c(2, 3))
  if(perconly) {
    return(plt2)
  } else {
    return(gridExtra::grid.arrange(plt1, plt2, plt3, nrow = 2, layout_matrix = lay))
  }
}

minimodel_plot_list <- function(dat, var, numbins = 50, span = 0.5, default_flag = 'dumdef1', label = NULL) {
  if(is.null(label) || is.na(label)) {
    label1 <- var
  } else {
    label1 <- label
  }
  dflt_col <- default_flag
  dff <- minimodel(dat = dat, var = var, numbins = numbins, default_flag = dflt_col)
  setDT(dff)
  vrs <- dff[, median(var), group]
  dff.unique <- dff[!duplicated(group)]
  cls <- copy(colnames(dff.unique))
  dff.unique[, var := NULL]
  dff.unique <- merge(dff.unique, vrs, all.x = TRUE, by = 'group')
  colnames(dff.unique)[ncol(dff.unique)] <- 'var'
  dff.unique <- dff.unique[, cls, with = FALSE]
  lspred <- loess(defrate ~ group, data = dff.unique, span = span)
  dff[, loesspd := predict(lspred, newdata = dff)]
  dff.unique[, loesspd := predict(lspred, newdata = dff.unique)]
  ar_var <- AR(dff, 'loesspd', default_flag = dflt_col)
  plt1 <- ggplot(dff, aes_string(x = 'var')) + geom_density() + ylab('Density') + xlab(label1) + ggtitle('Ratio distribution', paste0(label1, ', AR after transform = ', round(ar_var, 4))) + theme_minimal()
  dff <- dff.unique[, c(1, 4:6)]
  dff1 <- melt(dff, id.vars = 'group')
  plt2 <- ggplot(dff1, aes(x = group, y = value, colour = variable)) + geom_point(data = dff1[variable == 'defrate']) + geom_line(data = dff1[variable == 'loesspd']) + scale_x_continuous(labels = function(x) paste0(floor(100 / numbins * x), '%')) + ggtitle('Transform in percentile space') + xlab('Percentile') + ylab('Default rate') + theme_minimal() + theme(legend.position = 'none') + scale_y_continuous(labels = scales::percent)
  dff1 <- melt(dff[, c('var', 'defrate', 'loesspd')], id.vars = 'var')
  plt3 <- ggplot(dff1, aes(x = var, y = value, colour = variable)) + geom_point(data = dff1[variable == 'defrate']) + geom_line(data = dff1[variable == 'loesspd']) + ggtitle('Transform in ratio space - Median per group') + xlab(label1) + ylab(NULL) + theme_minimal() + theme(legend.position = 'none') + scale_y_continuous(labels = scales::percent)
  return(list(plt1, plt2, plt3))
}

compare_cap_plot_2 <- function(dat, var1, var2, default_flag = 'dumdef1', lbl = NULL) {
  stopifnot(is.data.table(dat))
  if(var1 == var2) {
    stop("'var1' and 'var2' have to differ.")
  }

  dflt_col <- default_flag

  ars <- ar_compare(dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2)))], var1, var2, default_flag = dflt_col)[, c(4:5, 7)]

  tmp1 <- copy(dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2))), c(var1, dflt_col), with = FALSE])
  #setorderv(tmp1, c(var1, dflt_col), order = c(-1, -1))
  setorderv(tmp1, c(var1), order = c(-1))
  tmp1[, perc_sample := .I/.N]
  tmp1[, perc_default := cumsum(get(eval(dflt_col)))/sum(get(eval(dflt_col)))]

  tmp2 <- copy(dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2))), c(var2, dflt_col), with = FALSE])
  #setorderv(tmp2, c(var2, dflt_col), order = c(-1, -1))
  setorderv(tmp2, c(var2), order = c(-1))
  tmp2[, perc_default := cumsum(get(eval(dflt_col)))/sum(get(eval(dflt_col)))]

  tmp <- tmp1[, c('perc_sample', 'perc_default')]
  tmp[, `Model 2` := tmp2$perc_default]

  rm(tmp1, tmp2)

  colnames(tmp)[1:2] <- c('PercSample', 'Model 1')
  tmp <- melt(tmp, id.vars = 1, variable.name = 'Model')
  tmp[, Model := data.table::fifelse(Model == 'Model 1', var1, var2)]
  tmp[, Model := factor(Model, levels = c(var1, var2))]

  if(is.null(lbl)) {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) + geom_line() + xlab('Percentage of Sample Excluded') + ylab('Percentage of Defaulters Excluded') + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + ggtitle('CAP Plots', paste0('AR 1 = ', round(ars[, 1], 3), ', AR 2 = ', round(ars[, 2], 3), '; p-value = ', round(ars[, 3], 4))) + theme_minimal()
  } else {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) + geom_line() + xlab('Percentage of Sample Excluded') + ylab('Percentage of Defaulters Excluded') + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + ggtitle(lbl, paste0('AR 1 = ', round(ars[, 1], 3), ', AR 2 = ', round(ars[, 2], 3), '; p-value = ', round(ars[, 3], 4))) + theme_minimal()
  }
}

compare_cap_plot_n <- function(dat, vars, default_flag = 'dumdef1', lbl = NULL) {
  stopifnot(is.data.table(dat))
  dflt_col <- default_flag
  ids <- complete.cases(dat[, vars, with = FALSE])

  ars <- paste0(paste0(vars, ': ', round(unlist(dat[, lapply(vars, function(x) pdutils::AR_vec(get(eval(x)), get(eval(dflt_col))))]), 3)), collapse = ', ')

  lst <- list()
  for(i in seq_along(vars)) {
    lst[[i]] <- copy(dat[ids, c(vars[i], dflt_col), with = FALSE])
    setorderv(lst[[i]], vars[i], order = -1)
    lst[[i]][, perc_sample := .I/.N]
    lst[[i]][, perc_default := cumsum(get(eval(dflt_col)))/sum(get(eval(dflt_col)))]
  }

  .tmp <- lst[[1]][, .(perc_sample, `Model 1` = perc_default)]
  colnames(.tmp)[2] <- vars[1]
  tmp1 <- Reduce(cbind, lapply(lst[-1], function(z) z[, 'perc_default']))
  colnames(tmp1) <- vars[-1]
  .tmp <- cbind(.tmp, tmp1)

  .tmp <- melt(.tmp, id.vars = 1, variable.name = 'Model')
  .tmp[, Model := factor(Model, levels = vars)]

  if(is.null(lbl)) {
    ggplot(.tmp, aes(x = perc_sample, y = value, colour = Model, group = Model)) + geom_line() + xlab('Percentage of sample excluded') + ylab('Percentage of defaulters excluded') + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + ggtitle('CAP plots', paste0('ARs: ', ars)) + theme_minimal()
  } else {
    ggplot(.tmp, aes(x = perc_sample, y = value, colour = Model, group = Model)) + geom_line() + xlab('Percentage of sample excluded') + ylab('Percentage of defaulters excluded') + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + ggtitle(lbl, paste0('ARs: ', ars)) + theme_minimal()
  }
}

cap_plot_data <- function(dat, var1, default_flag = 'dumdef1') {
  stopifnot(is.data.table(dat))
  dflt_col <- default_flag

  ars <- AR(dff = dat, pd_name = var1, default_flag = dflt_col)

  tmp1 <- copy(dat[!is.na(get(eval(var1))), c(var1, dflt_col), with = FALSE])
  setorderv(tmp1, c(var1), order = c(-1))
  tmp1[, perc_sample := .I/.N]
  tmp1[, perc_default := cumsum(get(eval(dflt_col)))/sum(get(eval(dflt_col)))]

  ndef <- tmp1[, sum(get(eval(dflt_col)))]

  tmp1 <- tmp1[, c('perc_sample', 'perc_default')]
  tmp1[, Perfect := .I / ndef]
  tmp1[, Perfect := pmin(Perfect, 1.0)]

  return(as.data.frame(tmp1))
}

cap_plot <- function(dat, var1, default_flag = 'dumdef1', col1 = NULL, col2 = NULL, lbl = NULL) {
  stopifnot(is.data.table(dat))

  if(is.null(col1)) {
    colour1 <- '#F8766D'
  } else {
    colour1 <- col1
  }

  if(is.null(col2)) {
    colour2 <- '#00BFC4'
  } else {
    colour2 <- col2
  }

  arv <- AR(dff = dat, pd_name = var1, default_flag = default_flag)
  tmp <- cap_plot_data(dat, var1 = var1, default_flag = default_flag)
  setDT(tmp)

  colnames(tmp)[1:2] <- c('PercSample', var1)
  tmp <- melt(tmp, id.vars = 1, variable.name = 'Model')
  tmp[, Model := factor(Model, levels = c(var1, 'Perfect'))]

  if(is.null(lbl)) {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) + scale_colour_manual(values = c(colour1, colour2)) + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'black') + geom_line() + xlab('Percentage of Sample Excluded') + ylab('Percentage of Defaulters Excluded') + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + ggtitle('CAP Plot', paste0('AR = ', round(arv, 3))) + theme_minimal()
  } else {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) + scale_colour_manual(values = c(colour1, colour2)) + geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'black') + geom_line() + xlab('Percentage of Sample Excluded') + ylab('Percentage of Defaulters Excluded') + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + ggtitle(lbl, paste0('AR = ', round(arv, 3))) + theme_minimal()
  }
}

create_dummy_sample <- function(size, mean_pd, ar_target) {
  tmp <- data.table(pdo = pnorm(rnorm(size, qnorm(mean_pd), 1)))
  tmp[, pd := LDPD::QMMRecalibrate(mean_pd, pdo, rep(1, size), AR.target = ar_target)$condPD.ac]
  tmp[, default := rbinom(size, 1, pd)]
  tmp[, pdo := NULL]
  tmp <- tmp[sample(seq_len(size), size, replace = FALSE)]
  tmp[, id := .I]
  return(as.data.frame(tmp[, c(3, 1, 2)]))
}

mean_ci <- function(x, conf_level = 0.95, na.rm = FALSE) {
  if(!na.rm & any(is.na(x))) {
    stop('There are NAs in the data and na.rm is FALSE.')
  }
  y <- x[!is.na(x)]
  tst <- qt(p = 1 - (1 - conf_level) / 2, df = length(y) - 1)
  return(c(mean(y) - tst * sd(y) / sqrt(length(y)), mean(y), mean(y) + tst * sd(y) / sqrt(length(y))))
}

weighted_mean_ci <- function(x, weights, conf_level = 0.95, na.rm = FALSE) {
  if(!na.rm & any(is.na(x))) {
    stop('There are NAs in the data and na.rm is FALSE.')
  }
  wts <- weights[!is.na(x)]
  y <- x[!is.na(x)]
  nx <- length(y)
  vx <- Hmisc::wtd.var(y, wts, normwt = TRUE)
  tstat <- weighted.mean(y, wts) / sqrt(vx / nx)
  cint <- qt(1 - (1 - conf_level)/2, nx - 1)
  cint <- tstat + c(- cint, cint)
  return(c(cint[1] * sqrt(vx / nx), weighted.mean(y, wts), cint[2] * sqrt(vx / nx)))
}

woe <- function(dff, var, default_flag = 'dumdef1') {
  stopifnot(is.data.table(dff))

  dflt_col <- default_flag
  tmp <- copy(dff[, c(var, dflt_col), with = FALSE])

  if(!all(sort(unique(tmp[, get(eval(dflt_col))])) == c(0, 1))) {
    stop('The default_flag column must only contain 0S and 1s.')
  }

  allobs <- nrow(tmp)
  defaults <- as.numeric(tmp[, sum(get(eval(dflt_col)))])

  tst <- setNames(tmp[, .N, by = list(get(eval(var)), get(eval(dflt_col)))], c('variable', dflt_col, 'count'))
  tst[, categ := data.table::fifelse(get(eval(dflt_col)) == 0, 'good', 'bad')]
  tst <- dcast(tst, variable ~ categ, value.var = 'count')

  tst[, perc_bad := bad / defaults]
  tst[, perc_good := good / (allobs - defaults)]

  tst[, eval(c('perc_good', 'perc_bad')) := lapply(.SD, function(z) data.table::fifelse(is.na(z), 0, z)), .SDcols = c('perc_good', 'perc_bad')]
  tst[, WOE := log(perc_bad) - log(perc_good)]
  tst[, WOE := data.table::fifelse(!is.finite(WOE), 0, WOE)]
  return(as.data.frame(setNames(tst[, .(variable, WOE)], c(var, 'woe'))))
}

info_value <- function(dff, var, default_flag = 'dumdef1') {
  stopifnot(is.data.table(dff))

  dflt_col <- default_flag
  tmp <- copy(dff[, c(var, dflt_col), with = FALSE])

  if(!all(sort(unique(tmp[, get(eval(dflt_col))])) == c(0, 1))) {
    stop('The default_flag column must only contain 0S and 1s.')
  }

  allobs <- nrow(tmp)
  defaults <- as.numeric(tmp[, sum(get(eval(dflt_col)))])

  tst <- setNames(tmp[, .N, by = list(get(eval(var)), get(eval(dflt_col)))], c('variable', dflt_col, 'count'))
  tst[, categ := data.table::fifelse(get(eval(dflt_col)) == 0, 'good', 'bad')]
  tst <- dcast(tst, variable ~ categ, value.var = 'count')

  tst[, perc_bad := bad / defaults]
  tst[, perc_good := good / (allobs - defaults)]

  tst[, eval(c('perc_good', 'perc_bad')) := lapply(.SD, function(z) data.table::fifelse(is.na(z), 0, z)), .SDcols = c('perc_good', 'perc_bad')]
  tst[, WOE := log(perc_bad) - log(perc_good)]
  tst[, WOE := data.table::fifelse(!is.finite(WOE), 0, WOE)]

  tst[, iv := (perc_bad - perc_good) * WOE]
  return(sum(tst$iv))
}

psi <- function(x1, x2) {
  tmp1 <- data.table(x = x1)
  tmp2 <- data.table(x = x2)

  tmp1 <- tmp1[, .(Count1 = .N), x]
  tmp2 <- tmp2[, .(Count2 = .N), x]

  tst <- merge(tmp1, tmp2, by = 'x', all = TRUE)

  if(any(is.na(tmp$Count1) | is.na(tmp$Count2))) {
    warning('The two sets of values do not contain the same unique elements.')
  }

  tst[, Perc1 := Count1 / sum(Count1, na.rm = TRUE)]
  tst[, Perc2 := Count2 / sum(Count2, na.rm = TRUE)]

  tst[, `Partial PSI` := (Perc1 - Perc2) * log(Perc1 / Perc2)]
  return(as.data.frame(tst))
}

triangle_plot <- function(dat, var1, var2, pd1, pd2, default_flag = 'dumdef1') {
  stopifnot(is.data.table(dat))

  dflt_col <- default_flag

  tmp <- copy(dat[, c(var1, var2, pd1, pd2, dflt_col), with = FALSE])
  tmp <- tmp[complete.cases(tmp)]

  var1_bin <- 'var1_bin'
  var2_bin <- 'var2_bin'

  tmp[, eval(var1_bin) := .bincode(get(eval(var1)), breaks = quantile(get(eval(var1)), probs = c(0, 0.25, 0.5, 0.75, 1)), right = TRUE, include.lowest = TRUE)]
  tmp[, eval(var2_bin) := .bincode(get(eval(var2)), breaks = quantile(get(eval(var2)), probs = c(0, 0.25, 0.5, 0.75, 1)), right = TRUE, include.lowest = TRUE)]

  sumstat <- tmp[, .(actDR = mean(get(eval(dflt_col))), pd2DR = mean(get(eval(pd2))), pd1DR = mean(get(eval(pd1)))), keyby = .(var1_bin, var2_bin)]
  sumstat <- melt(sumstat, measure.vars = 3:5)
  sumstat[, variable := factor(variable, levels = c('actDR', 'pd2DR', 'pd1DR'))]
  setorder(sumstat, var1_bin, var2_bin, variable)
  sumstat[, id := paste(var1_bin, var2_bin, variable, sep = '_')]
  sumstat <- sumstat[, 4:5]

  positions <- data.table(id = rep(sumstat$id, each = 3), x = c(rep(c(0, 0, 1, 0, 0.5, 1, 0.5, 1, 1), 4), rep(c(0, 0, 1, 0, 0.5, 1, 0.5, 1, 1) + 1, 4), rep(c(0, 0, 1, 0, 0.5, 1, 0.5, 1, 1) + 2, 4), rep(c(0, 0, 1, 0, 0.5, 1, 0.5, 1, 1) + 3, 4)), y = c(c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1), c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 1, c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 2, c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 3, c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1), c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 1, c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 2, c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 3))

  datapoly <- merge(sumstat, positions, by = 'id')
  colnames(datapoly)[2] <- 'PD'
  datapoly[, x := x/4]
  datapoly[, y := y/4]
  datapoly[, PD := PD * 100]
  datapoly[, x_c := sum(x)/3, by = id]
  datapoly[, y_c := sum(y)/3, by = id]

  return(ggplot(datapoly) + geom_polygon(aes(x = x, y = y, fill = PD, group = id)) + geom_text(aes(x = x_c, y = y_c, label = round(PD, 2))) + scale_y_continuous(labels = scales::percent) + scale_x_continuous(labels = scales::percent) + theme(legend.position = "none") + xlab(paste0('Percentiles of ', var1)) + ylab(paste0('Percentiles of ', var2)) + ggtitle('Actual Default Rate and Predicted PDs', paste0('Top left: actual, right-hand side: ', pd1, ', bottom: ', pd2)) + scale_fill_gradient(low = '#0099f7', high = '#f11712'))
}

crossover_deleq <- function(mat, newmat, cr) {
  return((1 - cr) * mat + cr * newmat)
}

isapprox <- function(x, y, atol = 0, rtol = sqrt(.Machine$double.eps)) {
  return(x == y | (is.finite(x) & is.finite(y) & abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)))))
}

deleq <- function(fun, ..., NP, boxbounds, Emat, constr, x0 = NULL, cr = 0.7, f.param = 0.9, maxgen = 100, show.progress = TRUE) {
  if(!is.matrix(boxbounds) | ncol(boxbounds) != 2) {
    stop("'boxbounds' must be a two-column matrix.")
  }
  if(!is.matrix(Emat)) {
    stop("'Emat' must be a matrix.")
  }
  if(ncol(Emat) != nrow(boxbounds)) {
    stop("The number of rows of 'boxbounds' must be the same as the number of columns of 'Emat'.")
  }
  if(!is.matrix(constr)) {
    stop("'constr' must be a matrix.")
  }
  if(nrow(constr) != nrow(Emat)) {
    stop("Number of rows of 'constr' must be the same as that of 'Emat'.")
  }
  if(maxgen < 1 | maxgen %% 1 != 0) {
    stop("'maxgen' must be an integer greater than 0.")
  }
  if(cr < 0 | cr > 1) {
    stop("'cr' must be between 0 and 1.")
  }
  if(f.param < 0 | f.param > 1) {
    stop("'f.param' must be between 0 and 1.")
  }
  if(qr(Emat %*% t(Emat))$rank < nrow(Emat) & is.null(x0)) {
    stop("'Emat' * transpose('Emat') is singular and you have not provided an 'x0'.")
  }
  if(NP < nrow(boxbounds) * 10) {
    warning("'NP' is recommended to be at least 10 times the number of parameters.")
  }
  if(!is.null(x0)) {
    if(!is.matrix(x0) | ncol(x0) > 1) {
      stop("'x0' must be a one-column matrix.")
    }
    if(any(abs(Emat %*% x0 - constr) > sqrt(.Machine$double.eps))) {
      stop("'x0' does not satisfy the equality constraint(s).")
    }
  }
  gen <- 1
  if(is.null(x0)) {
    mat <- gen_init_pop(NP = NP, boxbounds = boxbounds, Emat = Emat, constr = constr)
  } else {
    mat <- gen_init_pop_x0(NP = NP, boxbounds = boxbounds, Emat = Emat, constr = constr, x0 = x0)
  }
  funvals <- apply(X = mat, MARGIN = 1, FUN = fun, ...)
  rbest <- which.max(funvals)
  trugen <- maxgen
  while(gen <= trugen) {
    if(show.progress) {
      print(paste0('Generation ', gen))
    }
    pbest <- rbest
    pbest_ind <- mat[pbest, ]
    pbest_val <- fun(pbest_ind, ...)
    newmat <- mutate_deleq(mat = mat, boxbounds = boxbounds, fParam = f.param)
    newmat <- crossover_deleq(mat = mat, newmat = newmat, cr = cr)
    funvals1 <- apply(X = newmat, MARGIN = 1, FUN = fun, ...)
    mat[funvals1 > funvals, ] <- newmat[funvals1 > funvals, ]
    funvals <- data.table::fifelse(funvals1 > funvals, funvals1, funvals)
    rbest <- which.max(funvals)
    if(!any(isapprox(Emat %*% matrix(mat[rbest, ], ncol = 1), constr))) {
      if(show.progress) {
        print('* Unfeasible, projecting back.')
      }
      trugen <- trugen - 1
      if(gen <= trugen) {
        projmat <- project_population(mat = mat, Emat = Emat, constr = constr)
        funvals1 <- apply(X = projmat, MARGIN = 1, FUN = fun, ...)
        rbest <- which.max(funvals1)
        if(pbest_val > fun(projmat[rbest, ], ...)) {
          projmat[rbest, ] <- pbest_ind
          rbest <- pbest
        }
        mat <- projmat
      } else {
        mat[rbest, ] <- pbest_ind
      }
    }
    gen <- gen + 1
    if(show.progress) {
      print(paste0("** Value = ", funvals[rbest]))
    }
  }
  return(list(params = mat[rbest, ], value = fun(mat[rbest, ], ...), generations = trugen))
}
