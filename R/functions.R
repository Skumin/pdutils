library(data.table)
library(ggplot2)

binomial_test <- function(n, p, a) {
  return(head(which(pbinom(c(0, seq_len(n)), n, p) >= a) - 1, 1))
}

binomial_test_full <- function(n, p, a, defs) {
  return(data.table::between(defs, binomial_test(n, p, a), binomial_test(n, p, 1 - a)))
}

bucket <- function(x, bins, na.rm = FALSE) {
  return(.bincode(x, quantile(x, probs = 0:bins/bins, na.rm = na.rm), right = TRUE, include.lowest = TRUE))
}

mann_whitney <- function(dat, pd_name, default_flag = 'dumdef1') {
  tmp <- copy(dat)
  setDT(tmp)
  cls <- setdiff(colnames(tmp), c(pd_name, default_flag))
  tmp[, eval(cls) := NULL]
  all <- nrow(tmp)
  defaults <- as.numeric(tmp[, sum(get(eval(default_flag)))])
  tmp[, rank := frank(get(eval(pd_name)))]
  output <- (as.numeric(tmp[get(eval(default_flag)) == 1, sum(rank)]) - defaults * (defaults + 1) / 2) / defaults / (all - defaults)
  return(output)
}

ar_compare <- function(dat, pd_name1, pd_name2, default_flag = 'dumdef1') {
  tmp <- copy(dat)
  setDT(tmp)
  cls <- setdiff(colnames(tmp), c(pd_name1, pd_name2, default_flag))
  tmp[, eval(cls) := NULL]
  MW1 <- mann_whitney(tmp, pd_name1, default_flag)
  MW2 <- mann_whitney(tmp, pd_name2, default_flag)
  all <- nrow(tmp)
  defaults <- as.numeric(tmp[, sum(get(eval(default_flag)))])
  nodefaults <- all - defaults

  setorderv(tmp, pd_name1)
  tmp[, rank1 := .I]
  setorderv(tmp, pd_name2)
  tmp[, rank2 := .I]

  data.default <- copy(tmp[get(eval(default_flag)) == 1])
  data.nodefault <- copy(tmp[get(eval(default_flag)) == 0])

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
  output <- data.frame(Group = c("all"), N = all, Ndef = defaults, AR1, AR2, diff, p.value = p)
  return(output)
}

yearmon <- function(dates) {
  return(as.integer(ifelse(data.table::month(dates) <= 9, paste0(data.table::year(dates), '0', data.table::month(dates)), paste0(data.table::year(dates), data.table::month(dates)))))
}

AR <- function(dff, pd_name, default_flag = "dumdef1") {
  stopifnot(is.data.table(dff))
  temp <- copy(dff[, c(pd_name, default_flag), with = FALSE])
  temp <- temp[complete.cases(temp[, c(pd_name, default_flag), with = FALSE])]
  temp <- as.matrix(temp)
  temp <- temp[order(temp[, 1], temp[, 2], decreasing = TRUE), ]
  arv <- (pracma::trapz(seq(0, 1, length = nrow(temp) + 1), c(0, cumsum(temp[, 2])/sum(temp[, 2])))- .5)/(pracma::trapz(seq(0, 1, length = nrow(temp) + 1), c(0, cumsum(sort(temp[, 2], decreasing = TRUE))/sum(temp[, 2]))) - .5)
  return(arv)
}

minimodel <- function(dat, var, numbins = 50, default_flag = 'dumdef1') {
  dff <- copy(dat[, c(var, default_flag), with = FALSE])
  setDT(dff)
  dff <- dff[complete.cases(dff)]
  setkeyv(dff, var)
  dff[, id := .I]
  colnames(dff) <- c('var', default_flag, 'cumweights')
  dff[, group := .bincode(cumweights, breaks = seq(0, nrow(dff), length = numbins + 1), include.lowest = TRUE)]
  dff[, defrate := mean(get(eval(default_flag))), by = group]
  setDF(dff)
  return(dff)
}

minimodel_plot <- function(dat, var, numbins = 50, span = 0.5, default_flag = 'dumdef1', perconly = FALSE, label = NULL) {
  if(is.null(label) || is.na(label)) {
    label1 <- var
  } else {
    label1 <- label
  }
  dff <- minimodel(dat = dat, var = var, numbins = numbins, default_flag = default_flag)
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
  ar_var <- AR(dff, 'loesspd', default_flag = default_flag)
  plt1 <- ggplot(dff, aes_string(x = 'var')) + geom_density() + ylab('Density') + xlab(label1) + ggtitle('Ratio Distribution', paste0(label1, ', AR After Transform = ', round(ar_var, 4))) + theme_minimal()
  dff <- dff.unique[, c(1, 4:6)]
  dff1 <- melt(dff, id.vars = 'group')
  plt2 <- ggplot(dff1, aes(x = group, y = value, colour = variable)) + geom_point(data = dff1[variable == 'defrate']) + geom_line(data = dff1[variable == 'loesspd']) + scale_x_continuous(labels = function(x) paste0(floor(100 / numbins * x), '%')) + ggtitle('Transform in Percentile Space') + xlab('Percentile') + ylab('Default Rate (%)') + theme_minimal() + theme(legend.position = 'none') + scale_y_continuous(labels = scales::percent)
  dff1 <- melt(dff[, c('var', 'defrate', 'loesspd')], id.vars = 'var')
  plt3 <- ggplot(dff1, aes(x = var, y = value, colour = variable)) + geom_point(data = dff1[variable == 'defrate']) + geom_line(data = dff1[variable == 'loesspd']) + ggtitle('Transform in Ratio Space - Median per Group') + xlab(label1) + ylab(NULL) + theme_minimal() + theme(legend.position = 'none') + scale_y_continuous(labels = scales::percent)
  lay <- rbind(c(1, 1), c(2, 3))
  if(perconly) {
    return(plt2)
  } else {
    return(gridExtra::grid.arrange(plt1, plt2, plt3, nrow = 2, layout_matrix = lay))
  }
}

compare_cap_plot <- function(dat, var1, var2, default_flag = 'dumdef1', lbl = NULL) {
  stopifnot(is.data.table(dat))

  ars <- ar_compare(dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2)))], var1, var2, default_flag = default_flag)[, c(4:5, 7)]

  tmp1 <- copy(dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2))), c(var1, default_flag), with = FALSE])
  setorderv(tmp1, c(var1, default_flag), order = c(-1, -1))
  tmp1[, perc_sample := .I/.N]
  tmp1[, perc_default := cumsum(get(eval(default_flag)))/sum(get(eval(default_flag)))]

  tmp2 <- copy(dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2))), c(var2, default_flag), with = FALSE])
  setorderv(tmp2, c(var2, default_flag), order = c(-1, -1))
  tmp2[, perc_default := cumsum(get(eval(default_flag)))/sum(get(eval(default_flag)))]

  tmp <- tmp1[, c('perc_sample', 'perc_default')]
  tmp[, `Model 2` := tmp2$perc_default]

  rm(tmp1, tmp2)

  colnames(tmp)[1:2] <- c('PercSample', 'Model 1')
  tmp <- melt(tmp, id.vars = 1, variable.name = 'Model')
  tmp[, Model := ifelse(Model == 'Model 1', var1, var2)]
  tmp[, Model := factor(Model, levels = c(var1, var2))]

  if(is.null(lbl)) {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) + geom_line() + xlab('Percentage of Sample Excluded') + ylab('Percentage of Defaulters Excluded') + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + ggtitle('CAP Plots', paste0('AR 1 = ', round(ars[, 1], 3), ', AR 2 = ', round(ars[, 2], 3), '; p-value = ', round(ars[, 3], 4))) + theme_minimal()
  } else {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) + geom_line() + xlab('Percentage of Sample Excluded') + ylab('Percentage of Defaulters Excluded') + scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) + ggtitle(lbl, paste0('AR 1 = ', round(ars[, 1], 3), ', AR 2 = ', round(ars[, 2], 3), '; p-value = ', round(ars[, 3], 4))) + theme_minimal()
  }
}
