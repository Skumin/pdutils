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
  return(1 - (1 - cpds)^(1/tnrs))
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
  cpds <- 1 - (1 - apds)^tnrs
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
  ids <- complete.cases(temp[, c(pd_name, default_flag), with = FALSE])
  temp <- temp[ids]
  temp <- as.matrix(temp)
  temp <- temp[order(temp[, 1], temp[, 2], decreasing = TRUE), ]
  arv <- (pracma::trapz(seq(0, 1, length = nrow(temp) + 1), c(0, cumsum(temp[, 2])/sum(temp[, 2])))- .5)/(pracma::trapz(seq(0, 1, length = nrow(temp) + 1), c(0, cumsum(sort(temp[, 2], decreasing = TRUE))/sum(temp[, 2]))) - .5)
  return(arv)
}

AR_vec <- function(pd, default) {
  temp <- matrix(c(pd, default), ncol = 2, byrow = FALSE)
  temp <- temp[complete.cases(temp), ]
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

cap_plot_data <- function(dat, var1, default_flag = 'dumdef1') {
  stopifnot(is.data.table(dat))

  ars <- AR(dff = dat, pd_name = var1, default_flag = default_flag)

  tmp1 <- copy(dat[!is.na(get(eval(var1))), c(var1, default_flag), with = FALSE])
  setorderv(tmp1, c(var1, default_flag), order = c(-1, -1))
  tmp1[, perc_sample := .I/.N]
  tmp1[, perc_default := cumsum(get(eval(default_flag)))/sum(get(eval(default_flag)))]

  ndef <- tmp1[, sum(get(eval(default_flag)))]

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
    stop("The number of rows of 'boxbounds' must be the same as the number of columns of 'EMat'.")
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
    funvals <- ifelse(funvals1 > funvals, funvals1, funvals)
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
