# PD modelling functions
minimodel <- function(dat, var, numbins = 50, default_flag = "dumdef1") {
  dflt_col <- default_flag
  dff <- copy(dat[, c(var, dflt_col), with = FALSE])
  setDT(dff)
  dff <- dff[complete.cases(dff)]
  setkeyv(dff, var)
  dff[, id := .I]
  colnames(dff) <- c("var", dflt_col, "cumweights")
  dff[, group := bucket(cumweights, bins = numbins)]
  dff[, defrate := mean(get(eval(dflt_col))), by = group]
  setDF(dff)

  return(dff)
}

trans_var_vec <- function(var, default_flag, numbins = 50, span = 0.75, na.rm = FALSE) {
  if (any(!sort(unique(default_flag)) %in% c(0, 1))) {
    stop("The default_flag column must only contains 0s and 1s.")
  }

  if (length(var) != length(default_flag)) {
    stop("pds and default_flag must have the same length.")
  }

  tbl <- data.table(varx = var, dflt = default_flag)
  tbl[, id := .I]

  tbl[, group := bucket(varx, bins = numbins, na.rm = TRUE)]
  temp <- tbl[, .(defrate = mean(dflt, na.rm = TRUE)), by = group]
  temp <- temp[!is.na(group)]
  mdl <- loess(defrate ~ group, data = temp, span = span)
  temp[, predl := predict(mdl)]
  tbl <- merge(tbl, temp[, .(group, predl)], by = "group", all.x = TRUE)
  setorder(tbl, id)

  return(tbl[, predl])
}

trans_var <- function(dat, var, numbins = 50, default_flag = "dumdef1", span = 0.75, na.rm = FALSE) {
  dflt_col <- default_flag
  tbl <- copy(dat[, c(var, dflt_col), with = FALSE])
  setDT(tbl)
  colnames(tbl) <- c("varx", "dflt")

  if (any(!sort(unique(tbl[, dflt])) %in% c(0, 1))) {
    stop("The default_flag column must only contains 0s and 1s.")
  }

  tbl[, id := .I]
  tbl[, group := bucket(varx, bins = numbins, na.rm = TRUE)]
  temp <- tbl[, .(defrate = mean(dflt, na.rm = TRUE)), by = group]
  temp <- temp[!is.na(group)]
  mdl <- loess(defrate ~ group, data = temp, span = span)
  temp[, predl := predict(mdl)]
  tbl <- merge(tbl, temp[, .(group, predl)], by = "group", all.x = TRUE)
  setorder(tbl, id)

  return(tbl[, predl])
}

minimodel_plot <- function(
    dat, var, numbins = 50, span = 0.5, default_flag = "dumdef1", perconly = FALSE, label = NULL
  ) {
  if (is.null(label) || is.na(label)) {
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
  dff.unique <- merge(dff.unique, vrs, all.x = TRUE, by = "group")
  colnames(dff.unique)[ncol(dff.unique)] <- "var"
  dff.unique <- dff.unique[, cls, with = FALSE]
  lspred <- loess(defrate ~ group, data = dff.unique, span = span)
  dff[, loesspd := predict(lspred, newdata = dff)]
  dff.unique[, loesspd := predict(lspred, newdata = dff.unique)]
  ar_var <- AR(dff, "loesspd", default_flag = dflt_col)

  plt1 <- ggplot(dff, aes_string(x = "var")) +
    geom_density() +
    ylab("Density") +
    xlab(label1) +
    ggtitle("Ratio distribution", paste0(label1, ", AR after transform = ", round(ar_var, 4))) +
    theme_minimal()

  dff <- dff.unique[, c(1, 4:6)]
  dff1 <- melt(dff, id.vars = "group")

  plt2 <- ggplot(dff1, aes(x = group, y = value, colour = variable)) +
    geom_point(data = dff1[variable == "defrate"]) +
    geom_line(data = dff1[variable == "loesspd"]) +
    scale_x_continuous(labels = function(x) paste0(floor(100 / numbins * x), "%")) +
    ggtitle("Transform in percentile space") +
    xlab("Percentile") +
    ylab("Default rate") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::percent)

  dff1 <- melt(dff[, c("var", "defrate", "loesspd")], id.vars = "var")

  plt3 <- ggplot(dff1, aes(x = var, y = value, colour = variable)) +
    geom_point(data = dff1[variable == "defrate"]) +
    geom_line(data = dff1[variable == "loesspd"]) +
    ggtitle("Transform in ratio space - Median per group") +
    xlab(label1) +
    ylab(NULL) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::percent)

  lay <- rbind(c(1, 1), c(2, 3))

  if (perconly) {
    return(plt2)
  } else {
    return(gridExtra::grid.arrange(plt1, plt2, plt3, nrow = 2, layout_matrix = lay))
  }
}

minimodel_plot_list <- function(dat, var, numbins = 50, span = 0.5, default_flag = "dumdef1", label = NULL) {
  if (is.null(label) || is.na(label)) {
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
  dff.unique <- merge(dff.unique, vrs, all.x = TRUE, by = "group")
  colnames(dff.unique)[ncol(dff.unique)] <- "var"
  dff.unique <- dff.unique[, cls, with = FALSE]
  lspred <- loess(defrate ~ group, data = dff.unique, span = span)
  dff[, loesspd := predict(lspred, newdata = dff)]
  dff.unique[, loesspd := predict(lspred, newdata = dff.unique)]
  ar_var <- AR(dff, "loesspd", default_flag = dflt_col)

  plt1 <- ggplot(dff, aes_string(x = "var")) +
    geom_density() +
    ylab("Density") +
    xlab(label1) +
    ggtitle("Ratio distribution", paste0(label1, ", AR after transform = ", round(ar_var, 4))) +
    theme_minimal()

  dff <- dff.unique[, c(1, 4:6)]
  dff1 <- melt(dff, id.vars = "group")

  plt2 <- ggplot(dff1, aes(x = group, y = value, colour = variable)) +
    geom_point(data = dff1[variable == "defrate"]) +
    geom_line(data = dff1[variable == "loesspd"]) +
    scale_x_continuous(labels = function(x) paste0(floor(100 / numbins * x), "%")) +
    ggtitle("Transform in percentile space") +
    xlab("Percentile") +
    ylab("Default rate") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::percent)

  dff1 <- melt(dff[, c("var", "defrate", "loesspd")], id.vars = "var")

  plt3 <- ggplot(dff1, aes(x = var, y = value, colour = variable)) +
    geom_point(data = dff1[variable == "defrate"]) +
    geom_line(data = dff1[variable == "loesspd"]) +
    ggtitle("Transform in ratio space - Median per group") +
    xlab(label1) +
    ylab(NULL) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::percent)

  return(list(plt1, plt2, plt3))
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

psi <- function(x1, x2) {
  tmp1 <- data.table(x = x1)
  tmp2 <- data.table(x = x2)

  tmp1 <- tmp1[, .(Count1 = .N), x]
  tmp2 <- tmp2[, .(Count2 = .N), x]

  tst <- merge(tmp1, tmp2, by = "x", all = TRUE)

  if(any(is.na(tmp$Count1) | is.na(tmp$Count2))) {
    warning("The two sets of values do not contain the same unique elements.")
  }

  tst[, Perc1 := Count1 / sum(Count1, na.rm = TRUE)]
  tst[, Perc2 := Count2 / sum(Count2, na.rm = TRUE)]

  tst[, `Partial PSI` := (Perc1 - Perc2) * log(Perc1 / Perc2)]

  return(as.data.frame(tst))
}

create_default_flag <- function(dff, id, date_col, default_status, default_flag_name = "default_flag", horizon = 1, quiet = FALSE) {
  stopifnot(is.data.table(dff))

  tmp <- copy(dff[, c(id, date_col, default_status), with = FALSE])
  tmp[, num_col := .I]
  setorderv(tmp, cols = c(date_col, id))

  ## All cases (dates x ID combos)
  defs <- data.table::CJ(
    Date = sort(unique(tmp[, get(eval(date_col))])), identifier = sort(unique(tmp[, get(eval(id))]))
  )
  defs <- merge(
    tmp[, c(date_col, id, default_status), with = FALSE],
    defs,
    by.x = c(date_col, id), by.y = c("Date", "identifier"), all = TRUE
  )
  colnames(defs) <- c("Date", "identifier", "dflt")
  setorder(defs, Date, identifier)

  ## Look only at defaulted ones
  defs <- defs[identifier %in% defs[dflt == 1, identifier]]
  defs <- defs[!is.na(dflt)]
  defs <- split(defs, by = "identifier")

  if (!quiet) {
    pb <- txtProgressBar(max = length(defs), initial = 0, style = 3)
    setTxtProgressBar(pb, 0)
    for (i in seq_along(defs)) {
      defs[[i]][, eval(default_flag_name) := integer(nrow(defs[[i]]))]
      for (j in seq_len(nrow(defs[[i]]))) {
        date_row <- defs[[i]][j, Date]
        if (any(defs[[i]][(Date >= date_row) & (Date < (date_row + 370 * horizon)), dflt] == 1)) {
          defs[[i]][j, eval(default_flag_name) := 1]
        }
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
  } else {
    for (i in seq_along(defs)) {
      defs[[i]][, eval(default_flag_name) := integer(nrow(defs[[i]]))]
      for (j in seq_len(nrow(defs[[i]]))) {
        date_row <- defs[[i]][j, Date]
        if (any(defs[[i]][(Date >= date_row) & (Date < (date_row + 370 * horizon)), dflt] == 1)) {
          defs[[i]][j, eval(default_flag_name) := 1]
        }
      }
    }
  }

  ## Rbind and merge
  defs <- rbindlist(defs)
  tmp <- merge(
    tmp,
    defs[, c("Date", "identifier", default_flag_name), with = FALSE],
    by.x = c(date_col, id), by.y = c("Date", "identifier"), all.x = TRUE
  )
  rm(defs)

  tmp[is.na(get(eval(default_flag_name))), eval(default_flag_name) := 0]
  tmp <- merge(dff[, -c(default_status), with = FALSE], tmp, by = c(id, date_col))
  setorder(tmp, num_col)
  tmp[, num_col := NULL]

  return(as.data.frame(tmp))
}
