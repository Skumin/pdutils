# Various functions to measure predictive power of models
## The Mann-Whitney U test statistics, also known as the two-sample Wilcoxon test statistic. Equal to the AUC
mann_whitney_vec <- function(pds, default_flag, na.rm = FALSE) {
  if (any(!sort(unique(default_flag)) %in% c(0, 1))) {
    stop("The default_flag column must only contains 0s and 1s.")
  }

  allobs <- length(pds)

  if (allobs != length(default_flag)) {
    stop("pds and default_flag must have the same length.")
  }

  dflt_na <- is.na(default_flag)
  pds_na <- is.na(pds)

  if (any(pds_na) | any(dflt_na)) {
    if (!na.rm) {
      stop("There are NAs in the two vectors and na.rm is FALSE.")
    } else {
      ids <- which(pds_na | dflt_na)
      pds <- pds[-ids]
      default_flag <- default_flag[-ids]
    }
  }

  defaults <- sum(default_flag)
  ranker <- frank(pds, ties.method = "average")

  # This is "Method 2" from here: https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Calculations
  output <- (sum(ranker[default_flag == 1]) - defaults * (defaults + 1) / 2) / defaults / (allobs - defaults)

  return(output)
}

mann_whitney <- function(dat, pd_name, default_flag = "dumdef1", na.rm = FALSE) {
  stopifnot(is.data.frame(dat))

  nms <- names(dat)
  stopifnot(pd_name %chin% nms, default_flag %chin% nms)

  if (is.data.table(dat)) {
    pd_name_tmp <- pd_name
    default_flag_tmp <- default_flag

    pds <- dat[, get(eval(pd_name_tmp))]
    dflt <- dat[, get(eval(default_flag_tmp))]
  } else {
    pds <- dat[, pd_name]
    dflt <- dat[, default_flag]
  }

  return(mann_whitney_vec(pds, dflt, na.rm = na.rm))
}

mann_whitney_multiclass <- function(dat, pred, resp, na.rm = FALSE) {
  stopifnot(is.data.table(dat))

  if(any(unlist(lapply(dat[, c(pred, resp), with = FALSE], class)) != "factor")) {
    stop("The pred and resp columns must be factors.")
  }

  resp_col <- resp
  tmp <- copy(dat[, c(pred, resp_col), with = FALSE])

  if (sum(complete.cases(tmp)) != nrow(tmp)) {
    if (!na.rm) {
      stop("There are NAs in the two columns and na.rm is FALSE.")
    } else {
      tmp <- tmp[complete.cases(tmp)]
    }
  }

  actl_g <- levels(tmp[, get(eval(resp_col))])
  num_g <- data.table(lab = actl_g, lab_id = seq_along(actl_g))

  tmp <- merge(tmp, setNames(num_g, c(resp_col, "actual_id")), by = resp_col, all.x = TRUE)
  tmp <- merge(tmp, setNames(num_g, c(pred, "pred_id")), by = pred, all.x = TRUE)

  actuals <- tmp[, actual_id]
  predicted <- tmp[, pred_id]

  mat <- matrix(0, nrow = length(actl_g), ncol = length(actl_g))

  for (i in seq_along(actl_g)) {
    for (j in seq_along(actl_g)) {
      if (j >= i) {
        next
      } else {
        mat[i, j] <- mann_whitney_vec(
          predicted[actuals == i | actuals == j], as.integer(actuals[actuals == i | actuals == j] == i)
        )
      }
    }
  }

  return(mean(mat[lower.tri(mat)]))
}

## The accuracy ratio, also known as the Gini coefficient
AR <- function(dff, pd_name, default_flag = "dumdef1", na.rm = FALSE) {
  return(mann_whitney(dat = dff, pd_name = pd_name, default_flag = default_flag, na.rm = na.rm) * 2 - 1)
}

AR_vec <- function(pds, default_flag, na.rm = FALSE) {
  return(mann_whitney_vec(pds = pds, default_flag = default_flag, na.rm = na.rm) * 2 - 1)
}

## Compares the ARs of two models
ar_compare <- function(dat, pd_name1, pd_name2, default_flag = "dumdef1", conf_level = 0.95, na.rm = FALSE) {
  stopifnot(conf_level > 0 & conf_level < 1)
  stopifnot(is.data.table(dat))
  stopifnot(pd_name1 %in% names(dat), pd_name2 %in% names(dat), default_flag %in% names(dat))

  .dflt_col <- default_flag
  .tmp <- copy(dat[, c(pd_name1, pd_name2, .dflt_col), with = FALSE])

  fltr <- complete.cases(.tmp)

  if (nrow(.tmp[fltr]) < nrow(dat)) {
    if (!na.rm) {
      stop("There are missing values in the data but `na.rm` is FALSE.")
    } else {
      .tmp <- .tmp[fltr]
    }
  }

  test_result <- delong_paired_test(
    .tmp[, get(eval(pd_name1))], .tmp[, get(eval(pd_name2))], .tmp[, get(eval(default_flag))], conf_level = conf_level,
    test_what = "AR"
  )

  output <- data.frame(
    Group = c("all"), N = nrow(.tmp), Ndef = .tmp[, sum(get(eval(default_flag)))], AR1 = test_result$metrics[1],
    AR2 = test_result$metrics[2], diff = -diff(test_result$metrics), p.value = test_result$p_value
  )
  return(output)
}

## The confidence interval for the AR
ar_ci_vec <- function(pds, default_flag, conf_level = 0.95, na.rm = FALSE) {
  stopifnot(conf_level > 0 & conf_level < 1)

  dflt_na <- is.na(default_flag)
  pds_na <- is.na(pds)

  if (any(pds_na) | any(dflt_na)) {
    if (!na.rm) {
      stop("There are NAs in the two vectors and na.rm is FALSE.")
    } else {
      ids <- which(pds_na | dflt_na)
      pds <- pds[-ids]
      default_flag <- default_flag[-ids]
    }
  }

  if (length(unique(default_flag)) == 1L) {
    return(as.numeric(c(NA, NA)))
  }

  ## Get the DeLong variance; note that `delong_auc_variance` returns the variance of the AUC, not the AR,
  ## so we need to reflect that
  auc_var <- delong_auc_variance(pds, default_flag, na.rm = FALSE)
  ar_sd <- 2 * sqrt(auc_var)
  ar_mid <- AR_vec(pds, default_flag)

  vals <- qnorm(c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2), mean = ar_mid, sd = ar_sd)

  return(pmax(pmin(vals, 1), 0))
}

ar_ci <- function(dat, pd_name, default_flag = "dumdef1", conf_level = 0.95, na.rm = FALSE) {
  stopifnot(is.data.frame(dat))

  nms <- names(dat)
  stopifnot(pd_name %chin% nms, default_flag %chin% nms)

  if (is.data.table(dat)) {
    pd_name_tmp <- pd_name
    default_flag_tmp <- default_flag

    pds <- dat[, get(eval(pd_name_tmp))]
    dflt <- dat[, get(eval(default_flag_tmp))]
  } else {
    pds <- dat[, pd_name]
    dflt <- dat[, default_flag]
  }

  return(ar_ci_vec(pds, dflt, na.rm = na.rm))
}

## The KS test for the difference between two distributions (testing the difference between two ECDFs)
kolmogorov_smirnov <- function(dat, pd_name, default_flag = "dumdef1", na.rm = FALSE) {
  stopifnot(is.data.table(dat))
  dflt_col <- default_flag
  tmp <- copy(dat[, c(pd_name, dflt_col), with = FALSE])

  if (any(sort(unique(unlist(tmp[, 2]))) != c(0, 1))) {
    stop("The default_flag column must contain 0s and 1s.")
  }

  if (sum(complete.cases(tmp)) != nrow(tmp)) {
    if (!na.rm) {
      stop("There are NAs in the two columns and na.rm is FALSE.")
    } else {
      tmp <- tmp[complete.cases(tmp)]
    }
  }

  nondefs <- tmp[get(eval(dflt_col)) == 0, get(eval(pd_name))]
  defs <- tmp[get(eval(dflt_col)) == 1, get(eval(pd_name))]
  cdf_nondefs <- ecdf(nondefs)
  cdf_defs <- ecdf(defs)

  xs <- sort(union(defs, nondefs))

  if (length(xs) <= 50) {
    warning(
      paste0(
        pd_name,
        " is either discrete or the dataset is too small. The location of the maximum difference may be imprecise."
      )
    )
  }

  x_max <- xs[which.max(cdf_nondefs(xs) - cdf_defs(xs))]
  y_val <- c(cdf_nondefs(x_max), cdf_defs(x_max))

  return(
    data.frame(
      x_pos = x_max, y_min = y_val[2], y_max = y_val[1], ks_stat = y_val[1] - y_val[2],
      p_val = ks.test(defs, nondefs)$p.val
    )
  )
}

## CAP plot comparisons
compare_cap_plot_2 <- function(dat, var1, var2, default_flag = "dumdef1", lbl = NULL) {
  stopifnot(is.data.table(dat))
  if (var1 == var2) {
    stop("'var1' and 'var2' have to differ.")
  }

  dflt_col <- default_flag

  ars <- ar_compare(
    dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2)))], var1, var2, default_flag = dflt_col
  )[, c(4:5, 7)]

  tmp1 <- copy(dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2))), c(var1, dflt_col), with = FALSE])
  setorderv(tmp1, c(var1), order = c(-1))
  tmp1[, perc_sample := .I/.N]
  tmp1[, perc_default := cumsum(get(eval(dflt_col)))/sum(get(eval(dflt_col)))]

  tmp2 <- copy(dat[!is.na(get(eval(var1))) & !is.na(get(eval(var2))), c(var2, dflt_col), with = FALSE])
  setorderv(tmp2, c(var2), order = c(-1))
  tmp2[, perc_default := cumsum(get(eval(dflt_col)))/sum(get(eval(dflt_col)))]

  tmp <- tmp1[, c("perc_sample", "perc_default")]
  tmp[, `Model 2` := tmp2$perc_default]

  rm(tmp1, tmp2)

  colnames(tmp)[1:2] <- c("PercSample", "Model 1")
  tmp <- melt(tmp, id.vars = 1, variable.name = "Model")
  tmp[, Model := data.table::fifelse(Model == "Model 1", var1, var2)]
  tmp[, Model := factor(Model, levels = c(var1, var2))]

  if (is.null(lbl)) {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) +
      geom_line() +
      xlab("Percentage of Sample Excluded") +
      ylab("Percentage of Defaulters Excluded") +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      ggtitle(
        "CAP Plots",
        paste0("AR 1 = ", round(ars[, 1], 3), ", AR 2 = ", round(ars[, 2], 3), "; p-value = ", round(ars[, 3], 4))
      ) +
      theme_minimal()
  } else {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) +
      geom_line() +
      xlab("Percentage of Sample Excluded") +
      ylab("Percentage of Defaulters Excluded") +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      ggtitle(
        lbl,
        paste0("AR 1 = ", round(ars[, 1], 3), ", AR 2 = ", round(ars[, 2], 3), "; p-value = ", round(ars[, 3], 4))
      ) +
      theme_minimal()
  }
}

compare_cap_plot_n <- function(dat, vars, default_flag = "dumdef1", lbl = NULL) {
  stopifnot(is.data.table(dat))
  dflt_col <- default_flag
  ids <- complete.cases(dat[, vars, with = FALSE])

  ars <- paste0(
    paste0(
      vars, ": ",
      round(
        unlist(dat[, lapply(vars, function(x) pdutils::AR_vec(get(eval(x)), get(eval(dflt_col)), na.rm = TRUE))]), 3
      )
    ),
    collapse = ", "
  )

  lst <- list()
  for (i in seq_along(vars)) {
    lst[[i]] <- copy(dat[ids, c(vars[i], dflt_col), with = FALSE])
    setorderv(lst[[i]], vars[i], order = -1)
    lst[[i]][, perc_sample := .I/.N]
    lst[[i]][, perc_default := cumsum(get(eval(dflt_col)))/sum(get(eval(dflt_col)))]
  }

  .tmp <- lst[[1]][, .(perc_sample, `Model 1` = perc_default)]
  colnames(.tmp)[2] <- vars[1]
  tmp1 <- Reduce(cbind, lapply(lst[-1], function(z) z[, "perc_default"]))
  colnames(tmp1) <- vars[-1]
  .tmp <- cbind(.tmp, tmp1)

  .tmp <- melt(.tmp, id.vars = 1, variable.name = "Model")
  .tmp[, Model := factor(Model, levels = vars)]

  if (is.null(lbl)) {
    ggplot(.tmp, aes(x = perc_sample, y = value, colour = Model, group = Model)) +
      geom_line() +
      xlab("Percentage of sample excluded") +
      ylab("Percentage of defaulters excluded") +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      ggtitle("CAP plots", paste0("ARs: ", ars)) +
      theme_minimal()
  } else {
    ggplot(.tmp, aes(x = perc_sample, y = value, colour = Model, group = Model)) +
      geom_line() +
      xlab("Percentage of sample excluded") +
      ylab("Percentage of defaulters excluded") +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      ggtitle(lbl, paste0("ARs: ", ars)) +
      theme_minimal()
  }
}

cap_plot_data <- function(dat, var1, default_flag = "dumdef1") {
  stopifnot(is.data.table(dat))
  dflt_col <- default_flag

  ars <- AR(dff = dat, pd_name = var1, default_flag = dflt_col)

  tmp1 <- copy(dat[!is.na(get(eval(var1))), c(var1, dflt_col), with = FALSE])
  setorderv(tmp1, c(var1), order = c(-1))
  tmp1[, perc_sample := .I/.N]
  tmp1[, perc_default := cumsum(get(eval(dflt_col)))/sum(get(eval(dflt_col)))]

  ndef <- tmp1[, sum(get(eval(dflt_col)))]

  tmp1 <- tmp1[, c("perc_sample", "perc_default")]
  tmp1[, Perfect := .I / ndef]
  tmp1[, Perfect := pmin(Perfect, 1.0)]

  return(as.data.frame(tmp1))
}

cap_plot <- function(dat, var1, default_flag = "dumdef1", col1 = NULL, col2 = NULL, lbl = NULL) {
  stopifnot(is.data.table(dat))

  if (is.null(col1)) {
    colour1 <- "#F8766D"
  } else {
    colour1 <- col1
  }

  if (is.null(col2)) {
    colour2 <- "#00BFC4"
  } else {
    colour2 <- col2
  }

  arv <- AR(dff = dat, pd_name = var1, default_flag = default_flag)
  tmp <- cap_plot_data(dat, var1 = var1, default_flag = default_flag)
  setDT(tmp)

  colnames(tmp)[1:2] <- c("PercSample", var1)
  tmp <- melt(tmp, id.vars = 1, variable.name = "Model")
  tmp[, Model := factor(Model, levels = c(var1, "Perfect"))]

  if (is.null(lbl)) {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) +
      scale_colour_manual(values = c(colour1, colour2)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "black") +
      geom_line() +
      xlab("Percentage of Sample Excluded") +
      ylab("Percentage of Defaulters Excluded") +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      ggtitle("CAP Plot", paste0("AR = ", round(arv, 3))) +
      theme_minimal()
  } else {
    ggplot(tmp, aes(x = PercSample, y = value, colour = Model, group = Model)) +
      scale_colour_manual(values = c(colour1, colour2)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "black") +
      geom_line() +
      xlab("Percentage of Sample Excluded") +
      ylab("Percentage of Defaulters Excluded") +
      scale_x_continuous(labels = scales::percent) +
      scale_y_continuous(labels = scales::percent) +
      ggtitle(lbl, paste0("AR = ", round(arv, 3))) +
      theme_minimal()
  }
}

## Weight of evidence
woe <- function(dff, var, default_flag = "dumdef1") {
  stopifnot(is.data.table(dff))

  dflt_col <- default_flag
  tmp <- copy(dff[, c(var, dflt_col), with = FALSE])

  if (!all(sort(unique(tmp[, get(eval(dflt_col))])) == c(0, 1))) {
    stop("The default_flag column must only contain 0S and 1s.")
  }

  allobs <- nrow(tmp)
  defaults <- as.numeric(tmp[, sum(get(eval(dflt_col)))])

  tst <- setNames(tmp[, .N, by = list(get(eval(var)), get(eval(dflt_col)))], c("variable", dflt_col, "count"))
  tst[, categ := data.table::fifelse(get(eval(dflt_col)) == 0, "good", "bad")]
  tst <- dcast(tst, variable ~ categ, value.var = "count")

  tst[, perc_bad := bad / defaults]
  tst[, perc_good := good / (allobs - defaults)]

  tst[,
    eval(c("perc_good", "perc_bad")) := lapply(
      .SD, function(z) data.table::fifelse(is.na(z), 0, z)
    ), .SDcols = c("perc_good", "perc_bad")
  ]

  tst[, WOE := log(perc_bad) - log(perc_good)]
  tst[, WOE := data.table::fifelse(!is.finite(WOE), 0, WOE)]

  return(as.data.frame(setNames(tst[, .(variable, WOE)], c(var, "woe"))))
}

## Information value
info_value <- function(dff, var, default_flag = "dumdef1") {
  stopifnot(is.data.table(dff))

  dflt_col <- default_flag
  tmp <- copy(dff[, c(var, dflt_col), with = FALSE])

  if(!all(sort(unique(tmp[, get(eval(dflt_col))])) == c(0, 1))) {
    stop("The default_flag column must only contain 0S and 1s.")
  }

  allobs <- nrow(tmp)
  defaults <- as.numeric(tmp[, sum(get(eval(dflt_col)))])

  tst <- setNames(tmp[, .N, by = list(get(eval(var)), get(eval(dflt_col)))], c("variable", dflt_col, "count"))
  tst[, categ := data.table::fifelse(get(eval(dflt_col)) == 0, "good", "bad")]
  tst <- dcast(tst, variable ~ categ, value.var = "count")

  tst[, perc_bad := bad / defaults]
  tst[, perc_good := good / (allobs - defaults)]

  tst[,
    eval(c("perc_good", "perc_bad")) := lapply(
      .SD, function(z) data.table::fifelse(is.na(z), 0, z)
    ), .SDcols = c("perc_good", "perc_bad")
  ]

  tst[, WOE := log(perc_bad) - log(perc_good)]
  tst[, WOE := data.table::fifelse(!is.finite(WOE), 0, WOE)]

  tst[, iv := (perc_bad - perc_good) * WOE]

  return(sum(tst$iv))
}

triangle_plot <- function(dat, var1, var2, pd1, pd2, default_flag = "dumdef1") {
  stopifnot(is.data.table(dat))

  dflt_col <- default_flag

  tmp <- copy(dat[, c(var1, var2, pd1, pd2, dflt_col), with = FALSE])
  tmp <- tmp[complete.cases(tmp)]

  var1_bin <- "var1_bin"
  var2_bin <- "var2_bin"

  tmp[,
    eval(var1_bin) := .bincode(
      get(eval(var1)),
      breaks = quantile(get(eval(var1)), probs = c(0, 0.25, 0.5, 0.75, 1)),
      right = TRUE,
      include.lowest = TRUE
    )
  ]

  tmp[,
    eval(var2_bin) := .bincode(
      get(eval(var2)),
      breaks = quantile(get(eval(var2)), probs = c(0, 0.25, 0.5, 0.75, 1)),
      right = TRUE,
      include.lowest = TRUE
    )
  ]

  sumstat <- tmp[
    , .(
      actDR = mean(get(eval(dflt_col))),
      pd2DR = mean(get(eval(pd2))),
      pd1DR = mean(get(eval(pd1)))
    ),
    keyby = .(var1_bin, var2_bin)
  ]

  sumstat <- melt(sumstat, measure.vars = 3:5)
  sumstat[, variable := factor(variable, levels = c("actDR", "pd2DR", "pd1DR"))]
  setorder(sumstat, var1_bin, var2_bin, variable)
  sumstat[, id := paste(var1_bin, var2_bin, variable, sep = "_")]
  sumstat <- sumstat[, 4:5]

  positions <- data.table(
    id = rep(sumstat$id, each = 3),
    x = c(rep(c(
      0, 0, 1, 0, 0.5, 1, 0.5, 1, 1
    ), 4), rep(c(
      0, 0, 1, 0, 0.5, 1, 0.5, 1, 1
    ) + 1, 4), rep(c(
      0, 0, 1, 0, 0.5, 1, 0.5, 1, 1
    ) + 2, 4), rep(c(
      0, 0, 1, 0, 0.5, 1, 0.5, 1, 1
    ) + 3, 4)),
    y = c(
      c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1),
      c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 1,
      c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 2,
      c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 3,
      c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1),
      c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 1,
      c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 2,
      c(0, 1, 1, 0, 0.5, 0, 0.5, 0, 1) + 3
    )
  )

  datapoly <- merge(sumstat, positions, by = "id")
  colnames(datapoly)[2] <- "PD"
  datapoly[, x := x / 4]
  datapoly[, y := y / 4]
  datapoly[, PD := PD * 100]
  datapoly[, x_c := sum(x) / 3, by = id]
  datapoly[, y_c := sum(y) / 3, by = id]

  return(
    ggplot(datapoly) +
      geom_polygon(aes(x = x, y = y, fill = PD, group = id)) +
      geom_text(aes(x = x_c, y = y_c, label = round(PD, 2))) +
      scale_y_continuous(labels = scales::percent) +
      scale_x_continuous(labels = scales::percent) +
      theme(legend.position = "none") +
      xlab(paste0("Percentiles of ", var1)) +
      ylab(paste0("Percentiles of ", var2)) +
      ggtitle(
        "Actual Default Rate and Predicted PDs", paste0("Top left: actual, right-hand side: ", pd1, ", bottom: ", pd2)
      ) +
      scale_fill_gradient(low = "#0099f7", high = "#f11712")
  )
}
