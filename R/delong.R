## Calculates the variance of the AUC using DeLong's algorithm
## This specific version is based on https://doi.org/10.1109/LSP.2014.2337313
delong_auc_variance <- function(pds, default_flag, na.rm = FALSE) {
  stopifnot(class(pds) == "numeric")
  stopifnot(class(default_flag) == "integer" | class(default_flag) == "numeric")

  len <- length(default_flag)
  stopifnot((length(pds) == len) & (len > 1))

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

  if (any(!sort(unique(default_flag)) %in% c(0, 1))) {
    stop("default_flag must only contains 0s and 1s.")
  }

  n_defs <- sum(default_flag)
  n <- len - n_defs

  ordr <- order(default_flag, decreasing = TRUE)
  default_flag <- default_flag[ordr]
  pds <- pds[ordr]

  def_inds <- seq_len(n_defs)
  nondef_inds <- seq.int(from = n_defs + 1, to = len)

  tx <- frank(pds[def_inds], ties.method = "average")
  ty <- frank(pds[nondef_inds], ties.method = "average")
  tz <- frank(pds, ties.method = "average")

  v01 <- (tz[def_inds] - tx) / n
  v10 <- 1 - (tz[nondef_inds] - ty) / n_defs

  sx <- var(v01)
  sy <- var(v10)

  return(sx / n_defs + sy / n)
}

## Calculates the covariance of the AUCs of two or more models using DeLong's algorithm
## Doesn't work well on big datasets due to the explicit call to `cov`
delong_auc_covariance <- function(pds, default_flag) {
  stopifnot("matrix" %in% class(pds))
  stopifnot(class(default_flag) == "integer" | class(default_flag) == "numeric")

  len <- length(default_flag)
  stopifnot((nrow(pds) == len) & (len > 1L))

  if (any(!sort(unique(default_flag)) %in% c(0, 1))) {
    stop("default_flag must only contains 0s and 1s.")
  }

  n_models <- ncol(pds)
  stopifnot(n_models >= 2L)

  n_defs <- sum(default_flag)
  n <- len - n_defs

  temp_frame <- cbind(pds, default_flag)
  temp_frame <- temp_frame[order(temp_frame[, n_models + 1L], decreasing = TRUE), ]

  tx <- matrix(ncol = n_models, nrow = n_defs)
  ty <- matrix(ncol = n_models, nrow = n)
  tz <- matrix(ncol = n_models, nrow = n_defs + n)

  def_inds <- seq_len(n_defs)
  nondef_inds <- seq.int(from = n_defs + 1L, to = len)

  for (i in seq_len(n_models)) {
    tx[, i] <- as.numeric(frank(temp_frame[def_inds, i], ties.method = "average"))
    ty[, i] <- as.numeric(frank(temp_frame[nondef_inds, i], ties.method = "average"))
    tz[, i] <- as.numeric(frank(temp_frame[, i], ties.method = "average"))
  }

  v01 <- (tz[def_inds, ] - tx) / n
  v10 <- 1.0 - (tz[nondef_inds, ] - ty) / n_defs

  sx <- cov(v01)
  sy <- cov(v10)

  return(sx / n_defs + sy / n)
}

## Calculates the (normal) confidence interval around the AUC using DeLong's algorithm and returns all three values
delong_auc_ci <- function(pds, default_flag, conf_level = 0.95, na.rm = FALSE) {
  stopifnot(conf_level > 0 & conf_level < 1)

  auc <- mann_whitney_vec(pds, default_flag, na.rm = na.rm)
  auc_sd <- sqrt(delong_auc_variance(pds, default_flag, na.rm = na.rm))

  return(qnorm(c((1 - conf_level) / 2, 0.5, 1 - (1 - conf_level) / 2), mean = auc, sd = auc_sd))
}

## Calculates the DeLong covariance for exactly 2 classifiers
delong_paired_test_prep <- function(pds1, pds2, default_flag) {
  ## Get the "DeLong placements" first
  m <- sum(default_flag)
  n <- length(default_flag) - m

  def_inds <- default_flag == 1
  nondef_inds <- default_flag == 0

  vr <- delong_placements(pds1[def_inds], pds1[nondef_inds])
  vs <- delong_placements(pds2[def_inds], pds2[nondef_inds])

  ## Check that the AUCs align
  auc1 <- mann_whitney_vec(pds1, default_flag)
  auc2 <- mann_whitney_vec(pds2, default_flag)

  stopifnot(abs(auc1 - vr$theta) <= 1e-16)
  stopifnot(abs(auc2 - vs$theta) <= 1e-16)

  ## Get the covariance matrix
  sx <- matrix(ncol = 2, nrow = 2)
  sx[1, 1] <- sum((vr$X - vr$theta) * (vr$X - vr$theta)) / (m - 1)
  sx[1, 2] <- sum((vr$X - vr$theta) * (vs$X - vs$theta)) / (m - 1)
  sx[2, 1] <- sum((vs$X - vs$theta) * (vr$X - vr$theta)) / (m - 1)
  sx[2, 2] <- sum((vs$X - vs$theta) * (vs$X - vs$theta)) / (m - 1)

  sy <- matrix(ncol = 2, nrow = 2)
  sy[1, 1] <- sum((vr$Y - vr$theta) * (vr$Y - vr$theta)) / (n - 1)
  sy[1, 2] <- sum((vr$Y - vr$theta) * (vs$Y - vs$theta)) / (n - 1)
  sy[2, 1] <- sum((vs$Y - vs$theta) * (vr$Y - vr$theta)) / (n - 1)
  sy[2, 2] <- sum((vs$Y - vs$theta) * (vs$Y - vs$theta)) / (n - 1)

  ## S is the variance-covariance matrix; it is equal to running delong_auc_covariance(cbind(pds1, pds2), default_flag)
  ## directly but faster in big datasets
  S <- sx / m + sy / n
  d <- vr$theta - vs$theta
  return(list(auc1 = auc1, auc2 = auc2, auc_diff = d, S_mat = S))
}

delong_paired_test <- function(pds1, pds2, default_flag, conf_level = 0.95, test_what = "AUC", na.rm = FALSE) {
  stopifnot(conf_level > 0 & conf_level < 1)
  stopifnot(test_what %chin% c("AUC", "AR"))
  stopifnot(class(pds1) == "numeric", class(pds2) == "numeric")
  stopifnot(class(default_flag) == "integer" | class(default_flag) == "numeric")

  len <- length(default_flag)
  stopifnot((length(pds1) == len) & (len > 1), length(pds2) == len)

  dflt_na <- is.na(default_flag)
  pds1_na <- is.na(pds1)
  pds2_na <- is.na(pds2)

  if (any(pds1_na) | any(dflt_na) | any(pds2_na)) {
    if (!na.rm) {
      stop("There are NAs in the three vectors and na.rm is FALSE.")
    } else {
      ids <- which(pds1_na | dflt_na | pds2_na)
      pds1 <- pds1[-ids]
      pds2 <- pds2[-ids]
      default_flag <- default_flag[-ids]
    }
  }

  if (any(!sort(unique(default_flag)) %in% c(0, 1))) {
    stop("default_flag must only contains 0s and 1s.")
  }

  z_vals <- qnorm(c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2))
  prep_data <- delong_paired_test_prep(pds1, pds2, default_flag)

  if (test_what == "AUC") {
    L <- c(1, -1)
    sig <- sqrt(L %*% prep_data$S_mat %*% L)[1, 1]
    conf_int <- prep_data$auc_diff + z_vals * sig
    zscore <- prep_data$auc_diff / sig

    return(
      list(
        Z = zscore, p_value = 2 * pnorm(-abs(zscore)), conf_int = conf_int, metrics = c(prep_data$auc1, prep_data$auc2)
      )
    )

  } else {
    L <- c(2, -2)
    sig <- sqrt(L %*% prep_data$S_mat %*% L)[1, 1]
    conf_int <- 2 * prep_data$auc_diff + z_vals * sig
    zscore <- 2 * prep_data$auc_diff / sig

    return(
      list(
        Z = zscore, p_value = 2 * pnorm(-abs(zscore)), conf_int = conf_int,
        metrics = 2 * c(prep_data$auc1, prep_data$auc2) - 1
      )
    )
  }
}
