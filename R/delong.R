## Calculates the variance of the AUC using DeLong's algorithm
## This specific version is based on https://doi.org/10.1109/LSP.2014.2337313
delong_auc_variance <- function(pds, default_flag, na.rm = FALSE) {
  stopifnot(class(pds) == "numeric")
  stopifnot(class(default_flag) == "integer" | class(default_flag) == "numeric")

  dflt_na <- is.na(default_flag)
  pds_na <- is.na(pds)

  if (any(pds_na) | any(dflt_na)) {
    if (!na.rm) {
      stop("There are NAs in the two columns and na.rm is FALSE.")
    } else {
      ids <- which(pds_na | dflt_na)
      pds <- pds[-ids]
      default_flag <- default_flag[-ids]
    }
  }

  len <- length(default_flag)
  stopifnot((length(pds) == len) & (len > 1))

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

  tx <- data.table::frank(pds[def_inds], ties.method = "average")
  ty <- data.table::frank(pds[nondef_inds], ties.method = "average")
  tz <- data.table::frank(pds, ties.method = "average")

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
    tx[, i] <- as.numeric(data.table::frank(temp_frame[def_inds, i], ties.method = "average"))
    ty[, i] <- as.numeric(data.table::frank(temp_frame[nondef_inds, i], ties.method = "average"))
    tz[, i] <- as.numeric(data.table::frank(temp_frame[, i], ties.method = "average"))
  }

  v01 <- (tz[def_inds, ] - tx) / n
  v10 <- 1.0 - (tz[nondef_inds, ] - ty) / n_defs

  sx <- cov(v01)
  sy <- cov(v10)

  return(sx / n_defs + sy / n)
}

## Calculates the (normal) confidence interval around the AUC using DeLong's algorithm and returns all three values
delong_auc_ci <- function(pds, default_flag, conf_level = 0.95, na.rm = FALSE) {
  auc <- pdutils::mann_whitney_vec(pds, default_flag, na.rm = na.rm)
  auc_sd <- sqrt(delong_auc_variance(pds, default_flag, na.rm = na.rm))

  return(qnorm(c((1 - conf_level) / 2, 0.5, 1 - (1 - conf_level) / 2), mean = auc, sd = auc_sd))
}
