# Contains various PD functions
## Annualized cumulative PDs from cumulative PDs
apd <- function(cpds, tenors = NULL) {
  if (any(diff(cpds) <= 0)) {
    stop("Cumulative default probabilities must be increasing.")
  }
  if (any(cpds <= 0)) {
    stop("Cumulative default probabilities must be greater than zero.")
  }
  if (is.null(tenors)) {
    tnrs <- seq_along(cpds)
  } else {
    tnrs <- tenors
  }

  return(1 - (1 - cpds) ^ (1 / tnrs))
}

## Cumulative PDs from annualized cumulative PDs
cpd <- function(apds, tenors = NULL) {
  if (any(apds <= 0)) {
    stop("Annualized default probabilities must be greater than zero.")
  }
  if (is.null(tenors)) {
    tnrs <- seq_along(apds)
  } else {
    tnrs <- tenors
  }

  cpds <- 1 - (1 - apds) ^ tnrs

  if (any(diff(cpds) <= 0)) {
    stop("The annualized default probabilities result in non-increasing cumulative defaut probabilities.")
  } else {
    return(cpds)
  }
}

## Forward PDs from cumulative annualized PDs
fpd <- function(cpds) {
  if (any(diff(cpds) <= 0)) {
    stop("Cumulative default probabilities must be increasing.")
  }

  if (any(cpds <= 0)) {
    stop("Cumulative default probabilities must be greater than zero.")
  }

  fpds <- c(cpds[1], sapply(seq_along(cpds)[-length(cpds)], function(x) 1 - (1 - cpds[x + 1])/(1 - cpds[x])))

  return(fpds)
}
