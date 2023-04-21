bcr <- function(nobs, ndefs, rho, tau, periods, ci = 0.7, simulations = 1000) {
  if (!class(nobs) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(ndefs) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(rho) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(tau) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(periods) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(ci) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(simulations) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (
    any(c(length(nobs), length(ndefs), length(rho), length(tau), length(periods), length(ci), length(simulations)) != 1)
  ) {
    stop("All arguments must be of length 1.")
  }

  if (ndefs > nobs) {
    stop("nobs must be greater than ndefs.")
  }

  if (ndefs <= 0 | nobs <= 0) {
    stop("nobs and ndefs must be greater than 0.")
  }

  if (nobs %% 1 != 0 | ndefs %% 1 != 0) {
    stop("The number of defaults and the number of observations must be integers.")
  }

  if (periods %% 1 != 0 | simulations %% 1 != 0) {
    stop("The number of periods and the numbers of simulations must be integers.")
  }

  if (periods < 1 | simulations < 1) {
    stop("periods and simulations must be greater than or equal to 1.")
  }

  if (class(rho) != "numeric" | class(tau) != "numeric" | class(ci) != "numeric") {
    stop("rho, tau, and ci must be floats.")
  }

  if (!data.table::between(rho, 0., 1., incbounds = TRUE) | !data.table::between(tau, 0., 1., incbounds = TRUE)) {
    stop("rho and tau must be between 0 and 1.")
  }

  if (ci <= 0 | ci > 1) {
    stop("ci must be greater than 0 and less than or equal to 1.")
  }

  return(pt_multi_pd(nobs, ndefs, rho, tau, periods, ci, simulations))
}

pluto_tasche <- function(nobs, ndefs, rho, tau, periods, ci = 0.7, simulations = 1000) {
  # From the highest risk to the lowest risk
  if (length(nobs) != length(ndefs)) {
    stop("nobs and ndefs must have the same length.")
  }

  if (length(nobs) < 2 | length(ndefs) < 2) {
    stop("nobs and ndefs must be vectors.")
  }

  if (!class(nobs) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(ndefs) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(rho) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(tau) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(periods) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(ci) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (!class(simulations) %in% c("numeric", "integer")) {
    stop("All arguments must be numeric.")
  }

  if (any(c(length(rho), length(tau), length(periods), length(ci), length(simulations)) != 1)) {
    stop("All arguments must be of length 1.")
  }

  if (any(ndefs < 0) | any(nobs <= 0)) {
    stop("All elements of nobs must be greater than 0 and of ndefs greater than or equal to 0.")
  }

  if (sum(ndefs) < 1) {
    stop("There must be at least one default.")
  }

  if (any(nobs %% 1 != 0) | any(ndefs %% 1 != 0)) {
    stop("The numbers of defaults and the numbers of observations must be integers.")
  }

  if (periods %% 1 != 0 | simulations %% 1 != 0) {
    stop("The number of periods and the numbers of simulations must be integers.")
  }

  if (periods < 1 | simulations < 1) {
    stop("periods and simulations must be greater than or equal to 1.")
  }

  if (class(rho) != "numeric" | class(tau) != "numeric" | class(ci) != "numeric") {
    stop("rho, tau, and ci must be floats.")
  }

  if (!data.table::between(rho, 0., 1., incbounds = TRUE) | !data.table::between(tau, 0., 1., incbounds = TRUE)) {
    stop("rho and tau must be between 0 and 1.")
  }

  if (ci <= 0 | ci > 1) {
    stop("ci must be greater than 0 and less than or equal to 1.")
  }

  # A basic check to see if the default rate is decreasing on average. In future better checks may be needed
  if (mean(diff(ndefs/nobs)) > 0) {
    warning(
      "The default rate is, on average, increasing, when moving through grades. Make sure nobs and ndefs are sorted from the worst (i.e., highest risk) to the best grade."
    )
  }

  # Run the C++ function here, after all checks have been performed
  return(pt_multi_pd_full(nobs, ndefs, rho, tau, periods, ci, simulations))
}
