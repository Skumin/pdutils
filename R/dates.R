# Date functions
yearmon <- function(dates) {
  mnths <- data.table::month(dates)
  yrs <- data.table::year(dates)
  return(as.integer(data.table::fifelse(mnths <= 9, stringi::stri_c(yrs, "0", mnths), stringi::stri_c(yrs, mnths))))
}

end_of_month <- function(dates) {
  if (class(dates) != "Date") {
    stop("dates must be of class 'Date'.")
  }
  return(lubridate::ceiling_date(dates, "months") - 1)
}

start_of_month <- function(dates) {
  if (class(dates) != "Date") {
    stop("dates must be of class 'Date'.")
  }

  return(lubridate::floor_date(dates, "months"))
}

end_of_quarter <- function(dates) {
  if (class(dates) != "Date") {
    stop("dates must be of class 'Date'.")
  }
  yr <- data.table::year(dates)
  qrt <- data.table::quarter(dates)

  return(
    end_of_month(
      anytime::anydate(
        paste0(yr, "-", data.table::fcase(qrt == 1, "03", qrt == 2, "06", qrt == 3, "09", default = "12"), "-01")
      )
    )
  )
}
