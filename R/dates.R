# Date functions
yearmon <- function(dates) {
  mnths <- month(dates)
  yrs <- year(dates)
  return(as.integer(fifelse(mnths <= 9, stri_c(yrs, "0", mnths), stri_c(yrs, mnths))))
}

end_of_month <- function(dates) {
  if (class(dates) != "Date") {
    stop("dates must be of class 'Date'.")
  }
  return(ceiling_date(dates, "months") - 1)
}

start_of_month <- function(dates) {
  if (class(dates) != "Date") {
    stop("dates must be of class 'Date'.")
  }

  return(floor_date(dates, "months"))
}

end_of_quarter <- function(dates) {
  if (class(dates) != "Date") {
    stop("dates must be of class 'Date'.")
  }
  yr <- year(dates)
  qrt <- quarter(dates)

  return(
    end_of_month(
      anydate(
        paste0(yr, "-", fcase(qrt == 1, "03", qrt == 2, "06", qrt == 3, "09", default = "12"), "-01")
      )
    )
  )
}
