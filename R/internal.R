#'
set_NA_to <- function(x, value=0){
  x[is.na(x)] <- value
  x
}

survfit_to_survfun <- function(sf){
  res <- stepfun(sf$time, c(1, sf$surv))
  res$lower <- sf$lower
  res$upper <- sf$upper
  class(res) <- c("survfun", class(res))
  res
}



