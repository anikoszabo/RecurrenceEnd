#'
set_NA_to <- function(x, value=0){
  x[is.na(x)] <- value
  x
}

survfit_to_survfun <- function(sf){
  stepfun(sf$time, c(1, sf$surv))
}


survfitCI_to_survfun <- function(sf){
  list(lower = stepfun(sf$time, c(1, sf$lower)),
       upper = stepfun(sf$time, c(1, sf$upper)))
}

