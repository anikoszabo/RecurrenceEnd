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

# combine list of stepfuns into functions with pointwise CLs
#' @importFrom stats isoreg
get_limits <- function(funlist, times, conf.level){
  lwr <- numeric(length(times))
  upr <- numeric(length(times))
  ci_quants <- c((1-conf.level)/2, (1+conf.level)/2)
  for (idx in seq_along(times)){
    fun_values <- sapply(funlist, function(f)f(times[idx]))
    qs <- quantile(fun_values, ci_quants)
    lwr[idx] <- qs[1]
    upr[idx] <- qs[2]
  }
  lwr <- c(1, lwr)
  upr <- c(1, upr)
  # ensure limits are non-increasing
  if (any(diff(lwr) > 0)) {
    isolwr <- isoreg(times, -lwr)
    lwr <- -isolwr$y
  }
  if (any(diff(upr) > 0)) {
    isoupr <- isoreg(times, -upr)
    upr <- -isoupr$y
  }
  list(lower = stepfun(times, lwr),
       upper = stepfun(times, upr))

}

