#' Methods for 'survfun' objects
#'
#' List of methods that apply to \code{survfun} objects created by the
#' \code{\link{estimate_end}} function.
#'
#' @rdname survfun
#'
#' @param x object of class 'survfun'
#' @param object object of class 'survfun'
#' @param conf.int logical, whether confidence limits / intervals should be shown.
#' Defaults to \code{FALSE}.
#' @param conf.lty integer, linetype for plotted confindence limits
#' @param conf.col color of plotted confindence limits
#' @param ylim numeric vector of length 2 defining the y limits of the plot. Defaults to \code{c(0,1)}.
#' @param do.points logical, specifying whether points should be plotted at the closed end
#' of each step. Defaults to \code{FALSE}.
#' @param ... additional parameters, only used for the \code{plot} method, where they
#' are passed to \code{\link[stats]{plot.stepfun}}
#' @export
#' @importFrom graphics plot
plot.survfun <- function(x, conf.int = FALSE, ylim = c(0,1), conf.lty = 2,
                         conf.col = 1, do.points = FALSE, ...){
  plot(x$fit, do.points=do.points, ylim = ylim, ...)
  if (conf.int & !is.null(x$ci)){
    lines(x$ci$lower, do.points=do.points, col=conf.col, lty=conf.lty)
    lines(x$ci$upper, do.points=do.points, col=conf.col, lty=conf.lty)
  }
}

#' @rdname survfun
#' @importFrom graphics lines
#' @export
lines.survfun <- function(x, conf.int = FALSE, conf.lty = 2,
                         conf.col = 1, ...){
  lines(x$fit, do.points=FALSE,  ...)
  if (conf.int & !is.null(x$ci)){
    lines(x$ci$lower, do.points=FALSE, col=conf.col, lty=conf.lty)
    lines(x$ci$upper, do.points=FALSE, col=conf.col, lty=conf.lty)
  }
}

#' @rdname survfun
#' @param times numeric vector of times (x-values) for which prediction is requested.
#' Defaults to all jump-points of the function.
#' @importFrom stats predict
#' @export
predict.survfun <- function(object, times = NULL, conf.int = FALSE, ...){
  if (is.null(times)){
    times <- get("x", envir = environment(object$fit))
  }
  res <- list(time = times, pred = object$fit(times))
  if (conf.int & !is.null(object$ci)){
    ci <- list(lower = object$ci$lower(times), upper = object$ci$upper(times))
    res <- c(res, ci)
  }
  res
}

#' @rdname survfun
#' @param probs numeric vector of probabilities. Defaults to all quartiles.
#' @importFrom stats quantile
#' @export
quantile.survfun <- function(x, probs = c(0.25, 0.5, 0.75), conf.int=TRUE, ...){
  # quantile of a non-increasing step-function that goes from 1 to 0
  q.step <- function(stepf){
    xx <- get("x", envir = environment(stepf))
    yy <- get("y", envir = environment(stepf))
    # quantiles are of the original distribution, not of survival function
    idx <- findInterval(probs, c(0, 1-yy))
    # values beyond the right tail should be set to NA (not reached)
    qs <- c(xx, NA)[idx]
    names(qs) <- probs
    qs
  }
  res <- list(quantile = q.step(x$fit))
  if (conf.int & !is.null(x$ci)){
    ci <- list(lower = q.step(x$ci$lower), upper = q.step(x$ci$upper))
    res <- c(res, ci)
  }
  res
}

#' @rdname survfun
#' @param na.rm logical, required as part of the generic for \code{median}. Ignored.
#' @importFrom stats median
#' @export
median.survfun <- function(x, na.rm = FALSE, conf.int=TRUE, ...){
  quantile(x, probs = 0.5, conf.int = conf.int, ...)
}

