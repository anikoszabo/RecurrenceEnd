% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_end.R
\name{estimate_end}
\alias{estimate_end}
\title{Estimate the distribution of the ending event of a recurrent event process}
\usage{
estimate_end(
  formula,
  method = c("naive", "threshold", "quantile", "NPMLE"),
  threshold = 0,
  quantile = 0.95,
  data,
  subset,
  na.action,
  estimand = c("end", "last_event"),
  verbose = FALSE,
  known_recur = NULL,
  IPSW = FALSE,
  bootCI = FALSE,
  conf.level = 0.95,
  bootB = 100
)
}
\arguments{
\item{formula}{a two-sided formula with a \code{\link[reda]{Recur}} object specifying the
recurrent event process on the left-hand side, and the predictors of its intensity
on the right-hand side. Origin and terminal events in \code{Recur} are not respected.}

\item{method}{character string specifying the estimation method}

\item{threshold}{optional numeric value, specifies the threshold for the "threshold" method.
Defaults to 0, which is equivalent to the "naive" method.}

\item{quantile}{optional numeric value between 0 and 1 (exclusive) specifying the
quantile used for the "quantile" method.}

\item{data}{an optional data frame containing the variables in the formula}

\item{subset}{an optional vector specifying a subset of observations to be used}

\item{na.action}{a function which indicates what should happen when the data
contains NAs. The default is set by the na.action setting of \code{options},
which in turn defaults to \code{na.omit}.}

\item{estimand}{character, indicating whether the distribution of the unobserved ending
event should be estimated ('end'), or the distribution of the last recurrent event
before the ending event should be estimated ('last_event')}

\item{verbose}{logical value, if TRUE then information is displayed during computation}

\item{known_recur}{optional list specifying the known recurrence survival distribution
function for each subject. It should have two components: \code{S0} - the baseline survival
function, and \code{coefs} - the coefficients of the recurrence model specified in \code{formula}}

\item{IPSW}{logical, indicating whether inverse-probability of selection weights based on
no censoring before first event time should be used}

\item{bootCI}{logical, indicating whether pointwise bootstrap confidence intervals should
be computed. Use carefully for the "NPMLE" method, as it may take long.}

\item{conf.level}{numeric, confidence level for pointwise confidence interval}

\item{bootB}{integer, number of bootstrap replicates}
}
\value{
an object of class \code{\link{survfun}} with the estimated survival distribution
}
\description{
Estimate the distribution of the ending event of a recurrent event process
}
