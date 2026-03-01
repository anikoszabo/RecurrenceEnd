#' Construct a recurrent-event prediction engine
#'
#' @description
#' \code{recur_engine} builds a \emph{classed engine object} that encapsulates how predictions
#' (e.g., survival, cumulative hazard, hazard) are obtained for the
#' recurrent-event estimators. The class of the returned
#' object (e.g., \code{recur_engine_coxph}) controls S3 dispatch in downstream
#' generics like \code{\link{recur_fit}} and \code{\link{recur_predictfun}}.
#'
#' @param name Character scalar naming the engine (e.g., "coxph", "known").
#'   The object will be tagged with class \code{recur_engine_<name>}.
#' @param model Optional pre-fitted model object to embed in the engine.
#'   If supplied, \code{recur_fit()} may skip fitting and reuse this object.
#' @param ... Optional additional fields to store in the engine.
#'
#'
#' @return An S3 object of class \code{c(paste0("recur_engine_", name), "recur_engine")}
#'   containing at least the following elements:
#'   \describe{
#'     \item{name}{engine name}
#'     \item{model}{a pre-fitted model object or \code{NULL}}
#'   }
#'
#' @seealso \code{\link{recur_fit}} and \code{\link{recur_predictfun}} for methods that
#' need to be defined for the engine to work, and
#' \code{\link{estimate_end}} for intended use.
#'
#' @examples
#' eng <- recur_engine(name = "coxph", list(centered = FALSE))
#' class(eng)  # "recur_engine_coxph" "recur_engine"
#'
#' @export

recur_engine <- function(name, model = NULL, ...) {
  structure(
    list(name = name, model = model, ...),
    class = c(paste0("recur_engine_", name), "recur_engine")
  )
}


#' Fit a model for a recurrent-event prediction engine
#'
#' @description
#' \code{recur_fit} is a generic to fit the underlying model required for predictions. The first argument must be a classed engine
#' created by \code{recur_engine()}
#'
#' @param engine A classed engine object returned by \code{recur_engine()}; e.g.,
#'   class \code{recur_engine_coxph}.
#' @param formula A one-sided model formula with predictors (\code{~ covariates}).
#' @param data A data.frame containing variables referenced by \code{formula} and the following variables:
#' @param ... Additional arguments passed to the model.
#' \describe{
#'   \item{time1}{the beginning of time segments;}
#'   \item{time2}{the end of time segments;}
#'   \item{id}{subject identifiers;}
#'   \item{event}{event indicator;}
#'   \item{terminal}{indicator of terminal event}
#'   }
#' Typically, it will be created by a call to \code{\link[reda]{Recur}}. See
#' \code{\link[reda]{Recur-class}} for more information.
#'
#' @return The same \code{engine} object, typically with its \code{model}
#'  field filled by a fitted object.
#'
#' @details The default method of \code{recur_fit} adds a \code{model} field
#' containing only the \code{formula} argument passed to it.
#'
#' @examples
#' eng <- recur_engine(name = "coxph")
#' eng <- recur_fit(eng,  ~ Z.1 + Z.2, data = SimulatedData_recur)
#'
#' @seealso \code{\link{recur_engine}} for defining engines, and
#' \code{\link{recur_predictfun}} for an additional method that needs to be defined
#' @export
#' @rdname recur_fit
recur_fit <- function(engine, formula, data, ...) {
  UseMethod("recur_fit")
}

#' @exportS3Method
#' @rdname recur_fit
recur_fit.default <- function(engine, formula, data, ...) {
  res <- engine
  res$model <- list(formula = formula)
  res
}

#' Create function that makes gap-time predictions for recurrent-event models
#'
#' @description
#' \code{recur_predictfun} is a generic to produce time-indexed predictions required by
#' the estimators in \code{\link{estimate_end}}: the survival  \eqn{S(y|X, t)},
#'  cumulative hazard \eqn{H(y|X, t)}, or hazard \eqn{h(y|X, t)} for gap-time \eqn{y} after
#'  an event at time \eqn{t} for an observation with covariates \eqn{X}.
#'
#'
#' @param engine A classed engine object created by \code{\link{recur_engine}} and
#'   usually fitted by \code{\link{recur_fit}}.
#' @param newdata A data.frame with covariates for all subjects to be predicted
#' @param eventtimes Numeric vector with the given event time \eqn{t} for
#' each observation in \code{newdata}.
#' @param type Character; one of "survival", "cumhaz", or "hazard".
#' @param log Logical; whether the log-transformed value is required
#'
#' @return a function that will make predictions on \code{newdata} from the following inputs
#'  \describe{
#'   \item{index}{Integer; indicates the observation number from \code{newdata} for prediction.}
#'   \item{gaptimes}{Numeric vector of gap-times \eqn{y}.}
#'  }
#'
#' @examples
#' eng <- recur_engine(name = "coxph")
#' eng <- recur_fit(eng,  ~ Z.1 + Z.2, data = SimulatedData_recur)
#' predfun <- recur_predictfun(eng, newdata = SimulatedData_recur[1:5, ], eventtimes = 2,
#'                          type = "survival", log=TRUE)
#' predfun(2, gaptimes =seq(0, 2, by = 0.1))
#'
#' @seealso \code{\link{recur_engine}} for defining engines, and
#' \code{\link{recur_fit}} for an additional method that needs to be defined
#' @export
#' @rdname recur_predictfun
recur_predictfun <- function(
  engine,
  newdata,
  eventtimes = NULL,
  type = c("survival", "cumhaz", "hazard"),
  log = FALSE
) {
  UseMethod("recur_predictfun")
}
