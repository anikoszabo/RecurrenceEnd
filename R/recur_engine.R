#' Construct a recurrent-event prediction engine
#'
#' @description
#' Builds a **classed engine object** that encapsulates how predictions
#' (e.g., survival, cumulative hazard, hazard) are obtained for the
#' recurrent-event estimators in \pkg{RecurrenceEnd}. The class of the returned
#' object (e.g., `"recur_engine_coxph"`) controls S3 dispatch in downstream
#' generics like [recur_fit()] and [recur_predictfun()].
#'
#' @param name Character scalar naming the engine (e.g., `"coxph"`, `"known_cox"`).
#'   The object will be tagged with class `"recur_engine_<name>"`.
#' @param model Optional **pre-fitted model object** to embed in the engine.
#'   If supplied, [recur_fit()] may skip fitting and reuse this object.
#' @param ... Optional additional fields to store in the engine.
#'
#' @return An S3 object of class `c(paste0("recur_engine_", name), "recur_engine")`
#'   containing (at least) the elements:
#'   \describe{
#'     \item{`name`}{— engine name}
#'     \item{`model`}{— a pre-fitted model object or `NULL`}
#'   }
#'
#'
#' @examples
#' \dontrun{
#' eng <- recur_engine(name = "coxph", config = list(centered = FALSE))
#' class(eng)  # "recur_engine_coxph" "recur_engine"
#' }
#' @seealso [recur_fit()],[recur_predictfun()], [estimate_end()]
#

recur_engine <- function(name, model = NULL, ...) {
  structure(
    list(name = name, model = model, ...),
    class = c(paste0("recur_engine_", name), "recur_engine")
  )
}


#' Fit a model for a recurrent-event prediction engine
#'
#' @description
#' Generic to fit the underlying model required for predictions. The first argument must be a classed engine
#' created by [recur_engine()].
#'
#' @param engine A classed engine object returned by [recur_engine()]; e.g.,
#'   class `"recur_engine_coxph"`.
#' @param formula A one-sided model formula with predictors (` ~ covariates`).
#' @param data A `data.frame` containing variables referenced by `formula` and the following variables:
#' \describe{
#'   \item{time1:}{the beginning of time segments;}
#'   \item{time2:}{the end of time segments;}
#'   \item{id:}{subject identifiers;}
#'   \item{event:}{event indicator;}
#'   \item{terminal:}{indicator of terminal event}
#'   }
#' See \code{\link[reda]{Recur-class}} for more information
#' @param ... Additional arguments passed to the model
#'
#' @return The same `engine` object, typically with its `model` field filled
#'   by a fitted object.
#'
#' @details The default method adds a `model` field containing only the `formula`.
#'
#' @examples
#' \dontrun{
#' eng <- recur_engine("coxph")
#' eng <- recur_fit(eng,  ~ Z1 + Z2, data = SimulatedData)
#' }
#' @seealso [recur_predict()], [estimate_end()]
#' @export
#' @rdname recur_engine
#'
recur_fit <- function(engine, formula, data, ...) {
  UseMethod("recur_fit")
}

#' @exportS3Method
#' @rdname recur_engine
recur_fit.default <- function(engine, formula, data, ...) {
  res <- engine
  res$model <- list(formula = formula)
  res
}

#' Create function that makes gap-time predictions for recurrent-event models
#'
#' @description
#' Generic to produce time-indexed predictions required by the estimators: survival `S(y|X, t)`,
#'  cumulative hazard `H(y|X, t)`, or hazard `h(y|X, t)`, for gap-time `y` after
#'  an event at time `t` for an observation with covariates `X`.
#'
#'
#' @param engine A classed engine object created by [recur_engine()] and
#'   usually fitted by [recur_fit()].
#' @param newdata A `data.frame` with covariates for all subjects to be predicted
#' @param eventtimes Numeric vector with the given event time `t` for each observation in `newdata`
#' @param type Character; one of `"survival"`, `"cumhaz"`, or `"hazard"`.
#' @param log Logical; whether the log-transformed value is required
#'
#' @return a function that will make predictions on `newdata` from the following inputs
#'  \describe{
#'   \item{index}{Integer; indicates the observation number from `newdata` for prediction.}
#'   \item{gaptimes}{Numeric vector of gap-times `y`.}
#'  }
#'
#' @examples
#' \dontrun{
#' predfun <- recur_predict(eng, newdata = SimulatedData[1:5, ], eventtimes = 2,
#'                          type = "survival", log=TRUE)
#' predfun(2, gaptimes =seq(0, 2, by = 0.1))
#' }
#' @seealso [recur_fit()], [estimate_end()]
#' @export
#' @rdname recur_engine
recur_predictfun <- function(
  engine,
  newdata,
  eventtimes = NULL,
  type = c("survival", "cumhaz", "hazard"),
  log = FALSE
) {
  UseMethod("recur_predictfun")
}
