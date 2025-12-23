#' @title Create a gap-time prediction function for a known model engine
#' @inheritParams recur_predictfun
#' @return A function(index, gaptimes).
#' @details
#' The following additional items need to be present in the engine:
#' \describe{
#' \item{`H0fun`}{a function giving the baseline gap-time cumulative hazard function}
#' \item{`coefs`}{optional numeric vector of coefficients used for the Cox model}
#' }
#'
#' @examples
#' \dontrun{
#' ke <- recur_engine("known",  H0fun = function(x) 0.5*x)
#' survfun <- recur_predictfun(ke, newdata = SimulatedData[1:5, ], eventtime = 2,
#'            type = "survival")
#' survfun(2, gaptimes =seq(0, 2, by = 0.1))
#' }
#' @export
#' @rdname recur_engine_known

recur_predictfun.recur_engine_known <- function(engine, newdata,  eventtimes=NULL,
                                                type=c("survival", "cumhaz",  "hazard"), log=FALSE) {
  model <- engine$model

  type <- match.arg(type)

  # Create  model matrix
  fla <- model$formula
  if (length(fla) == 3) {fla <- fla[-2]} # drop lhs
  X <- stats::model.matrix(fla, data = newdata)
  X <- X[, -1, drop=FALSE]  # Drop intercept

  # linear predictor
  if (is.null(engine$coefs)){  # if no predictors
    lin_pred <- rep(0, nrow(X))
  } else {
    lin_pred <- X %*% engine$coefs
  }

  H0_fun <- engine$H0fun
  resfun <- function(index, gaptimes){}
  if (type == "cumhaz"){
    body(resfun) <- if (log){
      quote(H0_fun(gaptimes) * exp(lin_pred[index]))
      } else {
      quote(log(H0_fun(gaptimes)) + exp(lin_pred[index]))
      }
  } else if (type == "survival"){
    body(resfun) <- if (log){
      quote(-H0_fun(gaptimes) * exp(lin_pred[index]))
    } else {
      quote(exp(-H0_fun(gaptimes) * exp(lin_pred[index])))
    }
  }  else if (type == "hazard"){
    body(resfun) <- quote(error("Not implemented"))
  }
  environment(resfun) <- list2env(list(H0_fun=H0_fun, lin_pred=lin_pred))

  resfun
}

