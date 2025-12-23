
#' @title Fit a Cox PH with frailty engine
#' @description
#' S3 method for fitting a Cox proportional hazards model with frailty as the prediction
#' engine. Uses \pkg{survival} to fit the model and stores the
#' fitted object in `engine$model`.
#'
#' @inheritParams recur_fit
#' @param ties, robust, model, x, y Optional arguments forwarded to
#'   [survival::coxph()] (if used internally).
#'
#' @return The updated engine object of class `"recur_engine_coxph"` with
#'   a fitted model embedded and any engine-specific metadata recorded.
#'
#' @examples
#' \dontrun{
#' eng <- recur_engine("coxph")
#' eng <- recur_fit(eng, Recur(time = time, id = patient.id, event = indicator) ~ Z1 + Z2,
#'                  data = SimulatedData, ties = "efron")
#' }
#' @seealso [survival::coxph()], [recur_predict.recur_engine_coxph()]
##' @export
##' @method recur_fit recur_engine_coxph
recur_fit.recur_engine_coxph <- function(engine, formula, data, ...){

  # fit the model
  surv_fla <- stats::update(formula, survival::Surv(time2-time1, event) ~ . + frailty(id))
  environment(surv_fla) <- list2env(data)
  mod <- survival::coxph(surv_fla, data=data, ...)

  engine$model <- mod
  engine
}



#' @title Create a gap-time prediction funcition for a Cox PH with frailty engine
#' @inheritParams recur_predictfun
#' @return A function(index, gaptimes).
#' @examples
#' \dontrun{
#' survfun <- recur_predict(eng, newdata = SimulatedData[1:5, ], eventtime = 2,
#'            type = "survival")
#' survfun(2, gaptimes =seq(0, 2, by = 0.1), type = "survival")
#' }
#' @export

recur_predictfun.recur_engine_coxph <- function(engine, newdata,  eventtimes=NULL,
                                                type=c("survival", "cumhaz",  "hazard"), log=FALSE) {
  model <- engine$model

  type <- match.arg(type)

  # Create  model matrix
  fla <- model$formula
  if (length(fla) == 3) {fla <- fla[-2]} # drop lhs
  X <- stats::model.matrix(fla, data = newdata)
  X <- X[, -1, drop=FALSE]  # Drop intercept

  # Extract random effects
  frailty_term <- attr(model$terms, "specials")$frailty #includes response in counting
  frailty_col <- model$assign[[frailty_term-1]]
  random_effects <- model$frail[X[,frailty_col]]
  X <- X[, -frailty_col, drop=FALSE]

  # linear predictor
  beta <- model$coefficients
  beta[is.na(beta)] <- 0 # protect from weird NA coefficient from coxph
  if (is.null(beta)){  # if no predictors
    lin_pred <- random_effects
  } else {
    lin_pred <- X %*% beta + random_effects
  }

  # create step function for baseline survival
  cum.base <- survival::basehaz(model, centered = FALSE)
  H0_fun <- stats::stepfun(cum.base$time, c(0, cum.base$hazard))

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

