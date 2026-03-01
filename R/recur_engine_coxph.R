#' Fit a Cox frailty model for the coxph engine
#' @description
#' \code{recur_fit.recur_engine_coxph} is the S3 method for fitting a Cox proportional hazards model with frailty as the prediction
#' engine. It uses \pkg{survival} to fit the model and stores the
#' fitted object in \code{engine$model}.
#'
#' @examples
#' eng <- recur_engine("coxph")
#' eng <- recur_fit(eng, ~ Z.1 + Z.2,data = SimulatedData_recur, ties = "efron")
#' @export
#' @method recur_fit recur_engine_coxph
#' @rdname recur_fit
recur_fit.recur_engine_coxph <- function(engine, formula, data, ...) {
  # fit the model
  surv_fla <- stats::update(
    formula,
    survival::Surv(time2 - time1, event) ~ . + frailty(id, sparse=TRUE)
  )
  environment(surv_fla) <- list2env(data)
  mod <- survival::coxph(surv_fla, data = data, ...)

  engine$model <- mod
  engine
}

#' Create a gap-time prediction function for the Cox frailty model engine
#' @description \code{recur_predictfun.recur_engine_coxph} is the S3 method for
#' creating a prediction function based on a fitted Cox proportional hazards model
#' with frailty.
#'
#' @examples
#' eng <- recur_engine("coxph") |>
#'   recur_fit(formula= ~ Z.1 + Z.2,data = SimulatedData_recur)
#' nd <- data.frame(id = 1:4, Z.1 = c(0,0,1,1), Z.2 = c(0, 0.1, 0.2, 0.4))
#' survfun <- recur_predictfun(eng, newdata = nd, eventtimes = 2,
#'            type = "survival")
#' # predict for row 2 in 'nd'
#' survfun(2, gaptimes =seq(0, 2, by = 0.1))
#'
#' @export
#' @rdname recur_predictfun

recur_predictfun.recur_engine_coxph <- function(
  engine,
  newdata,
  eventtimes = NULL,
  type = c("survival", "cumhaz", "hazard"),
  log = FALSE
) {
  model <- engine$model

  type <- match.arg(type)

  # Create  model matrix
  fla <- model$formula
  if (length(fla) == 3) {
    fla <- fla[-2]
  } # drop lhs
  X <- stats::model.matrix(fla, data = newdata)
  X <- X[, -1, drop = FALSE] # Drop intercept

  # Extract random effects
  frailty_term <- attr(model$terms, "specials")$frailty #includes response in counting
  frailty_col <- model$assign[[frailty_term - 1]]
  random_effects <- model$frail[X[, frailty_col]]
  X <- X[, -frailty_col, drop = FALSE]

  # linear predictor
  beta <- model$coefficients
  beta[is.na(beta)] <- 0 # protect from weird NA coefficient from coxph
  if (is.null(beta)) {
    # if no predictors
    lin_pred <- random_effects
  } else {
    lin_pred <- X %*% beta + random_effects
  }

  # create step function for baseline survival
  cum.base <- survival::basehaz(model, centered = FALSE)
  H0_fun <- stats::stepfun(cum.base$time, c(0, cum.base$hazard))

  resfun <- function(index, gaptimes) {}
  if (type == "cumhaz") {
    body(resfun) <- if (log) {
      quote(H0_fun(gaptimes) * exp(lin_pred[index]))
    } else {
      quote(log(H0_fun(gaptimes)) + exp(lin_pred[index]))
    }
  } else if (type == "survival") {
    body(resfun) <- if (log) {
      quote(-H0_fun(gaptimes) * exp(lin_pred[index]))
    } else {
      quote(exp(-H0_fun(gaptimes) * exp(lin_pred[index])))
    }
  } else if (type == "hazard") {
    body(resfun) <- quote(error("Not implemented"))
  }
  environment(resfun) <- list2env(list(H0_fun = H0_fun, lin_pred = lin_pred))

  resfun
}
