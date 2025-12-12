#' Cox frailty model engine for fitting recurrent event data
#'
#' @return A list with elements:
#'  - fit(formula, data,...) - a function that fits the model.
#'      The data can be assumed to have columns 'id', 'time1', 'time2', 'event',
#'       and all predictors
#'      The formula will be one-sided, stating any predictors.
#'
#'  - predfun_logSurv(fit_obj, newdata, ....) - a function that sets up prediction
#'      for observations in newdata starting from 'time1'. newdata has the same
#'      structure as the dgitata argument above. fit_obj will be the output of a
#'      call to the corresponding fit function. The output should be a function of
#'      (i, times) that predicts the  gap-time log-survival based on the fitted model
#'      for subject i of newdata. Time 0 for 'times' corresponds to 'time1' in newdata.

coxf_engine <- function(){
  fit <- function(formula, data, ...){
    # fit the model
    surv_fla <- stats::update(formula, survival::Surv(time2-time1, event) ~ . + frailty(id))
    environment(surv_fla) <- list2env(data)
    mod <- survival::coxph(surv_fla, data=data)
    mod
  }

  predfun_logSurv <- function(fit_obj, newdata, ...){

    # Create  model matrix
    fla <- fit_obj$formula
    X <- stats::model.matrix(fla[-2], data = newdata)
    X <- X[, -1, drop=FALSE]  # Drop intercept

    # Extract random effects
    frailty_term <- attr(fit_obj$terms, "specials")$frailty #includes response in counting
    frailty_col <- fit_obj$assign[[frailty_term-1]]
    random_effects <- fit_obj$frail[X[,frailty_col]]
    X <- X[, -frailty_col, drop=FALSE]

    # linear predictor
    beta <- fit_obj$coefficients
    beta[is.na(beta)] <- 0 # protect from weird NA coefficient from coxph
    if (is.null(beta)){  # if no predictors
      lin_pred <- random_effects
    } else {
      lin_pred <- X %*% beta + random_effects
    }

    # create step function for baseline survival
    cum.base <- survival::basehaz(fit_obj, centered = FALSE)
    cum.base$lsurv <- -cum.base$hazard
    lS0_fun <- stats::stepfun(x = c(0,cum.base$time),
                             y = c(0, cum.base$lsurv, utils::tail(cum.base$lsurv,1)-10))

    resfun <- function(i, times, ...){
      lS0_fun(times) * exp(lin_pred[i])
    }
    environment(resfun) <- list2env(list(lS0_fun=lS0_fun, lin_pred=lin_pred))

    resfun
  }

  list(fit = fit, predfun_logSurv = predfun_logSurv)
}

#' Cox model with known baseline survival and coefficients engine
#'
#' @param lS0 function(times) that gives baseline log-survival at 'times'
#' @param coefs coefficients of the recurrence model, should match formula
knownS_engine <- function(lS0, coefs=NULL){
  fit <- function(formula, data, ...){
    list(formula = formula)
  }

  predfun_logSurv <- function(fit_obj, newdata, ...){

    # Create  model matrix
    formula <- fit_obj$formula
    X <- stats::model.matrix(formula, data = newdata)

    # linear predictor
    if (is.null(coefs)){  # if no predictors
      lin_pred <- rep(0, nrow(X))
    } else {
      lin_pred <- X %*% coefs
    }

    resfun <- function(i, times, ...){
      lS0(times) * exp(lin_pred[i])
    }
    environment(resfun) <- list2env(list(lS0=lS0, lin_pred=lin_pred))

    resfun
  }

  list(fit = fit, predfun_logSurv = predfun_logSurv)
}
