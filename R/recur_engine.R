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

  predfun_logSurv <- function(cox_model, newdata, ...){

    # Create  model matrix
    fla <- cox_model$formula
    X <- stats::model.matrix(surv_fla[-2], data = newdata)
    X <- X[, -1, drop=FALSE]  # Drop intercept

    # Extract random effects
    frailty_term <- attr(cox_model$terms, "specials")$frailty #includes response in counting
    frailty_col <- cox_model$assign[[frailty_term-1]]
    random_effects <- cox_model$frail[X[,frailty_col]]
    X <- X[, -frailty_col, drop=FALSE]

    # linear predictor
    beta <- cox_model$coefficients
    beta[is.na(beta)] <- 0 # protect from weird NA coefficient from coxph
    if (is.null(beta)){  # if no predictors
      lin_pred <- random_effects
    } else {
      lin_pred <- X %*% beta + random_effects
    }

    # create step function for baseline survival
    cum.base <- survival::basehaz(cox_model, centered = FALSE)
    cum.base$surv <- exp(-cum.base$hazard)
    s0_fun <- stats::stepfun(x = c(0,cum.base$time),
                             y = c(1, cum.base$surv, utils::tail(cum.base$surv,1)*1e-10))

    resfun <- function(i, times, ...){
      log(S0_fun(times)) * exp(lin_pred[i])
    }

    resfun
  }


}
