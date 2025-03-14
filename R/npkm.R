# Define object of class 'npkm' for our unobserved end-time model
# assumes a row per patient with Tmax and C
# and a proportional hazards model so S_i = S0_i ^ (exp(lin_pred))
# and S0 - a function that generates gap-time survival probs for this patient
#' @importFrom stats stepfun model.matrix quantile
#' @importFrom survival basehaz
#' @importFrom utils tail

#' @export
npkm <- function(time, censor, lin_pred, S0, weights = NULL){
  if ((length(time) != length(censor)) | length(time) != length(lin_pred)) {
    stop("Arguments should have equal length")
  }
  if (!is.function(S0)){
    stop("S0 should be a function")
  }

  if (is.null(weights)){
    weights <- rep(1, length(time))
  } else {
    if (length(time) != length(weights)) {
      stop("Arguments should have equal length")
    }
  }
  ak <- sort(unique(time))

  res <- list(time = time, censor = censor, lin_pred = lin_pred, S0fun = S0,
              ak = ak, weights = weights)

  class(res) <- "npkm"
  res
}

npkm_from_mod <- function(trail_dat, cox_model, weights = NULL){
  # compute baseline hazard and linear predictor for trailing gap
  # Create  model matrix
  fla <- cox_model$formula
  X <- stats::model.matrix(fla[-2], data = trail_dat)
  X <- X[, -1, drop=FALSE]  # Drop intercept

  # Extract random effects - assumes Cox model has frailty term
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

  # Create 'npkm' object
  mod_npkm <- npkm(time=trail_dat$time1, censor=trail_dat$time2,
                   lin_pred=drop(lin_pred), s0_fun, weights = weights)

  mod_npkm
}

npkm_known_S <- function(trail_dat, formula, S0, coefs, weights = NULL){
  # compute linear predictor for trailing gap
  # Create  model matrix
  X <- stats::model.matrix(formula[-2], data = trail_dat)
  X <- X[, -1, drop=FALSE]  # Drop intercept

  # linear predictor
  if (is.null(coefs)){  # if no predictors
    lin_pred <- 0
  } else {
    lin_pred <- X %*% coefs
  }

  # Create 'npkm' object
  mod_npkm <- npkm(time=trail_dat$time1, censor=trail_dat$time2,
                   lin_pred=drop(lin_pred), S0, weights = weights)

  mod_npkm
}

#' @exportS3Method nspmix::length npkm
length.npkm <- function(x) length(x$time)

#' @exportS3Method nspmix::weight npkm
weight.npkm <- function(x, beta) x$weights


# lower and upper bounds on the support points
#' @exportS3Method nspmix::suppspace npkm
suppspace.npkm = function(x, beta) c(0, Inf)

#' @exportS3Method nspmix::gridpoints npkm
gridpoints.npkm = function(x, beta, grid=100) {
  tms <- x$time
  pp = unique(quantile(tms, seq(0, 1, len=grid),
                       type = 1))
  pp
}


# mix - discrete density on
#' @exportS3Method nspmix::initial npkm
initial.npkm = function(x, beta=NULL, mix=NULL, kmax=NULL) {
  tms <- x$time
  if(is.null(mix) || is.null(mix$pt)) {
    mi = min(tms)
    ma = max(tms)
    if (is.null(kmax)){
      pp = unique(tms)
    } else {
      if (kmax == 1)
        return(list(beta=beta,  mix=nspmix::disc(ma)))

      pp = unique(stats::quantile(tms, seq(0, 1, len=min(20, kmax)),
                           type=1))
    }

    mix = nspmix::disc(pp)
  }
  list(beta=beta, mix=mix)
}


#' @exportS3Method nspmix::logd npkm
logd.npkm = function(x, beta, pt, which=c(1,0,0)) {
  dl = vector("list", 3)
  names(dl) = c("ld","db","dt")
  n = length(x$time)
  k = length(pt)
  res <- matrix(NA, nrow=n, ncol = k)

  eval_pts <- outer(x$cens, pt, pmin) - matrix(x$time, nrow=n, ncol=k)
  for (i in 1:n){
    neg <- eval_pts[i,] < 0
    finite <- !neg

    res[i, !finite] <- -Inf
    res[i, finite] <- log(x$S0(eval_pts[i,finite])) * exp(x$lin_pred[i])
  }

  if(which[1] == 1) {
    dl$ld = res
  }
  if(which[3] == 1) {
    # not implemented, hopefully not needed
    dl$dt = matrix(NA, nrow=n, ncol = k)
  }
  dl
}
