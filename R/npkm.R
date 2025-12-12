# Define object of class 'npkm' for our unobserved end-time model
# assumes a row per patient with Tmax and C
# lS_fun(i,times) should give log-survival for subject i
#' @importFrom stats stepfun model.matrix quantile
#' @importFrom survival basehaz
#' @importFrom utils tail

npkm <- function(time, censor, terminal, lS_fun, weights = NULL){
  if ((length(time) != length(censor)) | length(time) != length(terminal)) {
    stop("Arguments should have equal length")
  }
  if (!is.function(predfun)){
    stop("predfun should be a function")
  }

  if (is.null(weights)){
    weights <- rep(1, length(time))
  } else {
    if (length(time) != length(weights)) {
      stop("Arguments should have equal length")
    }
  }
  ak <- sort(unique(time))

  res <- list(time = time, censor = censor, terminal = terminal, lS_fun=lS_fun,
              ak = ak,  weights = weights)

  class(res) <- "npkm"
  res
}

npkm_from_engine <- function(data, pred_formula, rows_to_predict, engine, weights = NULL){

  fit_data <- data[-rows_to_predict, ,drop=FALSE]
  mod <- engine$fit(formula = pred_formula, data=fit_data)

  pred_data <- data[rows_to_predict, ,drop=FALSE]
  predfun <- engine$predfun_logSurv(fit_obj = mod, newdata = pred_data)
  # Create 'npkm' object
  mod_npkm <- npkm(time=pred_dat$time1, censor=pred_dat$time2,
                   terminal=pred_dat$terminal,
                   lS_fun = predfun,
                   weights = weights)

  mod_npkm
}

# temporary function that makes fewer changes than npkm_from_engine
npkm_from_mod <- function(pred_data, model, engine, weights = NULL){
  # assumes model was fitted by engine

  predfun <- engine$predfun_logSurv(fit_obj = model, newdata = pred_data)
  # Create 'npkm' object
  mod_npkm <- npkm(time=pred_data$time1, censor=pred_data$time2,
                   terminal=pred_data$terminal,
                   lS_fun = predfun,
                   weights = weights)

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
  tms <- x$ak
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
    after_terminal <- x$terminal[i] & (pt >= x$cens[i])
    finite <- !neg &  !after_terminal

    res[i, !finite] <- -Inf
    res[i, finite] <- x$lS_fun(i, times=eval_pts[i,finite])
    if (sum(finite) == 0){# no pt is x$time to x$cens interval and x$terminal==1
       res[i, after_terminal] <- -100
    }

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
