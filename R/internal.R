#'
set_NA_to <- function(x, value=0){
  x[is.na(x)] <- value
  x
}

survfit_to_survfun <- function(sf){
  stepfun(sf$time, c(1, sf$surv))
}


survfitCI_to_survfun <- function(sf){
  list(lower = stepfun(sf$time, c(1, sf$lower)),
       upper = stepfun(sf$time, c(1, sf$upper)))
}

# combine list of stepfuns into functions with pointwise CLs
#' @importFrom stats isoreg
get_limits <- function(funlist, times, conf.level){
  lwr <- numeric(length(times))
  upr <- numeric(length(times))
  ci_quants <- c((1-conf.level)/2, (1+conf.level)/2)
  for (idx in seq_along(times)){
    fun_values <- sapply(funlist, function(f)f(times[idx]))
    qs <- quantile(fun_values, ci_quants)
    lwr[idx] <- qs[1]
    upr[idx] <- qs[2]
  }
  lwr <- c(1, lwr)
  upr <- c(1, upr)
  # ensure limits are non-increasing
  if (any(diff(lwr) > 0)) {
    isolwr <- isoreg(times, -lwr)
    lwr <- -isolwr$y
  }
  if (any(diff(upr) > 0)) {
    isoupr <- isoreg(times, -upr)
    upr <- -isoupr$y
  }
  list(lower = stepfun(times, lwr),
       upper = stepfun(times, upr))

}


# estimate last-event distribution from end distribution
post_process <- function(npkm_fit, npkm){

  n <- length(npkm)
  D_pts <- npkm_fit$mix$pt  # distribution atoms for D
  D_pr <- npkm_fit$mix$pr  # distribution weights

  # subject-specific distribution of D
  Sik <- exp(logd.npkm(x = npkm, beta=NA, pt = D_pts)$ld)
  numerator <- Sik * matrix(D_pr, nrow=n, ncol=length(D_pts), byrow=TRUE)
  ss_D_pr <- numerator/rowSums(numerator)

  # self-consistency:
  # all.equal(colMeans(ss_D_pr), pr)

  # subject-specific distribution of T*
  T_pts <- npkm$ak  # atoms for the distribution of T*
  ss_T_pr <- matrix(0, nrow=n, ncol=length(T_pts))
  for (i in 1:n){
    # P(T*_i=T_{i,m_i} | D_i = a_k)
    pk <- numeric(length(D_pts))
    pk[D_pts >= np$time[i] & D_pts <= np$censor[i]] <- 1
    eval_pts <- D_pts - np$time[i]
    beyond_censor <- (D_pts > np$censor[i]) & !np$terminal[i]
    pk[beyond_censor] <- np$S0(eval_pts[beyond_censor]) ^ exp(np$lin_pred[i])

    idx_Ti <- which.min(abs(T_pts - np$time[i])) # index of T_{i,m_i} in T_pts

    cond_pr <- matrix(0, nrow = length(T_pts), ncol=length(D_pts))
    cond_pr[idx_Ti, ] <- pk
    if (max(D_pts) > np$censor[i]){
      idx_D_Ci <- min(which(D_pts > np$censor[i])) # index of first observation after C_i in D_pts
      for (k in idx_D_Ci : length(D_pts)){
        idx <- which(T_pts <= D_pts[k] & T_pts > np$censor[i])
        cond_pr[idx, k] <- (1 - pk[k]) / length(idx)
      }
    }
    ss_T_pr[i, ] <- cond_pr %*% ss_D_pr[i,]
  }

  list(pt=T_pts, pr = colMeans(ss_T_pr))
}
