library(devtools)
library(nspmix)
load_all()

mix_to_stepfun <- function(mix){
  res <- stats::stepfun(mix$pt, 1-cumsum(c(0, mix$pr)))
  res
}

tilt_mix <- function(mix, beta){
  tmix <- mix
  tmix$pr <- exp(beta * mix$pt) * mix$pr
  tmix$pr <- tmix$pr / sum(tmix$pr)
  tmix
}

nptilt <- function(time, censor, terminal, covs_end,
                   lin_pred_recur, S0_recur, weights = NULL){
  if ((length(time) != length(censor)) | length(time) != length(lin_pred_recur) |
      length(time) != length(lin_pred_recur) | length(time) != nrow(covs_end)) {
    stop("Arguments should have equal length")
  }
  if (!is.function(S0_recur)){
    stop("S0_recur should be a function")
  }

  if (is.null(weights)){
    weights <- rep(1, length(time))
  } else {
    if (length(time) != length(weights)) {
      stop("Arguments should have equal length")
    }
  }
  ak <- sort(unique(time))

  res <- list(time = time, censor = censor, terminal = terminal, covs_end = covs_end,
              lin_pred_recur = lin_pred_recur,
              S0fun = S0_recur, ak = ak,  weights = weights)

  class(res) <- "nptilt"
  res
}

nptilt_known_S <- function(trail_dat, formula_end,
                           formula_recur, S0_recur,
                           coefs_recur, weights = NULL){
  # compute linear predictor for trailing gap
  # Create  model matrix for recurrent events
  X <- stats::model.matrix(formula_recur[-2], data = trail_dat)
  X <- X[, -1, drop=FALSE]  # Drop intercept

  # linear predictor
  if (is.null(coefs_recur)){  # if no predictors
    lin_pred_recur <- rep(0, nrow(X))
  } else {
    lin_pred_recur <- X %*% coefs
  }

  # Create  model matrix for ending time
  X_end <- stats::model.matrix(formula_end[-2], data = trail_dat)
  X_end <- X_end[, -1, drop=FALSE]  # Drop intercept

  # Create 'nptilt' object
  res <- nptilt(time=trail_dat$time1, censor=trail_dat$time2,
                terminal=trail_dat$terminal, covs_end = Z_end,
                lin_pred_recur=drop(lin_pred), S0_recur = S0_recur,
                weights = weights)

  mod_npkm
}

# simulate example data with two groups
ld1 <- 1; ld2 <- 2
lr0 <- 5
ss <- 0.1
n1 <- n2 <- 150
aa1 <- sim_recur_end(n=n1, lambda = ld1, lambda_r = lr0, sigma2=ss,
                    lambda_c = 1, C_min = 1, censored.ok = FALSE)
aa1$Z <- 0
aa2 <- sim_recur_end(n=n2, lambda = ld2, lambda_r = lr0, sigma2=ss,
                     lambda_c = 1, C_min = 1, censored.ok = FALSE)
aa2$patient.id <- n1+aa2$patient.id
aa2$Z <- 1
aa <- rbind(aa1, aa2)

a <- as.data.frame(with(aa, Recur(time=time, id=patient.id, event=indicator))@.Data)
a$Z <- aa$Z

trail_a <- a[a$event == 0, ]

npkm_a <- npkm_known_S(trail_a,
                        formula= Recur(time=time, id=patient.id, event=indicator)~ 1,
                        S0 = function(x)pexp(x, rate=lr0, lower=FALSE),
                        coefs = NULL)
npkm_a1 <- npkm_known_S(trail_a[trail_a$Z==0,],
                        formula= Recur(time=time, id=patient.id, event=indicator)~ 1,
                       S0 = function(x)pexp(x, rate=lr0, lower=FALSE),
                       coefs = NULL)
npkm_a2 <- npkm_known_S(trail_a[trail_a$Z==1,],
                        formula= Recur(time=time, id=patient.id, event=indicator)~ 1,
                        S0 = function(x)pexp(x, rate=lr0, lower=FALSE),
                        coefs = NULL)
mix0 <- disc(npkm_a$ak)
fit_a <- cnm(npkm_a, init=list(mix=mix0), model="proportions", plot="null", verbose=FALSE)
mix1 <- fit_a$mix

mix01<- disc(npkm_a1$ak)
fit_a1 <- cnm(npkm_a1, init=list(mix=mix01), model="proportions", plot="null", verbose=FALSE)
mix11 <- fit_a1$mix

mix02<- disc(npkm_a2$ak)
fit_a2 <- cnm(npkm_a2, init=list(mix=mix02), model="proportions", plot="null", verbose=FALSE)
mix12 <- fit_a2$mix

loglik(mix1, npkm_a)
loglik(mix12, npkm_a2) +  loglik(mix11, npkm_a1)

plot(mix_to_stepfun(mix1))
lines(mix_to_stepfun(mix11), col=2)
lines(mix_to_stepfun(mix12), col=3)

plot(mix_to_stepfun(mix1))
lines(mix_to_stepfun(tilt_mix(mix1, beta=0.5)), col=3)
lines(mix_to_stepfun(tilt_mix(mix1, beta=-0.5)), col=4)
