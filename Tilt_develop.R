library(devtools)
library(nspmix)
document()
load_all()

mix_to_stepfun <- function(mix) {
  res <- stats::stepfun(mix$pt, 1 - cumsum(c(0, mix$pr)))
  res
}


tilt_mix <- function(mix, omega) {
  tmix <- mix
  tmix$pr <- exp(omega * mix$pt) * mix$pr
  tmix
}


tilt_mix_free <- function(mix, beta0, omega) {
  tmix <- mix
  tmix$pr <- exp(beta0 + omega * mix$pt) * mix$pr
  tmix
}


logd_nptilt = function(x, beta, mix, which = c(1, 0, 0)) {
  # needs mixture
  pt <- mix$pt
  dl = vector("list", 3)
  names(dl) = c("ld", "db", "dt")
  n = length(x$time)
  k = length(pt)
  res <- matrix(NA, nrow = n, ncol = k)

  omega <- c(x$covs_end %*% beta)

  eval_pts <- outer(x$censor, pt, pmin) - matrix(x$time, nrow = n, ncol = k)
  for (i in 1:n) {
    tilt_sum <- sum(mix$pr * exp(omega[i] * pt))
    neg <- eval_pts[i, ] < 0
    after_terminal <- x$terminal[i] & (pt >= x$censor[i])
    finite <- !neg & !after_terminal

    res[i, !finite] <- -Inf
    res[i, finite] <-
      x$lS_fun(i, gaptimes = eval_pts[i, finite]) +
      pt[finite] * omega[i] -
      log(tilt_sum)

    if (sum(finite) == 0) {
      # no pt is in the x$time to x$cens interval and x$terminal==1
      res[i, after_terminal] <- -100
    }
  }

  if (which[1] == 1) {
    dl$ld = res
  }
  if (which[2] == 1) {
    # not implemented yet
    dl$db = matrix(NA, nrow = n, ncol = k)
  }
  if (which[3] == 1) {
    # not needed
    dl$dt = matrix(NA, nrow = n, ncol = k)
  }
  dl
}

loglik_nptilt = function(mix, x, beta = NULL, attr = FALSE) {
  ld = logd_nptilt(x, beta, mix, which = c(1, 0, 0))$ld
  ma = matMaxs(ld)
  dmix = drop(exp(ld - ma) %*% mix$pr) + 1e-100
  logd = log(dmix) + ma
  ll = sum(weight(x, beta) * logd)

  if (attr) {
    attr(ll, "dmix") = dmix
    attr(ll, "logd") = logd # log(mixture density)
  }
  ll
}

# simulate example data with two groups
ld1 <- 1
ld2 <- 2
lr0 <- 5
ss <- 0.1
n1 <- n2 <- 150
aa1 <- sim_recur_end(
  n = n1,
  lambda = ld1,
  lambda_r = lr0,
  sigma2 = ss,
  lambda_c = 1,
  C_min = 1,
  censored.ok = FALSE
)
aa1$Z <- 0
aa2 <- sim_recur_end(
  n = n2,
  lambda = ld2,
  lambda_r = lr0,
  sigma2 = ss,
  lambda_c = 1,
  C_min = 1,
  censored.ok = FALSE
)
aa2$patient.id <- n1 + aa2$patient.id
aa2$Z <- 1
aa <- rbind(aa1, aa2)

a <- as.data.frame(
  with(aa, Recur(time = time, id = patient.id, event = indicator))@.Data
)
a$Z <- aa$Z
a$Z2 <- rnorm(nrow(a)) # an extra covariate with coef = 0

trail_a <- a[a$event == 0, ]

# fit 1-sample curve
ks <- recur_engine("known", H0fun = function(x) lr0 * x) |>
  recur_fit(~1, data = subset(a, event == 1))
npkm_a <- npkm_from_engine(pred_data = trail_a, engine = ks)

mix0 <- disc(npkm_a$ak)
fit_a <- cnm(
  npkm_a,
  init = list(mix = mix0),
  model = "proportions",
  plot = "null",
  verbose = FALSE
)
mix1 <- fit_a$mix


# two-sample explorations
nptilt_a <- nptilt_from_engine(
  pred_data = trail_a,
  engine = ks,
  formula_end = ~Z
)

init1 <- initial(nptilt_a)
mix1 <- init1$mix
loglik_nptilt(mix = mix1, x = nptilt_a, beta = c(0))
# same as
loglik(mix1, npkm_a)

loglik_nptilt(mix = mix1, x = nptilt_a, beta = c(-0.2))
curve(
  Vectorize(loglik_nptilt, vectorize.args = "beta")(
    mix = mix1,
    x = nptilt_a,
    beta = x
  ),
  from = -1.5,
  to = 1
)

# should be same with llex approach
loglik(mix = mix1, x = nptilt_a, beta = c(-0.2))

# check derivative wrt beta
nspmix:::dll0(x = nptilt_a, mix = mix1, beta = -0.2, which = c(1, 0, 0, 1))
eps <- 0.001
(loglik(mix = mix1, x = nptilt_a, beta = c(-0.2 + eps)) -
  loglik(mix = mix1, x = nptilt_a, beta = c(-0.2 - eps))) /
  (2 * eps)
# extra term
llexdb(mix = mix1, x = nptilt_a, beta = c(-0.2))
(llex(mix = mix1, x = nptilt_a, beta = c(-0.2 + eps)) -
  llex(mix = mix1, x = nptilt_a, beta = c(-0.2 - eps))) /
  (2 * eps)


# second derivative
U <- nspmix:::dll0(x = nptilt_a, mix = mix1, beta = 0, which = c(0, 0, 0, 1))$db
I <- -(nspmix:::dll0(
  x = nptilt_a,
  mix = mix1,
  beta = -0 + eps,
  which = c(0, 0, 0, 1)
)$db -
  nspmix:::dll0(
    x = nptilt_a,
    mix = mix1,
    beta = -0 - eps,
    which = c(0, 0, 0, 1)
  )$db) /
  (2 * eps)
U / sqrt(I)


## Two predictors
a$Z2 <- rnorm(nrow(a)) # an extra covariate with coef = 0

nptilt_b <- nptilt_from_engine(
  pred_data = trail_a,
  engine = ks,
  formula_end = ~ Z + Z2
)

loglik(mix = mix1, x = nptilt_b, beta = c(-0.2, 0))
# same as
loglik(mix = mix1, x = nptilt_a, beta = c(-0.2))

# check derivative wrt beta2
nspmix:::dll0(
  x = nptilt_b,
  mix = mix1,
  beta = c(-0.2, 0),
  which = c(1, 0, 0, 1)
)
eps <- 0.001
(loglik(mix = mix1, x = nptilt_b, beta = c(-0.2, eps)) -
  loglik(mix = mix1, x = nptilt_b, beta = c(-0.2, -eps))) /
  (2 * eps)

## Second derivative?

###########
# try optimization
# profile likelihood? probably not, as llex is ignored?
r0 <- nspmix:::pll(
  mix0 = mix1,
  x = nptilt_a,
  beta = c(-0),
  model = "proportions"
)
r1 <- nspmix:::pll(
  mix0 = mix1,
  x = nptilt_a,
  beta = c(-0.7),
  model = "proportions"
)
r2 <- nspmix:::pll(
  mix0 = mix1,
  x = nptilt_a,
  beta = c(0.7),
  model = "proportions"
)


beta_vec <- seq(-0.7, -0.6, by = 0.005)
pll_vec <- sapply(beta_vec, function(b) {
  nspmix:::pll(mix = mix1, x = nptilt_a, beta = b, model = "proportions")$ll
})
plot(beta_vec, pll_vec, type = "l")

### BFGS?
# very unstable, b/c too many dimensions? but sometimes gives good result
# maybe modifying matrix D would help?
# works well if only beta is modified
b1 <- nspmix:::bfgs(mix = mix1, beta = 0, x = nptilt_a, which = c(0, 0, 1))
b2 <- nspmix:::bfgs(mix = mix1, beta = -0.1, x = nptilt_a, which = c(1, 0, 1))
b3 <- nspmix:::bfgs(mix = mix1, beta = 0.1, x = nptilt_a, which = c(0, 0, 1))
b4 <- nspmix:::bfgs(
  mix = mix1,
  beta = -0.2,
  x = nptilt_a,
  which = c(1, 0, 1),
  D = diag(c(rep(-1, length(mix1$pt)), -1))
)
c(b1$beta, b2$beta, b3$beta, b4$beta)
c(b1$ll, b2$ll, b3$ll, b4$ll)

##### fit non-parametric model to both groups
npkm_a <- npkm_from_engine(ks, trail_a)

ks1 <- recur_engine("known", H0fun = function(x) lr0 * x) |>
  recur_fit(~1, data = subset(a, event == 1 & Z == 0))
npkm_a1 <- npkm_from_engine(ks1, trail_a[trail_a$Z == 0, ])

ks2 <- recur_engine("known", H0fun = function(x) lr0 * x) |>
  recur_fit(~1, data = subset(a, event == 1 & Z == 1))
npkm_a2 <- npkm_from_engine(ks1, trail_a[trail_a$Z == 1, ])


mix0 <- disc(npkm_a$ak)
fit_a <- cnm(
  npkm_a,
  init = list(mix = mix0),
  model = "proportions",
  plot = "null",
  verbose = FALSE
)
mix1 <- fit_a$mix

mix01 <- disc(npkm_a1$ak)
fit_a1 <- cnm(
  npkm_a1,
  init = list(mix = mix01),
  model = "proportions",
  plot = "null",
  verbose = FALSE
)
mix11 <- fit_a1$mix

mix02 <- disc(npkm_a2$ak)
fit_a2 <- cnm(
  npkm_a2,
  init = list(mix = mix02),
  model = "proportions",
  plot = "null",
  verbose = FALSE
)
mix12 <- fit_a2$mix

loglik(mix1, npkm_a)
loglik(mix12, npkm_a2) + loglik(mix11, npkm_a1)

plot(mix_to_stepfun(mix1))
lines(mix_to_stepfun(mix11), col = 2)
lines(mix_to_stepfun(mix12), col = 3)

plot(mix_to_stepfun(mix1))
lines(mix_to_stepfun(tilt_mix(mix1, omega = 0.5)), col = 3)
lines(mix_to_stepfun(tilt_mix(mix1, omega = -0.5)), col = 4)
