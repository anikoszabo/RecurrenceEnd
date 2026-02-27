# Define object of class 'nptilt' for exponential tilting of
# the unobserved end-time model
# assumes a row per patient with Tmax and C
# lS_fun(i,times) should give log-survival for subject i
# covs_end gives the covariate matrix for the tilting model
#' @importFrom stats stepfun model.matrix quantile
#' @importFrom survival basehaz
#' @importFrom utils tail

nptilt <- function(time, censor, terminal, lS_fun, covs_end, weights = NULL) {
  if (
    (length(time) != length(censor)) |
      length(time) != length(terminal) |
      length(time) != nrow(covs_end)
  ) {
    stop("Arguments should have equal length")
  }
  if (!is.function(lS_fun)) {
    stop("lS_fun should be a function")
  }

  if (is.null(weights)) {
    weights <- rep(1, length(time))
  } else {
    if (length(time) != length(weights)) {
      stop("Arguments should have equal length")
    }
  }
  ak <- sort(unique(time))

  res <- list(
    time = time,
    censor = censor,
    terminal = terminal,
    lS_fun = lS_fun,
    covs_end = covs_end,
    ak = ak,
    weights = weights
  )

  class(res) <- "nptilt"
  res
}


nptilt_from_engine <- function(engine, pred_data, formula_end, weights = NULL) {
  lsurv_fun <- recur_predictfun(
    engine,
    newdata = pred_data,
    eventtimes = pred_data$time1,
    type = "survival",
    log = TRUE
  )
  # Create  model matrix for ending time
  X_end <- stats::model.matrix(formula_end, data = pred_data)
  X_end <- X_end[, -1, drop = FALSE] # Drop intercept

  # Create 'nptilt' object
  mod_nptilt <- nptilt(
    time = pred_data$time1,
    censor = pred_data$time2,
    terminal = pred_data$terminal,
    lS_fun = lsurv_fun,
    covs_end = X_end,
    weights = weights
  )

  mod_nptilt
}

#' @exportS3Method nspmix::length nptilt
length.nptilt <- function(x) length(x$time)

#' @exportS3Method nspmix::weight nptilt
weight.nptilt <- function(x, beta) x$weights

#' @exportS3Method nspmix::suppspace nptilt
suppspace.nptilt = function(x, beta) c(-Inf, Inf)

#' @exportS3Method nspmix::gridpoints nptilt
gridpoints.nptilt = function(x, beta, grid = 100) {
  tms <- x$time
  pp = unique(quantile(tms, seq(0, 1, len = grid), type = 1))
  pp
}


#' @exportS3Method nspmix::initial nptilt
initial.nptilt = function(x, beta = NULL, mix = NULL, kmax = NULL) {
  if (is.null(beta)) {
    beta <- rep(0, ncol(x$covs_end))
  }
  tms <- x$ak
  if (is.null(mix) || is.null(mix$pt)) {
    mi = min(tms)
    ma = max(tms)
    if (is.null(kmax)) {
      pp = unique(tms)
    } else {
      if (kmax == 1) {
        return(list(beta = beta, mix = nspmix::disc(ma)))
      }

      pp = unique(stats::quantile(
        tms,
        seq(0, 1, len = min(20, kmax)),
        type = 1
      ))
    }

    mix0 = nspmix::disc(pp)
    class(x) <- "npkm"
    fit = cnm(
      x,
      init = list(mix = mix0),
      model = "proportions",
      plot = "null",
      verbose = FALSE
    )
    mix = fit$mix
  }
  list(beta = beta, mix = mix)
}


#' @exportS3Method nspmix::logd nptilt
# logd with tilt-sum denominator separated into llex
# does not need the mixture
logd.nptilt = function(x, beta, pt, which = c(1, 0, 0)) {
  dl = vector("list", 3)
  names(dl) = c("ld", "db", "dt")
  n = length(x$time)
  k = length(pt)
  res <- matrix(NA, nrow = n, ncol = k)

  omega <- c(x$covs_end %*% beta)

  eval_pts <- outer(x$censor, pt, pmin) - matrix(x$time, nrow = n, ncol = k)
  for (i in 1:n) {
    neg <- eval_pts[i, ] < 0
    after_terminal <- x$terminal[i] & (pt >= x$censor[i])
    finite <- !neg & !after_terminal

    res[i, !finite] <- -Inf
    res[i, finite] <- x$lS_fun(i, gaptimes = eval_pts[i, finite]) +
      pt[finite] * omega[i]
    if (sum(finite) == 0) {
      # no pt is x$time to x$cens interval and x$terminal==1
      res[i, after_terminal] <- -100
    }
  }

  if (which[1] == 1) {
    dl$ld = res
  }
  if (which[2] == 1) {
    # derivative of log-density term wrt beta
    dbres = array(dim = c(n, k, length(beta)))
    for (j in 1:length(beta)) {
      dbres[,, j] <- outer(x$covs_end[, j], pt)
    }
    dl$db = dbres
  }
  if (which[3] == 1) {
    # not needed
    dl$dt = matrix(NA, nrow = n, ncol = k)
  }
  dl
}

#' @exportS3Method nspmix::llex nptilt
# extra term for log-likelihood with tilt-sum denominators
llex.nptilt = function(x, beta, mix) {
  omega <- c(x$covs_end %*% beta)
  sum_tilt <- drop(exp(outer(omega, mix$pt)) %*% mix$pr)
  res <- -sum(weight(x, beta) * log(sum_tilt))
  res
}

#' @exportS3Method nspmix::llexdb nptilt
# derivative of llex w/respect to beta
llexdb.nptilt = function(x, beta, mix) {
  omega <- c(x$covs_end %*% beta)
  terms <- outer(omega, mix$pt)
  sum_tilt <- drop(exp(terms) %*% mix$pr)
  db_sum_tilt = matrix(NA, nrow = length(x), ncol = length(beta))
  for (j in 1:length(beta)) {
    db_sum_tilt[, j] <- (exp(terms) * outer(x$covs_end[, j], mix$pt)) %*% mix$pr
  }

  res <- -colSums(weight(x, beta) * db_sum_tilt / sum_tilt)
  res
}
