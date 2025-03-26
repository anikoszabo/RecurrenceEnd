
#' Estimate the distribution of the ending event of a recurrent event process
#'
#' @param formula a two-sided formula with a \code{\link[reda]{Recur}} object specifying the
#' recurrent event process on the left-hand side, and the predictors of its intensity
#' on the right-hand side. Origin and terminal events in \code{Recur} are not respected.
#' @param method character string specifying the estimation method
#' @param threshold optional numeric value, specifies the threshold for the "threshold" method.
#' Defaults to 0, which is equivalent to the "naive" method.
#' @param quantile optional numeric value between 0 and 1 (exclusive) specifying the
#' quantile used for the "quantile" method.
#' @param data an optional data frame containing the variables in the formula
#' @param subset an optional vector specifying a subset of observations to be used
#' @param na.action a function which indicates what should happen when the data
#' contains NAs. The default is set by the na.action setting of \code{options},
#' which in turn defaults to \code{na.omit}.
#' @param verbose logical value, if TRUE then information is displayed during computation
#' @param known_recur optional list specifying the known recurrence survival distribution
#' function for each subject. It should have two components: \code{S0} - the baseline survival
#' function, and \code{coefs} - the coefficients of the recurrence model specified in \code{formula}
#' @param IPSW logical, indicating whether inverse-probability of selection weights based on
#' no censoring before first event time should be used
#' @param bootCI logical, indicating whether pointwise bootstrap confidence intervals should
#' be computed. Use carefully for the "NPMLE" method, as it may take long.
#' @param conf.level numeric, confidence level for pointwise confidence interval
#' @param bootB integer, number of bootstrap replicates
#' @return a list with components \code{fit} - an object of class \link[stats]{stepfun} with the
#' estimated survival distribution, \code{method} - the method used, and additional
#' information based on the method used
#' @export
#' @importFrom stats stepfun model.frame model.extract update quantile
#' @importFrom survival survfit coxph Surv

estimate_end <- function(formula,
                    method=c("naive", "threshold","quantile", "NPMLE"),
                    threshold = 0, quantile=0.95, data, subset, na.action,
                    verbose=FALSE, known_recur=NULL, IPSW = FALSE,
                    bootCI = FALSE, conf.level = 0.95, bootB = 100){
  method <- match.arg(method)

  # based on reda::rateReg
  if (missing(formula))
    stop("Argument 'formula' is required.")
  if (missing(data))
    data <- environment(formula)
  if (!missing(subset)) {
    sSubset <- substitute(subset)
    subIdx <- eval(sSubset, data, parent.frame())
    if (!is.logical(subIdx))
      stop("'subset' must be logical")
    subIdx <- subIdx & !is.na(subIdx)
    data <- data[subIdx, ]
  }
  mcall <- match.call(expand.dots = FALSE)
  mmcall <- match(c("formula", "data", "na.action"), names(mcall),
                  0L)
  mcall <- mcall[c(1L, mmcall)]
  mcall$data <- data
  mcall$drop.unused.levels <- TRUE
  mcall[[1L]] <- quote(stats::model.frame)
  mf <- eval(mcall, parent.frame())
  resp <- stats::model.extract(mf, "response")
  if (!reda::is.Recur(resp))
    stop("Response in the formula must be a 'Recur' object.")

  if (length(na.action <- attr(mf, "na.action"))) {
    if (verbose)
      message("Observations with missing value in covariates ",
              "are removed.\nChecking the new dataset again...\n",
              appendLF = FALSE)
    resp <- reda::check_Recur(resp, check = "hard")
  }

  Dat <- as.data.frame(resp@.Data)
  # arrange data into nice order
  ord <- resp@ord
  Dat <- Dat[ord,]

  # create indicators for largest event time, latest time-point (censoring time)
  # allow for no-event intervals within the recurrent event process
  events <- Dat$event==1
  nevents <- tapply(events, Dat["id"], sum)
  if (any(nevents == 0)){
    stop("All subjects must to have at least one observed recurrent event")
  }
  last_event_idx <- tapply(seq_len(nrow(Dat)), Dat["id"],
                           function(x)max(x[events[x]]))
  last_event <- Dat$time2[last_event_idx]
  last_time <- Dat$time2[resp@last_idx]  # censoring times

  if (IPSW){
    first_event_idx <- tapply(seq_len(nrow(Dat)), Dat["id"],
                             function(x)min(x[events[x]]))
    first_event <- Dat$time2[first_event_idx]
    G_cens <- function(t){mean(t <= last_time)}
    p_select <- sapply(first_event, G_cens)
    weights <- 1/p_select
    Dat$.weights <- weights[Dat[["id"]]]
  } else {
    weights <- NULL
    Dat$.weights <- 1
  }

  if (method == "naive"){
    res <- survival::survfit(survival::Surv(last_event) ~ 1, weights = weights,
                             conf.int = conf.level)
    output <- list(fit = survfit_to_survfun(res),
                ci = survfitCI_to_survfun(res),
                method = method)
  }
  if (method == "threshold"){
    if (!isTRUE(threshold >= 0)) stop("'threshold' should be a positive number")
    res <- survival::survfit(survival::Surv(last_event, last_time-last_event >= threshold ) ~ 1,
                             weights = weights,  conf.int = conf.level)
    output <- list(fit = survfit_to_survfun(res),
                ci = survfitCI_to_survfun(res),
                method = method,
                threshold = threshold)
  }
  if (method == "quantile"){
    if (!isTRUE(quantile > 0 & quantile < 1))
      stop("'quantile' should be between 0 and 1 (exclusive)")
    # estimate gap-time distribution to recurrent event part
    gap_fit <- survival::survfit(survival::Surv(time2 - time1, event) ~ 1,
                                 data = Dat[-resp@last_idx,],
                                 weights=Dat[-resp@last_idx,]$.weights)
    thresh <- stats::quantile(gap_fit, probs = quantile, conf.int = FALSE)
    res <- survival::survfit(survival::Surv(last_event, last_time-last_event >= thresh ) ~ 1,
                             weights = weights,  conf.int = conf.level)

    output <- list(fit = survfit_to_survfun(res),
                ci = survfitCI_to_survfun(res),
                method = method,
                quantile = quantile,
                threshold = thresh)
  }

  if (method == "NPMLE"){
    # add predictors to Dat - everything on the RHS of 'formula'
    covars <- all.vars(formula[-2])
    if (any(names(Dat) %in% covars))
      stop(paste("Please avoid using '", toString(names(Dat)), "' as predictor names in 'formula'"))

    if (length(na.action)>0 & inherits(na.action, "omit")){
      Dat <- cbind(Dat, data[-na.action, ][ord,covars,drop=FALSE])
    } else {
     Dat <- cbind(Dat, data[ord,covars,drop=FALSE])
    }

   # check for only one 'trailing' event
    next_to_last <- (resp@last_idx - 1)[resp@last_idx > resp@first_idx]
    if (!all(events[next_to_last]))
      stop("There should be only one censored interval after the last event")

    tmax_unique <- sort(unique(last_event))

    # data from 'trailing gaps'
    trailDat <- Dat[resp@last_idx,]

    # Fit gap-time model to recurrent event part
    surv_fla <- stats::update(formula, survival::Surv(time2-time1, event) ~ . + frailty(id))
    environment(surv_fla) <- list2env(Dat)

    if (!is.null(known_recur)){
      mod_npkm <- npkm_known_S(trail_dat = trailDat, formula = formula,
                               S0 = known_recur$S0, coefs = known_recur$coefs,
                               weights = weights)
      mod <- NULL
    } else {
      # ignore warning that is due to approx zero estimate for frailty variance
      suppressWarnings(
        mod <- survival::coxph(surv_fla, data=Dat[-resp@last_idx,])
        #mod <- my_coxph(surv_fla, data=Dat[-resp@last_idx,], weights=.weights)
      )

      # Create 'npkm' object
      mod_npkm <- npkm_from_mod(trail_dat = trailDat, cox_model = mod,
                                weights = weights)
    }
    # fit NPMLE
    npmix_fit <- nspmix::cnm(mod_npkm, init=list(mix=nspmix::disc(mod_npkm$ak)),
                     model="proportions", plot="null", verbose=verbose)

    res <- stats::stepfun(npmix_fit$mix$pt, 1-cumsum(c(0, npmix_fit$mix$pr)))
    output <- list(fit = res,
                method = method,
                 model = mod)
  }

  if (bootCI){
    if (method %in% c("threshold", "naive", "quantile")){
      thresh <- switch(method, naive = 0, threshold = threshold,
                       quantile = NA)
      boot_times <- res$time # step-points of original fit
      bootres <- list()
      for (b in 1:bootB){
        idx <- sample.int(length(last_event), replace = TRUE)
        if (method == "quantile"){
          last_event_id <- Dat$id[last_event_idx]
          Dat_events_b <- merge(Dat[-resp@last_idx,],
                                data.frame(id = last_event_id[idx]))
          # re-estimate quantile
          gap_fit_b <- survival::survfit(survival::Surv(time2 - time1, event) ~ 1,
                                     data = Dat_events_b,
                                     weights=Dat_events_b$.weights)
          thresh <- stats::quantile(gap_fit_b, probs = quantile, conf.int = FALSE)
        }

        res_b <- survival::survfit(
          survival::Surv(last_event[idx], last_time[idx]-last_event[idx] >= thresh ) ~ 1,
                         weights = weights[idx],  conf.int = conf.level)

        bootres <- c(bootres, list(survfit_to_survfun(res_b)))
      }
    }
    if (method == "NPMLE"){
      boot_times <- npmix_fit$mix$pt # step-points of original fit
      Dat_events <- Dat[-resp@last_idx,]
      bootres <- list()
      for (b in 1:bootB){
        idx <- sample.int(length(last_event), replace = TRUE)
        # select bootstrapped subjects from Dat_events & trailDat, assigning new IDs
        last_event_id <- Dat$id[last_event_idx]
        newid <- seq_along(last_event_id)
        id_map <- data.frame(.origid = last_event_id[idx], id = newid)
        Dat_events_b <- Dat_events
        names(Dat_events_b)[names(Dat_events_b) == "id"] <- ".origid"
        Dat_events_b <- merge(Dat_events_b, id_map, by = ".origid", all.x = FALSE, all.y = TRUE)
        trailDat_b <- trailDat
        names(trailDat_b)[names(trailDat_b) == "id"] <- ".origid"
        trailDat_b <- merge(trailDat_b, id_map, by = ".origid", all.x = FALSE, all.y = TRUE)

        if (!is.null(known_recur)){
          mod_npkm_b <- npkm_known_S(trail_dat = trailDat_b, formula = formula,
                                   S0 = known_recur$S0, coefs = known_recur$coefs,
                                   weights = weights[idx])
        } else {
          # ignore warning that is due to approx zero estimate for frailty variance
          environment(surv_fla) <- list2env(Dat_events_b)
          suppressWarnings(
            mod_b <- survival::coxph(surv_fla, data=Dat_events_b)
          )

          # Create 'npkm' object
          mod_npkm_b <- npkm_from_mod(trail_dat = trailDat_b, cox_model = mod_b,
                                      weights = weights[idx])
        }
        # fit NPMLE
        npmix_fit_b <- nspmix::cnm(mod_npkm_b, init=list(mix=nspmix::disc(mod_npkm_b$ak)),
                                 model="proportions", plot="null", verbose=verbose)

        res_b <- stats::stepfun(npmix_fit_b$mix$pt, 1-cumsum(c(0, npmix_fit_b$mix$pr)))
        bootres <- c(bootres, list(res_b))
      }
    }

    boot_ci <- get_limits(bootres, times = boot_times, conf.level = conf.level)
    output$ci <- boot_ci
  }

  output
}


