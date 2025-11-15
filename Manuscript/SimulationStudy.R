# Simulation Study

library(RecurrenceEnd)
library(doParallel)
library(foreach)
library(dplyr)
library(tidyr)

# Helper function for Mean Squared Error (MSE)
compute_mse <- function(bias) {
  return(mean(bias^2))
}

# Helper function for the new bias calculation
compute_bias <- function(quantile_exp, km_function, target_probs) {
  # Evaluate the KM step function at the exponential quantiles
  survival_probs <- sapply(quantile_exp, km_function)

  # Compute the bias as the difference between the KM survival probabilities and target probabilities
  bias <- survival_probs - target_probs
  names(bias) <- paste0(target_probs * 100, "%")
  return(bias)
}

# Main simulation function
# method_list - named list with method name, its parameters, and RHS of formula to use
run_setting <- function(nsim, method_list, quantiles, n, mu_d, EN, p_c, C_min,
                        sigma2, cov_fun=NULL, beta=NULL){
  lambda_d <- log(2)/mu_d
  lambda_r <- EN*lambda_d
  lambda_c <- (lambda_d*p_c)/(exp(-lambda_d*C_min) - p_c)
  true_quantiles <- qexp(quantiles, rate = lambda_d, lower.tail = FALSE)


  res_list <- list()
  for (i in 1: nsim){
    # simulate a dataset without fully censored observations
    observed_data <- sim_recur_end(n = n, lambda_d = lambda_d, lambda_r = lambda_r,
                                sigma2 = sigma2, cov_fun = cov_fun, beta = beta,
                                lambda_c = lambda_c, C_min=C_min, censored.ok = FALSE)

    resmat <- matrix(NA, nrow=length(method_list)+1, ncol = length(quantiles)+2)
    # create version for full-data estimate
    observed_data2 <- observed_data |>
      mutate(minimum = pmin(recurrence_end, C),
             min_indicator = ifelse(recurrence_end <= C, 1, 0)) |>
      select(patient.id, minimum, min_indicator) |>
      unique()

    full_data_km <- survfit(Surv(minimum, min_indicator) ~ 1, data=observed_data2)
    full_data_est <- RecurrenceEnd:::survfit_to_survfun(full_data_km)
    full_data_res <- compute_bias(quantile_exp = true_quantiles,
                                km_function = full_data_est,
                                target_probs = quantiles)
    resmat[1, ] <- c(i, 0, full_data_res)

    # apply all other methods
    for (idx in seq_along(method_list)){
      mm <- method_list[[idx]]
      mm$formula <- update(Recur(time = time, id = patient.id, event = indicator) ~ .,
                           mm$formula)
      est <- do.call("estimate_end", c(mm, list(data=observed_data)))
      res <- compute_bias(quantile_exp = true_quantiles,
                          km_function = est$fit,
                          target_probs = quantiles)
      resmat[idx+1,] <- c(i, idx, res)
    }

    res_list <- c(res_list, list(resmat))
  }
  browser

  resdat <- do.call(rbind, res_list)
  colnames(resdat) <- c("Rep", "Method", sprintf("Q%d%%", round(quantiles*100)))
  resdat <- as_tibble(resdat) |>
    mutate(Method = c("Full data", names(method_list))[Method + 1])
  resdat
}

# Loop through settings with parallel processing
loop_settings <- function(nsim, method_list, quantiles, param_grid, ...){

  # Set up parallel backend
  num_cores <- detectCores() - 1  # Use all but one core
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  loop <- foreach(i = 1:nrow(param_grid), .combine = bind_rows,
                  .multicombine = TRUE,
                  .packages = c("dplyr", "survival", "RecurrenceEnd"),
                  .export = c("run_setting", "compute_bias", "compute_mse"))

  res <- loop %do% {
    params <- c(list(nsim=nsim, method_list = method_list, quantiles=quantiles),
                param_grid[i,], ...)
    res1 <- do.call(run_setting, params)
    bind_cols(res1, param_grid[i,])
  }
}


######################
# Simulation study 1 #
######################
set.seed(111020251)
quantiles0 <- c(0.75, 0.50, 0.25)
nsim0 <- 500

# Define parameter grid
param_grid1 <- expand.grid(mu_d = c(1, 2, 5),
                          p_c = c(0.15, 0.3),
                          EN = c(5, 10, 20),
                          C_min = 1,
                          sigma2 = 0.25,
                          n = 300)
method_list1 <- list(
  "Naive" = list(method="naive", formula = ~1),
  "Threshold" = list(method="quantile", quantile = 0.9, formula = ~1),
  "NPMLE" = list(method="NPMLE", formula = ~1)
)

# run simulation 1
system.time(
 bias_res1 <- loop_settings(nsim=nsim0, method_list = method_list1,
                           quantiles=quantiles0, param_grid = param_grid1)
)

# aggregate results
avebias_res1 <- bias_res1 |>
  group_by(Method, mu_d, p_c, EN) |>
  summarize(across(starts_with("Q"), ~mean(.x, na.rm=TRUE)),
            .groups="drop") |>
  mutate(setting = as.integer(factor(interaction(mu_d, p_c, EN))))

save(avebias_res1, file = "Manuscript/Data/Simulation1.RData")

######################
# Simulation study 2 #
######################
set.seed(111020252)

cov_fun0 <- function(n, ...) {
  matrix(rbinom(n, 1, 0.5), ncol = 1)
}
# Define parameter grid
param_grid2 <- expand.grid(mu_d =  2,
                           p_c = 0.15,
                           EN = 10,
                           C_min = 1,
                           sigma2 = 0.25,
                           n = 300,
                           beta = 1)

method_list2 <- list(
  "Naive" = list(method="naive", formula = ~1),
  "Threshold" = list(method="quantile", quantile = 0.9, formula = ~1),
  "NPMLE unadjusted" = list(method="NPMLE", formula = ~1),
  "NPMLE adjusted" = list(method="NPMLE", formula = ~Z)
)

system.time(
  bias_res2 <- loop_settings(nsim=nsim0, method_list = method_list2, quantiles=quantiles0,
                          param_grid = param_grid2, cov_fun = cov_fun0)
)

save(bias_res2, file = "Manuscript/Data/Simulation2.RData")
