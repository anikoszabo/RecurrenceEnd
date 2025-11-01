library(RecurrenceEnd)
library(doParallel)
library(foreach)
library(dplyr)
library(tidyr)
# Simulation Study

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



# Set up parallel backend
num_cores <- detectCores() - 1  # Use all but one core
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define parameter grid
param_grid <- expand.grid(mu_d = c(1, 2, 5),
                          p_c = c(0.15, 0.3),
                          EN = c(5, 10, 20))
cov_fun <- function(n, ...) {
  matrix(rbinom(n, 1, 0.5), ncol = 1)
}
c_min <- 1
beta <- 1
sigma2 <- 0.25
n = 300


iter = 500
quantiles <- c(0.75, 0.50, 0.25)

bias_results <- list()
final_bias_data <- foreach(i = 1:nrow(param_grid), .combine = rbind,
                           .packages = c("dplyr", "survival", "RecurrenceEnd")) %dopar% {

                             # Extract parameter values for this iteration
                             mu_d <- param_grid$mu_d[i]
                             p_c <- param_grid$p_c[i]
                             EN <- param_grid$EN[i]
                             lambda_d <- log(2)/mu_d
                             lambda_r <- EN*lambda_d
                             lambda_c <- (lambda_d*p_c)/(exp(-lambda_d*c_min) - p_c)


                             # Storage for this parameter combination
                             bias_npmle_adjusted <- as.data.frame(matrix(NA, ncol = 3, nrow = iter))
                             colnames(bias_npmle_adjusted) <- c("25%", "50%", "75%")
                             bias_naive <- bias_full_data <- bias_data_driven <- bias_npmle_unadjusted <- bias_npmle_adjusted
                             true_quantiles <- qexp(quantiles, rate = lambda_d, lower.tail = FALSE)


                             for (k in 1:iter) {
                               observed_data <- sim_recur_end(n = n, lambda_d = lambda_d, lambda_r = lambda_r, sigma2 = sigma2,
                                                              cov_fun = cov_fun, beta = beta, lambda_c = lambda_c, C_min=c_min)

                               observed_data <- subset(observed_data, nevents > 0)

                               observed_data2 <- observed_data |>
                                 mutate(minimum = pmin(recurrence_end, C),
                                        min_indicator = ifelse(recurrence_end <= C, 1, 0)) |>
                                 select(patient.id, minimum, min_indicator) |>
                                 unique()

                               # # Apply survival analysis methods
                               naive_km <- estimate_end(Recur(time = time, id = patient.id, event = indicator) ~ Z,
                                                        method = "naive", threshold = 0, data = observed_data)

                               data_driven_km <- estimate_end(Recur(time = time, id = patient.id, event = indicator) ~ Z,
                                                              method = "quantile", quantile = 0.9, data = observed_data)

                               npmle_adjusted <- estimate_end(Recur(time = time, id = patient.id, event = indicator) ~ Z,
                                                              method = "NPMLE", data = observed_data)

                               npmle_unadjusted <- estimate_end(Recur(time = time, id = patient.id, event = indicator) ~ 1,
                                                                method = "NPMLE", data = observed_data)

                               full_data_km <- survfit(Surv(minimum, min_indicator) ~ 1, data=observed_data2)

                               full_data <- RecurrenceEnd:::survfit_to_survfun(full_data_km)

                               # # Compute bias
                               bias_npmle_adjusted[k, ] <- compute_bias(quantile_exp = true_quantiles,
                                                                        km_function = npmle_adjusted$fit,
                                                                        target_probs = quantiles)
                               bias_npmle_unadjusted[k, ] <- compute_bias(quantile_exp = true_quantiles,
                                                                          km_function = npmle_unadjusted$fit,
                                                                          target_probs = quantiles)
                               bias_naive[k, ] <- compute_bias(quantile_exp = true_quantiles,
                                                               km_function = naive_km$fit,
                                                               target_probs = quantiles)
                               bias_full_data[k, ] <- compute_bias(quantile_exp = true_quantiles,
                                                                   km_function = full_data,
                                                                   target_probs = quantiles)
                               bias_data_driven[k, ] <- compute_bias(quantile_exp = true_quantiles,
                                                                     km_function = data_driven_km$fit,
                                                                     target_probs = quantiles)
                             }


                             # # Add method labels
                             bias_npmle_adjusted$Method <- "NPMLE_adjusted"
                             bias_npmle_unadjusted$Method <- "NPMLE_unadjusted"
                             bias_naive$Method <- "Naive"
                             bias_full_data$Method <- "full_data"
                             bias_data_driven$Method <- "data_driven"

                             # Combine results for this parameter set
                             bias_combined <- rbind(bias_npmle_adjusted, bias_npmle_unadjusted, bias_naive,
                                                    bias_full_data, bias_data_driven)

                             # Store results with parameter values
                             bias_combined$mu_d <- mu_d
                             bias_combined$p_c <- p_c
                             bias_combined$EN <- EN
                             bias_combined$n <- nrow(observed_data2)
                             bias_results[[i]] <- bias_combined

                           }
# Convert list to dataframe
final_bias_data <- bind_rows(bias_list)
