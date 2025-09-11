
set.seed(2453)
ld0 <- 1
lr0 <- 5
ss <- 0.1
b0 <- c(-0.3,0.5)
aa <- sim_recur_end(n=100, lambda = ld0, lambda_r = lr0, sigma2=ss,
                    cov_fun = function(n)cbind(rbinom(n, size=1, prob = 0.5), rnorm(n)),
                    beta = b0, lambda_c = 1, C_min = 1)
SimulatedData <- subset(aa, nevents > 0) |>
  transform(recurrence_end = NULL,
            C = NULL,
            xi = NULL,
            nevents = NULL)


usethis::use_data(SimulatedData, overwrite = TRUE)
