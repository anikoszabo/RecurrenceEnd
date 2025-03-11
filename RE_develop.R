library(devtools)

use_package("nspmix")
use_package("reda")
use_package("survival")
use_build_ignore("RE_develop.R")
use_gpl3_license()
use_github()

re <- as.package(".")
load_all(re)
document(re)

check(re)

set.seed(245)
a <- sim_recur_end(n=10, lambda_d = 0.1, lambda_r = 3, sigma2=0.1,
                   cov_fun = function(n)rbinom(n, size=1, prob = 0.5),
                   beta = 1, lambda_c = 0.2, C_min = 1)

library(reda)
library(survival)
library(nspmix)
res_n <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z,
                     method="naive", data=a)
res_np <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z,
                      method="NPMLE", data=a)
res_np2 <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z,
                       method="NPMLE", data=a,
                       known_recur = list(S0 = function(x)pexp(x, rate=3, lower=FALSE),
                                          coefs = 1))
