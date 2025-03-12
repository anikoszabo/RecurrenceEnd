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


library(reda)
library(survival)
library(nspmix)
set.seed(2453)
a <- sim_recur_end(n=1000, lambda_d = 0.2, lambda_r = 1, sigma2=0.1,
                   cov_fun = function(n)cbind(rbinom(n, size=1, prob = 0.5), rnorm(n)),
                   beta = c(1,-0.2), lambda_c = 0.2, C_min = 1)
res_n <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ 1,
                     method="naive", data=a)
res_np0 <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ 1,
                      method="NPMLE", data=a)
res_np <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z.1 + Z.2,
                       method="NPMLE", data=a)
res_np2 <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z.1 + Z.2 + log(xi),
                       method="NPMLE", data=a,
                       known_recur = list(S0 = function(x)pexp(x, rate=1, lower=FALSE),
                                          coefs = c(1,-1,1)))
a0 <- subset(a, indicator==0)
a0$obs <- pmin(a0$disease_onset, a0$C)
a0$delta <- ifelse(a0$disease_onset <= a0$C, 1, 0)
km <- survfit(Surv(obs, delta) ~ 1, data=a0)

curve(pexp(x, rate=0.2, lower=FALSE), from=0, to=max(a$time), ylim=c(0,1), ylab="")
lines(res_n$fit, do.points=FALSE, col="red")
lines(res_np0$fit, do.points=FALSE, col="magenta")
lines(res_np$fit, do.points=FALSE, col="blue")
lines(res_np2$fit, do.points=FALSE, col="green")
lines(km, conf.int=FALSE, col="brown")
#legend("topright", legend=c("Naive", "NPMLE", "NPMLE_known"), col=2:4, lty=1)
