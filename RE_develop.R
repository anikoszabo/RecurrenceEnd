library(devtools)

use_package("nspmix")
use_package("reda", type = "Depends")
use_package("survival", type = "Depends")
use_build_ignore("RE_develop.R")
use_gpl3_license()
use_github()
desc::desc_add_author("Aniko", "Szabo", "aszabo@mcw.edu", role=c("aut","cre"))
desc::desc_add_author("Juntian", "Wang", "juwang@mcw.edu", role=c("aut"))
desc::desc_add_author("Duo", "Yu", "duoyu@mcw.edu", role=c("aut"))

use_readme_rmd()

usethis::use_data_raw(name="SimulatedData")

re <- as.package(".")
load_all(re)
document(re)
run_examples()

build_readme()
check(re)

set.seed(2453)
ld0 <- 1
lr0 <- 5
ss <- 0.1
b0 <- c(0,0)
aa <- sim_recur_end(n=300, lambda = ld0, lambda_r = lr0, sigma2=ss,
                   cov_fun = function(n)cbind(rbinom(n, size=1, prob = 0.5), rnorm(n)),
                   beta = b0, lambda_c = 1, C_min = 1)
a <- subset(aa, nevents > 0)
res_n <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ 1,
                     method="naive", data=a)
res_q <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ 1,
                      method="quantile", data=a, quantile=0.99)
res_np0 <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ 1,
                      method="NPMLE", data=a, IPSW=FALSE)
res_np0b <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ 1,
                        method="NPMLE", data=a, IPSW=TRUE)
res_np <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z.1 + Z.2,
                       method="NPMLE", data=a, IPSW=FALSE)
res_npb <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z.1 + Z.2,
                       method="NPMLE", data=a, IPSW=TRUE)
res_np2 <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z.1 + Z.2 + log(xi),
                       method="NPMLE", data=a,
                       known_recur = list(S0 = function(x)pexp(x, rate=lr0, lower=FALSE),
                                          coefs = c(b0,1)), IPSW=FALSE)
res_np2b <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z.1 + Z.2 + log(xi),
                        method="NPMLE", data=a,
                        known_recur = list(S0 = function(x)pexp(x, rate=1, lower=FALSE),
                                           coefs = c(1,-1,1)), IPSW=TRUE)

a0 <- subset(a, indicator==0)
a0$obs <- pmin(a0$disease_onset, a0$C)
a0$delta <- ifelse(a0$disease_onset <= a0$C, 1, 0)
km <- survfit(Surv(obs, delta) ~ 1, data=a0)

curve(pexp(x, rate=ld0, lower=FALSE), from=0, to=max(a$time[a$indicator==1]),
      ylim=c(0,1), ylab="")
lines(res_n$fit, do.points=FALSE, col="red")
lines(res_q$fit, do.points=FALSE, col="pink")
lines(res_np0$fit, do.points=FALSE, col="magenta")
lines(res_np0b$fit, do.points=FALSE, col="magenta", lty=2)
lines(res_np$fit, do.points=FALSE, col="blue")
lines(res_npb$fit, do.points=FALSE, col="blue", lty=2)
lines(res_np2$fit, do.points=FALSE, col="green")
lines(res_np2b$fit, do.points=FALSE, col="green", lty=2)
lines(km, conf.int=FALSE, col="brown")
#legend("topright", legend=c("Naive", "NPMLE", "NPMLE_known"), col=2:4, lty=1)

# bootstrap
res_qboot <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ 1,
                       method="quantile", data=a, quantile=0.99, bootCI = TRUE,
                       bootB = 100)
res_npboot <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z.1 + Z.2,
                       method="NPMLE", data=a, bootCI = TRUE,
                       bootB = 100)
curve(pexp(x, rate=ld0, lower=FALSE), from=0, to=max(a$time[a$indicator==1]),
      ylim=c(0,1), ylab="")
lines(res_npboot$fit, do.points=FALSE, col="blue")
lines(res_npboot$ci$lower, do.points=FALSE, col="blue", lty=2)
lines(res_npboot$ci$upper, do.points=FALSE, col="blue", lty=2)
lines(res_qboot$fit, do.points=FALSE, col="pink")
lines(res_qboot$ci$lower, do.points=FALSE, col="pink", lty=2)
lines(res_qboot$ci$upper, do.points=FALSE, col="pink", lty=2)

# terminal events
# randomly select 20% of censoring indicators that are larger than the true event time
# to be terminal events
set.seed(3456)
all_pts <- unique(a$patient.id)
possible_pts <- unique(a$patient.id[a$disease_onset < a$C])
term <- ifelse(all_pts %in% possible_pts,
               rbinom(n = length(possible_pts), size=1, prob = 0.2), 0)

tres_q <- estimate_end(Recur(time=time, id=patient.id, event=indicator, terminal = term) ~ 1,
                      method="quantile", data=a, quantile=0.95)
tres_np <- estimate_end(Recur(time=time, id=patient.id, event=indicator, terminal = term) ~ Z.1 + Z.2,
                        method="NPMLE", data=a)


curve(pexp(x, rate=ld0, lower=FALSE), from=0, to=max(a$time[a$indicator==1]),
      ylim=c(0,1), ylab="")
lines(res_q$fit, do.points=FALSE, col="red")
lines(tres_q$fit, do.points=FALSE, col="pink")
lines(res_np$fit, do.points=FALSE, col="blue")
lines(tres_np$fit, do.points=FALSE, col="lightblue")
lines(tres_np2$fit, do.points=FALSE, col="lightblue", lty=2)

