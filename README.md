
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RecurrenceEnd

<!-- badges: start -->

<!-- badges: end -->

The goal of RecurrenceEnd is to implement several non-parametric methods
for estimating the distribution of the unobservable terminating event of
a recurrentevent process.

## Installation

You can install the development version of RecurrenceEnd from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("anikoszabo/RecurrenceEnd")
```

## Example

``` r
library(RecurrenceEnd)
#> Loading required package: reda
#> Loading required package: survival
#> Warning: package 'survival' was built under R version 4.3.3
```

The `SimulatedData` data set contains simulated recurrent events for 97
subjects. The recurrent events continue until an unobserved end-time,
when they stop. The goal is to estimate the distribution of the
unobserved end-time across the subjects.

The data has one row for each occurrence of the recurrent event within a
subject, and a last row for the follow-up time during which no
additional events have occurred. The following table shows data from 3
subjects with varying numbers of recurrent events. The ‘time’ variable
records the event/last follow-up time, the ‘indicator’ distinguishes
between events (1) and last follow-up (0), while Z.1 and Z.2 are
covariates that affect the rate of the recurrent event process.

``` r
SimulatedData[SimulatedData$patient.id %in% c(5,6,7),]
#>    patient.id       time indicator Z.1        Z.2
#> 25          5 0.05351576         1   1 0.04779898
#> 26          5 5.01451459         0   1 0.04779898
#> 27          6 0.03134527         1   0 0.98247007
#> 28          6 0.13323338         1   0 0.98247007
#> 29          6 0.23782361         1   0 0.98247007
#> 30          6 0.29015648         1   0 0.98247007
#> 31          6 0.33696037         1   0 0.98247007
#> 32          6 0.41951277         1   0 0.98247007
#> 33          6 0.43688110         1   0 0.98247007
#> 34          6 0.44779014         1   0 0.98247007
#> 35          6 0.85480934         1   0 0.98247007
#> 36          6 0.89867815         1   0 0.98247007
#> 37          6 1.77427821         0   0 0.98247007
#> 38          7 0.24106241         1   0 1.01133452
#> 39          7 0.31224971         1   0 1.01133452
#> 40          7 0.32755919         1   0 1.01133452
#> 41          7 0.43923456         1   0 1.01133452
#> 42          7 1.47966919         0   0 1.01133452
```

We can visualize the data to get a better understanding of its
structure. In the plot, each gray line represents the follow-up of a
patients, with the blue points showing the occurrences of the recurrent
event in that patient. We can see that for some subjects the recurrent
events seem to keep going until the last follow-up time (and likely
beyond it as well), while for others there is a long stretch of
follow-up time without any recurrent events, indicating that the
unobserved ending event has probably happened.

``` r
library(reReg) # for plotting recurrent event data
#> Warning: package 'reReg' was built under R version 4.3.3
plotEvents(Recur(time=time, id=patient.id, event=indicator) ~ 1,
           data=SimulatedData, 
           control=list(cex=1, width=0.2, recurrent.color="navy")) 
```

<img src="man/figures/README-plot-1.png" width="100%" />

This package implements four methods for estimating the distribution of
the unobserved end-time. The NPMLE method is recommended, the others are
ad-hoc methods implemented for comparison.

``` r
res_np <- estimate_end(Recur(time=time, id=patient.id, event=indicator) ~ Z.1+Z.2, 
                       data = SimulatedData, bootCI = TRUE, bootB = 50, 
                       method="NPMLE")

plot(res_np, col="blue", conf.int = TRUE, conf.lty=3, conf.col="blue")
```

<img src="man/figures/README-FitModel-1.png" width="100%" />

The values at specific time-points can be obtained using the `predict`
function, while quantiles can be extracted using the `quantile` and
`median` functions.

``` r
# probability of ongoing recurrent events at time = 1
predict(res_np, times=1, conf.int = TRUE)
#> $time
#> [1] 1
#> 
#> $pred
#> [1] 0.3844874
#> 
#> $lower
#> [1] 0.2829808
#> 
#> $upper
#> [1] 0.4774988

# median ending time
median(res_np, conf.int=TRUE)
#> $quantile
#>        50 
#> 0.7829867 
#> 
#> $lower
#>        50 
#> 0.5809653 
#> 
#> $upper
#>        50 
#> 0.9056862
```
