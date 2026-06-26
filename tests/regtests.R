library("hhh4ZI")

## fit and simulate ZI without lags
fit <- hhh4ZI(measles[1:156,10], list(zi = list(lag = NULL)))
summary(fit)
sim <- simulate(fit, seed = 20260626, y.start = 0)
## hhh4ZI 0.3.2 gave Error : missing value where TRUE/FALSE needed
stopifnot(mean(observed(sim) == 0) > 0.9)  # cf. ZI prob. fit$gamma[1]
