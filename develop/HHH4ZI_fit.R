# library(hhh4ZI)
data("measles", package = "hhh4ZI")
measles <- aggregate(measles, by = "time", nfreq = 26)

stsObj <- measles
neW <- neighbourhood(measles)


f.1 <- addSeason2formula(f = ~-1 + fe(1, unitSpecific = TRUE),
                           S = 1, period = 26)
f.2 <- addSeason2formula(f = ~-1 + fe(1, unitSpecific = FALSE),
                             S = 1, period = 26)
f.3 <- addSeason2formula(f = ~-1 + fe(1, unitSpecific = FALSE))
# f.ar <- ~ -1 + fe(1, unitSpecific = TRUE)
# f.ne <- ~-1 + fe(1, unitSpecific = TRUE)
fitZi <- hhh4ZI(stsObj,
                control = list(
                  ar = list(f = f.3,        # a formula "exp(x'lamba)*y_t-lag" (ToDo: matrix)
                            offset = 1,      # multiplicative offset
                            lag = 1),        # autoregression on y_i,t-lag
                  ne = list(f = ~1,        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
                            offset = 1,      # multiplicative offset
                            lag = 1,         # regression on y_j,t-lag
                            weights = neW == 1,  # weights w_ji
                            scale = NULL,    # such that w_ji = scale * weights
                            normalize = TRUE), # w_ji -> w_ji / rowSums(w_ji), after scaling
                  end = list(f = f.1,        # a formula "exp(x'nu) * n_it"
                             offset = 1),    # optional multiplicative offset e_it
                  zi = list(f = f.1,
                            lag = 1        # can be a scalar or vector
                  ),
                  family = c("NegBin1"), # or a factor of length nUnit for Negbin
                  #subset = 2:nrow(stsObj),   # epidemic components require Y_{t-lag}
                  optimizer = list(stop = list(tol = 1e-5, niter = 100), # control arguments
                                   regression = list(method = "nlminb"), # for optimization
                                   variance = list(method = "Nelder-Mead")),  # <- or "Nelder-Mead"
                  verbose = TRUE,           # level of reporting during optimization
                  start = list(fixed = NULL, # list of start values, replacing initial
                               random = NULL,  # values from fe() and ri() in 'f'ormulae
                               sd.corr = NULL),
                  data = list(t = stsObj@epoch - min(stsObj@epoch)), # named list of covariates
                  keep.terms = TRUE  # whether to keep interpretControl(control, stsObj)
                )
)
summary(fitZi)
fitZi2 <- update(fitZi,subset.upper=370)
summary(fitZi2)
pred <- oneStepAhead(fitZi, 318, type = "rolling", which.start = "current")
x <- hhh4ZI:::quantile.oneStepAhead_hhh4ZI(pred)
hhh4ZI:::plot.oneStepAhead_hhh4ZI(pred, unit = 11)
scores(pred)
scores(fitZi2)
pit(pred, units = 16)

fitH <- hhh4(stsObj,
             control = list(
               ar = list(f = ~-1 + ri(type = "iid", corr = "all"),        # a formula "exp(x'lamba)*y_t-lag" (ToDo: matrix)
                         offset = 1,      # multiplicative offset
                         lag = 1),        # autoregression on y_i,t-lag
               ne = list(f = ~-1 + ri(type = "iid", corr = "all"),        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
                         offset = 1,      # multiplicative offset
                         lag = 1,         # regression on y_j,t-lag
                         weights = neW == 1,  # weights w_ji
                         scale = NULL,    # such that w_ji = scale * weights
                         normalize = TRUE), # w_ji -> w_ji / rowSums(w_ji), after scaling
               end = list(f = ~-1 + ri(type = "iid", corr = "all"),        # a formula "exp(x'nu) * n_it"
                          offset = 1),    # optional multiplicative offset e_it
               family = c("NegBin1"), # or a factor of length nUnit for Negbin
               #subset = 2:nrow(stsObj),   # epidemic components require Y_{t-lag}
               optimizer = list(stop = list(tol = 1e-5, niter = 200), # control arguments
                                regression = list(method = "nlminb"), # for optimization
                                variance = list(method = "Nelder-Mead")),  # <- or "Nelder-Mead"
               verbose = TRUE,           # level of reporting during optimization
               start = list(fixed = NULL, # list of start values, replacing initial
                            random = NULL,  # values from fe() and ri() in 'f'ormulae
                            sd.corr = NULL),
               data = list(t = stsObj@epoch - min(stsObj@epoch)), # named list of covariates
               keep.terms = TRUE  # whether to keep interpretControl(control, stsObj)
             ), check.analyticals = FALSE
)
summary(fitH)
fitH_ZI <- hhh4ZI(fitH, control = list(f =  ~-1 + ri(type = "iid", corr = "all"),
                                       lag = 1        # can be a scalar or vector
))

fitH_ZI <- hhh4ZI(fitH,
                  control = list(
                    ar = list(f = ~1,        # a formula "exp(x'lamba)*y_t-lag" (ToDo: matrix)
                              offset = 1,      # multiplicative offset
                              lag = 1),        # autoregression on y_i,t-lag
                    ne = list(f = ~-1 + ri(type = "iid", corr = "all"),        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
                              offset = 1,      # multiplicative offset
                              lag = 1,         # regression on y_j,t-lag
                              weights = neW == 1,  # weights w_ji
                              scale = NULL,    # such that w_ji = scale * weights
                              normalize = TRUE), # w_ji -> w_ji / rowSums(w_ji), after scaling
                    end = list(f = ~-1 + ri(type = "iid", corr = "all"),        # a formula "exp(x'nu) * n_it"
                               offset = 1),    # optional multiplicative offset e_it
                    zi = list(f =  ~-1 + ri(type = "iid", corr = "all"),
                           lag = 1),       # can be a scalar or vector
                    family = c("NegBin1"), # or a factor of length nUnit for Negbin
                    #subset = 2:nrow(stsObj),   # epidemic components require Y_{t-lag}
                    optimizer = list(stop = list(tol = 1e-5, niter = 200), # control arguments
                                     regression = list(method = "nlminb"), # for optimization
                                     variance = list(method = "Nelder-Mead")),  # <- or "Nelder-Mead"
                    verbose = TRUE,           # level of reporting during optimization
                    start = list(fixed = NULL, # list of start values, replacing initial
                                 random = NULL,  # values from fe() and ri() in 'f'ormulae
                                 sd.corr = NULL),
                    data = list(t = stsObj@epoch - min(stsObj@epoch)), # named list of covariates
                    keep.terms = TRUE  # whether to keep interpretControl(control, stsObj)
                  )
)
surveillance:::plot(fitH)
#plotHHH4ZI_fitted(fitH_ZI, decompose = TRUE)


plotHHH4ZI_season(fitZi)
plotHHH4ZI_neweights(fitH_ZI)
plot(fitH_ZI, type = "ri", component = "zi.ri(iid)")
plot(fitH_ZI, type = "neweights")
plot(fitH_ZI, type = "fitted", unit = 3)
plot(fitH_ZI, type = "map")
