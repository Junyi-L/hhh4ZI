################################################################################
### Demo of hhh4ZI() modelling of measles in Germany - data("measles")
### based on
###
### Lu and Meyer (2023): A zero-inflated endemic-epidemic model
### with an application to measles time series in Germany.
###
### RUNNING THE WHOLE SCRIPT TAKES about 1 hour!
################################################################################

library(hhh4ZI)
options(digits = 5)

# for reproducibility, as ri() parameters are initialized randomly
RNGversion("3.6.3")
set.seed(10052021)

ptm <- proc.time()


##################
# Data preparation
##################

# compute adjusted proportion of vaccinated school starters by year and state
# using the same assumptions as in Herzog et al. (2011)
PassRate <- VacPass/TotalKids
Dosis1Rate <- PassRate * Dosis1/100 + (1 - PassRate) * 0.5 * Dosis1/100

# Impute two missing values in Hamburg via LOCF
Dosis1Rate[c("2017", "2018"), "Hamburg"] <- Dosis1Rate["2016", "Hamburg"]
stopifnot(!anyNA(Dosis1Rate))

# setup covariate matrix for dim(measles), i.e., replicate each row 26 times
Dosis1Rate <- as.matrix(Dosis1Rate)
Coverage <- Dosis1Rate[rep(1 : nrow(Dosis1Rate), each = measles@freq),]
rm(list = c("PassRate", "Dosis1Rate"))

# covariate data for the model: time index and log(proportion unvaccinated)
covarlist <- list(t = epoch(measles),
                  logVac0 = log(1 - Coverage))
Vac0 <-  (1 - 0.92 * Coverage)
# higher-than-default niter is needed for rolling forecasts
control.optimizer <-  list(stop = list(tol = 1e-5, niter = 200),
                           regression = list(method = "nlminb"),
                           variance = list(method = "Nelder-Mead"))

##################################################
# Fit the models
##################################################

#############################
####### Model P0 ############
#############################

# endemic component: Intercept + sine/cosine terms
f.end <- addSeason2formula(f = ~ 1, S = 1, period = 26)
# autoregressive component: Intercept + proportion unvaccinated
f.ar <- ~ 1 + logVac0

model.P0 <- list(ar = list(f = f.ar),
                 end = list(f = f.end, offset = population(measles)),
                 data = covarlist,
                 family = "Poisson",
                 optimizer = control.optimizer)
# fit the model
result.P0 <- hhh4(measles, model.P0)
summary(result.P0, amplitudeShift = TRUE, maxEV = TRUE)


#############################
####### Model NB1 ############
#############################

f.end <- addSeason2formula(f = ~ 1, S = 2, period = 52)
f.ar <- addSeason2formula(f = ~ 1, S = 2, period = 52)

model.NB1 <- list(ar = list(f = f.ar, offset =  Vac0 ),
                  end = list(f = f.end, offset = population(measles) * Vac0),
                  data = covarlist,
                  family = "NegBinM",
                  optimizer = control.optimizer)

result.NB1 <- hhh4(measles, model.NB1)
summary(result.NB1, amplitudeShift = TRUE, maxEV = TRUE)

#############################
####### Model NB2 ############
#############################

f.end <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "none"),
                           S = 2, period = 52)
f.ar <- addSeason2formula(f = ~ -1 + ri(type = "iid", corr = "none"),
                          S = 2, period = 52)

model.NB2 <- list(ar = list(f = f.ar, offset =  Vac0),
                  end = list(f = f.end, offset = population(measles) * Vac0),
                  data = covarlist,
                  family = "NegBinM",
                  optimizer = control.optimizer)

result.NB2 <- hhh4(measles, model.NB2)
summary(result.NB2, amplitudeShift = TRUE, maxEV = TRUE)

#############################
####### Model NB3############
#############################

f.end <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "all"),
                           S = 2, period = 52)
f.ar <- addSeason2formula(f = ~ -1 + ri(type = "iid", corr = "all"),
                          S = 2, period = 52)

model.NB3 <- list(ar = list(f = f.ar, offset =  Vac0 ),
                  end = list(f = f.end, offset = population(measles) * Vac0 ),
                  data = covarlist,
                  family = "NegBinM",
                  optimizer = control.optimizer)

result.NB3 <- hhh4(measles, model.NB3)
summary(result.NB3, amplitudeShift = TRUE, maxEV = TRUE)


#############################
####### Model ZI1 ############
#############################

f.end <- addSeason2formula(f = ~ 1, S = 2, period = 52)
f.ar <- addSeason2formula(f = ~ 1, S = 2, period = 52)
f.zi <- ~ 1

model.ZI1 <- list(ar = list(f = f.ar, offset =  Vac0 ),
                  end = list(f = f.end, offset = population(measles) * Vac0 ),
                  zi = list(f = f.zi),
                  data = covarlist,
                  family = "NegBinM",
                  optimizer = control.optimizer)

result.ZI1 <- hhh4ZI(measles, model.ZI1)
summary(result.ZI1, amplitudeShift = TRUE, maxEV = TRUE)

#############################
####### Model ZI2 ############
#############################
f.end <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "none"),
                           S = 2, period = 52)
f.ar <- addSeason2formula(f = ~-1  +  ri(type = "iid", corr = "none"),
                          S = 2, period = 52)
f.zi <- ~-1  +  ri(type = "iid", corr = "none")

model.ZI2 <- list(ar = list(f = f.ar, offset =  Vac0 ),
                  end = list(f = f.end, offset = population(measles) * Vac0 ),
                  zi = list(f = f.zi),
                  data = covarlist,
                  family = "NegBinM",
                  optimizer = control.optimizer)

result.ZI2 <- hhh4ZI(measles, model.ZI2)
summary(result.ZI2, amplitudeShift = TRUE, maxEV = TRUE)

#############################
####### Model ZI3 ############
#############################

f.end <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "all") ,
                           S = 2, period = 52)
f.ar <- addSeason2formula(f = ~-1  +  ri(type = "iid", corr = "all"),
                          S = 2, period = 52)
f.zi <- ~-1  +  ri(type = "iid", corr = "all")

model.ZI3 <- list(ar = list(f = f.ar, offset =  Vac0 ),
                  end = list(f = f.end, offset = population(measles) * Vac0 ),
                  zi = list(f = f.zi),
                  data = covarlist,
                  family = "NegBinM",
                  optimizer = control.optimizer)

result.ZI3 <- hhh4ZI(measles, model.ZI3)
summary(result.ZI3, amplitudeShift = TRUE, maxEV = TRUE)

#############################
####### Model ZI4 ############
#############################

f.end <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "none"),
                           S = 2, period = 52)
f.ar <- addSeason2formula(f = ~-1  +  ri(type = "iid", corr = "none"),
                          S = 2, period = 52)
f.zi <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "none"),
                          S = 2, period = 52)

model.ZI4 <- list(ar = list(f = f.ar, offset =  Vac0 ),
                  end = list(f = f.end, offset = population(measles) * Vac0 ),
                  zi = list(f = f.zi),
                  data = covarlist,
                  family = "NegBinM",
                  optimizer = control.optimizer)

result.ZI4 <- hhh4ZI(measles, model.ZI4)
summary(result.ZI4, amplitudeShift = TRUE, maxEV = TRUE)

#############################
####### Model ZI5 ############
#############################

f.end <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "all"),
                           S = 2, period = 52)
f.ar <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "all"),
                          S = 2, period = 52)
f.zi <- addSeason2formula(f = ~-1 +  ri(type = "iid", corr = "all"),
                          S = 2, period = 52)

model.ZI5 <- list(ar = list(f = f.ar, offset =  Vac0 ),
                  end = list(f = f.end, offset = population(measles) * Vac0 ),
                  zi = list(f = f.zi),
                  data = covarlist,
                  family = "NegBinM",
                  optimizer = control.optimizer,
                  verbose = TRUE)

result.ZI5 <- hhh4ZI(measles, model.ZI5)
summary(result.ZI5, amplitudeShift = TRUE, maxEV = TRUE)


###########################################################
## Exemplary summary of model ZI3
###########################################################
plotHHH4ZI_season(result.ZI3, components = c("en", "ar"),
                  period = 52, xlab = "biweek")
plotHHH4ZI_maxEV(result.ZI3)
plotHHH4ZI_fitted(units = 1:8,  result.ZI3, par.settings = list(mfrow=c(4,2)))
plotHHH4ZI_fitted(units = 9:16, result.ZI3, par.settings = list(mfrow=c(4,2)))
plotHHH4ZI_maps(result.ZI3, prop = TRUE,
                which = c("mean", "endemic", "epi.own", "zi"),
                main = c("fitted counts", "endemic proportion",
                         "autoregressive proportion", "ZI probability"),
                par.settings = list(add.line = list(col = "white")))
plotHHH4ZI_ri(result.ZI3, component = "ar.ri(iid)", main = "autoregressive", exp = TRUE)
plotHHH4ZI_ri(result.ZI3, component = "end.ri(iid)", main = "endemic", exp = TRUE)
plotHHH4ZI_ri(result.ZI3, component = "zi.ri(iid)", main = "zero inflation", exp = TRUE)

######################################################################
# Compare the predictive performance of the models by computing
# one-step-ahead predictions to be assessed by proper scoring rules
######################################################################
tp <- 261
i <- 1
type <- "rolling"
model_names <- c("P0",
                 "NB1", "NB2", "NB3",
                 "ZI1", "ZI2", "ZI3", "ZI4", "ZI5")
score_table <- data.frame()
score_list <- list()

# forecast
for (model_name in model_names){
  modeli <- get(paste0("result.", model_name))
  cat("\n--- forecasting with", model_name, "---\n")
  osai <- oneStepAhead(modeli, tp = tp, type = type, which.start = "final")
  score_list[[i]] <-
    scorei <- scores(osai)
  score_table <- rbind(score_table,
                       c(colMeans(scorei), sapply(as.data.frame(scorei), max)))
  i <- i + 1
}

ref_score <- score_list[[7]][,1]
Test_score <- sapply(score_list, FUN = function(x){
  permutationTest(x[,1], ref_score, nPermutation = 9999)$pVal.permut
})
Test_score[7] <- NA
colnames(score_table) <- c(colnames(scorei)[1:4], colnames(scorei)[1:4])

printformat <- function(x) {
  paste0(formatC(round(x, digits = 2), format='f', digits=2 ), " (",rank(x), ")")
}
printformat2 <- function(x) {
  round(x, digits = 2)
}

score_table2 <- data.frame(model_names, printformat(score_table[,1]),
                           sapply(Test_score, printformat2),
                           printformat(score_table[,5]),
                           printformat(score_table[,3]),
                           printformat(score_table[,2]),
                           printformat(score_table[,4]))
colnames(score_table2) <- c("Model", "LS", "p-value","maxLS", "DSS","RPS", "SES")
print(score_table2)

cat("Time elapsed: ", proc.time() - ptm, "\n")
