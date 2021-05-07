library("hhh4ZI")

## source this script with varying 'i' to find that hhh4ZI fits more reliably
if (!exists("i")) i <- 1
measlesi <- suppressMessages(measles[,i])

## endemic-only zero-inflated model
control <- list(end = list(f = addSeason2formula(~1, period=26)),
                zi  = list(f = ~1, lag = 1),
                family = "NegBin1")

fit <- hhh4ZI(measlesi, control)
summary(fit)

## we can fit this simple model using pscl::zeroinfl for comparison
df <- as.data.frame(measlesi, tidy = TRUE)
df$t <- df$epoch - 1
df$Ylag <- c(NA, head(df$observed, -1))
psclfit <- pscl::zeroinfl(
    observed ~ 1 + sin(2*pi*t/26) + cos(2*pi*t/26) | 1 + Ylag,
    data = df, subset = epoch > 1, dist = "negbin")
summary(psclfit)

## compare logLik
stopifnot(all.equal(logLik(fit), logLik(psclfit), tolerance = 1e-6))

## compare coefficients
comp_coef <- cbind(hhh4ZI = coef(fit),
                   zeroinfl = c(coef(psclfit), 1/psclfit$theta))
print(comp_coef)
stopifnot(all.equal(comp_coef[,1], comp_coef[,2], tolerance = 1e-4))

## compare estimated standard errors
comp_se <- cbind(hhh4ZI = head(sqrt(diag(vcov(fit))), -1),
                 zeroinfl = sqrt(diag(vcov(psclfit))))
print(comp_se)
stopifnot(all.equal(comp_se[,1], comp_se[,2], tolerance = 1e-4))

## compare whole estimated variance-covariance matrix
stopifnot(all.equal(vcov(psclfit), vcov(fit)[1:5, 1:5],
                    check.attributes = FALSE, tolerance = 1e-3))
