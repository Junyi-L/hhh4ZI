################################################################################
### The following are modified versions of functions from the surveillance package
### and wrappers around them.
### See below the original copyright declaration.
################################################################################

################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Scoring rules as discussed in:
### Predictive model assessment for count data
### Czado, C., Gneiting, T. & Held, L. (2009)
### Biometrics 65:1254-1261
###
### Copyright (C) 2010-2012 Michaela Paul, 2014-2015,2017-2019 Sebastian Meyer
################################################################################


## logarithmic score
## logs(P,x) = -log(P(X=x))

.logs <- function (px)
  -log(px)

logs <- function (x, mu, size=NULL, gamma = NULL) {
  if (is.null(size)) {
    - dpois(x, lambda=mu, log=TRUE)
  } else if (is.null(gamma)) {
    - dnbinom(x, mu=mu, size=size, log=TRUE)
  } else {
    - VGAM::dzinegbin(x, munb = mu, size = size, pstr0 = gamma, log = TRUE)
  }
}


## squared error score
## ses(P,x) = (x-mu_p)^2

ses <- function (x, mu, size=NULL, gamma = NULL) {
  if(is.null(gamma)) (x-mu)^2 else
    (x- (1 - gamma) * mu)^2
}


## normalized squared error score (IMPROPER)
## nses(P,x) = ((x-mu_p)/sigma_p)^2

nses <- function (x, mu, size=NULL, gamma = NULL) {
  sigma2 <- if (is.null(size)) mu else
    if(is.null(gamma)) mu * (1 + mu/size) else
      (1 - gamma) * (1 + mu/size + gamma * mu) * mu
  if(is.null(gamma)) ((x-mu)^2) / sigma2 else
    ((x- (1 - gamma) * mu)^2) / sigma2
}


## Dawid-Sebastiani score
## dss(P,x) = ((x-mu_p)/sigma_p)^2 + 2*log(sigma_p)

.dss <- function (meanP, varP, x)
  (x-meanP)^2 / varP + log(varP)

dss <- function (x, mu, size=NULL, gamma = NULL)
  .dss(meanP = if(is.null(gamma)) mu else (1 - gamma) * mu,
       varP = if (is.null(size)) mu else
         if(is.null(gamma)) mu * (1 + mu/size) else
           (1 - gamma) * (1 + mu/size + gamma * mu) * mu,
       x = x)


## ranked probability score
## rps(P,x) = sum_0^Kmax {P(X<=k) - 1(x<=k)}^2

## for a single prediction (general formulation)
.rps <- function (P, ..., x, kmax, tolerance = sqrt(.Machine$double.eps))
{
  ## compute P(X<=k)
  k <- 0:kmax
  Pk <- P(k, ...)

  ## check precision
  if ((1 - Pk[length(Pk)])^2 > tolerance)
    warning("finite sum approximation error larger than tolerance=",
            format(tolerance))

  ## compute the RPS
  sum((Pk - (x <= k))^2)
}

## for a single Poisson prediction
rps_1P <- function (x, mu, k=40, tolerance=sqrt(.Machine$double.eps)) {
  ## return NA for non-convergent fits (where mu=NA)
  if (is.na(mu)) return(NA_real_)
  ## determine the maximum number of summands as Kmax=mean+k*sd
  kmax <- ceiling(mu + k*sqrt(mu))
  ## compute the RPS
  .rps(P = ppois, lambda = mu, x = x,
       kmax = kmax, tolerance = tolerance)
}

## for a single NegBin prediction
rps_1NB <- function (x, mu, size, k=40, tolerance=sqrt(.Machine$double.eps)) {
  ## return NA for non-convergent fits (where mu=NA)
  if (is.na(mu)) return(NA_real_)
  ## determine the maximum number of summands as Kmax=mean+k*sd
  sigma2 <- mu * (1 + mu/size)
  kmax <- ceiling(mu + k*sqrt(sigma2))
  ## compute the RPS
  .rps(P = pnbinom, mu = mu, size = size, x = x,
       kmax = kmax, tolerance = tolerance)
}

## for a single zero inflated NegBin prediction
rps_1ZINB <- function (x, mu, size, gamma,
                       k=40, tolerance=sqrt(.Machine$double.eps)) {
  ## return NA for non-convergent fits (where mu=NA)
  if (is.na(mu)) return(NA_real_)
  ## determine the maximum number of summands as Kmax=mean+k*sd
  mean <- (1 - gamma) * mu
  sigma2 <- (1 - gamma) * (1 + mu/size + gamma * mu) * mu
  kmax <- ceiling(mean + k * sqrt(sigma2))
  ## compute the RPS
  .rps(P = VGAM::pzinegbin, munb = mu, size = size, pstr0 = gamma, x = x,
       kmax = kmax, tolerance = tolerance)
}

## vectorized version
rps <- function (x, mu, size=NULL, gamma = NULL, k=40, tolerance=sqrt(.Machine$double.eps)) {
  res <- if (is.null(size)) {
    mapply(rps_1P, x=x, mu=mu,
           MoreArgs=list(k=k, tolerance=tolerance),
           SIMPLIFY=TRUE, USE.NAMES=FALSE)
  } else if(is.null(gamma)){
    mapply(rps_1NB, x=x, mu=mu, size=size,
           MoreArgs=list(k=k, tolerance=tolerance),
           SIMPLIFY=TRUE, USE.NAMES=FALSE)
  } else {
    mapply(rps_1ZINB, x=x, mu=mu, size=size, gamma = gamma,
           MoreArgs=list(k=k, tolerance=tolerance),
           SIMPLIFY=TRUE, USE.NAMES=FALSE)

  }
  attributes(res) <- attributes(x)  # set dim and dimnames
  res
}


#' @title Proper Scoring Rules for \code{hhh4ZI} Models
#' @description The following scores are implemented:
#' logarithmic score (logs), ranked probability score (rps),
#' Dawid-Sebastiani score (dss), squared error score (ses).
#' These are extended versions of the corresponding frunctions in \pkg{surveillance}
#' to handle zero-inflated negative binomial predictions,
#' which use an additional zero inflation parameter \code{gamma}.
#' @name scores
NULL

## adapted from surveillance:::scores.default()
.scores <- function(x, mu, size = NULL, gamma = NULL,
                           which = c("logs", "rps", "dss", "ses"),
                           sign = FALSE)
{
  ## compute individual scores (these have the same dimensions as x)
  scorelist <- lapply(X = setNames(nm = which), FUN = do.call,
                      args = alist(x = x, mu = mu, size = size, gamma = gamma),
                      envir = environment())

  ## append sign of x-mu
  if (sign)
    scorelist <- c(scorelist,
      list("sign" = if(is.null(gamma)) sign(x-mu) else sign(x - (1-gamma)*mu)))

  ## gather scores in an array
  simplify2array(scorelist, higher = TRUE)
}


### apply scoring rules to a set of oneStepAhead() forecasts
#' @inheritParams surveillance::scores.oneStepAhead
#' @param x an object of class \code{"oneStepAhead_hhh4ZI"} or \code{"hhh4ZI"},
#' respectively.
#' @param units integer or character vector indexing the units for which
#' to compute the scores. By default, all units are considered.
#' @param ... unused (argument of the generic).
#' @rdname scores
#' @importFrom surveillance scores
#' @export
scores.oneStepAhead_hhh4ZI <- function (x, which = c("logs","rps","dss","ses"),
                                        units = NULL, sign = FALSE, individual = FALSE,
                                        ...)
{
  chkDots(...)
  y <- x$observed  # observed counts during the prediction window
  mu <- x$mu
  ## transform overdispersion to dnbinom() parameterization
  size <- psi2size.oneStepAhead(x) # -> NULL or full dim(y) matrix
  gamma <- x$gamma
  ## select units
  if (!is.null(units)) {
    y <- y[,units,drop=FALSE]
    mu <- mu[,units,drop=FALSE]
    gamma <- gamma[,units,drop=FALSE]
    size <- size[,units,drop=FALSE] # works with size = NULL
  }
  nUnits <- ncol(y)
  if (nUnits == 1L)
    individual <- TRUE  # no need to apply rowMeans() below

  result <- .scores(x = y, mu = mu, size = size, gamma = gamma,
                           which = which, sign = sign)

  ## average over units if requested
  if (individual) {
    drop(result)
  } else {
    apply(X=result, MARGIN=3L, FUN=rowMeans)
    ## this gives a nrow(y) x (5L+sign) matrix (or a vector in case nrow(y)=1)
  }
}


## calculate scores with respect to fitted values
#' @param subset subset of time points for which to compute the scores.
#' @rdname scores
#' @importFrom surveillance scores
#' @export
scores.hhh4ZI <- function (x, which = c("logs","rps","dss","ses"),
                         subset = x$control$subset, units = seq_len(x$nUnit),
                         sign = FALSE, ...)
{
  chkDots(...)
  result <- .scores(
    x = x$stsObj@observed[subset, units, drop = FALSE],
    mu = x$mu[match(subset, x$control$subset), units, drop = FALSE],
    gamma = x$gamma[match(subset, x$control$subset), units, drop = FALSE],
    size = surveillance:::psi2size.hhh4(x, subset, units),
    which = which, sign = sign)
  rownames(result) <- subset
  drop(result)
}
