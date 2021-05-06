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
### Non-randomized version of the PIT histogram as discussed in:
### Predictive model assessment for count data
### Czado, C., Gneiting, T. & Held, L. (2009)
### Biometrics 65:1254-1261
###
### Copyright (C) 2010-2012 Michaela Paul, 2013-2015,2017,2019 Sebastian Meyer
### $Revision$
### $Date$
################################################################################


## a convenient wrapper for Poisson / NegBin / ZI-NegBin predictions
.pit <- function (x, mu, size = NULL,gamma = NULL, ...)
{
  if (is.null(size)) {
    surveillance:::pit.default(x = x, pdistr = "ppois", lambda = mu, ...)
  } else if(is.null(gamma)){
    surveillance:::pit.default(x = x, pdistr = "pnbinom", mu = mu, size = size, ...)
  } else {
    surveillance:::pit.default(x = x, pdistr = VGAM::pzinegbin, munb = mu, size = size, pstr0 = gamma, ...)
  }
}


## pit-methods for oneStepAhead() predictions and "hhh4" fits
## (similar to the scores-methods)
#' @title Non-Randomized Version of the PIT Histogram (for Count Data)
#' @rdname pit
#' @importFrom surveillance pit
#' @export
pit.oneStepAhead_hhh4ZI <- function (x, units = NULL, ...)
{
  if (is.null(units)) {
    .pit(x = x$observed,
         mu = x$mu,
         gamma = x$gamma,
         size = psi2size.oneStepAhead(x), ...)
  } else if(is.null(gamma)){
    .pit(x = x$observed[, units, drop = FALSE],
         mu = x$pred[, units, drop = FALSE],
         size = psi2size.oneStepAhead(x)[, units, drop = FALSE],
         ...)
  } else{
    .pit(x = x$observed[, units, drop = FALSE],
         mu = x$mu[, units, drop = FALSE],
         gamma = x$gamma[, units, drop = FALSE],
         size = psi2size.oneStepAhead(x)[, units, drop = FALSE],
         ...)

  }
}
#' @rdname pit
#' @importFrom surveillance pit
#' @export
pit.hhh4ZI <- function (x, subset = x$control$subset, units = seq_len(x$nUnit), ...)
{
  .pit(x = x$stsObj@observed[subset, units, drop = FALSE],
       mu = x$fitted.values[match(subset, x$control$subset), units, drop = FALSE],
       size = surveillance:::psi2size.hhh4(x, subset, units),
       gamma = x$gamma[match(subset, x$control$subset), units, drop = FALSE],
       ...)
}
