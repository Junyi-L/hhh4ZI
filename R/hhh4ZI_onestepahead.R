################################################################################
### The following are modified versions of functions from the surveillance package
### and wrappers around them.
### See below the original copyright declaration.
################################################################################

################################################################################
### Compute one-step-ahead predictions at a series of time points
###
### Copyright (C) 2011-2012 Michaela Paul, 2012-2018 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

#' @title Predictive Model Assessment for hhh4ZI Models
#' @description  The function oneStepAhead computes successive one-step-ahead
#'  predictions for a (random effects) HHH4ZI model fitted by \code{hhh4ZI.}
#' This function is the equivalent of \code{surveillance::oneStepAhead}
#' for model fits of class
#' \code{hhh4ZI}, obtained from \code{hhh4ZI}. The arguments are the
#' same as in \code{surveillance::oneStepAhead}.
#' Predictions can be inspected using the \code{quantile},
#' \code{confint} or \code{plot} methods.
#' @export
oneStepAhead <- function (result, tp, ...) UseMethod("oneStepAhead")

#' @rdname oneStepAhead
#' @export
oneStepAhead.hhh4 <- surveillance::oneStepAhead
formals(oneStepAhead.hhh4)$... <- substitute()  # consistent with the generic

#' @rdname oneStepAhead
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
oneStepAhead.hhh4ZI <- function(result, # hhh4ZI-object (i.e. a hhh4ZI model fit)
                                tp,     # scalar: one-step-ahead predictions for time
                                # points (tp+1):nrow(stsObj), or tp=c(from, to)
                                type = c("rolling", "first", "final"),
                                which.start = c("current", "final"), #if type="rolling"
                                keep.estimates = FALSE,
                                verbose = TRUE, # verbose-1 is used as verbose setting
                                # for sequentially refitted hhh4 models
                                cores = 1, # if which.start="final", the predictions
                                           # can be computed in parallel
                                ...)
{
  chkDots(...)  # there shouldn't be any
  stopifnot(inherits(result, "hhh4ZI"))
  type <- match.arg(type)
  if (type == "rolling" && !is.list(which.start)) {
    ## new in surveillance 1.10-0: if 'which.start' is a list, it is
    ## directly used as the 'start' argument for hhh4() in all time steps
    which.start <- match.arg(which.start)
    if (cores > 1 && which.start == "current")
      stop("no parallelization for 'type=\"rolling\"' ",
           "if 'which.start=\"current\"'")
  }

  ## get model terms
  model <- result[["terms"]]
  if (is.null(model))
    model <- result$terms <- terms(result)
  nTime <- model$nTime   # = nrow(result$stsObj)
  nUnits <- model$nUnits # = ncol(result$stsObj)
  dimPsi <- model$nOverdisp
  withPsi <- dimPsi > 0L
  psiIdx <- model$nFE + model$nd + seq_len(dimPsi)

  ## check that tp is within the time period of the data
  stopifnot(length(tp) %in% 1:2, tp >= 0)
  tpRange <- c(model$subset[1L], nTime-1L) # supported range
  if (any(tp > tpRange[2L]) || (type != "final" && any(tp < tpRange[1L]))) {
    stop("the time range defined by 'tp' must be a subset of ",
         tpRange[1L], ":", tpRange[2L])
  }
  if (length(tp) == 1) {
    tp <- c(tp, max(model$subset)-1L)  # historical default
    if (tp[1L] > tp[2L])  # probably unintended
      stop("'tp' larger than the default upper limit (", tp[2L], ")")
  }
  tps <- tp[1L]:tp[2L]  # this function actually works if tp[1] > tp[2]
  ntps <- length(tps)
  observed <- model$response[tps+1,,drop=FALSE]
  rownames(observed) <- tps+1

  ## adjust verbosity for model refitting
  verbose <- as.integer(verbose)
  result$control$verbose <- max(0, verbose - (ntps>1))
  if (type != "rolling" && verbose > 1L) verbose <- 1L
  do_pb <- verbose == 1L && interactive()

  ## initial fit
  fit <- if (type == "first") {
    if (do_pb)
      cat("\nRefitting model at first time point t =", tps[1L], "...\n")
    update.hhh4ZI(result, subset.upper = tps[1L], use.estimates = TRUE,
                  keep.terms = TRUE) # need "model" -> $terms
  } else result
  if (!fit$convergence) stop("initial fit did not converge")

  ## result templates (named and filled with NA's)
  pred <- matrix(NA_real_, nrow=ntps, ncol=nUnits,
                 dimnames=list(tps+1, colnames(observed)))
  gamma <- matrix(NA_real_, nrow=ntps, ncol=nUnits,
                  dimnames=list(tps+1, colnames(observed)))
  mu <- matrix(NA_real_, nrow=ntps, ncol=nUnits,
               dimnames=list(tps+1, colnames(observed)))

  if (withPsi)
    psi <- matrix(NA_real_, nrow=ntps, ncol=dimPsi,
                  dimnames=list(tps, names(model$initialTheta)[psiIdx]))
  if (keep.estimates) {
    coefficients <- matrix(NA_real_,
                           nrow=ntps, ncol=length(model$initialTheta),
                           dimnames=list(tps, names(model$initialTheta)))
    Sigma.orig <- matrix(NA_real_, nrow=ntps, ncol=model$nSigma,
                         dimnames=list(tps, names(result$Sigma.orig)))
    logliks <- matrix(NA_real_, nrow=ntps, ncol=2L,
                      dimnames=list(tps, c("loglikelihood", "margll")))
  }

  ## extract predictions and stuff for specific tp from fit
  getPreds <- function (fit, tp) {
    coefs <- unname(fit$coefficients)
    mu <- surveillance::meanHHH(coefs, fit$terms, subset=tp+1L, total.only=TRUE)
    gamma <- gammaZero(coefs, fit$terms, subset=tp+1L, d = 0)
    mean <- (1 - gamma) * mu

    c(list(pred = as.vector(mean),
           gamma = as.vector(gamma),
           mu = as.vector(mu)),
      if (withPsi) list(psi = coefs[psiIdx]),
      if (keep.estimates) list(
        coefficients=coefs,
        Sigma.orig=unname(fit$Sigma.orig),
        logliks=c(fit$loglikelihood, fit$margll))
    )
  }

  ## compute the predictions and save
  ## pred, psi, coefficients, Sigma.orig, and logliks
  if (cores > 1L) {

    ## return value template (unnamed NA vectors)
    resTemplate <- lapply(getPreds(fit, tps[1L]), "is.na<-", TRUE)

    ## run parallel
    res <- parallel::mclapply(tps, function (tp) {
      if (verbose)
        cat("One-step-ahead prediction @ t =", tp, "...\n")
      if (type == "rolling") { # update fit
        fit <- update.hhh4ZI(result, subset.upper=tp, use.estimates=TRUE,
                             start=if (is.list(which.start)) which.start,
                             verbose=FALSE, # chaotic in parallel
                             keep.terms=TRUE) # need "model" -> $terms
        if (!fit$convergence) {
          cat("WARNING: No convergence @ t =", tp, "!\n")
          return(resTemplate)
        }
      }
      getPreds(fit, tp)
    }, mc.preschedule=TRUE, mc.cores=cores)

    ## gather results
    .extractFromList <- function (what)
      t(vapply(res, "[[", resTemplate[[what]], what, USE.NAMES=FALSE))
    pred[] <- .extractFromList("pred")
    if (withPsi)
      psi[] <- .extractFromList("psi")
    if (keep.estimates) {
      coefficients[] <- .extractFromList("coefficients")
      Sigma.orig[] <- .extractFromList("Sigma.orig")
      logliks[] <- .extractFromList("logliks")
    }

  } else { ## sequential one-step ahead predictions

    if (do_pb) pb <- txtProgressBar(min=0, max=ntps, initial=0, style=3)
    for(i in seq_along(tps)) {
      if (verbose > 1L) {
        cat("\nOne-step-ahead prediction @ t =", tps[i], "...\n")
      } else if (do_pb) setTxtProgressBar(pb, i)

      if (type == "rolling") { # update fit
        fit.old <- fit # backup
        start <- if (is.list(which.start)) {
          which.start
        } else if (which.start == "current") surveillance:::hhh4coef2start(fit)
        ## else NULL
        fit <- update.hhh4ZI(result, subset.upper=tps[i],
                             start=start, # takes precedence
                             use.estimates=TRUE,
                             keep.terms=TRUE) # need "model" -> $terms
        if (!fit$convergence) {
          if (do_pb) cat("\n")
          cat("WARNING: No convergence @ t =", tps[i], "!\n")
          cat("WARNING: Last fit is used\n")
          ## FIXME: do a grid search ?
          fit <- fit.old
          #next
        }
      }

      res <- getPreds(fit, tps[i])

      ## gather results
      pred[i,] <- res$pred
      gamma[i,] <- res$gamma
      mu[i,] <- res$mu
      if (withPsi)
        psi[i,] <- res$psi
      if (keep.estimates) {
        coefficients[i,] <- res$coefficients
        Sigma.orig[i,] <- res$Sigma.orig
        logliks[i,] <- res$logliks
      }
    }
    if (do_pb) close(pb)

  }

  ## with shared overdispersion parameters we need to expand psi to ncol(pred)
  if (dimPsi > 1L && dimPsi != nUnits) {
    psi <- psi[,model$indexPsi,drop=FALSE]
  }

  ## done
  res <- c(list(pred = pred, mu = mu, gamma = gamma,
                observed = observed,
                psi = if (withPsi) psi else NULL,
                allConverged = all(!is.na(pred))),
           if (keep.estimates) list(coefficients = coefficients,
                                    Sigma.orig = Sigma.orig,
                                    logliks = logliks)
  )
  class(res) <- "oneStepAhead_hhh4ZI"
  res
}


## extract estimated overdispersion in dnbinom() parametrization, as full matrix
psi2size.oneStepAhead <- surveillance:::psi2size.oneStepAhead

# Quantiles of the one-step-ahead forecasts
#' @rdname oneStepAhead
#' @import stats
#' @export
quantile.oneStepAhead_hhh4ZI <- function (x, probs = c(2.5, 10, 50, 90, 97.5)/100, ...)
{
  stopifnot(is.vector(probs, mode = "numeric"), probs >= 0, probs <= 1,
            (np <- length(probs)) > 0)
  names(probs) <- paste(format(100*probs, trim=TRUE, scientific=FALSE, digits=3), "%")

  size <- psi2size.oneStepAhead(x)

  qs <- vapply(X = probs, FUN = VGAM::qzinegbin, FUN.VALUE = x$pred,
         munb = x$pred, size = size, pstr0 = x$gamma)

  ## one tp, one unit -> qs is a vector of length np
  ## otherwise, 'qs' has dimensions ntps x nUnit x np
  ## if nUnit==1, we return an ntps x np matrix, otherwise an array
  if (is.vector(qs)) {
    qs <- t(qs)
    rownames(qs) <- rownames(x$pred)
    qs
  } else if (dim(qs)[2L] == 1L) {
    matrix(qs, dim(qs)[1L], dim(qs)[3L], dimnames = dimnames(qs)[c(1L,3L)])
  } else qs
}

# Confidence intervals for one-step-ahead predictions
#' @rdname oneStepAhead
#' @export
confint.oneStepAhead_hhh4ZI <- function (object, parm, level = 0.95, ...)
{
  quantile.oneStepAhead_hhh4ZI(object, (1+c(-1,1)*level)/2, ...)
}

#  plot of one-step-ahead forecasts
#' @rdname oneStepAhead
#' @import graphics
#' @export
plot.oneStepAhead_hhh4ZI <- function (x, unit = 1, probs = 1:99/100,
                                      start = NULL, ...)
{
  stopifnot(length(unit) == 1, length(probs) > 1)

  ## select unit
  obs <- x$observed[,unit]
  ms <- x$pred[,unit]
  qs <- quantile.oneStepAhead_hhh4ZI(x, probs = probs)
  if (!is.matrix(qs))  # multi-unit predictions
    qs <- matrix(qs[,unit,], dim(qs)[1L], dim(qs)[3L],
                 dimnames = dimnames(qs)[c(1L,3L)])

  ## produce fanplot
  if (is.null(start))
    start <- as.integer(rownames(qs)[1L])
  surveillance::fanplot(quantiles = qs, probs = probs, means = ms,
          observed = obs, start = start, ...)
}

## pit
