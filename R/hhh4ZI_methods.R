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
### Standard methods for hhh4-fits
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2020 Sebastian Meyer
### $Revision$
### $Date$
################################################################################


#' @title Refit a \code{hhh4ZI} Model
#' @description Re-fit a \code{\link{hhh4ZI}} model with a modified control list,
#' equivalently to \code{\link[surveillance]{update.hhh4}}.
#' @param object a fitted \code{"\link{hhh4ZI}"} model.
#' @param ... elements modifying the original control list.
#' @param S a named list of numeric vectors (with names \code{"ar"},
#' \code{"end"}, etc) to adjust the number of harmonics in the model components
#' via \code{\link{addSeason2formula}}, or \code{NULL} (meaning no modification
#' of seasonal terms).
#' @param subset.upper refit on a subset of the data up to that time point
#' (used by \code{\link{oneStepAhead}}).
#' @param use.estimates logical indicating if \code{coef(object)} should be
#' used as starting values for the new fit.
#' @param evaluate logical indicating if the updated \code{hhh4ZI} call should
#' be evaluated or returned.
#' @importFrom utils modifyList
#' @export

update.hhh4ZI <- function (object, ..., S = NULL, subset.upper = NULL,
                           use.estimates = object$convergence, evaluate = TRUE)
{
  control <- object$control

  ## first modify the control list according to the components in ...
  extras <- list(...)
  control <- modifyList(control, extras)

  ## adjust start values
  control$start <- if (use.estimates) { # use parameter estimates
    surveillance:::hhh4coef2start(object)
  } else local({ # re-use previous 'start' specification
    ## for pre-1.8-2 "hhh4" objects,
    ## object$control$start is not necessarily a complete list:
    template <- eval(formals(hhh4ZI)$control$start)
    template[] <- object$control$start[names(template)]
    template
  })
  ## and update according to an extra 'start' argument
  if (!is.null(extras[["start"]])) {
    if (!is.list(extras$start) || is.null(names(extras$start))) {
      stop("'start' must be a named list, see 'help(\"hhh4\")'")
    }
    control$start[] <- mapply(
      FUN = function (now, extra) {
        if (is.null(names(extra))) {
          extra
        } else { # can retain non-extra values
          now[names(extra)] <- extra
          now
        }
      },
      control$start, extras$start[names(control$start)],
      SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
  }
  ## update initial values of parametric weight function
  if (use.estimates && length(coefW <- surveillance::coefW(object)) &&
      ! "weights" %in% names(extras$ne)) { # only if function is unchanged
    control$ne$weights$initial <- coefW
  }

  ## adjust seasonality
  if (!is.null(S)) {
    stopifnot(is.list(S), !is.null(names(S)),
              names(S) %in% c("ar", "ne", "end", "zi"))
    control[names(S)] <- mapply(function (comp, S) {
      comp$f <- surveillance::addSeason2formula(
        surveillance:::removeSeasonFromFormula(comp$f),
        period = object$stsObj@freq, S = S)
      comp
    }, control[names(S)], S, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  }

  ## use a different time range of the data (only changing the end)
  ## Note: surveillance < 1.15.0 disallowed subset.upper > max(control$subset)
  if (surveillance:::isScalar(subset.upper)) {
    if (subset.upper < control$subset[1L])
      stop("'subset.upper' is smaller than the lower bound of 'subset'")
    control$subset <- control$subset[1L]:subset.upper
  }

  ## fit the updated model or just return the modified control list
  if (evaluate) {
    hhh4ZI(object$stsObj, control)
  } else {
    control
  }
}

#' @export
summary.hhh4ZI <- function (object, maxEV = FALSE, ...)
{
  s <- NextMethod(maxEV = FALSE)
  if (maxEV && inherits(s, "summary.hhh4"))
      ## use hhh4ZI-specific createLambda()
      s$maxEV_range <- unique(range(getMaxEV(object), na.rm = TRUE))
  s
}
