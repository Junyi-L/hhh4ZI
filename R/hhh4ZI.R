################################################################################
### The following are modified versions of functions from the surveillance package
### and wrappers around them.
### See below the original copyright declaration.
################################################################################

################################################################################
### Copyright declaration from the surveillance package:
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### hhh4 is an extended version of algo.hhh for the sts-class
### The function allows the incorporation of random effects and covariates.
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2016,2019-2021 Sebastian Meyer
### $Revision$
### $Date$
################################################################################
ADVICEONERROR <- "\n  Try different starting values, more iterations, or another optimizer.\n"

#' @title Fitting zero-inflated HHH Models with Random Effects
#' and Neighbourhood Structure
#'
#' @description Fits a zero-inflated autoregressive negative binomial (\code{hhh4}) model
#' to a univariate or multivariate time series of counts.
#' The characteristic feature of \code{hhh4} models is the additive
#' decomposition of the conditional mean into \emph{epidemic} and
#' \emph{endemic} components (Held et al, 2005).
#' The inflated parameter is a logit-linear predictor and can have autoregressive terms.
#'
#' @param stsObj object of class \code{"\linkS4class{sts}"} containing the (multivariate)
#' count data time series.
#' @param control a list containing the model specification and control arguments,
#'  the parts relating to \code{hhh4} model is the same as in \code{surveillance::hhh4}:
#' \itemize{
#' \item{\code{ar}}{
#' Model for the autoregressive component given as
#' list with the following components:
#' \itemize{
#' \item{f = ~ -1}{
#' a formula specifying
#' \eqn{\log(\lambda_{it})}{log(\lambda_it)}.}
#' \item{offset = 1}{
#' optional multiplicative offset, either 1 or
#' a matrix of the same dimension as \code{observed(stsObj)}.}
#' \item{lag = 1}{
#' a positive integer meaning autoregression on
#' \eqn{y_{i,t-lag}}.}
#' }
#' }
#' \item{\code{ne}}{
#' Model for the neighbour-driven component given as
#'list with the following components:
#'  \itemize{
#'   \item{f = ~ -1}{
#'   a formula specifying \eqn{\log(\phi_{it})}{log(\phi_it)}}.
#'   \item{offset = 1}{
#'   optional multiplicative offset, either 1 or
#'      a matrix of the same dimension as \code{observed(stsObj)}.}
#'    \item{lag = 1}{
#'    a non-negative integer meaning dependency on
#'     \eqn{y_{j,t-lag}}.}
#'   \item{weights = neighbourhood(stsObj) == 1}{
#'      neighbourhood weights \eqn{w_{ji}}{w_ji}. The default
#'      corresponds to the original formulation by Held et al
#'      (2005), i.e., the spatio-temporal component incorporates an
#'      unweighted sum over the lagged cases of the first-order
#'      neighbours. See Paul et al (2008) and Meyer and Held (2014)
#'      for alternative specifications, e.g.,
#'      \code{\link{W_powerlaw}}.
#'      Time-varying weights are possible by specifying an
#'      array of \code{dim()} \code{c(nUnits, nUnits, nTime)}, where
#'      \code{nUnits=ncol(stsObj)} and \code{nTime=nrow(stsObj)}.}
#'    \item{scale = NULL}{
#'      optional matrix of the same dimensions as \code{weights} (or
#'     a vector of length \code{ncol(stsObj)}) to scale the
#'      \code{weights} to \code{scale * weights}.
#'    }
#'    \item{normalize = FALSE}{
#'     logical indicating if the (scaled) \code{weights} should be
#'      normalized such that each row sums to 1.
#'    }
#'  }
#'  }
#' \item{\code{end}}{
#' Model for the endemic component given as list
#' with the following components
#' \itemize{
#' \item{f = ~ 1 }{
#' a formula specifying \eqn{\log(\nu_{it})}{log(\nu_it)}}.
#' \item{offset = 1 }{
#' optional multiplicative offset \eqn{e_{it}}{e_it},
#' either 1 or a matrix of the same dimension as \code{observed(stsObj)}.}
#' }}
#' \item{\code{zi}}{
#' Model for the zero inflation component given as
#' list with the following components:
#' \itemize{
#' \item{f = ~ -1}{
#' a formula specifying
#' \eqn{\logit(\gamma_{it})}{logit(\gamma_it)}.}
#' \item{lag = 1}{
#' a positive integer or vector meaning autoregression on
#' \eqn{y_{i,t-lag}}}
#' \item{lag.unitSpecific}{ logical indicating if the autoregressive parameter
#' in the zero inflation part is unit specific. }
#' }
#' }
#'\item{\code{family}}{
#'Distributional family --
#' the Negative Binomial distribution. The
#' overdispersion parameter can be assumed to be the same for all
#' units (\code{"NegBin1"}), to vary freely over all units
#'(\code{"NegBinM"}), or to be shared by some units (specified by
#' a factor of length \code{ncol(stsObj)} such that its number of
#' levels determines the number of overdispersion parameters).
#' Note that \code{"NegBinM"} is equivalent to
#' \code{factor(colnames(stsObj), levels = colnames(stsObj))}.
#' }
#' \item{\code{subset}}{
#' Typically \code{2:nrow(obs)} if model contains
#' autoregression.}
#' \item{\code{optimizer}}{
#' a list of three lists of control arguments.
#'
#' The \code{"stop"} list specifies two criteria for the outer
#' optimization of regression and variance parameters: the relative
#' \code{tol}erance for parameter change using the criterion
#' \code{max(abs(x[i+1]-x[i])) / max(abs(x[i]))},
#' and the maximum number \code{niter} of outer iterations.
#'
#' Control arguments for the single optimizers are specified in the
#' lists named \code{"regression"} and \code{"variance"}.
#' \code{method="nlminb"} is the default optimizer for regression update,
#'  however, the \code{method}s from \code{\link{optim}} may also be specified
#' (as well as \code{"\link{nlm}"} but that one is not recommended here).
#' For the variance updates, only Nelder-Mead optimization
#' (\code{method="Nelder-Mead"}) is provided.
#' All other elements of these two lists are passed as
#' \code{control} arguments to the chosen \code{method}, e.g., if
#' \code{method="nlminb"} adding \code{iter.max=50} increases the
#' maximum number of inner iterations from 20 (default) to 50.
#' }
#'
#' \item{\code{verbose}}{
#' non-negative integer (usually in the range
#' \code{0:3}) specifying the amount of tracing information to be
#' output during optimization.}
#'
#' \item{\code{start}}{
#' a list of initial parameter values replacing
#' initial values set via \code{\link{fe}} and \code{\link{ri}}.
#' Since \pkg{surveillance} 1.8-2, named vectors are matched
#' against the coefficient names in the model (where unmatched
#' start values are silently ignored), and need not be complete,
#' e.g., \code{start = list(fixed = c("-log(overdisp)" = 0.5))}
#' (default: 2) for a \code{family = "NegBin1"} model.
#' In contrast, an unnamed start vector must specify the full set
#' of parameters as used by the model.}
#'
#' \item{\code{data}}{
#' a named list of covariates that are to be
#' included as fixed effects (see \code{\link{fe}}) in any of the 3
#' component formulae.
#' By default, the time variable \code{t} is available and used for
#' seasonal effects created by \code{\link{addSeason2formula}}.
#' In general, covariates in this list can be either vectors of
#' length \code{nrow(stsObj)} interpreted as time-varying but
#' common across all units, or matrices of the same dimension as
#' the disease counts \code{observed(stsObj)}.}
#' \item{\code{keep.terms}}{ logical indicating if the terms object
#' used in the fit is to be kept as part of the returned object.
#' This is usually not necessary, since the terms object is
#' reconstructed by the \code{\link{terms}}-method for class
#' \code{"hhh4ZI"} if necessary (based on \code{stsObj} and
#' \code{control}, which are both part of the returned
#' \code{"hhh4ZI"} object).}
#' }
#' @param check.analyticals {logical (or a subset of
#' \code{c("numDeriv", "maxLik")}), indicating if (how) the implemented
#' analytical score vector and Fisher information matrix should be
#' checked against numerical derivatives at the parameter starting values,
#' using the packages \pkg{numDeriv} and/or \pkg{maxLik}. If activated,
#' \code{hhh4} will return a list containing the analytical and numerical
#' derivatives for comparison (no ML estimation will be performed).
#' This is mainly intended for internal use by the package developers.}
#'
#' @return \code{hhh4ZI} returns an object of class \code{"hhh4ZI"},
#' which inherits from class \code{"hhh4"}, and
#' is a list containing the following components:
#' \itemize{
#' \item{coefficients}{
#' named vector with estimated (regression) parameters of the model}
#' \item{se}{
#' estimated standard errors (for regression parameters)}
#' \item{cov}{
#' covariance matrix (for regression parameters)}
#' \item{Sigma}{
#' estimated variance-covariance matrix of random effects}
#' \item{Sigma.orig}{
#' estimated variance parameters on internal scale used for optimization}
#' \item{call}{ the matched call }
#' \item{dim}{ vector with number of fixed and random effects in the model }
#' \item{loglikelihood}{ (penalized) loglikelihood evaluated at the MLE}
#' \item{margll}{ (approximate) log marginal likelihood should the model contain random effects  }
#' \item{convergence}{ logical. Did optimizer converge?}
#' \item{mu}{ The fitted mean values in hhh4 model part}
#' \item{fitted.values}{ fitted mean values in zero-inflated model}
#' \item{gamma}{ fitted zero inflation parameter}
#' \item{control}{ control object of the fit}
#' \item{terms}{ the terms object used in the fit if \code{keep.terms = TRUE}
#'   and \code{NULL} otherwise}
#' \item{stsObj}{ the supplied \code{stsObj} }
#' \item{lags}{ named integer vector of length three containing the lags
#'   used for the epidemic components \code{"ar"}, \code{"ne"} and \code{"zi"}
#'   respectively. The corresponding lag is \code{NA} if the component
#'   was not included in the model.}
#' \item{nObs}{ number of observations used for fitting the model}
#' \item{nTime}{ number of time points used for fitting the model }
#' \item{nUnit}{ number of units (e.g. areas) used for fitting the model}
#' \item{runtime}{ the \code{\link{proc.time}}-queried time taken
#'   to fit the model, i.e., a named numeric vector of length 5 of class
#'   \code{"proc_time"}}
#' }
#' @examples
#' data("measles", package = "hhh4ZI")
#' library(surveillance)
#' measles <- aggregate(measles, by = "time", nfreq = 26)
#' adjmat <- poly2adjmat(measles@map)
#' neW1 <- adjmat/colSums(adjmat)
#' fit <- hhh4ZI(measles,
#' control = list(
#'   ar = list(f = ~1,
#'             offset = 1,
#'             lag = 1),
#'   ne = list(f = ~1,
#'             offset = 1,
#'             lag = 1,
#'            weights = neW1,
#'            scale = NULL,
#'            normalize = TRUE),
#'  end = list(f = ~1,
#'              offset = 1),
#'   zi = list(f = ~1,
#'             lag = 1
#'  ),
#'   family = "NegBin1",
#'   optimizer = list(stop = list(tol = 1e-5, niter = 100),
#'                    regression = list(method = "nlminb")),
#'   verbose = TRUE,
#'   keep.terms = TRUE
#' )
#' )
#' summary(fit)
#' sim_data <- simulate(fit, simplify = FALSE)
#' @export
hhh4ZI <- function (object, control, ...) UseMethod("hhh4ZI")

#' @rdname hhh4ZI
#' @export
hhh4ZI.sts <- function(stsObj,
                       control = list(
                         ar = list(f = ~ -1,        # a formula "exp(x'lamba)*y_t-lag" (ToDo: matrix)
                                   offset = 1,      # multiplicative offset
                                   lag = 1),        # autoregression on y_i,t-lag
                         ne = list(f = ~ -1,        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
                                   offset = 1,      # multiplicative offset
                                   lag = 1,         # regression on y_j,t-lag
                                   weights = neighbourhood(stsObj) == 1,  # weights w_ji
                                   scale = NULL,    # such that w_ji = scale * weights
                                   normalize = FALSE), # w_ji -> w_ji / rowSums(w_ji), after scaling
                         end = list(f = ~ 1,        # a formula "exp(x'nu) * n_it"
                                    offset = 1),    # optional multiplicative offset e_it
                         zi = list(f = ~-1,
                                   lag = 1,
                                   lag.unitSpecific = FALSE# can be a scalar or vector
                         ),
                         family = c("NegBin1", "NegBinM"), # or a factor of length nUnit for Negbin
                         subset = 2:nrow(stsObj),   # epidemic components require Y_{t-lag}
                         optimizer = list(stop = list(tol = 1e-5, niter = 100), # control arguments
                                          regression = list(method = "nlminb"), # for optimization
                                          variance = list(method = "Nelder-Mead")),
                         verbose = FALSE,           # level of reporting during optimization
                         start = list(fixed = NULL, # list of start values, replacing initial
                                      random = NULL,  # values from fe() and ri() in 'f'ormulae
                                      sd.corr = NULL),
                         data = list(t = stsObj@epoch - min(stsObj@epoch)), # named list of covariates
                         keep.terms = FALSE  # whether to keep interpretControl(control, stsObj)
                       ), check.analyticals = FALSE,
                       ... # unused (argument of the generic)
){
  ptm <- proc.time()
  ## check control and set default values (for missing arguments)
  control <- setControl(control, stsObj)

  ## get model terms
  model <- interpretControl(control, stsObj)
  dimFixedEffects <- model$nFE + model$nd + model$nOverdisp
  dimRandomEffects <- model$nRE

  ## starting values
  #* -> better default values possible
  theta.start <- model$initialTheta
  Sigma.start <- model$initialSigma

  # check if initial values are valid
  # CAVE: there might be NA's in mu if there are missing values in Y
  mu <- surveillance::meanHHH(theta.start, model, total.only=TRUE)
  if(any(mu==0, na.rm=TRUE) || any(is.infinite(mu)))
    stop("some mean is degenerate (0 or Inf) at initial values")

  gamma <- gammaZero(theta.start, model)
  if(any(gamma == 0, na.rm=TRUE) || any(gamma == 1, na.rm=TRUE))
    stop("some gamma is degenerate (0 or 1) at initial values")

  ## check score vector and fisher information at starting values
  check.analyticals <- if (isTRUE(check.analyticals)) {
    if (length(theta.start) > 50) "maxLik" else "numDeriv"
  } else if (is.character(check.analyticals)) {
    match.arg(check.analyticals, c("numDeriv", "maxLik"), several.ok=TRUE)
  } else NULL
  if (length(check.analyticals) > 0L) {
    resCheck <- checkAnalyticals(model, theta.start, Sigma.start,
                                 methods=check.analyticals)
    return(resCheck)
  }
  #-------------------------------------------------------------------
  ## maximize loglikelihood (penalized and marginal)
  ## maximize loglikelihood (penalized and marginal)
  myoptim <- fitHHH4ZI(theta=theta.start,sd.corr=Sigma.start, model=model,
                       cntrl.stop       = control$optimizer$stop,
                       cntrl.regression = control$optimizer$regression,
                       cntrl.variance   = control$optimizer$variance,
                       verbose=control$verbose)


  ## extract parameter estimates
  convergence <- myoptim$convergence == 0
  thetahat <- myoptim$theta

  # fitted value
  mu <- surveillance::meanHHH(thetahat, model, total.only=TRUE)
  gamma <- gammaZero(thetahat, model, subset = model$subset, d = 0)
  mean <- (1 - gamma) * mu

  if (dimRandomEffects>0) {
    Sigma.orig <- myoptim$sd.corr
    Sigma.trans <- getSigmai(head(Sigma.orig,model$nVar),
                             tail(Sigma.orig,model$nCorr),
                             model$nVar)
    dimnames(Sigma.trans) <-
      rep.int(list(sub("^sd\\.", "",
                       names(Sigma.orig)[seq_len(model$nVar)])), 2L)
  } else {
    Sigma.orig <- Sigma.trans <- NULL
  }

  # compute covariance matrices of regression and variance parameters
  cov <- try(solve(myoptim$fisher), silent=TRUE)
  Sigma.cov <- if(dimRandomEffects>0) try(solve(myoptim$fisherVar), silent=TRUE)

  # check for degenerate fisher info
  if(inherits(cov, "try-error")){ # fisher info is singular
    if (control$verbose)
      cat("WARNING: Final Fisher information matrix is singular!\n")
    convergence <- FALSE
  } else if(any(!is.finite(diag(cov))) || any(diag(cov)<0)){
    if (control$verbose)
      cat("WARNING: non-finite or negative covariance of regression parameters!\n")
    convergence <- FALSE
  }
  if (!convergence) {
    if (control$verbose) {
      cat("Penalized loglikelihood =", myoptim$loglik, "\n")
      thetastring <- paste(round(thetahat,2), collapse=", ")
      thetastring <- strwrap(thetastring, exdent=10, prefix="\n", initial="")
      cat("theta = (", thetastring, ")\n")
    }
    warning("Results are not reliable!",
            if (any(surveillance:::splitParams(thetahat, model)$overdisp > 10)) { # FALSE for Poisson
              "\n  Overdispersion parameter close to zero; maybe try a Poisson model.\n"
            } else ADVICEONERROR)
  }

  ## gather results in a list -> "hhh4" object
  result <- list(coefficients = thetahat,
                 se=if (convergence) sqrt(diag(cov)), cov=cov,
                 Sigma=Sigma.trans,     # estimated covariance matrix of ri's
                 Sigma.orig=Sigma.orig, # variance parameters on original scale
                 Sigma.cov=Sigma.cov,   # covariance matrix of Sigma.orig
                 call=match.call(),
                 dim=c(fixed=dimFixedEffects,random=dimRandomEffects),
                 loglikelihood=myoptim$loglik, margll=myoptim$margll,
                 convergence=convergence,
                 mu = mu,
                 fitted.values = mean,
                 gamma = gamma,
                 control=control,
                 terms=if(control$keep.terms) model else NULL,
                 stsObj=stsObj,
                 lags=sapply(control[c("ar","ne")], function (comp)
                   if (comp$inModel) comp$lag else NA_integer_),
                 nObs=sum(!model$isNA[control$subset,]),
                 nTime=length(model$subset), nUnit=ncol(stsObj),
                 ## CAVE: nTime is not nrow(stsObj) as usual!
                 runtime=proc.time()-ptm)
  if (!convergence) {
    ## add (singular) Fisher information for further investigation
    result[c("fisher","fisherVar")] <- myoptim[c("fisher","fisherVar")]
  }
  class(result) <- c("hhh4ZI","hhh4")
  return(result)

}

setControl <- function (control, stsObj)
{
  stopifnot(is.list(control))
  nTime <- nrow(stsObj)
  nUnit <- ncol(stsObj)
  if(nTime <= 2) stop("too few observations")

  ## arguments in 'control' override any corresponding default arguments
  defaultControl <- eval(formals(hhh4ZI.sts)$control)
  environment(defaultControl$ar$f) <- environment(defaultControl$ne$f) <-
    environment(defaultControl$end$f) <-
    environment(defaultControl$zi$f)<- .GlobalEnv
  control <- modifyList(defaultControl, control)


  ## check that component specifications are list objects
  for (comp in c("ar", "ne", "end", "zi")) {
    if(!is.list(control[[comp]])) stop("'control$", comp, "' must be a list")
  }

  ## check lags in "ar" and "ne" components
  for (comp in c("ar", "ne")) {
    if (!surveillance:::isScalar(control[[comp]]$lag) || control[[comp]]$lag < (comp=="ar"))
      stop("'control$", comp, "$lag' must be a ",
           if (comp=="ar") "positive" else "non-negative", " integer")
    control[[comp]]$lag <- as.integer(control[[comp]]$lag)
  }

  # chek lags in "zi" component, it can be a vector
  if (!is.vector(control$zi$lag, mode = "numeric") & !is.null(control$zi$lag))
    stop("'control$zi$lag' must be a non-negative integer")
  control[["zi"]]$lag <- as.integer(control[["zi"]]$lag)

  ### check AutoRegressive component

  if (control$ar$isMatrix <- is.matrix(control$ar$f)) {
    ## this form is not implemented -> will stop() in interpretControl()
    if (any(dim(control$ar$f) != nUnit))
      stop("'control$ar$f' must be a square matrix of size ", nUnit)
    if (is.null(control$ar$weights)) { # use identity matrix
      control$ar$weights <- diag(nrow=nUnit)
    } else if (!is.matrix(control$ar$weights) ||
               any(dim(control$ar$weights) != nUnit)) {
      stop("'control$ar$weights' must be a square matrix of size ", nUnit)
    }
    control$ar$inModel <- TRUE
  } else if (inherits(control$ar$f, "formula")) {
    if (!is.null(control$ar$weights)) {
      warning("argument 'control$ar$weights' is not used")
      control$ar$weights <- NULL
    }
    # check if formula is valid
    control$ar$inModel <- surveillance:::isInModel(control$ar$f)
  } else {
    stop("'control$ar$f' must be either a formula or a matrix")
  }


  ### check NEighbourhood component

  if (!inherits(control$ne$f, "formula"))
    stop("'control$ne$f' must be a formula")
  control$ne$inModel <- surveillance:::isInModel(control$ne$f)

  if (control$ne$inModel) {
    if (nUnit == 1)
      stop("\"ne\" component requires a multivariate 'stsObj'")
    ## if ar$f is a matrix it includes neighbouring units => no "ne" component
    if (control$ar$isMatrix)
      stop("there must not be an extra \"ne\" component ",
           "if 'control$ar$f' is a matrix")
    ## check ne$weights specification
    surveillance:::checkWeights(control$ne$weights, nUnit, nTime,
                                neighbourhood(stsObj), control$data,
                                check0diag = control$ar$inModel)
    ## check optional scaling of weights
    if (!is.null(control$ne$scale)) {
      stopifnot(is.numeric(control$ne$scale))
      if (is.vector(control$ne$scale)) {
        stopifnot(length(control$ne$scale) == 1L ||
                    length(control$ne$scale) %% nUnit == 0,
                  !is.na(control$ne$scale))
      } else {
        surveillance:::checkWeightsArray(control$ne$scale, nUnit, nTime)
      }
    }
  } else {
    control$ne[c("weights", "scale", "normalize")] <- list(NULL, NULL, FALSE)
  }


  ### check ENDemic component

  if (!inherits(control$end$f, "formula"))
    stop("'control$end$f' must be a formula")
  control$end$inModel <- surveillance:::isInModel(control$end$f)


  ### check zero component
  if (!inherits(control$zi$f, "formula"))
    stop("'control$zi$f' must be a formula")
  control$zi$inModel <- surveillance:::isInModel(control$zi$f)


  ### check offsets

  for (comp in c("ar", "ne", "end")) {
    if (is.matrix(control[[comp]]$offset) && is.numeric(control[[comp]]$offset)){
      if (!identical(dim(control[[comp]]$offset), dim(stsObj)))
        stop("'control$",comp,"$offset' must be a numeric matrix of size ",
             nTime, "x", nUnit)
      if (any(is.na(control[[comp]]$offset)))
        stop("'control$",comp,"$offset' must not contain NA values")
    } else if (!identical(as.numeric(control[[comp]]$offset), 1)) {
      stop("'control$",comp,"$offset' must either be 1 or a numeric ",
           nTime, "x", nUnit, " matrix")
    }
  }


  ### stop if no component is included in the model

  if (length(comps <- componentsHHH4ZI(list(control=control))) == 0L)
    stop("none of the components 'ar', 'ne', 'end', 'zi' is included in the model")


  ### check remaining components of the control list

  if (is.factor(control$family)) {
    stopifnot(length(control$family) == nUnit)
    ## guard against misuse as family = factor("Poisson"), e.g., if taken
    ## from a data.frame of control options with "stringsAsFactors"
    if (nUnit == 1 && as.character(control$family) %in% defaultControl$family) {
      control$family <- as.character(control$family)
      warning("'family = factor(\"", control$family, "\")' is interpreted ",
              "as 'family = \"", control$family, "\"'")
    } else {
      control$family <- droplevels(control$family)
      names(control$family) <- colnames(stsObj)
    }
  } else {
    control$family <- match.arg(control$family, defaultControl$family)
    if(!(control$family %in% c("NegBin1", "NegBinM")))
      stop("control$family must be either NegBin1, NegBinM, or a factor")
  }

  if (!is.vector(control$subset, mode="numeric") ||
      !all(control$subset %in% seq_len(nTime)))
    stop("'control$subset' must be %in% 1:", nTime)
  lags <- c(ar = control$ar$lag, ne = control$ne$lag,
            zi = suppressWarnings(max(control$zi$lag)))
  maxlag <- suppressWarnings(max(lags[names(lags) %in% comps])) # could be -Inf
  if (control$subset[1L] <= maxlag) {
    warning("'control$subset' should be > ", maxlag, " due to epidemic lags")
  }

  if (!is.list(control$optimizer) ||
      any(! sapply(c("stop", "regression", "variance"),
                   function(x) is.list(control$optimizer[[x]]))))
    stop("'control$optimizer' must be a list of lists")

  control$verbose <- as.integer(control$verbose)
  if (length(control$verbose) != 1L || control$verbose < 0)
    stop("'control$verbose' must be a logical or non-negative numeric value")

  stopifnot(is.list(control$start))
  control$start <- local({
    defaultControl$start[] <- control$start[names(defaultControl$start)]
    defaultControl$start
  })
  if (!all(vapply(X = control$start,
                  FUN = function(x) is.null(x) || is.vector(x, mode="numeric"),
                  FUN.VALUE = TRUE, USE.NAMES = FALSE)))
    stop("'control$start' must be a list of numeric start values")

  stopifnot(length(control$keep.terms) == 1L, is.logical(control$keep.terms))

  ## Done
  return(control)
}

componentsHHH4ZI <- function (object)
  names(which(sapply(object$control[c("ar", "ne", "end","zi")], "[[", "inModel")))

AR <- function(lag){
  stsObj <- get("stsObj", envir = parent.frame(1), inherits = TRUE)
  Y <- surveillance::observed(stsObj)
  rbind(matrix(NA_real_, lag, ncol(Y)),
        Y[seq_len(nrow(Y) - lag),, drop = FALSE])
}
# interpret and check the specifications of each component
# control must contain all arguments, i.e. setControl was used
interpretControl <- function (control, stsObj)
{
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)

  Y <- observed(stsObj)


  ##########################################################################
  ##  get the model specifications for each of the three components
  ##########################################################################
  ar <- control$ar
  ne <- control$ne
  end <- control$end
  zi <- control$zi
  ziTrue <- !is.null(zi)

  ## for backwards compatibility with surveillance < 1.8-0, where the ar and ne
  ## components of the control object did not have an offset
  if (is.null(ar$offset)) ar$offset <- 1
  if (is.null(ne$offset)) ne$offset <- 1
  ## for backward compatibility with surveillance < 1.9-0
  if (is.null(ne$normalize)) ne$normalize <- FALSE

  ## create list of offsets of the three components
  Ym1 <- rbind(matrix(NA_integer_, ar$lag, nUnits), head(Y, nTime-ar$lag))
  Ym1.ne <- surveillance:::neOffsetFUN(Y, ne$weights, ne$scale, ne$normalize,
                                       neighbourhood(stsObj), control$data, ne$lag, ne$offset)
  offsets <- list(ar=ar$offset * Ym1, ne = Ym1.ne,
                  end = end$offset)
  ## -> offsets$zi and offsets$ar are offset * Y_t-lag
  ## -> offsets$ne is a function of the parameter vector 'd', which returns a
  ##    nTime x nUnits matrix -- or 0 (scalar) if there is no NE component
  ## -> offsets$end might just be 1 (scalar)

  ## Initial parameter vector 'd' of the neighbourhood weight function
  initial.d <- if (is.list(ne$weights)) ne$weights$initial else numeric(0L)
  dim.d <- length(initial.d)
  names.d <- if (dim.d == 0L) character(0L) else {
    paste0("neweights.", if (is.null(names(initial.d))) {
      if (dim.d==1L) "d" else paste0("d", seq_len(dim.d))
    } else names(initial.d))
  }

  ## determine all NA's
  isNA <- is.na(Y)
  if (ar$inModel)
    isNA <- isNA | is.na(offsets[[1L]])
  if (ne$inModel)
    isNA <- isNA | is.na(offsets[[2L]](initial.d))

  ## get terms for all components
  all.term <- NULL
  if(ar$isMatrix) stop("matrix-form of 'control$ar$f' is not implemented")
  if(ar$inModel) # ar$f is a formula
    all.term <- cbind(all.term,
                      checkFormula(ar$f, 1, control$data, stsObj))
  if(ne$inModel)
    all.term <- cbind(all.term,
                      checkFormula(ne$f, 2, control$data, stsObj))
  if(end$inModel)
    all.term <- cbind(all.term,
                      checkFormula(end$f,3, control$data, stsObj))


  zi.formula <- if(zi$lag[1] >0){

    update.formula(zi$f,
                   as.formula(paste("~ . +",
                                    paste0("fe(AR(", zi$lag, ")",
                                           ", unitSpecific =", zi$lag.unitSpecific, ")",
                                           collapse = " + "),
                                    collapse = "+"))
    )} else zi$f

  if(zi$inModel)
    all.term <- cbind(all.term,
                      checkFormula(zi.formula, 4, control$data, stsObj))

  dim.fe <- sum(unlist(all.term["dim.fe",]))
  dim.re.group <- unlist(all.term["dim.re",], use.names=FALSE)
  dim.re <- sum(dim.re.group)
  dim.var <- sum(unlist(all.term["dim.var",]))
  dim.corr <- sum(unlist(all.term["corr",]))

  if(dim.corr>0){
    if(dim.var!=dim.corr) stop("Use corr=\'all\' or corr=\'none\' ")
    dim.corr <- switch(dim.corr,0,1,3,6)
  }

  # the vector with dims of the random effects must be equal if they are correlated
  if(length(unique(dim.re.group[dim.re.group>0]))!=1 & dim.corr>0){
    stop("Correlated effects must have same penalty")
  }

  n <- c("ar","ne","end", "zi")[unlist(all.term["offsetComp",])]
  names.fe <- names.var <- names.re <- character(0L)
  for(i in seq_along(n)){
    .name <- all.term["name",i][[1]]
    names.fe <- c(names.fe, paste(n[i], .name, sep="."))
    if(all.term["random",i][[1]]) {
      names.var <- c(names.var, paste("sd", n[i], .name, sep="."))
      names.re <- c(names.re, paste(n[i], .name, if (.name == "ri(iid)") {
        colnames(stsObj)
      } else {
        seq_len(all.term["dim.re",i][[1]])
      }, sep = "."))
    }
  }
  index.fe <- rep(1:ncol(all.term), times=unlist(all.term["dim.fe",]))
  index.re <- rep(1:ncol(all.term), times=unlist(all.term["dim.re",]))

  # poisson or negbin model
  if(identical(control$family, "Poisson")){
    ddistr <- function(y,mu,size){
      dpois(y, lambda=mu, log=TRUE)
    }
    dim.overdisp <- 0L
    index.overdisp <- names.overdisp <- NULL
  } else { # NegBin
    ddistr <- function(y,mu,size){
      dnbinom(y, mu=mu, size=size, log=TRUE)
    }
    ## version that can handle size = Inf (i.e. the Poisson special case):
    ## ddistr <- function (y,mu,size) {
    ##     poisidx <- is.infinite(size)
    ##     res <- y
    ##     res[poisidx] <- dpois(y[poisidx], lambda=mu[poisidx], log=TRUE)
    ##     res[!poisidx] <- dnbinom(y[!poisidx], mu=mu[!poisidx],
    ##                              size=size[!poisidx], log=TRUE)
    ##     res
    ## }
    index.overdisp <- if (is.factor(control$family)) {
      control$family
    } else if (control$family == "NegBinM") {
      factor(colnames(stsObj), levels = colnames(stsObj))
      ## do not sort levels (for consistency with unitSpecific effects)
    } else { # "NegBin1"
      factor(character(nUnits))
    }
    names(index.overdisp) <- colnames(stsObj)
    dim.overdisp <- nlevels(index.overdisp)
    names.overdisp <- if (dim.overdisp == 1L) {
      "-log(overdisp)"
    } else {
      paste0("-log(", paste("overdisp", levels(index.overdisp), sep = "."), ")")
    }
  }

  environment(ddistr) <- getNamespace("stats")  # function is self-contained

  # parameter start values from fe() and ri() calls via checkFormula()
  initial <- list(
    fixed = c(unlist(all.term["initial.fe",]),
              initial.d,
              rep.int(2, dim.overdisp)),
    random = as.numeric(unlist(all.term["initial.re",])), # NULL -> numeric(0)
    sd.corr = c(unlist(all.term["initial.var",]),
                rep.int(0, dim.corr))
  )
  # set names of parameter vectors
  names(initial$fixed) <- c(names.fe, names.d, names.overdisp)
  names(initial$random) <- names.re
  names(initial$sd.corr) <- c(names.var, head(paste("corr",1:4,sep="."), dim.corr))

  # modify initial values according to the supplied 'start' values
  initial[] <- mapply(
    FUN = function (initial, start, name) {
      if (is.null(start))
        return(initial)
      if (is.null(names(initial)) || is.null(names(start))) {
        if (length(start) == length(initial)) {
          initial[] <- start
        } else {
          stop("initial values in 'control$start$", name,
               "' must be of length ", length(initial))
        }
      } else {
        ## we match by name and silently ignore additional start values
        start <- start[names(start) %in% names(initial)]
        initial[names(start)] <- start
      }
      return(initial)
    },
    initial, control$start[names(initial)], names(initial),
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )

  # Done
  result <- list(response = Y,
                 terms = all.term,
                 nTime = nTime,
                 nUnits = nUnits,
                 nFE = dim.fe,
                 nd = dim.d,
                 nOverdisp = dim.overdisp,
                 nRE = dim.re,
                 rankRE = dim.re.group,
                 nVar = dim.var,
                 nCorr = dim.corr,
                 nSigma = dim.var+dim.corr,
                 nGroups = ncol(all.term),
                 namesFE = names.fe,
                 indexFE = index.fe,
                 indexRE = index.re,
                 initialTheta = c(initial$fixed, initial$random),
                 initialSigma = initial$sd.corr,
                 offset = offsets,
                 family = ddistr,
                 indexPsi = index.overdisp,
                 subset = control$subset,
                 isNA = isNA,
                 ziTrue = ziTrue,
                 zi.lag.unitSpecific = control$zi$lag.unitSpecific
  )
  return(result)
}

# used to incorporate covariates and unit-specific effects
fe <- function (x, unitSpecific = FALSE, which = NULL, initial = NULL)
{
  stsObj <- get("stsObj", envir = parent.frame(1), inherits = TRUE)
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)
  if (!is.numeric(x)) {
    stop("Covariate '", deparse(substitute(x)), "' is not numeric\n")
  }
  lengthX <- length(x)
  if (lengthX == 1) {
    terms <- matrix(x, nTime, nUnits, byrow = FALSE)
    mult <- "*"
  }
  else if (lengthX == nTime) {
    terms <- matrix(x, nTime, nUnits, byrow = FALSE)
    mult <- "*"
  }
  else if (lengthX == nTime * nUnits) {
    if (!is.matrix(x)) {
      stop("Covariate '", deparse(substitute(x)), "' is not a matrix\n")
    }
    if ((ncol(x) != nUnits) | (nrow(x) != nTime)) {
      stop("Dimension of covariate '", deparse(substitute(x)),
           "' is not suitably specified\n")
    }
    terms <- x
    mult <- "*"
  }
  else {
    stop("Covariate '", deparse(substitute(x)), "' is not suitably specified\n")
  }
  intercept <- all(terms == 1)
  unitSpecific <- unitSpecific || !is.null(which)
  if (unitSpecific) {
    if (is.null(which)) {
      which <- rep.int(TRUE, nUnits)
    }
    else {
      stopifnot(is.vector(which, mode = "logical"), length(which) ==
                  nUnits)
    }
    terms[, !which] <- 0
  }
  dim.fe <- if (unitSpecific)
    sum(which)
  else 1
  if (is.null(initial)) {
    initial <- rep.int(0, dim.fe)
  }
  else if (length(initial) != dim.fe) {
    stop("initial values for '", deparse(substitute(x)),
         "' must be of length ", dim.fe)
  }
  summ <- if (unitSpecific)
    "colSums"
  else "sum"
  name <- deparse(substitute(x))
  if (unitSpecific)
    name <- paste(name, colnames(stsObj)[which], sep = ".")
  result <- list(terms = terms, name = name, Z.intercept = NULL,
                 which = which, dim.fe = dim.fe, initial.fe = initial,
                 dim.re = 0, dim.var = 0, initial.var = NULL, initial.re = NULL,
                 intercept = intercept, unitSpecific = unitSpecific, random = FALSE,
                 corr = FALSE, summ = summ, mult = mult)
  return(result)
}
# random intercepts
ri <- function (type = c("iid", "car"), corr = c("none", "all"), initial.fe = 0,
                initial.var = -0.5, initial.re = NULL)
{
  stsObj <- get("stsObj", envir = parent.frame(1), inherits = TRUE)
  if (ncol(stsObj) == 1)
    stop("random intercepts require a multivariate 'stsObj'")
  type <- match.arg(type)
  corr <- match.arg(corr)
  corr <- switch(corr, none = FALSE, all = TRUE)
  if (type == "iid") {
    Z <- 1
    dim.re <- ncol(stsObj)
    mult <- "*"
  }
  else if (type == "car") {
    K <- neighbourhood(stsObj)
    checkNeighbourhood(K)
    K <- K == 1
    ne <- colSums(K)
    K <- -1 * K
    diag(K) <- ne
    dimK <- nrow(K)
    if (qr(K)$rank != dimK - 1)
      stop("neighbourhood matrix contains islands")
    svdK <- svd(K)
    L <- svdK$u[, -dimK] %*% diag(sqrt(svdK$d[-dimK]))
    Z <- L %*% solve(t(L) %*% L)
    dim.re <- dimK - 1L
    mult <- "%*%"
  }
  stopifnot(length(initial.fe) == 1, length(initial.var) ==
              1)
  if (is.null(initial.re)) {
    initial.re <- rnorm(dim.re, 0, sd = sqrt(0.001))
  }
  else if (length(initial.re) != dim.re) {
    stop("'initial.re' must be of length ", dim.re)
  }
  result <- list(terms = 1, name = paste("ri(", type, ")",
                                         sep = ""), Z.intercept = Z, which = NULL, dim.fe = 1,
                 initial.fe = initial.fe, dim.re = dim.re, dim.var = 1,
                 initial.var = initial.var, initial.re = initial.re, intercept = TRUE,
                 unitSpecific = FALSE, random = TRUE, corr = corr, summ = "colSums",
                 mult = mult)
  return(result)
}

checkFormula <- function (f, component, data, stsObj)
{
  term <- terms.formula(f, specials = c("fe", "ri"))
  intercept.all <- attr(term, "intercept") == 1
  vars <- as.list(attr(term, "variables"))[-1]
  nVars <- length(vars)
  res <- if (intercept.all) {
    c(fe(1), list(offsetComp = component))
  }
  else {
    if (nVars == 0)
      stop("formula ", deparse(substitute(f)), " contains no variables")
    NULL
  }
  fe.raw <- setdiff(seq_len(nVars), unlist(attr(term, "specials")))
  for (i in fe.raw) res <- cbind(res, c(eval(substitute(fe(x),
                                                        list(x = vars[[i]])), envir = data), list(offsetComp = component)))
  for (i in attr(term, "specials")$fe) res <- cbind(res, c(eval(vars[[i]],
                                                                envir = data), list(offsetComp = component)))
  res <- cbind(res, deparse.level = 0)
  RI <- attr(term, "specials")$ri
  if (sum(unlist(res["intercept", ])) + length(RI) > 1)
    stop("There can only be one intercept in the formula ",
         deparse(substitute(f)))
  for (i in RI) res <- cbind(res, c(eval(vars[[i]], envir = data),
                                    list(offsetComp = component)))
  return(res)
}
gammaZero <- function(theta, model, subset = model$subset, d = 0, .ar = TRUE)
{
  ## unpack theta
  pars <- surveillance:::splitParams(theta, model)
  fixed <- pars$fixed
  random <- pars$random

  ## unpack model
  term <- model$terms
  nGroups <- model$nGroups

  comp <- unlist(term["offsetComp",])
  idxFE <- model$indexFE
  idxRE <- model$indexRE

  toMatrix <- function (x, r=model$nTime, c=model$nUnits)
    matrix(x, r, c, byrow=TRUE)

  unitNames <- dimnames(model$response)[[2L]]

  setColnames <- if (is.null(unitNames)) identity else
    function(x) "dimnames<-"(x, list(NULL, unitNames))

  pred <- nullMatrix <- toMatrix(0)

  for(i in seq_len(nGroups)[comp==4]){
    #browser()
    fe <- fixed[idxFE==i]
    if(!.ar & grepl("AR", term["name",i])) next
    if(term["unitSpecific",i][[1]]){
      fe <- nullMatrix
      which <- term["which",i][[1]]
      fe[,which] <- toMatrix(fixed[idxFE==i],c=sum(which))
    }
    if(term["random",i][[1]]){
      re <- random[idxRE==i]
      "%m%" <- get(term["mult",i][[1]])
      Z.re <- toMatrix(term["Z.intercept",i][[1]] %m% re)
    } else {
      Z.re <- 0
    }
    X <- term["terms",i][[1]]
    pred <- pred + X*fe + Z.re
  }
  x <- pred[subset,,drop=FALSE]
  pred <- setColnames(plogis(x))
  res <- if(!.ar) x else
    if(d == 1) gprime(x) else
      if(d == 2) gprime2(x) else
        pred
  return(res)

}
gprime <- function(x) exp(x)/(exp(x) + 1)^2
gprime2 <- function(x) exp(x) * (1 - exp(x))/(exp(x) + 1)^3
# for hhh4 and zero-inflated model
penLogLik <- function(theta, sd.corr, model, attributes=FALSE, individual = FALSE)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE

  ## unpack random effects
  if (dimRE > 0) {
    pars <- surveillance:::splitParams(theta, model)
    randomEffects <- pars$random
    sd   <- head(sd.corr, model$nVar)
    corr <- tail(sd.corr, model$nCorr)
    dimBlock <- model$rankRE[model$rankRE>0]
    Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }

  ############################################################

  ## evaluate dispersion
  psi <- surveillance::sizeHHH(theta, model,
                                subset = if (dimPsi > 1L) subset) # else scalar or NULL

  #psi might be numerically equal to 0 or Inf in which cases dnbinom (in meanHHH)
  #would return NaN (with a warning). The case size=Inf rarely happens and
  #corresponds to a Poisson distribution. Currently this case is not handled
  #in order to have the usual non-degenerate case operate faster.
  #For size=0, log(dnbinom) equals -Inf for positive x or if (x=0 and mu=0), and
  #zero if x=0 and mu>0 and mu<Inf. Thus, if there is at least one case in Y
  #(x>0, which is always true), we have that sum(ll.units) = -Inf, hence:
  if (any(psi == 0)) return(-Inf)

  ## evaluate mean
  mu <- surveillance::meanHHH(theta, model, total.only=TRUE)
  # if, numerically, mu=Inf, log(dnbinom) or log(dpois) both equal -Inf, hence:
  #if (any(is.infinite(mu))) return(-Inf)
  # however, since mu=Inf does not produce warnings below and this is a rare
  # case, it is faster to not include this conditional expression

  # evaluate gamma
  gamma <- if(model$ziTrue) gammaZero(theta, model) else 0

  ## penalization term for random effects
  lpen <- if (dimRE==0) 0 else { # there are random effects
    ##-.5*(t(randomEffects)%*%Sigma.inv%*%randomEffects)
    ## the following implementation takes ~85% less computing time !
    -0.5 * c(crossprod(randomEffects, Sigma.inv) %*% randomEffects)
  }

  ## log-likelihood
  ll.individ <- log(gamma * (Y == 0) + (1 - gamma) * exp(model$family(Y,mu,psi)))
  # in model$family, log = TRUE
  if(individual) return(ll.individ)
  ll.units <- .colSums(ll.individ,
                       length(subset), model$nUnits, na.rm=TRUE)

  ## penalized log-likelihood
  ll <- sum(ll.units) + lpen

  ## Done
  if (attributes) {
    attr(ll, "loglik") <- ll.units
    attr(ll, "logpen") <- lpen
  }
  return(ll)
}

###------------------------------------------------------------------------
# only for zero-inflated model
penScore <- function(theta, sd.corr, model, individual = FALSE)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE
  term <- model$terms
  nGroups <- model$nGroups
  dimd <- model$nd

  ## unpack parameters
  pars <- surveillance:::splitParams(theta, model)
  if (dimRE > 0) {
    randomEffects <- pars$random
    sd   <- head(sd.corr, model$nVar)
    corr <- tail(sd.corr, model$nCorr)
    dimBlock <- model$rankRE[model$rankRE>0]
    Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }

  ## evaluate dispersion
  psi <- surveillance::sizeHHH(theta, model,
                                subset = if (dimPsi > 1L) subset) # else scalar or NULL

  ## evaluate mean
  mu <- surveillance::meanHHH(theta, model)
  meanTotal <- mu$mean

  ## evaluate gamma
  gamma <- gammaZero(theta, model)

  # evaluate logLik
  llZi <- penLogLik(theta, sd.corr, model, attributes=FALSE, individual = TRUE)

  # hhh part density
  llHHH <- model$family(Y,meanTotal,psi)
  # in model$family, log = TRUE
  ############################################################

  ## helper function for derivatives, hhh part
  derivHHH.factor <- if(dimPsi > 0L){ # NegBin
    psiPlusMu <- psi + meanTotal    # also used below for calculation of grPsi
    psiYpsiMu <- (psi+Y) / psiPlusMu
    Y/meanTotal - psiYpsiMu
  } else { # Poisson
    Y/meanTotal - 1
  }
  derivHHH <- function (dmu) derivHHH.factor * dmu

  ## helper function for zero part
  derivZero.factor <- exp(-llZi + llHHH) * (1 - gamma)
  derivZero <- function(dllHHH) derivZero.factor * dllHHH

  derivGamma.factor <- exp(-llZi) * ((Y == 0) - exp(llHHH))
  derivGamma <-  function(dgammaZero) derivGamma.factor * dgammaZero

  score_list <- list()
  score_list$terms <- vector(mode = "list", length = nGroups)
  # individual scores for gamma is not needed in calculation for fisher
  gamma1 <- gammaZero(theta, model, subset = model$subset, d = 1)
  mean.comp <- mu[c("epi.own","epi.neighbours","endemic")]
  ## go through groups of parameters and compute the gradient of each component

  grad.fe <- numeric(0L)
  grad.re <- numeric(0L)

  for(i in seq_len(nGroups)){
    comp <- term["offsetComp",i][[1]]
    Xit<- term["terms",i][[1]] # eiter 1 or a matrix with values
    if(is.matrix(Xit)){
      Xit <- Xit[subset,,drop=FALSE]
    }
    if(comp != 4){
      dTheta_HHH <- derivHHH(mean.comp[[comp]]*Xit) # time * region
      dTheta_HHH[isNA] <- 0   # dTheta must not contain NA's (set NA's to 0)
      dTheta <- derivZero(dTheta_HHH)
      score_list$terms[[i]] <- dTheta_HHH
      if(term["unitSpecific",i][[1]]){
        which <- term["which",i][[1]]
        dimi <- sum(which)
        if(dimi < model$nUnits)
          dTheta <- dTheta[,which,drop=FALSE]
        dTheta <- .colSums(dTheta, length(subset), dimi)
        grad.fe <- c(grad.fe,dTheta)

      }  else if(term["random",i][[1]]){
        Z <- term["Z.intercept",i][[1]]
        "%m%" <- get(term["mult",i][[1]])
        dThetamZ <- dTheta %m% Z
        dRTheta <- .colSums(dThetamZ, length(subset), term["dim.re",i][[1]])
        # region specific
        grad.re <- c(grad.re, dRTheta)
        grad.fe <- c(grad.fe, sum(dTheta))
        # intercept
      } else{
        grad.fe <- c(grad.fe, sum(dTheta))
      }

    } else{
      dgammaZero <- gamma1 * Xit
      dgammaZero[isNA] <- 0   # dTheta must not contain NA's (set NA's to 0)
      dsgamma <- derivGamma(dgammaZero)

      if(term["unitSpecific",i][[1]]){
        which <- term["which",i][[1]]
        dimi <- sum(which)
        if(dimi < model$nUnits)
          dsgamma <- dsgamma[,which,drop=FALSE]
        dsgamma <- .colSums(dsgamma, length(subset), dimi)
        grad.fe <- c(grad.fe,dsgamma)
      }  else if(term["random",i][[1]]){
        Z <- term["Z.intercept",i][[1]]
        "%m%" <- get(term["mult",i][[1]])
        dsgammamZ <- dsgamma %m% Z
        dRsgamma <- .colSums(dsgammamZ, length(subset), term["dim.re",i][[1]])
        grad.re <- c(grad.re, dRsgamma)
        grad.fe <- c(grad.fe, sum(dsgamma))
      } else{
        grad.fe <- c(grad.fe, sum(dsgamma))
      }
    }
  } # for loop
  gradients <- list(fe=grad.fe, re=grad.re)



  # mu["epi.own"] = lambda_rt * y_r,t-1
  # mu["epi.neighbours"] = phi_rt * sum(omega_qr * y_q,r-1)
  # mu["endemic"] = nu_rt

  ## gradient for parameter vector of the neighbourhood weights
  # grd <- if (dimd > 0L) {
  #   dneOffset <- model$offset[[2L]](pars$d, type="gradient")
  #   ##<- this is always a list (of length dimd) of matrices
  #   onescore.d <- function (dneoff) {
  #     dmudd <- mu$ne.exppred * dneoff[subset,,drop=FALSE]
  #     grd.terms <- derivHHH(dmudd)
  #     grd.termsZi <- derivZero(grd.terms)
  #     if(individual) grd.termsZi else sum(grd.termsZi, na.rm=TRUE)
  #   }
  #   if(individual){
  #     score_list$d <- clapply(dneOffset, onescore.d)
  #     numeric(0L)
  #   }else unlist(clapply(dneOffset, onescore.d), recursive=FALSE, use.names=FALSE)
  # } else numeric(0L)

  grd <- if (dimd > 0L & !individual) {
    dneOffset <- model$offset[[2L]](pars$d, type="gradient")
    ##<- this is always a list (of length dimd) of matrices
    onescore.d <- function (dneoff) {
      dmudd <- mu$ne.exppred * dneoff[subset,,drop=FALSE]
      grd.terms <- derivHHH(dmudd)
      grd.termsZi <- derivZero(grd.terms)
      sum(grd.termsZi, na.rm=TRUE)
    }
    unlist(clapply(dneOffset, onescore.d), recursive=FALSE, use.names=FALSE)
  }else numeric(0L)
  if(individual & dimd > 0L){
    dneOffset <- model$offset[[2L]](pars$d, type="gradient")
    onescore.d.HHH <- function (dneoff) {
      dmudd <- mu$ne.exppred * dneoff[subset,,drop=FALSE]
      grd.terms <- derivHHH(dmudd)
      grd.terms[is.na(grd.terms)] <- 0
      grd.terms
    }
    score_list$d <- clapply(dneOffset, onescore.d.HHH)
  }

  ## gradient for overdispersion parameter psi
  grPsi <- if(dimPsi > 0L){
    dPsiMat <- psi * (digamma(Y+psi) - digamma(psi) + log(psi) + 1
                      - log(psiPlusMu) - psiYpsiMu)
    dPsiMatZi <- derivZero(dPsiMat)
    score_list$Psi <- dPsiMat
    surveillance:::.colSumsGrouped(dPsiMatZi, model$indexPsi)
  } else numeric(0L)

  ## add penalty to random effects gradient
  s.pen <- if(dimRE > 0) c(Sigma.inv %*% randomEffects) else numeric(0L)
  if(length(gradients$re) != length(s.pen))
    stop("oops... lengths of s(b) and Sigma.inv %*% b do not match")
  grRandom <- c(gradients$re - s.pen)

  res <- c(gradients$fe, grd, grPsi, grRandom)

  ## Done
  if(individual){
    return(score_list) # individual scores are not penalized
  } else return(res)

}
##------------------------------------------------------------------------

penFisher <- function(theta, sd.corr, model, attributes=FALSE)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE
  term <- model$terms
  nGroups <- model$nGroups
  dimd  <- model$nd
  dimFE <- model$nFE
  idxFE <- model$indexFE
  idxRE <- model$indexRE
  indexPsi <- model$indexPsi

  ## unpack parameters
  pars <- surveillance:::splitParams(theta, model)
  if (dimRE > 0) {
    randomEffects <- pars$random
    sd   <- head(sd.corr, model$nVar)
    corr <- tail(sd.corr, model$nCorr)
    dimBlock <- model$rankRE[model$rankRE>0]
    Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }

  ## evaluate dispersion
  psi <- surveillance::sizeHHH(theta, model,
                                subset = if (dimPsi > 1L) subset) # else scalar or NULL

  ## evaluate mean
  mu <- surveillance::meanHHH(theta, model)
  meanTotal <- mu$mean

  # evaluate individual score
  score_list <- penScore(theta, sd.corr, model, individual = TRUE)

  ## evaluate gamma
  gamma <- gammaZero(theta, model)
  gamma1 <- gammaZero(theta, model, subset = model$subset, d = 1)
  gamma2 <- gammaZero(theta, model, subset = model$subset, d = 2)
  # evaluate logLik
  llZi <- penLogLik(theta, sd.corr, model, attributes=FALSE, individual = TRUE)

  # hhh part density
  llHHH <- model$family(Y,meanTotal,psi)
  # in model$family, log = TRUE

  ############################################################

  ## helper functions for derivatives:
  if (dimPsi > 0L) { # negbin
    psiPlusY <- psi + Y
    psiPlusMu <- psi + meanTotal
    psiPlusMu2 <- psiPlusMu^2
    psiYpsiMu <- psiPlusY / psiPlusMu
    psiYpsiMu2 <- psiPlusY / psiPlusMu2
    deriv2HHH.fac1 <- psiYpsiMu2 - Y / (meanTotal^2)
    deriv2HHH.fac2 <- Y / meanTotal - psiYpsiMu
    ## psi-related derivatives
    dThetadPsi.fac <- psi * (psiYpsiMu2 - 1/psiPlusMu)
    dThetadPsi <- function(dTheta){
      dThetadPsi.fac * dTheta
    }
    dPsiMat <- psi * (digamma(psiPlusY) - digamma(psi) + log(psi) + 1
                      - log(psiPlusMu) - psiYpsiMu)  # as in penScore()
    dPsidPsiMat <- psi^2 * (
      trigamma(psiPlusY) - trigamma(psi) + 1/psi - 1/psiPlusMu -
        (meanTotal-Y)/psiPlusMu2) + dPsiMat
  } else { # poisson
    deriv2HHH.fac1 <- -Y / (meanTotal^2)
    deriv2HHH.fac2 <- Y / meanTotal - 1
  }
  deriv2HHH <- function(dTheta_l, dTheta_k, dTheta_lk){
    dTheta_l * dTheta_k * deriv2HHH.fac1 + dTheta_lk * deriv2HHH.fac2
  }
  deriv2HHHZi <- function(score_i, score_j, deriv2_ij){
    # (1 - gamma) * exp(llHHH) *
    #   (- exp(-2 * llZi + llHHH) * (1 - gamma) * score_j * score_i +
    #      exp(-llZi) * (score_j * score_i + deriv2_ij))
    (1 - gamma) * exp(llHHH -llZi) *
      (- exp(- llZi + llHHH) * (1 - gamma) * score_j * score_i +
         (score_j * score_i + deriv2_ij))
  }

  dThetadGamma <- function(score_Theta, dGamma){
    # exp(llHHH) * score_Theta *
    #   (- exp(-2 * llZi) * ((Y == 0) - exp(llHHH)) * dGamma * (1 - gamma)
    #    - exp(-llZi) * dGamma
    #    )

    exp(llHHH) * score_Theta *
      (exp(-llZi) *( -exp(-llZi) *((Y == 0) - exp(llHHH)) *
                       dGamma * (1 - gamma) - dGamma)
      )
  }

  d2Gamma <- function(dGamma1, dGamma2, dGamma_2){
    # ((Y == 0) - exp(llHHH)) *
    #   ((- exp(-2 * llZi) * ((Y == 0) - exp(llHHH)) * dGamma2) * dGamma1 +
    #      exp(-llZi) * dGamma_2)
    ((Y == 0) - exp(llHHH)) *
      (exp(-llZi) *(-exp(-llZi) * ((Y == 0) - exp(llHHH)) * dGamma2 * dGamma1
                    +   dGamma_2 ))
  }
  ## go through groups of parameters and compute the hessian of each component
  computeFisher <- function(mean.comp){
    # initialize hessian
    hessian.FE.FE <- matrix(0,dimFE,dimFE)
    hessian.FE.RE <- matrix(0,dimFE,dimRE)
    hessian.RE.RE <- matrix(0,dimRE,dimRE)

    hessian.FE.Psi <- matrix(0,dimFE,dimPsi)
    hessian.Psi.RE <- matrix(0,dimPsi,dimPsi+dimRE) # CAVE: contains PsiPsi and PsiRE

    hessian.FE.d <- matrix(0,dimFE,dimd)
    hessian.d.d <- matrix(0,dimd,dimd)
    hessian.d.Psi <- matrix(0,dimd,dimPsi)
    hessian.d.RE <- matrix(0,dimd,dimRE)

    ## derivatives wrt neighbourhood weight parameters d
    if (dimd > 0L) {
      phi.doff <- function (dneoff) {
        mu$ne.exppred * dneoff[subset,,drop=FALSE]
      }
      ## for type %in% c("gradient", "hessian"), model$offset[[2L]] always
      ## returns a list of matrices. It has length(pars$d) elements for the
      ## gradient and length(pars$d)*(length(pars$d)+1)/2 for the hessian.
      dneOffset <- model$offset[[2L]](pars$d, type="gradient")
      dmudd <- lapply(dneOffset, phi.doff)
      d2neOffset <- model$offset[[2L]](pars$d, type="hessian")
      d2mudddd <- lapply(d2neOffset, phi.doff)
      ## d l(theta,x) /dd dd (fill only upper triangle, BY ROW)
      ij <- 0L
      for (i in seq_len(dimd)) {
        for (j in i:dimd) {
          ij <- ij + 1L  #= dimd*(i-1) + j - (i-1)*i/2  # for j >= i
          ## d2mudddd contains upper triangle by row (=lowertri by column)
          d2ij <- deriv2HHH(dmudd[[i]], dmudd[[j]], d2mudddd[[ij]])

          score_di <- score_list$d[[i]]
          score_dj <- score_list$d[[j]]
          d2ij <- deriv2HHHZi(score_di, score_dj, d2ij)

          hessian.d.d[i,j] <- sum(d2ij, na.rm=TRUE)
        }
      }
    }

    if (dimPsi > 0L) {
      ## d l(theta,x) /dpsi dpsi
      score_Psi <- score_list$Psi
      dPsidPsi <- deriv2HHHZi(score_Psi, score_Psi, dPsidPsiMat)

      dPsidPsi <- surveillance:::.colSumsGrouped(dPsidPsi, indexPsi)
      hessian.Psi.RE[,seq_len(dimPsi)] <- if (dimPsi == 1L) {
        dPsidPsi
      } else {
        diag(dPsidPsi)
      }
      ## d l(theta) / dd dpsi
      for (i in seq_len(dimd)) {      # will not be run if dimd==0
        ## dPsi.i <- colSums(dThetadPsi(dmudd[[i]]),na.rm=TRUE)
        ## hessian.d.Psi[i,] <- if(dimPsi==1L) sum(dPsi.i) else dPsi.i[order(indexPsi)]
        dddPsi <- dThetadPsi(dmudd[[i]])
        score_di <- score_list$d[[i]]
        score_Psi <- score_list$Psi
        dddPsi <- deriv2HHHZi(score_di, score_Psi, dddPsi)

        hessian.d.Psi[i,] <- surveillance:::.colSumsGrouped(dddPsi, indexPsi)
      }
    }

    ##
    i.fixed <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]]
        "%mj%" <- get(term["mult",j][[1]])
        hessian.FE.RE[idxFE==i,idxRE==j] <<- colSums(didj %mj% Z.j)
        ##<- didj must not contain NA's (all NA's set to 0)
        dIJ <- sum(didj,na.rm=TRUE)     # fixed on 24/09/2012
      } else if(unitSpecific.j){
        dIJ <- colSums(didj,na.rm=TRUE)[ which.j ]
      } else {
        dIJ <- sum(didj,na.rm=TRUE)
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ
    }
    ##
    i.unit <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]]
        "%mj%" <- get(term["mult",j][[1]])
        dIJ <- colSums(didj %mj% Z.j)   # didj must not contain NA's (all NA's set to 0)
        hessian.FE.RE[idxFE==i,idxRE==j] <<- diag(dIJ)[ which.i, ] # FIXME: does not work if type="car"
        dIJ <- dIJ[ which.i ]           # added which.i subsetting in r432
      } else if(unitSpecific.j){
        dIJ <- diag(colSums(didj))[ which.i, which.j ]
      } else {
        dIJ <- colSums(didj)[ which.i ]
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ
    }
    ##
    i.random <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]]
        "%mj%" <- get(term["mult",j][[1]])
        hessian.FE.RE[idxFE==i,idxRE==j] <<- colSums(didj %mj% Z.j)
        if (j != i)  # otherwise redundant (duplicate)
          hessian.FE.RE[idxFE==j,idxRE==i] <<- colSums(didj %m% Z.i)

        if(length(Z.j)==1 & length(Z.i)==1){ # both iid
          Z <- Z.i*Z.j
          hessian.RE.RE[which(idxRE==i),idxRE==j] <<- diag(colSums( didj %m% Z))
        } else if(length(Z.j)==1 & length(Z.i)>1){         #*
          Z.j <- diag(nrow=model$nUnits)
          for(k in seq_len(ncol(Z.j))){
            Z <- Z.i*Z.j[,k]
            hessian.RE.RE[idxRE==i,which(idxRE==j)[k]] <<- colSums( didj %m% Z)
          }
        } else if(length(Z.j)>1 & length(Z.i)==1){         #*
          Z.i <- diag(nrow=model$nUnits)
          for(k in seq_len(ncol(Z.i))){
            Z <- Z.i[,k]*Z.j
            hessian.RE.RE[which(idxRE==i)[k],idxRE==j] <<- colSums( didj %mj% Z)
          }
        } else { # both CAR
          for(k in seq_len(ncol(Z.j))){
            Z <- Z.i*Z.j[,k]
            hessian.RE.RE[which(idxRE==i)[k],idxRE==j] <<- colSums( didj %m% Z)
          }
        }
        dIJ <- sum(didj)
      } else if(unitSpecific.j){
        dIJ <- colSums(didj %m% Z.i)
        hessian.FE.RE[idxFE==j,idxRE==i] <<- diag(dIJ)[ which.j, ]
        dIJ <- dIJ[ which.j ]
      } else {
        hessian.FE.RE[idxFE==j,idxRE==i] <<- colSums(didj %m% Z.i)
        dIJ <- sum(didj)
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ
    }
    ##----------------------------------------------

    for(i in seq_len(nGroups)){ #go through rows of hessian
      comp.i <- term["offsetComp",i][[1]]

      if(comp.i != 4){
        # parameter group belongs to which components
        # get covariate value
        Xit <- term["terms",i][[1]] # eiter 1 or a matrix with values
        if(is.matrix(Xit)){
          Xit <- Xit[subset,,drop=FALSE]
        }
        m.Xit <- mean.comp[[comp.i]] * Xit

        random.i <- term["random",i][[1]]
        unitSpecific.i <- term["unitSpecific",i][[1]]

        ## fill psi-related entries and select fillHess function
        if (random.i) {
          Z.i <- term["Z.intercept",i][[1]]   # Z.i and %m% (of i) determined here
          "%m%" <- get(term["mult",i][[1]])   # will also be used in j's for loop
          fillHess <- i.random
          if (dimPsi > 0L) {
            dThetadPsiMat <- dThetadPsi(m.Xit)

            score_Thetai <- score_list$terms[[i]]
            score_Psi <- score_list$Psi
            dThetadPsi <- deriv2HHHZi(score_Thetai, score_Psi, dThetadPsiMat)

            hessian.FE.Psi[idxFE==i,] <- surveillance:::.colSumsGrouped(dThetadPsi, indexPsi)
            dThetadPsi.i <- .colSums(dThetadPsi %m% Z.i, length(subset), term["dim.re",i][[1]], na.rm=TRUE)
            if (dimPsi==1L) {
              hessian.Psi.RE[,dimPsi + which(idxRE==i)] <- dThetadPsi.i
            } else {
              hessian.Psi.RE[cbind(indexPsi,dimPsi + which(idxRE==i))] <- dThetadPsi.i
              ## FIXME: does not work with type="car"
            }
          }
        } else if (unitSpecific.i) {
          which.i <- term["which",i][[1]]
          fillHess <- i.unit
          if (dimPsi > 0L) {
            dThetadPsiMat <- dThetadPsi(m.Xit)

            score_Thetai <- score_list$terms[[i]]
            score_Psi <- score_list$Psi
            dThetadPsi <- deriv2HHHZi(score_Thetai, score_Psi, dThetadPsiMat)

            dThetadPsi.i <- .colSums(dThetadPsi, length(subset), model$nUnits, na.rm=TRUE)
            if (dimPsi==1L) {
              hessian.FE.Psi[idxFE==i,] <- dThetadPsi.i[which.i]
            } else {
              hessian.FE.Psi[cbind(which(idxFE==i),indexPsi[which.i])] <-
                dThetadPsi.i[which.i]
            }
          }
        } else {
          fillHess <- i.fixed
          if (dimPsi > 0L) {
            ## dPsi <- colSums(dThetadPsi(m.Xit),na.rm=TRUE)
            ## hessian.FE.Psi[idxFE==i,] <- if (dimPsi==1L) sum(dPsi) else dPsi[order(indexPsi)]
            dThetadPsiMat <- dThetadPsi(m.Xit)

            score_Thetai <- score_list$terms[[i]]
            score_Psi <- score_list$Psi
            dThetadPsi <- deriv2HHHZi(score_Thetai, score_Psi, dThetadPsiMat)


            hessian.FE.Psi[idxFE==i,] <- surveillance:::.colSumsGrouped(dThetadPsi, indexPsi)
          }
        }

        ## fill pars$d-related entries
        for (j in seq_len(dimd)) {      # will not be run if dimd==0
          didd <- deriv2HHH(dTheta_l = m.Xit, dTheta_k = dmudd[[j]],
                            dTheta_lk = if (comp.i == 2) dmudd[[j]] * Xit else 0)
          didd[isNA] <- 0
          score_dj <- score_list$d[[j]]
          score_Thetai <- score_list$terms[[i]]
          didd <- deriv2HHHZi(score_Thetai, score_dj, didd)

          hessian.FE.d[idxFE==i,j] <- if (unitSpecific.i) {
            colSums(didd,na.rm=TRUE)[which.i]
          } else sum(didd)
          if (random.i) hessian.d.RE[j,idxRE==i] <- colSums(didd %m% Z.i)
        }

        ## fill other (non-psi, non-d) entries (only upper triangle, j >= i!)
        for(j in i:nGroups){
          comp.j <- term["offsetComp",j][[1]]

          Xjt <- term["terms",j][[1]] # eiter 1 or a matrix with values
          if(is.matrix(Xjt)){
            Xjt <- Xjt[subset,,drop=FALSE]
          }
          # if param i and j do not belong to the same component, d(i)d(j)=0
          m.Xit.Xjt <- if (comp.i != comp.j) 0 else m.Xit * Xjt
          if(comp.j !=4){
            didj <- deriv2HHH(dTheta_l = m.Xit, dTheta_k = mean.comp[[comp.j]]*Xjt,
                              dTheta_lk = m.Xit.Xjt)
            didj[isNA]<-0

            score_di <- score_list$terms[[i]]
            score_dj <- score_list$terms[[j]]
            didj <- deriv2HHHZi(score_di, score_dj, didj)


          }else{
            # didj <- deriv2HHH(dTheta_l = m.Xit, dTheta_k = mean.comp[[comp.j]]*Xjt,
            #                   dTheta_lk = m.Xit.Xjt)
            dGamma <- gamma1 * Xjt
            score_Theta <- score_list$terms[[i]]
            didj <- dThetadGamma(score_Theta, dGamma)
            didj[isNA]<-0

          }

          random.j <- term["random",j][[1]]
          unitSpecific.j <- term["unitSpecific",j][[1]]
          which.j <- term["which",j][[1]]
          #browser()
          fillHess()
        }

      }else{
        # parameter group belongs to which components
        # get covariate value
        Xit <- term["terms",i][[1]] # eiter 1 or a matrix with values
        if(is.matrix(Xit)){
          Xit <- Xit[subset,,drop=FALSE]
        }
        #m.Xit <- mean.comp[[comp.i]] * Xit

        random.i <- term["random",i][[1]]
        unitSpecific.i <- term["unitSpecific",i][[1]]

        ## fill psi-related entries and select fillHess function
        if (random.i) {
          Z.i <- term["Z.intercept",i][[1]]   # Z.i and %m% (of i) determined here
          "%m%" <- get(term["mult",i][[1]])   # will also be used in j's for loop
          fillHess <- i.random
          if (dimPsi > 0L) {
            score_Psi <- score_list$Psi
            dGammai <- gamma1 * Xit
            dThetadPsi <-  dThetadGamma(score_Psi, dGammai)

            hessian.FE.Psi[idxFE==i,] <- surveillance:::.colSumsGrouped(dThetadPsi, indexPsi)
            dThetadPsi.i <- .colSums(dThetadPsi %m% Z.i, length(subset), term["dim.re",i][[1]], na.rm=TRUE)
            if (dimPsi==1L) {
              hessian.Psi.RE[,dimPsi + which(idxRE==i)] <- dThetadPsi.i
            } else {
              hessian.Psi.RE[cbind(indexPsi,dimPsi + which(idxRE==i))] <- dThetadPsi.i
              ## FIXME: does not work with type="car"
            }
          }
        } else if (unitSpecific.i) {
          which.i <- term["which",i][[1]]
          fillHess <- i.unit
          if (dimPsi > 0L) {
            score_Psi <- score_list$Psi
            dGammai <- gamma1 * Xit
            dThetadPsi <-  dThetadGamma(score_Psi, dGammai)


            dThetadPsi.i <- .colSums(dThetadPsi, length(subset), model$nUnits, na.rm=TRUE)
            if (dimPsi==1L) {
              hessian.FE.Psi[idxFE==i,] <- dThetadPsi.i[which.i]
            } else {
              hessian.FE.Psi[cbind(which(idxFE==i),indexPsi[which.i])] <-
                dThetadPsi.i[which.i]
            }
          }
        } else {
          fillHess <- i.fixed
          if (dimPsi > 0L) {
            ## dPsi <- colSums(dThetadPsi(m.Xit),na.rm=TRUE)
            ## hessian.FE.Psi[idxFE==i,] <- if (dimPsi==1L) sum(dPsi) else dPsi[order(indexPsi)]
            score_Psi <- score_list$Psi
            dGammai <- gamma1 * Xit
            dThetadPsi <-  dThetadGamma(score_Psi, dGammai)

            hessian.FE.Psi[idxFE==i,] <- surveillance:::.colSumsGrouped(dThetadPsi, indexPsi)
          }
        }

        ## fill pars$d-related entries
        for (j in seq_len(dimd)) {      # will not be run if dimd==0

          dGammai <- gamma1 * Xit
          score_dj <- score_list$d[[j]]
          didd <- dThetadGamma(score_dj, dGammai)
          didd[isNA] <- 0

          hessian.FE.d[idxFE==i,j] <- if (unitSpecific.i) {
            colSums(didd,na.rm=TRUE)[which.i]
          } else sum(didd)
          if (random.i) hessian.d.RE[j,idxRE==i] <- colSums(didd %m% Z.i)
        }

        ## fill other (non-psi, non-d) entries (only upper triangle, j >= i!)
        for(j in i:nGroups){
          comp.j <- term["offsetComp",j][[1]]

          Xjt <- term["terms",j][[1]] # eiter 1 or a matrix with values
          if(is.matrix(Xjt)){
            Xjt <- Xjt[subset,,drop=FALSE]
          }
          # # if param i and j do not belong to the same component, d(i)d(j)=0
          # m.Xit.Xjt <- if (comp.i != comp.j) 0 else m.Xit * Xjt
          if(comp.j !=4){
            dGammai <- gamma1 * Xit
            score_dj <- score_list$terms[[j]]
            didd <- dThetadGamma(score_dj, dGammai)
            didd[isNA] <- 0

          }else{
            # didj <- deriv2HHH(dTheta_l = m.Xit, dTheta_k = mean.comp[[comp.j]]*Xjt,
            #                   dTheta_lk = m.Xit.Xjt)
            dGammai <- gamma1 * Xit
            dGammaj <- gamma1 * Xjt
            dGamma2 <- gamma2 * Xit * Xjt
            didj <- d2Gamma(dGammai, dGammaj, dGamma2)
            didj[isNA]<-0

          }

          random.j <- term["random",j][[1]]
          unitSpecific.j <- term["unitSpecific",j][[1]]
          which.j <- term["which",j][[1]]

          fillHess()
        }
      }
    }

    #########################################################
    ## fill lower triangle of hessians and combine them
    ########################################################
    hessian <- rbind(cbind(hessian.FE.FE,hessian.FE.d,hessian.FE.Psi,hessian.FE.RE),
                     cbind(matrix(0,dimd,dimFE),hessian.d.d,hessian.d.Psi,hessian.d.RE),
                     cbind(matrix(0,dimPsi,dimFE+dimd),hessian.Psi.RE),
                     cbind(matrix(0,dimRE,dimFE+dimd+dimPsi),hessian.RE.RE))

    hessian[lower.tri(hessian)] <- 0  # CAR blocks in hessian.RE.RE were fully filled
    diagHessian <- diag(hessian)
    fisher <- -(hessian + t(hessian))
    diag(fisher) <- -diagHessian

    return(fisher)
  }

  fisher <- computeFisher(mu[c("epi.own","epi.neighbours","endemic")])
  #browser()
  ## add penalty for random effects
  pen <- matrix(0, length(theta), length(theta))
  #browser()
  Fpen <- if(dimRE > 0){
    thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
    pen[thetaIdxRE,thetaIdxRE] <- Sigma.inv
    fisher + pen
  } else fisher
  #browser()
  ## Done
  if(attributes){
    attr(Fpen, "fisher") <- fisher
    attr(Fpen, "pen") <- pen
  }
  Fpen
}

#-----------------------------------
# sigma for dim = 4
getSigmai <- function(sd,    # vector of length dim with log-stdev's
                      correlation, # vector of length dim with correlation
                      # parameters, 0-length if uncorrelated
                      dim
){
  if(dim==0) return(NULL)

  Sigma.i <- if (length(correlation) == 0L) diag(exp(2*sd), dim) else {
    D <- diag(exp(sd), dim)
    L <- diag(nrow=dim)
    L[2,1:2] <- c(correlation[1],1)/surveillance:::sqrtOf1pr2(correlation[1])
    if (dim >=3) {
      L[3,1:3] <- c(correlation[2:3],1)/surveillance:::sqrtOf1pr2(correlation[2])
      L[3,2:3] <- L[3,2:3]/surveillance:::sqrtOf1pr2(correlation[3])
    }
    if(dim == 4){
      L[4,1:4] <- c(correlation[4:6],1)/surveillance:::sqrtOf1pr2(correlation[4])
      L[4,2:4] <- L[4,2:4]/surveillance:::sqrtOf1pr2(correlation[5])
      L[4,3:4] <- L[4,3:4]/surveillance:::sqrtOf1pr2(correlation[6])
    }
    D %*% tcrossprod(L) %*% D  # ~75% quicker than D %*% L %*% t(L) %*% D
  }
  return(Sigma.i)
}

# sigmainv for dim = 4
getSigmaiInv <- function(sd,    # vector of length dim with log-stdev's
                         correlation, # vector of length dim with correlation
                         # parameters, 0-length if uncorrelated
                         dim
){

  if(dim==0) return(NULL)

  Sigma.i.inv <- if (length(correlation) == 0L) diag(exp(-2*sd), dim) else {
    r <- correlation
    Dinv <- diag(exp(-sd), dim)
    L <- diag(nrow=dim)
    L[2,1:2] <- c(-r[1],surveillance:::sqrtOf1pr2(r[1]))
    if(dim>=3){
      L[3,1] <- r[1]*r[3]-r[2]*surveillance:::sqrtOf1pr2(r[3])
      L[3,2] <- -L[2,2]*r[3]
      L[3,3] <- surveillance:::sqrtOf1pr2(r[2])*surveillance:::sqrtOf1pr2(r[3])
    }
    if(dim == 4){
      L[4,1] <- -r[4]*surveillance:::sqrtOf1pr2(r[5])*
        surveillance:::sqrtOf1pr2(r[6]) +
        r[1]*r[5]*surveillance:::sqrtOf1pr2(r[6]) -
        r[6] * L[3,1]
      L[4,2] <- -r[5]*surveillance:::sqrtOf1pr2(r[1])*
        surveillance:::sqrtOf1pr2(r[6]) - r[6]*L[3,2]
      L[4,3] <- -r[6]* L[3,3]
      L[4,4] <- surveillance:::sqrtOf1pr2(r[4])*
        surveillance:::sqrtOf1pr2(r[5])*
        surveillance:::sqrtOf1pr2(r[6])
    }
    Dinv %*% crossprod(L) %*% Dinv  # ~75% quicker than Dinv %*% t(L) %*% L %*% Dinv
  }

  return(Sigma.i.inv)
}

getSigmaInv <- function(sd, correlation, dimSigma, dimBlocks, SigmaInvi=NULL){
  if(is.null(SigmaInvi)){
    SigmaInvi <- getSigmaiInv(sd,correlation,dimSigma)
  }
  if(length(unique(dimBlocks))==1){  # kronecker product formulation possible
    kronecker(SigmaInvi,diag(nrow=dimBlocks[1]))
    # the result is a symmetric matrix if SigmaInvi is symmetric
  } else { # kronecker product not possible -> correlation=0
    diag(rep.int(diag(SigmaInvi),dimBlocks))
  }
}

getSigma <- function(sd, correlation, dimSigma, dimBlocks, Sigmai=NULL){
  if(is.null(Sigmai)){
    Sigmai <- getSigmai(sd,correlation,dimSigma)
  }
  if(length(unique(dimBlocks))==1){  # kronecker product formulation possible
    kronecker(Sigmai,diag(nrow=dimBlocks[1]))
    # the result is a symmetric matrix if Sigmai is symmetric
  } else { # kronecker product not possible -> correlation=0
    diag(rep.int(diag(Sigmai),dimBlocks))
  }
}

## Approximate marginal likelihood for variance components
## Parameter and model unpacking at the beginning (up to the ###...-line) is
## identical in marScore() and marFisher()
marLogLik <- function(sd.corr, theta, model, fisher.unpen=NULL, verbose=FALSE){

  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma

  if(dimSigma == 0){
    return(-Inf)
  }

  if(any(is.na(sd.corr))){
    # in order to avoid nlminb from running into an infinite loop (cf. bug
    # report #15052), we have to emergency stop() in this case.
    # As of R 2.15.2, nlminb() throws an error if it receives NA from
    # any of the supplied functions.
    stop("NAs in variance parameters.", ADVICEONERROR)
  }

  sd   <- head(sd.corr,dimVar)
  corr <- tail(sd.corr,dimCorr)

  pars <- surveillance:::splitParams(theta,model)
  randomEffects <- pars$random
  dimRE <- model$nRE

  dimBlocks <- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)

  # if not given, calculate unpenalized part of fisher info
  if(is.null(fisher.unpen)){
    fisher.unpen <- attr(penFisher(theta, sd.corr, model,attributes=TRUE), "fisher")
  }

  # add penalty to fisher
  fisher <- fisher.unpen
  thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
  fisher[thetaIdxRE,thetaIdxRE] <- fisher[thetaIdxRE,thetaIdxRE] + Sigma.inv

  ############################################################

  # penalized part of likelihood
  # compute -0.5*log(|Sigma|) - 0.5*RE' %*% Sigma.inv %*% RE
  # where -0.5*log(|Sigma|) = -dim(RE_i)*[Sum(sd_i) -0.5*log(1+corr_i^2)]
  ##lpen <- -0.5*(t(randomEffects)%*%Sigma.inv%*%randomEffects)
  ## the following implementation takes ~85% less computing time !
  lpen <- -0.5 * c(crossprod(randomEffects, Sigma.inv) %*% randomEffects)
  loglik.pen <- sum(-dimBlocks*sd) + lpen
  if(dimCorr >0){
    loglik.pen <- loglik.pen + 0.5*dimBlocks[1]*sum(log(1+corr^2))
  }

  ## approximate marginal likelihood
  logdetfisher <- determinant(fisher,logarithm=TRUE)$modulus
  lmarg <- loglik.pen -0.5*c(logdetfisher)

  return(lmarg)
}
updateParams_nlminb <- surveillance:::updateParams_nlminb
updateParams_nlm <- surveillance:::updateParams_nlm
updateParams_optim <- function (start, ll, sc = NULL, fi = NULL, ..., control)
{
  ## Note: "fi" is not used in optim
  method <- control[["method"]]; control$method <- NULL
  lower <- control[["lower"]]; control$lower <- NULL
  upper <- control[["upper"]]; control$upper <- NULL
  res <- optim(start, ll, sc, ...,    # Note: control$fnscale is negative
               method=method, lower=lower, upper=upper, control=control)
  if (any(is.finite(c(lower, upper)))) checkParBounds(res$par, lower, upper)
  ## Done
  list(par=res$par, ll=res$value,
       rel.tol= surveillance:::getRelDiff(res$par, start),
       convergence=res$convergence, message=res$message)
}
##------------------------------------------------------------------------
## fitHHH is the main workhorse where the iterative optimization is performed
fitHHH4ZI <- function(theta, sd.corr, model,
                      cntrl.stop=list(tol=1e-5, niter=100),
                      cntrl.regression=list(method="nlminb"),
                      cntrl.variance=list(method="Nelder-Mead"),
                      verbose=0, shrinkage=FALSE)
{
  dimFE.d.O <- model$nFE + model$nd + model$nOverdisp
  dimRE <- model$nRE

  getUpdater <- function (cntrl, start, ...) {
    method <- cntrl$method; cntrl$method <- NULL
    if (length(start) == 1 && method == "Nelder-Mead") {
      method <- "Brent"
      message("Switched optimizer from \"Nelder-Mead\" to \"Brent\"",
              " (dim(", deparse(substitute(start)), ")=1)")
    }
    list(paste("updateParams",
               if (method %in% c("nlminb", "nlm", "nr"))
                 method else "optim", sep="_"),
         control = surveillance:::setOptimControl(method, cntrl, ...))
  }

  ## ## artificial lower bound on intercepts of epidemic components
  ## reg.lower <- rep.int(-Inf, length(theta))
  ## reg.lower[grep("^(ar|ne)\\.(1|ri)", model$namesFE)] <- -20

  ## set optimizer for regression parameters
  updateRegressionControl <- getUpdater(cntrl.regression, theta,
                                        ## lower=reg.lower,
                                        iter.max=if(dimRE==0) 100,
                                        verbose=verbose+(dimRE==0))
  updateRegression <- function (theta, sd.corr)
    do.call(updateRegressionControl[[1]],
            alist(theta, penLogLik, penScore, penFisher,
                  sd.corr=sd.corr, model=model,
                  control=updateRegressionControl[[2]]))

  ## set optimizer for variance parameters
  if(cntrl.variance$method != "Nelder-Mead")
    stop("only method Nelder-Mead is available for variance part")
  updateVarianceControl <- getUpdater(cntrl.variance, sd.corr,
                                      lower=-5, upper=5, verbose=verbose)
  updateVariance <- function (sd.corr, theta, fisher.unpen)
    do.call(updateVarianceControl[[1]],
            alist(sd.corr, marLogLik,
                  theta=theta, model=model,
                  fisher.unpen=fisher.unpen, verbose=verbose>1,
                  control=updateVarianceControl[[2]]))

  ## Let's go
  if (verbose>0) {
    cat(as.character(Sys.time()), ":",
        if (dimRE == 0) "Optimization of regression parameters" else
          "Iterative optimization of regression & variance parameters", "\n")
  }

  if (dimRE == 0) { # optimization of regression coefficients only
    parReg <- updateRegression(theta, sd.corr)
    theta <- parReg$par
    if ((convergence <- parReg$convergence) != 0 && !is.null(parReg$message))
      cat("! Non-convergence message from optimizer:", parReg$message, "\n")
  } else { # swing between updateRegression & updateVariance
    convergence <- 99
    i <- 0
    while(convergence != 0 && (i < cntrl.stop$niter)){
      i <- i+1
      if (verbose>0) cat("\n")

      ## update regression coefficients
      parReg <- updateRegression(theta, sd.corr)
      theta <- parReg$par
      fisher.unpen <- attr(penFisher(theta, sd.corr, model, attributes=TRUE),
                           "fisher")

      if(verbose>0)
        cat("Update of regression parameters: ",
            "max|x_0 - x_1| / max|x_0| =", parReg$rel.tol, "\n")

      if(parReg$convergence != 0) {
        if (!is.null(parReg$message))
          cat("! Non-convergence message from optimizer:",
              parReg$message, "\n")
        cat("Update of regression coefficients in iteration ",
            i, " unreliable\n")
      }

      if(parReg$convergence > 20 && shrinkage){
        cat("\n\n***************************************\nshrinkage",
            0.1*theta[abs(theta)>10],"\n")
        theta[abs(theta)>10] <- 0.1*theta[abs(theta)>10]
        diag(fisher.unpen) <- diag(fisher.unpen)+1e-2
      }

      ## update variance parameters
      parVar <- updateVariance(sd.corr, theta, fisher.unpen)

      if(verbose>0)
        cat("Update of variance parameters:  max|x_0 - x_1| / max|x_0| =",
            parVar$rel.tol, "\n")

      if(parVar$convergence!=0) {
        if (!is.null(parVar$message)) print(parVar$message)
        cat("Update of variance parameters in iteration ", i, " unreliable\n")
      }

      ## NA values in sd.corr cause a stop() already in marLogLik()
      ## if(any(is.na(parVar$par))){
      ##     updateVarianceControl[[1]] <- "updateParams_optim"
      ##     updateVarianceControl[[2]]$method <-
      ##         if (length(sd.corr) == 1L) "Brent" else "Nelder-Mead"
      ##     cat("  WARNING: at least one updated variance parameter is not a number\n",
      ##         "\t-> NO UPDATE of variance\n",
      ##         "\t-> SWITCHING to robust", dQuote(updateVarianceControl[[2]]$method),
      ##         "for variance updates\n")
      ## } else
      sd.corr <- parVar$par

      ## overall convergence ?
      if( (parReg$rel.tol < cntrl.stop$tol) && (parVar$rel.tol < cntrl.stop$tol)
          && (parReg$convergence==0) && (parVar$convergence==0) )
        convergence <- 0

      ## exit loop if no more change in parameters (maybe false convergence)
      if (parReg$rel.tol == 0 && parVar$rel.tol == 0)
        break
    }
  }

  if(verbose > 0) {
    cat("\n")
    cat(as.character(Sys.time()), ":", if (convergence==0)
      "Optimization converged" else "Optimization DID NOT CONVERGE", "\n\n")
  }

  ll <- penLogLik(theta=theta,sd.corr=sd.corr,model=model)
  fisher <- penFisher(theta=theta,sd.corr=sd.corr,model=model)
  dimnames(fisher) <- list(names(theta), names(theta))
  margll <- marLogLik(sd.corr=sd.corr, theta=theta, model=model)

  list(theta=theta, sd.corr=sd.corr,
       loglik=ll, margll=margll,fisher=fisher,
       convergence=convergence, dim=c(fixed=dimFE.d.O,random=dimRE))

}

## check analytical score functions and Fisher informations for
## a given model (the result of interpretControl(control, stsObj))
## and given parameters theta (regression par.) and sd.corr (variance par.).
## This is a wrapper around functionality of the numDeriv and maxLik packages.
checkAnalyticals <- function (model,
                              theta = model$initialTheta,
                              sd.corr = model$initialSigma,
                              methods = c("numDeriv","maxLik"))
{
  cat("\nPenalized log-likelihood:\n")
  resCheckPen <- sapply(methods, function(derivMethod) {
    if (requireNamespace(derivMethod)) {
      do.call(paste("checkDerivatives", derivMethod, sep="."),
              args=alist(penLogLik, penScore, penFisher, theta,
                         sd.corr=sd.corr, model=model))
    }
  }, simplify=FALSE, USE.NAMES=TRUE)
  if (length(resCheckPen) == 1L) resCheckPen <- resCheckPen[[1L]]

  resCheckMar <- if (length(sd.corr) == 0L) list() else {
    cat("\nMarginal log-likelihood:\n")
    fisher.unpen <- attr(penFisher(theta, sd.corr, model, attributes=TRUE),
                         "fisher")
    resCheckMar <- sapply(methods, function(derivMethod) {
      if (requireNamespace(derivMethod)) {
        do.call(paste("checkDerivatives", derivMethod, sep="."),
                args=alist(marLogLik, marScore, marFisher, sd.corr,
                           theta=theta, model=model, fisher.unpen=fisher.unpen))
      }
    }, simplify=FALSE, USE.NAMES=TRUE)
    if (length(resCheckMar) == 1L) resCheckMar[[1L]] else resCheckMar
  }

  list(pen = resCheckPen, mar = resCheckMar)
}
checkDerivatives.numDeriv <- surveillance:::checkDerivatives.numDeriv
checkDerivatives.maxLik <- surveillance:::checkDerivatives.maxLik
