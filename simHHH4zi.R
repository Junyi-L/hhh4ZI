simulate <- function (stsObj, control, ...) UseMethod("simulate")
simulate.hhh4ZI <- function (object, # result from a call to hhh4ZI
                           nsim=1, # number of replicates to simulate
                           seed=NULL,
                           y.start=NULL, # initial counts for epidemic components
                           subset=1:nrow(object$stsObj),
                           coefs=coef(object), # coefficients used for simulation
                           components=c("ar","ne","end"), # which comp to include
                           simplify=nsim>1, # counts array only (no full sts)
                           ...)
{ 
  ## Determine seed (this part is copied from stats:::simulate.lm with
  ## Copyright (C) 1995-2012 The R Core Team)
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)                     # initialize the RNG if necessary
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ## END seed
  
  cl <- match.call()
  theta <- if (missing(coefs)) coefs else surveillance:::checkCoefs(object, coefs)
  stopifnot(subset >= 1, subset <= nrow(object$stsObj))
  
  ## lags
  lag.ar <- object$control$ar$lag
  lag.ne <- object$control$ne$lag
  lag.gamma <- object$control$zi$lag
  maxlag <- max(lag.ar, lag.ne, lag.gamma)
  
  ## initial counts
  nUnits <- object$nUnit
  if (is.null(y.start)) { # set starting value to mean observed (in subset!)
    y.means <- ceiling(colMeans(observed(object$stsObj)[subset,,drop=FALSE]))
    y.start <- matrix(y.means, maxlag, nUnits, byrow=TRUE)
  } else {
    if (is.vector(y.start)) y.start <- t(y.start)
    if (ncol(y.start) != nUnits)
      stop(sQuote("y.start"), " must have nUnits=", nUnits, " columns")
    if (nrow(y.start) < maxlag)
      stop("need 'y.start' values for lag=", maxlag, " initial time points")
  }
  
  ## store model terms in the hhh4 object because we request them repeatedly
  ## (within get_exppreds_with_offsets() and directly afterwards)
  ## CAVE: for an ri()-model, building the terms affects the .Random.seed,
  ## so doing that twice would yield different simulations than pre-1.16.2
  if (is.null(object$terms))
    object$terms <- terms.hhh4ZI(object)
  
  ## get fitted exppreds nu_it, phi_it, lambda_it (incl. offsets, t in subset)
  exppreds <- surveillance:::get_exppreds_with_offsets(object, subset = subset, theta = theta)
  
  ## extract overdispersion parameters (simHHH4 assumes psi->0 means Poisson)
  model <- terms.hhh4ZI(object)
  psi <- surveillance:::splitParams(theta,model)$overdisp
  if (length(psi) > 1) # "NegBinM" or shared overdispersion parameters
    psi <- psi[model$indexPsi]
  
  ## weight matrix/array of the ne component
  neweights <- getNEweights(object, coefW(theta))
  
  ## set predictor to zero if not included ('components' argument)
  stopifnot(length(components) > 0, components %in% c("ar", "ne", "end"))
  getComp <- function (comp) {
    exppred <- exppreds[[comp]]
    if (comp %in% components) exppred else "[<-"(exppred, value = 0)
  }
  ar <- getComp("ar")
  ne <- getComp("ne")
  end <- getComp("end")

  gamma_end <- gammaZero(theta, model, subset, .ar = FALSE)
  gamma_ar <- theta[grepl("zi.AR", names(theta))]
  gamma_ar <- if(! model$zi.lag.unitSpecific) matrix(rep(gamma_ar, nUnits), ncol = nUnits) else
    matrix(gamma_ar, ncol = nUnits, byrow = TRUE)
  
  ## simulate
  simcall <- quote(
    simHHH4ZI(ar, ne, end, psi,gamma_end, gamma_ar,
              neweights, y.start, lag.ar, lag.ne, lag.gamma)
  )
  
  if (!simplify) {
    ## result template
    res0 <- object$stsObj[subset,]
    setObserved <- function (observed) {
      res0@observed[] <- observed
      res0
    }
    simcall <- call("setObserved", simcall)
  }
  res <- if (nsim==1) eval(simcall) else
    replicate(nsim, eval(simcall),
              simplify=if (simplify) "array" else FALSE)
  if (simplify) {
    dimnames(res)[1:2] <- list(subset, colnames(model$response))
    attr(res, "initial") <- y.start
    attr(res, "stsObserved") <- object$stsObj[subset,]
    class(res) <- c("hhh4ZIsims","hhh4sims")
  }
  
  ## Done
  attr(res, "call") <- cl
  attr(res, "seed") <- RNGstate
  res
}

terms.hhh4ZI <- function (x, ...)
{
  if (is.null(x$terms))
    interpretControl(x$control,x$stsObj) else x$terms
}

simHHH4ZI <- function(ar,     # lambda_it (nTime x nUnits matrix)
                      ne,     # phi_it (nTime x nUnits matrix)
                      end,    # nu_it (nTime x nUnits matrix, offset included)
                      psi,    # overdisp param(s) or numeric(0) (psi->0 = Poisson)
                      gamma_end,  # evaluated logit-scale linear predictor of gamma_it except lagged obs (nTime x nUnits matrix)
                      gamma_ar, #  coefficient of lagged obs in gamma_it
                      #(lag.gamma (in the order of lag.gamma) x nUnits matrix)
                      neW,    # weight matrix/array for neighbourhood component
                      start,  # starting counts (vector of length nUnits, or
                      # matrix with nUnits columns if lag > 1)
                      lag.ar = 1,
                      lag.ne = lag.ar,
                      lag.gamma = lag.ar){
  
  
  nTime <- nrow(end)
  nUnits <- ncol(end)
  
  ## check and invert psi since rnbinom() uses different parametrization
  size <- if (length(psi) == 0 ||
              isTRUE(all.equal(psi, 0, check.attributes=FALSE))) {
    stop("psi must be provided and not equal to zero")
  } else {
    if (!length(psi) %in% c(1, nUnits))
      stop("'length(psi)' must be ",
           paste(unique(c(1, nUnits)), collapse = " or "),
           " (number of units)")
    1/psi
  }
  
  # simulate from NegBin model
  ## unit-specific 'mean's and variance = mean + psi*mean^2
  ## where 'size'=1/psi and length(psi) == 1 or length(mean)
  rdistr <- function(n, mean) rnbinom(n, mu = mean, size = size)
  
  ## weighted sum of counts of other (neighbouring) regions
  ## params: y - vector with (lagged) counts of regions
  ##         W - nUnits x nUnits adjacency/weight matrix (0=no neighbour)
  wSumNE <- if (is.null(neW) || all(neW == 0)) { # includes the case nUnits==1
    function (y, W) numeric(nUnits)
  } else function (y, W) .colSums(W * y, nUnits, nUnits)
  
  ## W * y = (W[1, ]* y[1]; W[2, ]* y[2]; ...)
  
  ## initialize matrices for means mu_i,t and simulated data y_i,t
  mu <- y <- omega <- gamma <- matrix(0, nTime, nUnits)
  y <- rbind(start, y)
  nStart <- nrow(y) - nrow(mu)        # usually just 1 for lag=1
  timeDependentWeights <- length(dim(neW)) == 3
  if (!timeDependentWeights) neWt <- neW
  for(t in seq_len(nTime)){
    #browser()
    if (timeDependentWeights) neWt <- neW[,,t]
    gamma[t,] <- plogis(if(lag.gamma[1] > 0) {
      gamma_end[t,] +  (if(length(lag.gamma) >1) 
        colSums(y[nStart+t-lag.gamma,] * gamma_ar) else 
          y[nStart+t-lag.gamma,] * gamma_ar)
      } else gamma_end[t,])
    omega[t,] <- ifelse(runif(nUnits) < gamma[t,], 1, 0)
    
    #browser()
    ## mean mu_i,t = lambda*y_i,t-1 + phi*sum_j wji*y_j,t-1 + nu_i,t
    mu[t,] <-
      ar[t,] * y[nStart+t-lag.ar,] +
      ne[t,] * wSumNE(y[nStart+t-lag.ne,], neWt) +
      end[t,]
    ## Sample from NegBin or zero with that mean
    y[nStart+t,] <-
      ifelse(omega[t,] == 1, rep(0, nUnits),
             rdistr(nUnits, mu[t,]))
    
  }
  
  ## return simulated data without initial counts
  # list(y = y[-seq_len(nStart),,drop=FALSE], mu = mu, omega = omega)
  y[-seq_len(nStart),,drop=FALSE]
}

