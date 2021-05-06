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
### Plot-method(s) for fitted hhh4() models
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2020 Sebastian Meyer
### $Revision$
### $Date$
################################################################################
#' @title lots for Fitted hhh4ZI-models
#'
#' @description This function is the equivalent of \code{surveillance::plot.hhh4} for model fits of class
#' \code{hhh4ZI}, obtained from \code{hhh4ZI}. The arguments are the
#' same as in \code{surveillance::plot.hhh4}, with one more component option for zero inflation part..
#'
#' @import graphics
#' @export
# no maxEV
plot.hhh4ZI <- function (x,
                         type = c("fitted", "season", "maps", "ri", "neweights"),
                         ...)
{
  stopifnot(x$convergence)
  cl <- sys.call()  # not match.call() because plotHHH4_season() has no 'x'
  ## remove the type argument from the call
  if (is.null(names(cl)) && nargs() > 1L) { # unnamed call plot(x, type)
    cl[[3L]] <- NULL  # remove the second argument
  } else {
    cl$type <- NULL
  }
  cl[[1L]] <- as.name(paste("plotHHH4ZI", match.arg(type), sep="_"))
  eval(cl, envir = parent.frame())
}


###
### Time series of fitted component means and observed counts for selected units
###
#' @rdname plot.hhh4ZI
#' @importFrom grDevices n2mfrow
#' @export
plotHHH4ZI_fitted <- function (x, units = 1, names = NULL,
                               col = c("grey85", "blue", "orange"),
                               pch = 19, pt.cex = 0.6, pt.col = 1,
                               par.settings = list(),
                               legend = TRUE, legend.args = list(),
                               legend.observed = FALSE,
                               decompose = NULL, total = FALSE, meanHHH = NULL, ...)
{
  if (total) {
    units <- "Overall"  # only used as a label
  } else if (is.null(units)) {
    units <- seq_len(x$nUnit)
  }
  if (!is.null(names)) stopifnot(length(units) == length(names))
  if (isTRUE(decompose)) decompose <- colnames(x$stsObj)

  ## get decomposed mean => no need to compute it in each plotHHH4_fitted1()
  gamma <- x$gamma
  if (is.null(meanHHH)) {
    meanHHH <- if (is.null(decompose)) {
      mean_HHH <- surveillance::meanHHH(x$coefficients, surveillance:::terms.hhh4(x))
      lapply(mean_HHH[1:5],"*", 1 - gamma)
    } else {
      mean_HHH <- surveillance::decompose.hhh4(x)
      x2 <- array(dim = dim(mean_HHH), dimnames = dimnames(mean_HHH))
      for (i in 1 : (dim(mean_HHH)[3]))
        x2[,,i] <- mean_HHH[,,i] * (1-gamma)
      x2
    }
  }


  ## check color vector
  col <- if (is.null(decompose) && length(col) == 4) {
    ## compatibility with surveillance < 1.10-0
    pt.col <- col[4L]
    rev(col[-4L])
  } else {
    surveillance:::plotHHH4_fitted_check_col_decompose(col, decompose)
  }

  ## setup graphical parameters
  if (is.list(par.settings)) {
    par.defaults <- list(mfrow = sort(n2mfrow(length(units))),
                         mar = c(4,4,2,0.5)+.1, las = 1)
    par.settings <- modifyList(par.defaults, par.settings)
    opar <- do.call("par", par.settings)
    on.exit(par(opar))
  }

  ## legend options
  if (is.logical(legend)) legend <- which(legend)
  if (!is.list(legend.args)) {
    if (length(legend) > 0)
      warning("ignored 'legend' since 'legend.args' is not a list")
    legend <- integer(0L)
  }
  if (length(legend) > 0) {
    legendidx <- 1L + c(
      if (legend.observed && !is.na(pch)) 0L,
      if (is.null(decompose)) {
        which(c("ne","ar","end") %in% componentsHHH4ZI(x))
      } else seq_along(col))
    default.args <- list(
      x="topright", col=c(pt.col,rev(col))[legendidx], lwd=6,
      lty=c(NA,rep.int(1,length(col)))[legendidx],
      pch=c(pch,rep.int(NA,length(col)))[legendidx],
      pt.cex=pt.cex, pt.lwd=1, bty="n", inset=0.02,
      legend=if (is.null(decompose)) {
        c("observed","spatiotemporal","autoregressive","endemic")[legendidx]
      } else c("observed", rev(decompose), "endemic")[legendidx]
    )
    legend.args <- modifyList(default.args, legend.args)
  }

  ## plot fitted values region by region
  meanHHHunits <- vector(mode="list", length=length(units))
  names(meanHHHunits) <- if (is.character(units)) units else colnames(x$stsObj)[units]
  for(i in seq_along(units)) {
    meanHHHunits[[i]] <- plotHHH4ZI_fitted1(x, unit=units[i], main=names[i],
                                            col=col, pch=pch, pt.cex=pt.cex, pt.col=pt.col,
                                            decompose=decompose, total=total, meanHHH=meanHHH, ...)
    if (i %in% legend) do.call("legend", args=legend.args)
  }
  invisible(meanHHHunits)
}

### plot estimated component means for a single region
#' @rdname plot.hhh4ZI
#' @export
plotHHH4ZI_fitted1 <- function(x, unit=1, main=NULL,
                               col=c("grey85", "blue", "orange"),
                               pch=19, pt.cex=0.6, pt.col=1, border=col,
                               start=x$stsObj@start, end=NULL, xaxis=NULL,
                               xlim=NULL, ylim=NULL, xlab="", ylab="No. infected",
                               hide0s=FALSE, decompose=NULL, total=FALSE, meanHHH=NULL)
{
  stsObj <- x$stsObj
  if (!total && is.character(unit) &&
      is.na(unit <- match(.unit <- unit, colnames(stsObj))))
    stop("region '", .unit, "' does not exist")
  if (is.null(main))
    main <- if (total) "Overall" else colnames(stsObj)[unit]
  if (isTRUE(decompose))
    decompose <- colnames(stsObj)

  ## get observed counts
  obs <- if (total) rowSums(observed(stsObj)) else observed(stsObj)[,unit]

  ## time range for plotting
  start0 <- surveillance:::yearepoch2point(stsObj@start, stsObj@freq, toleft=TRUE)
  start <- surveillance:::yearepoch2point(start, stsObj@freq)
  tp <- start0 + seq_along(obs)/stsObj@freq # all observation time points
  if (start < start0 || start > tp[length(tp)])
    stop("'start' is not within the time range of 'x$stsObj'")
  end <- if(is.null(end)) tp[length(tp)] else surveillance:::yearepoch2point(end,stsObj@freq)
  stopifnot(start < end)
  tpInRange <- which(tp >= start & tp <= end)            # plot only those
  tpInSubset <- intersect(x$control$subset, tpInRange)   # fitted time points

  ## use time indexes as x-values for use of addFormattedXAxis()
  if (is.list(xaxis)) {
    tp <- seq_along(obs)
    start <- tpInRange[1L]
    end <- tpInRange[length(tpInRange)]
  }

  ## get fitted component means
  if (is.null(meanHHH)) {
    meanHHH <- if (is.null(decompose)) {
      mean_HHH <- surveillance::meanHHH(x$coefficients, surveillance:::terms.hhh4(x))
      lapply(mean_HHH[1:5],"*", 1 - gamma)
    } else {
      mean_HHH <- surveillance::decompose.hhh4(x)
      x2 <- array(dim = dim(mean_HHH), dimnames = dimnames(mean_HHH))
      for (i in 1 : (dim(mean_HHH)[3]))
        x2[,,i] <- mean_HHH[,,i] * (1-gamma)
      x2

    }
  }
  meanHHHunit <- if (is.null(decompose)) {
    if (total) {
      sapply(meanHHH, rowSums)
    } else {
      sapply(meanHHH, "[", i=TRUE, j=unit)
    }
  } else {
    if (!setequal(decompose, dimnames(meanHHH)[[3L]][-1L]))
      stop("'decompose' must be (a permutation of) the fitted units")
    if (total) {
      apply(meanHHH[,,c("endemic",decompose)], c(1L, 3L), sum)
    } else {
      meanHHH[,unit,c("endemic",decompose)]
    }
  }
  stopifnot(is.matrix(meanHHHunit), !is.null(colnames(meanHHHunit)),
            nrow(meanHHHunit) == length(x$control$subset))
  meanHHHunit <- meanHHHunit[x$control$subset %in% tpInRange,,drop=FALSE]
  if (any(is.na(meanHHHunit))) { # -> polygon() would be wrong
    ## could be due to wrong x$control$subset wrt the epidemic lags
    ## a workaround is then to set 'start' to a later time point
    stop("predicted mean contains missing values")
  }

  ## check color vector
  col <- if (is.null(decompose) && length(col) == 4L) {
    ## compatibility with surveillance < 1.10-0
    pt.col <- col[4L]
    rev(col[-4L])
  } else {
    surveillance:::plotHHH4_fitted_check_col_decompose(col, decompose)
  }

  ## establish basic plot window
  if (is.null(ylim)) ylim <- c(0, max(obs[tpInRange],na.rm=TRUE))
  plot(c(start,end), ylim, xlim=xlim, xlab=xlab, ylab=ylab, type="n",
       xaxt = if (is.list(xaxis)) "n" else "s")
  if (is.list(xaxis)) do.call("addFormattedXAxis", c(list(x = stsObj), xaxis))
  title(main=main, line=0.5)

  ## draw polygons
  if (is.null(decompose)) {
    non0 <- which(c("end", "ar", "ne") %in% componentsHHH4ZI(x))
    surveillance:::plotComponentPolygons(
      x = tp[tpInSubset],
      y = meanHHHunit[,c("endemic", "epi.own", "epi.neighbours")[non0],drop=FALSE],
      col = col[non0], border = border[non0], add = TRUE)
  } else {
    non0 <- apply(X = meanHHHunit > 0, MARGIN = 2L, FUN = any)
    surveillance:::plotComponentPolygons(x = tp[tpInSubset], y = meanHHHunit[, non0, drop = FALSE],
                                         col = col[non0], border = border[non0], add = TRUE)
  }

  ## add observed counts within [start;end]
  ptidx <- if (hide0s) intersect(tpInRange, which(obs > 0)) else tpInRange
  points(tp[ptidx], obs[ptidx], col=pt.col, pch=pch, cex=pt.cex)

  ## invisibly return the fitted component means for the selected region
  invisible(meanHHHunit)
}

###
### Maps of the fitted mean components averaged over time
###
#' @rdname plot.hhh4ZI
#' @importFrom grDevices n2mfrow
#' @importFrom sp spplot
#' @import methods
#' @export
plotHHH4ZI_maps <- function (x,
                             which = c("mean", "endemic",
                                       "epi.own", "epi.neighbours", "zi"),
                             prop = FALSE, main = which, zmax = NULL, col.regions = NULL,
                             labels = FALSE, sp.layout = NULL, ...,
                             map = x$stsObj@map, meanHHH = NULL)
{
  which <- match.arg(which, several.ok = TRUE)
  if (is.null(col.regions))
    col.regions <- surveillance:::.hcl.colors(10)

    ## extract district-specific mean components
  if (is.null(meanHHH)) {
    gamma <- x$gamma
    meanHHH <- surveillance::meanHHH(x$coefficients, surveillance:::terms.hhh4(x))
    meanHHH <- lapply(meanHHH[1:5],"*", 1 - gamma)
  }

  element <- meanHHH[c("mean", "endemic", "epi.own", "epi.neighbours")]
  element$zi <- gamma
  ## select relevant components and convert to an array
  meanHHH <- simplify2array(
    element,
    higher = TRUE)

  ## convert to proportions
  if (prop) {
    meanHHH[,,-1L] <- meanHHH[,,-1L,drop=FALSE] / c(meanHHH[,,1L])
  }

  ## select only 'which' components
  meanHHH <- meanHHH[,,which,drop=FALSE]

  ## check map
  map <- as(map, "SpatialPolygonsDataFrame")
  if (!all(dimnames(meanHHH)[[2L]] %in% row.names(map))) {
    stop("'row.names(map)' do not cover all fitted districts")
  }

  ## average over time
  comps <- as.data.frame(colMeans(meanHHH, dims = 1))

  ## attach to map data
  map@data <- cbind(map@data, comps[row.names(map),,drop=FALSE])

  ## color key range
  if (is.null(zmax)) {
    zmax <- if (prop) {
      ceiling(10*sapply(comps, max))/10
    } else ceiling(sapply(comps, max))
    ## sub-components should have the same color range
    .idxsub <- setdiff(seq_along(zmax), match("mean", names(zmax)))
    zmax[.idxsub] <- suppressWarnings(max(zmax[.idxsub]))
  }

  ## add sp.layout item for district labels
  if (!is.null(layout.labels <- layout.labels(map, labels))) {
    sp.layout <- c(sp.layout, list(layout.labels))
  }

  ## produce maps
  grobs <- mapply(
    FUN = function (zcol, main, zmax)
      if (is.na(zmax)) { # automatic color breaks over range of values
        if(zcol == "zi") zmax = 1
        spplot(map, zcol = zcol, main = main,
               cuts = length(col.regions) - 1L,
               col.regions = col.regions, sp.layout = sp.layout, ...)
      } else { # breakpoints from 0 to zmax
        if(zcol == "zi") zmax = 1
        spplot(map, zcol = zcol, main = main,
               at = seq(0, zmax, length.out = length(col.regions) + 1L),
               col.regions = col.regions, sp.layout = sp.layout, ...)
      },
    zcol = names(comps), main = main, zmax = zmax,
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  if (length(grobs) == 1L) {
    grobs[[1L]]
  } else {
    mfrow <- sort(n2mfrow(length(grobs)))
    gridExtra::grid.arrange(grobs = grobs, nrow = mfrow[1L], ncol = mfrow[2L])
  }
}


###
### Map of estimated random intercepts of a specific component
###
#' @rdname plot.hhh4ZI
#' @importFrom grDevices cm.colors
#' @importFrom sp spplot
#' @import methods
#' @export
plotHHH4ZI_ri <- function (x, component, exp = FALSE,
                           at = list(n = 10), col.regions = cm.colors(100),
                           colorkey = TRUE, labels = FALSE, sp.layout = NULL,
                           gpar.missing = list(col="darkgrey", lty=2, lwd=2),
                           ...)
{
  ranefmatrix <- ranef.hhh4ZI(x, tomatrix=TRUE)
  if (is.null(ranefmatrix)) stop("model has no random effects")
  stopifnot(length(component) == 1L)
  if (is.na(comp <- pmatch(component, colnames(ranefmatrix))))
    stop("'component' must (partially) match one of ",
         paste(dQuote(colnames(ranefmatrix)), collapse=", "))

  map <- as(x$stsObj@map, "SpatialPolygonsDataFrame")
  if (length(map) == 0L) stop("'x$stsObj' has no map")
  map$ranef <- ranefmatrix[,comp][row.names(map)]
  .range <- c(-1, 1) * max(abs(map$ranef), na.rm = TRUE)  # 0-centered
  if (exp) {
    map$ranef <- exp(map$ranef)
    .range <- exp(.range)
  }

  if (is.list(at)) {
    at <- modifyList(list(n = 10, range = .range), at)
    at <- if (exp) {
      stopifnot(at$range[1] > 0)
      scales::log_breaks(n = at$n)(at$range)
    } else {
      seq(at$range[1L], at$range[2L], length.out = at$n)
    }
    if (exp && isTRUE(colorkey))
      colorkey <- list(at = log(at),
                       labels = list(at = log(at), labels = at))
  }

  if (is.list(gpar.missing) && any(is.na(map$ranef))) {
    sp.layout <- c(sp.layout,
                   c(list("sp.polygons", map[is.na(map$ranef),]),
                     gpar.missing))
  }
  if (!is.null(layout.labels <- layout.labels(map, labels))) {
    sp.layout <- c(sp.layout, list(layout.labels))
  }

  spplot(map[!is.na(map$ranef),], zcol = "ranef",
         sp.layout = sp.layout, col.regions = col.regions,
         at = at, colorkey = colorkey, ...)
}

ranef.hhh4ZI <- function (object, tomatrix = FALSE, intercept = FALSE, ...)
{
  if (object$dim[2L] > 0){
    ranefvec <- tail(surveillance:::coef.hhh4(object, ...), object$dim[2L])
  } else return(NULL)
  if (intercept) tomatrix <- TRUE
  if (!tomatrix) return(ranefvec)

  ## transform to a nUnits x c matrix (c %in% 1:3)
  model <- terms.hhh4ZI(object)
  idxRE <- model$indexRE
  idxs <- unique(idxRE)
  mat <- vapply(X = idxs, FUN = function (idx) {
    RE <- ranefvec[idxRE==idx]
    Z <- model$terms["Z.intercept",][[idx]]
    "%m%" <- get(model$terms["mult",][[idx]])
    c(Z %m% RE)
  }, FUN.VALUE = numeric(model$nUnits), USE.NAMES = FALSE)
  dimnames(mat) <- list(
    colnames(model$response),
    model$namesFE[match(idxs, model$indexFE)]
  )

  if (intercept) {
    FE <- object$coefficients[colnames(mat)]
    mat <- t(t(mat) + FE)
  }

  return(mat)
}

###
### Plot estimated seasonality (sine-cosine terms) of one or several hhh4-fits
### either as multiplicative effect on the 'components' (intercept=FALSE)
### or with intercept=TRUE, which only makes sense if there are no further
### non-centered covariates and offsets.
###
#' @rdname plot.hhh4ZI
#' @importFrom grDevices n2mfrow
#' @export
plotHHH4ZI_season <- function (...,
                               components = NULL, intercept = FALSE,
                               xlim = NULL, ylim = NULL,
                               xlab = NULL, ylab = "", main = NULL,
                               par.settings = list(), matplot.args = list(),
                               legend = NULL, legend.args = list(),
                               refline.args = list(), unit = 1)
{
  objnams <- unlist(lapply(match.call(expand.dots=FALSE)$..., deparse))
  objects <- surveillance:::getHHH4list(..., .names = objnams)
  freq <- attr(objects, "freq")
  components <- if (is.null(components)) {
    intersect(c("end", "ar", "ne", "zi"), unique(unlist(
      lapply(objects, componentsHHH4ZI), use.names = FALSE)))
  } else {
    match.arg(components, choices = c("ar", "ne", "end", "zi","maxEV"),
              several.ok = TRUE)
  }

  ## x-axis
  if (is.null(xlim))
    xlim <- c(1,freq)
  if (is.null(xlab))
    xlab <- if (freq==52) "week" else if (freq==12) "month" else "time"

  ## auxiliary function for an argument list "x" with named "defaults" list
  withDefaults <- function(x, defaults)
  {
    if (is.null(x)) defaults else if (is.list(x)) {
      if (is.null(names(x))) {    # x must be complete
        stopifnot(length(x) == length(defaults))
        setNames(x, names(defaults))
      } else modifyList(defaults, x) # x might be a subset of parameters
    } else if (is.atomic(x)) {
      setNames(rep(list(x), length(defaults)), names(defaults))
    } else stop("'", deparse(substitute(x)), "' is not suitably specified")
  }

  ## component-specific arguments
  ylim <- withDefaults(ylim,
                       list(ar=NULL, ne=NULL, end=NULL, zi = NULL, maxEV=NULL))
  ylab <- withDefaults(ylab,
                       list(ar=expression(hat(lambda)),
                            ne=expression(hat(phi)),
                            end=expression(hat(nu)),
                            zi = expression(hat(gamma)),
                            maxEV="dominant eigenvalue"))
  main <- withDefaults(main,
                       list(ar="autoregressive component",
                            ne="spatiotemporal component",
                            end="endemic component",
                            zi = "zero-inflation component",
                            maxEV="dominant eigenvalue"))
  anyMain <- any(unlist(lapply(main, nchar),
                        recursive=FALSE, use.names=FALSE) > 0)

  ## basic graphical settings
  if (is.list(par.settings)) {
    par.defaults <- list(mfrow=sort(n2mfrow(length(components))),
                         mar=c(4,5,if(anyMain) 2 else 1,1)+.1, las=1)
    par.settings <- modifyList(par.defaults, par.settings)
    opar <- do.call("par", par.settings)
    on.exit(par(opar))
  }

  ## line style
  matplot.args <- modifyList(list(type="l", col=c(1,2,6,3), lty=c(1,3,2,4),
                                  lwd=1.7, cex=1, pch=NULL),
                             matplot.args)

  ## legend options
  if (is.null(legend)) legend <- length(objects) > 1
  if (is.logical(legend)) legend <- which(legend)
  if (!is.list(legend.args)) {
    if (length(legend) > 0)
      warning("ignored 'legend' since 'legend.args' is not a list")
    legend <- integer(0L)
  }
  if (length(legend) > 0) {
    default.args <- c(
      list(x="topright", inset=0.02, legend=names(objects), bty="n"),
      matplot.args[c("col", "lwd", "lty", "pch")],
      with(matplot.args, list(pt.cex=cex, text.col=col))
    )
    legend.args <- modifyList(default.args, legend.args)
  }

  ## plot seasonality in individual model components
  seasons <- list()
  for(comp in setdiff(components, "maxEV")){
    s2 <- lapply(objects, getSeason, component = comp, unit = unit)
    seasons[[comp]] <- exp(vapply(s2, FUN = if (intercept) {
      function (intseas) do.call("+", intseas)
    } else {
      function (intseas) intseas$season  # disregard intercept
    }, FUN.VALUE = numeric(freq), USE.NAMES = TRUE))
    do.call("matplot",              # x defaults to 1:freq
            c(list(seasons[[comp]], xlim=xlim, ylim=ylim[[comp]],
                   xlab=xlab, ylab=ylab[[comp]], main=main[[comp]]),
              matplot.args))
    if (is.list(refline.args) && !intercept && any(seasons[[comp]] != 1))
      do.call("abline", modifyList(list(h=1, lty=3, col="grey"),
                                   refline.args))
    if (match(comp, components) %in% legend)
      do.call("legend", legend.args)
  }

  ## plot seasonality of dominant eigenvalue
  # if ("maxEV" %in% components) {
  #   seasons[["maxEV"]] <- vapply(objects, FUN = function (obj) {
  #     getMaxEV_season(obj)$maxEV.season
  #   }, FUN.VALUE = numeric(freq), USE.NAMES = TRUE)
  #   do.call("matplot",
  #           c(list(seasons[["maxEV"]], xlim=xlim,
  #                  ylim=if (is.null(ylim[["maxEV"]]))
  #                    c(0,max(2,seasons[["maxEV"]])) else ylim[["maxEV"]],
  #                  xlab=xlab, ylab=ylab[["maxEV"]],
  #                  main=main[["maxEV"]]), matplot.args))
  #   if (is.list(refline.args))
  #     do.call("abline", modifyList(list(h=1, lty=3, col="grey"),
  #                                  refline.args))
  #   if (4 %in% legend) do.call("legend", legend.args)
  # }

  ## invisibly return the data that has been plotted
  invisible(seasons)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get estimated intercept and seasonal pattern in the different components
# CAVE: other covariates and offsets are ignored
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getSeason <- function(x, component = c("end", "ar", "ne", "zi"), unit = 1)
{
  stopifnot(inherits(x, "hhh4ZI"))
  component <- match.arg(component)
  startseason <- surveillance:::getSeasonStart(x)
  freq <- x$stsObj@freq
  if (is.character(unit)) unit <- match(unit, colnames(x$stsObj))

  ## return -Inf is component is not in the model (-> exp(-Inf) = 0)
  if (!component %in% componentsHHH4ZI(x))
    return(list(intercept=-Inf, season=rep.int(-Inf, freq)))

  ## get the intercept
  est <- surveillance:::fixef.hhh4(x, reparamPsi=FALSE)
  intercept <- unname(est[grep(paste0("^", component, "\\.(1|ri)"), names(est))])
  if (length(intercept) == 0) {
    intercept <- 0 # no intercept (not standard)
  } else if (length(intercept) > 1) { # unit-specific intercepts
    if (length(intercept) != ncol(x$stsObj))
      stop(component,"-component has incomplete unit-specific intercepts")
    intercept <- intercept[unit]
    if (is.na(intercept)) stop("the specified 'unit' does not exist")
  }

  ## get seasonality terms (relying on sin(2*pi*t/52)-kind coefficient names)
  coefSinCos <- est[grep(paste0("^",component, "\\.(sin|cos)\\("), names(est))]
  if (unitspecific <- length(grep(").", names(coefSinCos), fixed=TRUE))) {
    if (unitspecific < length(coefSinCos))
      stop("cannot handle partially unit-specific seasonality")
    coefSinCos <- coefSinCos[grep(paste0(").",colnames(x$stsObj)[unit]),
                                  names(coefSinCos), fixed=TRUE)]
    ## drop .unitname-suffix since non-syntactic (cannot reformulate())
    names(coefSinCos) <- sub("\\)\\..+$", ")", names(coefSinCos))
  }
  if (length(coefSinCos)==0)
    return(list(intercept=intercept, season=rep.int(0,freq)))
  fSinCos <- reformulate(
    sub(paste0("^",component,"\\."), "", names(coefSinCos)),
    intercept=FALSE)
  mmSinCos <- model.matrix(fSinCos,
                           data=data.frame(t=startseason-1 + seq_len(freq)))

  ## Done
  list(intercept=intercept, season=as.vector(mmSinCos %*% coefSinCos))
}



###
### plot neighbourhood weight as a function of distance (neighbourhood order)
###
#' @rdname plot.hhh4ZI
#' @export
plotHHH4ZI_neweights <- surveillance::plotHHH4_neweights
