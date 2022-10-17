#' @title Refit a HHH4 model with a ZI Component
#' @description Refits a previous \code{\link{hhh4}} fit with a ZI component
#' using \code{\link{hhh4ZI}}.
#' @param object a \code{"\link{hhh4}"} fit.
#' @param control either a list of specifications for the \code{zi} component
#' of \code{\link{hhh4ZI}} or a \code{hhh4} control list with added \code{zi}
#' specification.
#' @param ... further (non-\code{zi}) control elements.
#' @importFrom utils upgrade
#' @export
upgrade.hhh4 <- function(object,
                         control = list(f = ~1, lag = 1, lag.unitSpecific = FALSE),
                         ...)
{
  stopifnot(is.list(control))

  ## use update.hhh4() to merge control arguments
  args <- if ("f" %in% names(control))  # zi control list
    list(zi = control, ...)
  else  # hhh4 control list with zi list included
    c(control, list(...))
  args$object <- quote(object)
  args$evaluate <- FALSE
  control <- do.call(update.hhh4ZI, args)

  # use fitted coefficients in hhh4 object as start value, dimension is different
  # when having random effect
  control$start$sd.corr <- NULL
  ## call main function
  hhh4ZI(object$stsObj, control)
}
