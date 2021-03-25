#' @title Refit a HHH4 model to a HHH4ZI model
#' @export
hhh4ZI.hhh4 <- function(object, # a HHH4 object
                        control = list(f = ~-1,
                                  lag = 1,
                                  lag.unitSpecific = FALSE
                        ),... # control list for zero model part
                        ){
  stopifnot(is.list(control))
  # all_control <- object$control
  # all_control$zi <- control

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
  hhh4ZI.sts(object$stsObj, control)
}
