#' @title Refit a HHH4 model to a HHH4ZI model
#' @export
hhh4ZI.hhh4 <- function(object, # a HHH4 object
                        zi = list(f = ~-1,
                                  lag = 1,
                                  lag.unitSpecific = FALSE
                        ) # control list for zero model part
                        ){
  stopifnot(is.list(zi))
  control <- object$control
  control$zi <- zi
  # use fitted coefficients in hhh4 object as start value
  control$start <- surveillance:::hhh4coef2start(object)
  control$start$sd.corr <- NULL
  ## call main function
  hhh4ZI.sts(object$stsObj, control)
}
