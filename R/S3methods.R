#' @method BIC lciv
#' @export
BIC.lciv <- function(object, ...){
  AIC(object, k = log(nrow(object$gradientObs)), ...)
}
