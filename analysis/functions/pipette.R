#' Simulate pipetting error
#' 
#' A wide variety of sources can contribute error, and each of these results in
#' errors of differing magnitudes. This is a coarse simulation of general error.
#' 
#' @param vol Volume to pipette in microliters
#' 
#' @example pipette(c(1, 22.5, 360))
#' 
#' @references Rainin's poster on 'good pipetting practice'
#' 
#' @export
pipette <- function(vol){
  if(any(vol <= 0)){
    stop('volume must be positive')
  }
  MEAN <- vol
  N <- length(vol)
  SCALE <- rep(0.01, N)
  SHAPE <- MEAN/SCALE
  return(rgamma(n = N, shape = SHAPE, scale = SCALE)) # use shape and rate (alpha and beta) with mean = alpha/beta
}
