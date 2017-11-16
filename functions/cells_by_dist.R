#' Compute an example density of cells in a sample
#' 
#' Computes a density (cells per liter) in a sample at a given distance from 
#' organisms of a given mass
#' 
#' @param dist distance from organism in meters
#' @param mass mass of organism in grams
#' 
#' @return 
#' 
#' @examples 
#' suppose we have 100 individuals, each 10 grams, all 10 meters from our sample
#' N_ind <- 100
#' masss <- rep(1, N_ind)
#' dists <- rep(10, N_ind)
#' cells <- cells_by_dist(dists, masss)
#' sum(cells)
#' 
#' @export
cells_by_dist <- function(dist, mass){
  if(length(dist) != length(mass)){
    stop('dist and mass must be vectors of equal length')
  }
  
  # formulation 1:
  LAMBDA <- (mass^0.5)/dist 
  
  # formulation 2:
  # SOME_CONSTANT <- 1
  # LAMBDA <- rexp(n = length(dist), rate = SOME_CONSTANT*mass/dist)
  
  cells <- rpois(n = length(dist), lambda = LAMBDA)
  
  return(cells)
}

# plot(dists, cells,
#      xlim = range(dists), ylim = c(0, 10),
#      pch = 1, lwd = 2, col = 'slateblue',
#      xlab = 'meters',
#      las = 1)
