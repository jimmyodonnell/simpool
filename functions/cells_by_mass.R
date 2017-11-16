#' Cells per gram of biomass
#' 
#' compute a number of cells shed per second per gram of biomass
#' eDNA shedding rate range: 165-3368 pg of DNA per hour per gram of biomass
#' 3.3e9 bp = 3.59 pg
#' mass of one base pair in grams = 1.0794e-21
#' number of mt genomes per cell: 220-1720 (DOI:10.1002/jcp.1041360316)
#' number of mt genomes per mitochondrion: 1
#' basepairs in mt genome: 16000
#' 
#' @param mass mass of an organism, in grams
#' 
#' @return (int) A single number of cells shed per second
#' 
#' @examples
#' cells_by_mass(20)
#' 
#' @references Sassoubre Environ. Sci. Technol., 2016, 50 (19), pp 10456–10464
#' 
#' @export
cells_by_mass <- function(mass)
{
  LAMBDA <- mass
  cells <- rpois(n = 1, lambda = LAMBDA)
  return(cells)
}