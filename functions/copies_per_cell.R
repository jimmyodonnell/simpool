#' Number of mt genome copies per cell
#' 
#' This function returns counts of mitochondrial genomes for each of N cells.
#' 
#' @param cells Number of cells 
#' 
#' @return number of copies for each of 'cells' cells
#' 
#' @references number of mt genomes per cell: 220-1720 (DOI:10.1002/jcp.1041360316)
#' 
#' @examples copies_per_cell(100)
#' 
#' @export
copies_per_cell <- function(cells)
{
  MEAN <- 1000
  VARIANCE <- 1.5e5
  SIZE <- (MEAN^2)/VARIANCE
  COPIES <- rnbinom(n = cells, mu = MEAN, size = SIZE)
  return(COPIES)
}

# boxplot(copies, ylim = c(0, 3000), las = 1)
# fivenum(copies)
# plot(copies, ylim = c(0,3000), las = 1)
# abline(h = fivenum(copies), col = hsv(1,0.5,1), lty = 2)
