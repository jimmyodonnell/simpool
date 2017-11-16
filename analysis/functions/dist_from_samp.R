#' Calculate distance in three dimensions between a point and a set
#' 
#' Uses the Pythagorean formula to compute the distance between a single point 
#' and each of a set of points in three dimensions. Where 
#' D = sqrt(sum((p[i] - q[i])^2))
#' for each of three (i) dimensions of points p and q. 
#' 
#' @param sample_coords matrix Coordinates of sample (nrow = 1, ncol = 3).
#' @param other_coords matrix Coordinates of other points (nrow = N, ncol = 3).
#' 
#' @example
#' test_samp <- c(0,0,0)
#' test_coords <- matrix(data = c(0,0,1,0,2,0,3,0,0), ncol = 3, byrow = TRUE)
#' dist_from_samp(test_samp, test_coords)
#' 
#' @export
dist_from_samp <- function(sample_coords, other_coords)
{
  result <- sqrt(
    (other_coords[,1] - sample_coords[1])^2 + 
    (other_coords[,2] - sample_coords[2])^2 + 
    (other_coords[,3] - sample_coords[3])^2
    )
  return(result)
}
