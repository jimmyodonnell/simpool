#' Place individuals in a volume
#' 
#' Simulate coordinates of a set of individuals within a 3-dimensional space by 
#'   drawing from the beta distribution.
#' 
#' @param N_ind (numeric, length=1) Number of individuals
#' @param dims (numeric, length=3) Dimensions of volume (X,Y,Z)
#' @param distr_x (character, length=1) Description of the distribution of 
#'   individuals in dimension 1 (X). If no value (NULL) is given, one will be 
#'   chosen at random. 
#'   options: 'even', 'center', 'lower', 'lower_x', 'upper', 'upper_x'
#' @param distr_y same as distr_x
#' @param distr_z same as distr_x
#' 
#' @return A matrix of 3 columns (1 per dimension) for each individual (rows).
#' 
#' @examples
#' place_inds(N_ind = 20, dims = c(50, 25, 3), 
#'   distr_x = 'center', distr_y = 'lower', distr_z = 'upper_x')
#' 
#' @export
place_inds <- function(
  N_ind, dims, distr_x = NULL, distr_y = NULL, distr_z = NULL)
{

  # these values describe the distribution of individuals in the volume
  pref_names <- c('even', 'center', 'lower', 'lower_x', 'upper', 'upper_x')
  pref_options <- data.frame(
    name = pref_names, 
    par1 = as.numeric(rep(NA, length(pref_names))), 
    par2 = as.numeric(rep(NA, length(pref_names)))
    )
  pref_options[pref_options$name == 'even',2:3]    <- c(1, 1)
  pref_options[pref_options$name == 'center',2:3]  <- c(10, 10)
  pref_options[pref_options$name == 'lower',2:3]   <- c(1, 2)
  pref_options[pref_options$name == 'lower_x',2:3] <- c(1, 10)
  pref_options[pref_options$name == 'upper',2:3]   <- c(2, 1)
  pref_options[pref_options$name == 'upper_x',2:3] <- c(10, 1)
  
  pref_opts_full <- paste(c(pref_names, 'NULL'), collapse = ", ")

  distr <- character()
  for(i in c('distr_x', 'distr_y', 'distr_z')){
    if(is.null(get(i))){
      distr[i] <- sample(pref_names, size = 1)
    } else {
      distr[i] <- get(i)
    }
  }

  for(i in seq_along(distr)){
    if(!(distr[i] %in% pref_names)){
      msg <- paste('invalid distribution option. Options: ', pref_opts_full)
      stop(msg)
    }
  }

  PARS <- matrix(NA, nrow = 3, ncol = 2)
  for(i in seq_along(distr)){
    PARS[i,] <- as.numeric(pref_options[pref_options$name == distr[i],2:3])
  }

  # output: matrix of coordinates for each individual
  coords <- matrix(data = NA, nrow = N_ind, ncol = 3)
  
  for(i in 1:3){
    coords[,i] <- dims[i] * 
      rbeta(n = N_ind, shape1 = PARS[i,1], shape2 = PARS[i,2])
  }

  colnames(coords) <- c('X', 'Y', 'Z')

  return(coords)
  
}
