#' Simulate templates for PCR
#' 
#' This function is to be used in concert with do_pcr. 
#' If 'stochastic' = TRUE, output is draws from multinomial distribution with 
#' 
#' @param N_sp Int > 0
#' 
#' @param N_out Int > 0
#' 
#' @param evenness String: 'even', 'skew.lin', 'skew.low', 'skew.med', 'skew.hi'
#' 
#' @param stochastic Logical. Should counts be stochastic?
#' 
#' @return Numeric vector of counts with sum = N_out. Note that outputs are not 
#' guaranteed to be integers if stochastic = FALSE. You are welcomed to wrap the
#' function in round() if you require integer outputs; however, this does not 
#' guarantee the total sum == N_out (but it will be close). 
#' When stochastic = TRUE, outputs are draws from a multinomial distribution 
#' where n = N_out, and the 'evenness' parameter controls the distribution of 
#' probabilities of each event (which in this context are 'species').
#' Note also that while the total length out is guaranteed to be = N_sp, there 
#' is no guarantee that all elements will be greater than zero. 
#' Aaaaand just ONE last thing... output is not sorted, so if the function is 
#' repeated multiple times, the order of the elements (species) are not 
#' guaranteed to coincide with their rank ordering.  
#' 
#' @examples sim_templates(4, 100, 'skew.lin', FALSE) # 25 25 25 25
#' 
#' @export
sim_templates <- function(N_sp, N_out, evenness, stochastic){
  
  even.opts <- c('even', 'skew.lin', 'skew.low', 'skew.med', 'skew.hi')
  
  if(!(evenness %in% even.opts)){
    these <- paste(even.opts, collapse = ", ")
    stop(paste('evenness must be one of:', these))
  }
  
  x <- 1:N_sp
  p <- function(Y){Y/sum(Y)}
  
  # even
  if(evenness == 'even'){
    y <- rep(1/N_sp, N_sp)
  }
  
  # Here is a linearly decreasing option
  if(evenness == 'skew.lin'){
    B0 <- (N_sp+1)/N_sp
    B1 <- -1/(2*N_sp)
    y <- B0 + B1*x
  }

  # skew.low
  if(evenness == 'skew.low'){
    Beta <- 1/(N_sp^0.5) # 0.8 = 0.125*10
    y <- 1/((N_sp*Beta) + x)
  }

  # skew.med
  if(evenness == 'skew.med'){
    Beta <- 1/(N_sp^1)# 0.01*length(x) # 0.1 = 0.01 * 10
    y <- 1/((N_sp*Beta) + x)
  }

  # skew.hi
  if(evenness == 'skew.hi'){
    Beta <- 0 # same as: y <- 1/(x)
    y <- 1/((N_sp*Beta) + x) # same as: y <- 1/(x)
  }
  
  pY <- p(y)

  if(stochastic){
    return(rmultinom(n = 1, size = N_out, prob = pY)[,1])
  }else{
    return(pY * N_out)
  }
  
}