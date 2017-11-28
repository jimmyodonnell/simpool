#' Generate primer efficiencies
#' 
#' Given an indication of mismatch variation ('mmv'), return some representative
#' values of primer efficiency. 
#' The idea here is to simulate primer efficiency for some number N of templates 
#' which may vary in how well they are amplified by a given primer set. 
#' A value of 1 indicates perfect efficiency, while 0 indicates the template 
#' would not be amplified at all (complete 'dropout'). 
#' 
#' The 'mmv' argument can accept a number of options. First, it will accept a 
#' string matching one of the preformatted options. 
#' These options control the two parameters of the beta distribution from 
#' which values are drawn. 
#' If the option 'custom' is used for mmv, the values of mmv.beta1 and mmv.beta2 
#' will be used as the parameters for a beta distribution from which values are 
#' drawn.
#' Otherwise, it requires a numeric vector. 
#' If a numeric vector is used for mmv, the output values are sampled from it. 
#' If the vector is shorter than N, values are sampled with replacement.
#' 
#' @param N Integer. Number of values to generate.
#' 
#' @param mmv String or numeric indicating distribution of mismatch variation. 
#' See details.
#'
#' @param mmv.beta Numeric; length = 2; >0. Parameters of beta distribution.
#'
#' @return Numeric vector of length N with values in {0,1}.
#'
#' @examples
#' primer_eff(10, 'low')
#' 
#' @export
primer_eff <- function(N, mmv, mmv.beta = c(1,1)){
  mmv.opts <- c('absent', 'low', 'medium', 'high', 'custom')
  if(mmv[1] %in% mmv.opts){
    if(mmv  == 'absent'){
      eff_out <- rep(1, length = N)
    }
    if(mmv  == 'low'){
      eff_out <- rbeta(n = N, shape1 = 10, shape2 = 1)
    }
    if(mmv  == 'medium'){
      eff_out <- rbeta(n = N, shape1 = 2, shape2 = 1)
    }
    if(mmv  == 'high'){
      eff_out <- rbeta(n = N, shape1 = 1, shape2 = 1)
    }
    if(mmv == 'custom'){
      eff_out <- rbeta(n = N, shape1 = mmv.beta[1], shape2 = mmv.beta[2])
    }
  }else{
    if(all((0 <= mmv) & (mmv <= 1))){
      tooshort <- length(mmv) < N
      eff_out <- sample(mmv, size = N, replace = tooshort)
    }else{
      mmv.opts.msg <- paste(mmv.opts, collapse = ', ')
      msg <- 'value of mmv must be numeric with values in {0,1}, or one of:'
      stop(paste(msg, mmv.opts.msg))
    }
  }
  return(eff_out)
}
