

#' Piecewise Exponential Product Partition Model
#' @export
#' @aliases peppm
#' @param time vector of observed failure times.
#' @param status vector of failure indicators
#' @param a_rates shape parameter of the gamma distribution (prior for failure rates).
#' @param b_rates scale parameter of the gamma distribution (prior for failure rates).
#' @param cohesion type of prior cohesion (1 to 4).
#' @param a_beta shape1 parameter of the beta distribution (prior for p - cohesion 4).
#' @param b_beta shape2 parameter of the beta distribution (prior for p - cohesion 4).
#' @param nburnin number of iterations to be discarded.
#' @param npost desired posterior sample size
#' @param nlag number of jumps to eliminate autocorrelation of the chain.
#' @return Posterior sample of the number of intervals, failure rates, the auxiliary vector U, and the logarithm of the prior predictive distribution (log data factor).
#' @examples 
#' # Small chain used here due to time constraints. 
#' data(telecom)
#' 
#' # Prior cohesion 1:
#' fit1 <- with(telecom, peppm(time, status, cohesion=1, nburnin = 0, nlag = 1, npost = 100))
#' # Prior cohesion 2:
#' fit2 <- with(telecom, peppm(time, status, cohesion=2, nburnin = 0, nlag = 1, npost = 100))
#' # Prior cohesion 3:
#' fit3 <- with(telecom, peppm(time, status, cohesion=3, nburnin = 0, nlag = 1, npost = 100))
#' # Prior cohesion 4:
#' fit4 <- with(telecom, peppm(time, status, cohesion=4, nburnin = 0, nlag = 1, npost = 100))  
#'    

peppm <- function(time, status, a_rates=1, b_rates=1, cohesion=1, 
                  a_beta=1, b_beta=1, nburnin=10000, npost=20000, nlag=10){
  
  nsim <- nburnin + nlag*npost
  ftgrid <- unique(time[status==1])
  m <- length(ftgrid)
  U0 <- rep(0, m-1)
  samp <- gibbs(U0=U0, ftgrid=ftgrid, time=time, status=status, 
                                       a_rates=a_rates, b_rates=b_rates, 
                                       cohesion=cohesion, a_beta=a_beta, 
                                       b_beta=b_beta, npost=npost, 
                                       nburnin=nburnin, nlag=nlag)
  output <- list(b=samp$b, rates=samp$rates, U=samp$U, lpred=samp$lpred, ftgrid=ftgrid)
  class(output) <- "peppm"
  return(output)
}


