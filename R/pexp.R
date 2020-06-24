

#' Probability function, distribution function, quantile function and random generation for the Piecewise Exponential (PE) distribution.
#' @name PExp
#' @aliases PExp
#' @param x vector of time points.
#' @param q	 vector of quantiles.
#' @param p	 vector of probabilities.
#' @param n	 number of random values to return.
#' @param tgrid vector of time grid knots.
#' @param rates vector of failure rates.
#' @param log,log.p	 logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail	logical; if TRUE (default), probabilities are \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#' @return dpexp gives the (log) probability function, ppexp gives the (log) distribution function, qpexp gives the quantile function, and rpexp generates random deviates.
#' @description Probability function, distribution function, quantile function and random generation for the Piecewise Exponential (PE) distribution.
#' 
#' @examples 
#' n <- 10
#' tgrid <- c(0, 1, 3, 7, Inf)
#' rates <- c(0.5, 4, 0.8, 0.1)
#' x <- sort(rpexp(n, tgrid=tgrid, rates=rates))
#' Fx <- ppexp(x, tgrid, rates)
#' y <- qpexp(Fx, tgrid, rates)
#' # checking:
#' x==y

#' @rdname PExp
#' @export
dpexp <- function(x, tgrid, rates, log=FALSE){
  Ht <- Hpexp(x, tgrid, rates)
  id <- as.numeric(cut(x, tgrid, include.lowest = TRUE))
  log_ft <- log(rates[id]) - Ht
  if(log==FALSE)
    return(exp(log_ft))
  else
    return(log_ft)
}

#' @rdname PExp
#' @export
#' 
ppexp <- function(q, tgrid, rates, lower.tail=TRUE, log.p=FALSE){
  Ht <- Hpexp(q, tgrid, rates)
  St <- exp(-Ht)

  if(lower.tail==TRUE)
  {
    if(log.p==FALSE)
      return(1-St)
    else
      return(log(1-St))
  }else{
    if(log.p==FALSE)
      return(St)
    else
      return(log(St))
  }
}


#' @rdname PExp
#' @export
#'
#'  
qpexp <- function(p, tgrid, rates, lower.tail=TRUE, log.p=FALSE){
  if(lower.tail==TRUE){
    u <- p
  }else{
    u <- 1-p
  }

  h <- diff(tgrid)
  n <- length(p)
  area <- c(0,cumsum(h*rates))
  Ft <- 1-exp(-area)
  id <- as.numeric(cut(u, unique(Ft), include.lowest = TRUE))

  if(n==1){
    if(id==1){
      q <- -log(1-u)/rates[1]
    }else{
      q <- tgrid[id] - (area[id] + log(1-u))/rates[id]
    }
  }else{
    q  <- rep(0, n)
    q[which(id==1)] <- -log(1-u[which(id==1)])/rates[1]
    aux <- which(id>1)
    q[aux] <- tgrid[id[aux]] - (area[id[aux]] + log(1-u[aux]))/rates[id[aux]]
  }
  return(q)
}


#' @rdname PExp
#' @export
#' 
rpexp <- function(n, tgrid, rates){
  u <- stats::runif(n)
  x <- qpexp(u, tgrid, rates)
  return(x)
}




#---------------------------------------------

#' Hazard and cumulative hazard functions of the PE distribution
#' @name pehaz
#' @param x vector of time points.
#' @param tgrid vector of time grid knots.
#' @param rates vector of failure rates.
#' @return hpexp gives the hazard function and Hpexp gives the cumulative hazard function of the PE distribution.
#' 

#' @rdname pehaz
#' @export
#' 
hpexp <- function(x, tgrid, rates){
  id <- as.numeric(cut(x, tgrid, include.lowest = TRUE))
  return(rates[id])
}


#' @rdname pehaz
#' @export
#'
Hpexp <- function(x, tgrid, rates){
  h <- diff(tgrid)
  n <- length(x)
  id <- as.numeric(cut(x, tgrid, include.lowest = TRUE))
  area <- c(0, cumsum(h*rates))
  Ht <- (x-tgrid[id])*rates[id] + area[id]
  return(Ht)
}


#----------------------------------

#' Time grid for the PE distribution
#' @aliases timeGrid
#' @export
#' @param time Vector of failure times
#' @param status Vector of failure indicators
#' @param n.int Optional. Number of intervals. If \code{NULL}, the number of intervals is
#' set to be equal to the number of distinct observed failure times.
#' @description This function make use of the observed times and failure indicators to create a time grid for the PE distribution.
#' @return the time grid needed to specify the PE distribution.
#' 
#' @examples 
#' data(telecom)
#' tgrid1 <- with(telecom, timeGrid(time, status))
#' tgrid1
#' tgrid2 <- with(telecom, timeGrid(time, status, n.int = 4))
#' tgrid2

timeGrid <- function(time, status, n.int=NULL){
  o <- order(time)
  time <- time[o]
  status <- status[o]
  time.aux <- unique(time[status==1])
  if(is.null(n.int)){
    n.int <- length(time.aux)
  }

  m <- length(time.aux)
  if(n.int > m){
    tgrid <- c(0,unique(time[status==1]))
    tgrid[length(tgrid)] <- Inf
  }
  else{
    b <- min(m,n.int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    tgrid <- c(0,time.aux[idf])
    tgrid[length(tgrid)] <- Inf
  }
  return(tgrid)
}


#' Function to identify the times' intervals
#' @export
#' @aliases findInt
#' @param time vector of times.
#' @param tgrid time grid of the PE distribution.
#' @return indicator of times's intervals
#' 
#' @examples 
#' data(telecom)
#' tgrid <- with(telecom, timeGrid(time, status))
#' tgrid
#' findInt(telecom$time, tgrid) 

findInt <- function(time, tgrid){
  id <- as.numeric(cut(time, tgrid, include.lowest = TRUE))
  return(id)
}
